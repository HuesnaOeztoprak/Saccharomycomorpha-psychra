!/usr/bin/env bash
# 06_filter_mastertable.sh
# Filter and functionally annotate the new transporter_master.tsv
# Columns in master table:
#  1:gene_id  2:best_TCDB_id  3:blast_UniProt  4:blast_evalue  5:blast_pident
#  6:eggnog_TC  7:eggnog_PFAM  8:hmmer_PFAM  9:nTM  10:n_evidence  11:confidence

set -euo pipefail

INPUT="transporter_master.tsv"
OUTPUT="transporter_final.tsv"

echo "=== Step 1: Filtering and annotating ==="

awk -F"\t" '
BEGIN {
    OFS="\t"
    print "gene_id","best_TCDB_id","blast_UniProt","blast_evalue","blast_pident",\
          "eggnog_TC","eggnog_PFAM","hmmer_PFAM","nTM","n_evidence","confidence",\
          "tc_class","tc_subclass","transport_mechanism","pfam_family_group"
}

NR==1 { next }  # skip header

{
    gene=$1; tcdb=$2; uniprot=$3; evalue=$4; pident=$5
    egg_tc=$6; egg_pf=$7; hmm_pf=$8; nTM=$9+0; n_ev=$10+0; conf=$11

    # ── Decide best TC source ──────────────────────────────────────
    tc = (tcdb != "-" && tcdb != "Pfam_only") ? tcdb : \
         (egg_tc != "-" ? egg_tc : "-")

    # ── Parse TC number hierarchy ──────────────────────────────────
    # TC format: class.subclass.family.subfamily.TC#
    # e.g. 2.A.1.1.1 → class=2, subclass=A, family=1
    tc_class   = "-"
    tc_subclass = "-"
    if (tc != "-") {
        n = split(tc, parts, ".")
        if (n >= 1) tc_class    = parts[1]
        if (n >= 2) tc_subclass = parts[1]"."parts[2]
    }

    # ── Transport mechanism from TC class ─────────────────────────
    mech = "unknown"
    if      (tc_class == "1") mech = "channel_pore"
    else if (tc_class == "2") mech = "secondary_carrier"
    else if (tc_class == "3") mech = "primary_active_pump"
    else if (tc_class == "4") mech = "group_translocator"
    else if (tc_class == "5") mech = "transmembrane_electron_carrier"
    else if (tc_class == "8") mech = "accessory_factor"
    else if (tc_class == "9") mech = "incompletely_characterised"

    # ── Pfam-based family annotation ──────────────────────────────
    # Check both EggNOG Pfam and HMMER Pfam columns
    combined_pf = egg_pf","hmm_pf

    family = "other"
    if      (combined_pf ~ /MFS/)                              family = "MFS"
    else if (combined_pf ~ /ABC_tran|ABC_membrane/)            family = "ABC"
    else if (combined_pf ~ /Sugar_tr/)                         family = "Sugar_transporter"
    else if (combined_pf ~ /Aa_trans|Trp_Tyr_perm|GABA_trans/) family = "Amino_acid_transporter"
    else if (combined_pf ~ /Ammonium_transp|Amt/)              family = "Ammonium_transporter"
    else if (combined_pf ~ /Ion_trans|Na_K-ATPase/)            family = "Ion_transporter"
    else if (combined_pf ~ /MatE/)                             family = "MATE_efflux"
    else if (combined_pf ~ /P2X_receptor/)                     family = "P2X_ion_channel"
    else if (combined_pf ~ /Mito_carr/)                        family = "Mitochondrial_carrier"
    else if (combined_pf ~ /HAD|ATPase/)                       family = "ATPase_pump"
    else if (combined_pf ~ /Rft-1|Rft/)                        family = "Dolichol_transporter"
    else if (combined_pf ~ /SRP/)                              family = "SRP_related"
    else if (tc_subclass != "-") family = "TC_"tc_subclass"_no_pfam"

    # ── Apply final quality filters ───────────────────────────────
    # Include: any protein with at least 1 transporter-specific evidence
    # (conf is already set correctly by 05_combine script)
    # Exclude: n_evidence == 1 AND nTM == 0 (purely in silico, no membrane evidence)
    if (n_ev == 1 && nTM == 0) next

    print gene, tcdb, uniprot, evalue, pident, \
          egg_tc, egg_pf, hmm_pf, nTM, n_ev, conf, \
          tc_class, tc_subclass, mech, family
}
' "$INPUT" > "$OUTPUT"

echo "Done."
echo ""
echo "=== Summary ==="
echo "Total proteins in final table: $(tail -n +2 "$OUTPUT" | wc -l)"
echo ""
echo "--- Confidence breakdown ---"
awk -F"\t" 'NR>1 {print $11}' "$OUTPUT" | sort | uniq -c | sort -rn

echo ""
echo "--- Transport mechanism (TC class) ---"
awk -F"\t" 'NR>1 {print $14}' "$OUTPUT" | sort | uniq -c | sort -rn

echo ""
echo "--- Pfam family group ---"
awk -F"\t" 'NR>1 {print $15}' "$OUTPUT" | sort | uniq -c | sort -rn

echo ""
echo "--- TC subclass distribution (top 20) ---"
awk -F"\t" 'NR>1 && $13!="-" {print $13}' "$OUTPUT" \
  | sort | uniq -c | sort -rn | head -20


# ── Optional: write per-class subset files ────────────────────────────────
echo ""
echo "=== Writing per-class subset files ==="

for class_num in 1 2 3 4 9; do
    outfile="transporters_TC${class_num}.tsv"
    head -1 "$OUTPUT" > "$outfile"
    awk -F"\t" -v c="$class_num" 'NR>1 && $12==c' "$OUTPUT" >> "$outfile"
    count=$(tail -n +2 "$outfile" | wc -l)
    echo "  TC class $class_num → $count proteins → $outfile"
done

echo ""
echo "All done. Main output: $OUTPUT"

