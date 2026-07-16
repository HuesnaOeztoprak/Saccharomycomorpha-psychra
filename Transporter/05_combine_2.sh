✳️
usr/bin/env bash
# patch_and_rerun.sh
# 1. Adds missing MFS Pfam names to shortname whitelist
# 2. Reruns 05_combine_1.sh to pick up the new hits

set -euo pipefail

# ── Step 1: Patch the shortname whitelist ─────────────────────────────────
# MFS_1, MFS_3, MFS_5 are in your pfam.tbl but were missing from the whitelist
# because the keyword scan missed their parent Pfam accessions
echo "Patching transporter_pfam_shortnames.txt..."

for name in MFS_1 MFS_3 MFS_5; do
    if ! grep -qx "$name" transporter_pfam_shortnames.txt; then
        echo "$name" >> transporter_pfam_shortnames.txt
        echo "  Added: $name"
    else
        echo "  Already present: $name"
    fi
done

echo "Updated whitelist size: $(wc -l < transporter_pfam_shortnames.txt) names"


# ── Step 2: Rerun combine to pick up new HMMER hits ───────────────────────
echo ""
echo "Rerunning 05_combine_final.sh..."
bash 05_combine_final.sh


# ── Step 3: Run updated 06 with expanded family annotation ────────────────
echo ""
echo "Running updated filter/annotation..."

INPUT="transporter_master.tsv"
OUTPUT="transporter_final.tsv"

awk -F"\t" '
BEGIN {
    OFS="\t"
    print "gene_id","best_TCDB_id","blast_UniProt","blast_evalue","blast_pident",\
          "eggnog_TC","eggnog_PFAM","hmmer_PFAM","nTM","n_evidence","confidence",\
          "tc_class","tc_subclass","transport_mechanism","pfam_family_group"
}

NR==1 { next }

{
    gene=$1; tcdb=$2; uniprot=$3; evalue=$4; pident=$5
    egg_tc=$6; egg_pf=$7; hmm_pf=$8; nTM=$9+0; n_ev=$10+0; conf=$11

    tc = (tcdb != "-" && tcdb != "Pfam_only") ? tcdb : (egg_tc != "-" ? egg_tc : "-")

    tc_class    = "-"
    tc_subclass = "-"
    if (tc != "-") {
        n = split(tc, parts, ".")
        if (n >= 1) tc_class    = parts[1]
        if (n >= 2) tc_subclass = parts[1]"."parts[2]
    }

    mech = "unknown"
    if      (tc_class == "1") mech = "channel_pore"
    else if (tc_class == "2") mech = "secondary_carrier"
    else if (tc_class == "3") mech = "primary_active_pump"
    else if (tc_class == "4") mech = "group_translocator"
    else if (tc_class == "5") mech = "transmembrane_electron_carrier"
    else if (tc_class == "8") mech = "accessory_factor"
    else if (tc_class == "9") mech = "incompletely_characterised"

    combined_pf = egg_pf","hmm_pf

    # ── Family annotation — ordered by specificity ──────────────────────
    family = "other"

    # ABC superfamily
    if      (combined_pf ~ /ABC_tran|ABC_membrane|ABC2_membrane/)
                                                    family = "ABC"

    # MFS superfamily — all Pfam domain names
    else if (combined_pf ~ /MFS_1|MFS_2|MFS_3|MFS_4|MFS_5|MFS_1_like|MFS_MOT1|MFS_Mycoplasma/)
                                                    family = "MFS"

    # Sugar transporters (many are MFS but annotated separately)
    else if (combined_pf ~ /Sugar_tr/)              family = "Sugar_transporter"

    # Amino acid transporters
    else if (combined_pf ~ /Aa_trans|AA_permease|AA_permease_2|AA_permease_C|Trp_Tyr_perm|GABA_trans|UAA/)
                                                    family = "Amino_acid_transporter"

    # Ammonium transporters
    else if (combined_pf ~ /Ammonium_transp|Amt/)   family = "Ammonium_transporter"

    # Mitochondrial carriers
    else if (combined_pf ~ /Mito_carr/)             family = "Mitochondrial_carrier"

    # Ion transporters / cation efflux
    else if (combined_pf ~ /Ion_trans|Na_K-ATPase|Cation_efflux|Zip|ZT_dimer/)
                                                    family = "Ion_transporter"

    # ATPase pumps
    else if (combined_pf ~ /HAD|ATPase/)            family = "ATPase_pump"

    # Nucleotide sugar transporters
    else if (combined_pf ~ /Nuc_sug_transp|SLC35/)  family = "Nucleotide_sugar_transporter"

    # Drug/metabolite exporters (EamA = drug/metabolite transporter superfamily)
    else if (combined_pf ~ /EamA|MatE/)             family = "Drug_metabolite_transporter"

    # Sulfate transporters
    else if (combined_pf ~ /Sulfate_transp/)        family = "Sulfate_transporter"

    # Choline/ethanolamine transporters
    else if (combined_pf ~ /Choline_transpo/)       family = "Choline_transporter"

    # Xanthine/uracil permeases (nucleobase transporters)
    else if (combined_pf ~ /Xan_ur_permease/)       family = "Nucleobase_transporter"

    # Zinc transporters (Zip family)
    else if (combined_pf ~ /SBF|SBF_like/)          family = "SBF_transporter"

    # Citrate/malate transporters
    else if (combined_pf ~ /CitMHS/)                family = "Citrate_MHS_transporter"

    # Mitochondrial pyruvate carrier
    else if (combined_pf ~ /MPC/)                   family = "Mitochondrial_pyruvate_carrier"

    # Plastid/chloroplast transporters (TPT = triose phosphate translocator)
    else if (combined_pf ~ /TPT/)                   family = "Plastid_transporter"

    # CRT-like (chloroquine resistance transporter family)
    else if (combined_pf ~ /CRT-like/)              family = "CRT_like_transporter"

    # MATE efflux
    else if (combined_pf ~ /MatE/)                  family = "MATE_efflux"

    # P2X ion channels
    else if (combined_pf ~ /P2X_receptor/)          family = "P2X_ion_channel"

    # Dolichol-related
    else if (combined_pf ~ /Rft-1|Rft/)             family = "Dolichol_transporter"

    # Fall back to TC subclass if no Pfam family resolved
    else if (tc_subclass != "-")                    family = "TC_"tc_subclass"_no_pfam"

    # Quality filter: exclude n_evidence==1 with no TM support
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
echo "--- Pfam family group (all) ---"
awk -F"\t" 'NR>1 {print $15}' "$OUTPUT" | sort | uniq -c | sort -rn

echo ""
echo "--- TC subclass distribution ---"
awk -F"\t" 'NR>1 && $13!="-" {print $13}' "$OUTPUT" \
  | sort | uniq -c | sort -rn | head -20

echo ""
echo "--- Remaining unresolved (TC_X.X_no_pfam) ---"
awk -F"\t" 'NR>1 && $15 ~ /no_pfam/ {print $15}' "$OUTPUT" \
  | sort | uniq -c | sort -rn

echo ""
echo "=== Writing per-class subset files ==="
for class_num in 1 2 3 4 9; do
    outfile="transporters_TC${class_num}.tsv"
    head -1 "$OUTPUT" > "$outfile"
    awk -F"\t" -v c="$class_num" 'NR>1 && $12==c' "$OUTPUT" >> "$outfile"
    echo "  TC class $class_num → $(tail -n +2 "$outfile" | wc -l) proteins → $outfile"
done

echo ""
echo "All done. Main output: $OUTPUT"
