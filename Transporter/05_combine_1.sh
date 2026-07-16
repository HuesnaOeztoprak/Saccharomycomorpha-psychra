#✳️
#!/usr/bin/env bash

set -euo pipefail

ACC_WHITELIST="pfam_accessions_only.txt"   # PF00005, PF00083 etc. — one per line

if [[ ! -f "$ACC_WHITELIST" ]]; then
    echo "ERROR: $ACC_WHITELIST not found. Run get_transporter_pfam.py first."
    exit 1
fi
echo "Accession whitelist: $(wc -l < $ACC_WHITELIST) entries"

# ── 0. Build short-name whitelist from actual data ────────────────────────
# Extract pfam_shortname <-> pfam_accession pairs from pfam.tbl,
# keep only those whose accession is in the transporter whitelist.
grep -v '^#' pfam.tbl \
  | awk '{name=$1; acc=$2; sub(/\.[0-9]+$/,"",acc); print name"\t"acc}' \
  | sort | uniq \
  | awk 'NR==FNR{wl[$1]=1; next} $2 in wl {print $1}' \
      "$ACC_WHITELIST" - \
  > transporter_pfam_shortnames.txt

echo "Short-name whitelist derived from data: $(wc -l < transporter_pfam_shortnames.txt) names"
echo "Example names:"; head -5 transporter_pfam_shortnames.txt


# ── 1. BLAST: best hit per query ──────────────────────────────────────────
grep -v '^#' tcdb_hits.tsv \
  | sort -k1,1 -k11,11g \
  | awk '!seen[$1]++ {
      split($2, a, "|")
      print $1"\t"a[4]"\t"a[3]"\t"$11"\t"$3
    }' \
  > blast_best.txt
echo "BLAST best hits: $(wc -l < blast_best.txt)"


# ── 2. EggNOG TC ──────────────────────────────────────────────────────────
grep -v '^#' eggnog_TC.txt \
  | awk '{split($2,a,","); print $1"\t"a[1]}' \
  > eggnog_TC_clean.txt
echo "EggNOG TC entries: $(wc -l < eggnog_TC_clean.txt)"


# ── 3. EggNOG Pfam: filter by short-name whitelist ────────────────────────
# EggNOG Pfam column uses short names (ABC_tran, MatE, MFS_1 etc.)
grep -v '^#' eggnog_TC+Pfam.txt \
  | awk '{print $1"\t"$2}' \
  | awk -v wnfile="transporter_pfam_shortnames.txt" '
      BEGIN { while ((getline line < wnfile) > 0) {
          gsub(/[ \t\r\n]+/, "", line); wl[line]=1
        }
      }
      {
        n = split($2, domains, ",")
        hits = ""
        for (i=1; i<=n; i++) {
          d = domains[i]
          gsub(/^[ \t]+|[ \t]+$/, "", d)
          if (d in wl) hits = (hits=="" ? d : hits","d)
        }
        if (hits != "") print $1"\t"hits
      }' \
  > eggnog_pfams_filtered.txt
echo "EggNOG Pfam filtered: $(wc -l < eggnog_pfams_filtered.txt) genes"


# ── 4. HMMER Pfam: col1=name, col2=accession(versioned), col3=gene, col5=evalue
grep -v '^#' pfam.tbl \
  | awk -v wnfile="transporter_pfam_shortnames.txt" '
      BEGIN { while ((getline line < wnfile) > 0) {
          gsub(/[ \t\r\n]+/, "", line); wl[line]=1
        }
      }
      {
        pfam_name = $1
        gene      = $3
        evalue    = $5+0
        if (evalue <= 1e-5 && pfam_name in wl) {
          print gene"\t"pfam_name
        }
      }' \
  | sort -k1,1 -k2,2 | uniq \
  | awk '{
      if ($1==prev) val = val","$2
      else { if (prev!="") print prev"\t"val; prev=$1; val=$2 }
    } END { if (prev!="") print prev"\t"val }' \
  > hmmer_pfams_filtered.txt
echo "HMMER Pfam filtered: $(wc -l < hmmer_pfams_filtered.txt) genes"


# ── 5. DeepTMHMM ──────────────────────────────────────────────────────────
cp deeptmhmm_tm.txt deeptmhmm_clean.txt


# ── 6. Merge into master table ────────────────────────────────────────────
awk -v OFS="\t" '
FNR==NR {
  blastTC[$1]=$2; blastUP[$1]=$3; blastEV[$1]=$4; blastID[$1]=$5
  blast[$1]=1; next
}
FILENAME=="eggnog_TC_clean.txt"       { eggTC[$1]=$2; next }
FILENAME=="eggnog_pfams_filtered.txt" { eggPF[$1]=$2; next }
FILENAME=="hmmer_pfams_filtered.txt"  { hmPF[$1]=$2;  next }
FILENAME=="deeptmhmm_clean.txt"       { tm[$1]=$2;    next }

END {
  for (q in blast) transporter[q]=1
  for (q in eggTC) transporter[q]=1
  for (q in eggPF) transporter[q]=1
  for (q in hmPF)  transporter[q]=1

  print "gene_id","best_TCDB_id","blast_UniProt","blast_evalue","blast_pident",\
        "eggnog_TC","eggnog_PFAM","hmmer_PFAM","nTM","n_evidence","confidence"

  for (q in transporter) {
    n = 0
    if (q in blast) n++
    if (q in eggTC) n++
    if (q in eggPF) n++
    if (q in hmPF)  n++
    tm_count = (q in tm ? tm[q]+0 : 0)
    if (tm_count >= 2 && n >= 1) n++

    if      (n >= 3) conf = "HIGH"
    else if (n == 2) conf = "MEDIUM"
    else             conf = "LOW"

    if (q in blast)      best_tc = blastTC[q]
    else if (q in eggTC) best_tc = eggTC[q]
    else                 best_tc = "Pfam_only"

    print q,
      best_tc,
      (q in blastUP ? blastUP[q] : "-"),
      (q in blastEV ? blastEV[q] : "-"),
      (q in blastID ? blastID[q] : "-"),
      (q in eggTC   ? eggTC[q]   : "-"),
      (q in eggPF   ? eggPF[q]   : "-"),
      (q in hmPF    ? hmPF[q]    : "-"),
      tm_count, n, conf
  }
}
' blast_best.txt \
  eggnog_TC_clean.txt \
  eggnog_pfams_filtered.txt \
  hmmer_pfams_filtered.txt \
  deeptmhmm_clean.txt \
| sort -t$'\t' -k11,11 -k1,1 \
> transporter_master.tsv

echo ""
echo "=== Done. Summary ==="
echo "  Total proteins in master table: $(tail -n +2 transporter_master.tsv | wc -l)"
echo "  HIGH confidence:   $(awk -F'\t' '$11=="HIGH"'   transporter_master.tsv | wc -l)"
echo "  MEDIUM confidence: $(awk -F'\t' '$11=="MEDIUM"' transporter_master.tsv | wc -l)"
echo "  LOW confidence:    $(awk -F'\t' '$11=="LOW"'    transporter_master.tsv | wc -l)"
echo ""
echo "  Pfam filtering sanity check:"
echo "    eggnog_pfams_filtered : $(wc -l < eggnog_pfams_filtered.txt) genes"
echo "    hmmer_pfams_filtered  : $(wc -l < hmmer_pfams_filtered.txt) genes"
