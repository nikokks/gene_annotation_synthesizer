# Split features into several tabular files from the xml parser output
# One file will contain genes metadata, one will contain the
# sequence, one the full names and one the terms.

infile="$1"
outdir="$2"
mkdir -p "${outdir}"

# Write features
echo -e "id\tevidence\tcitations\ttaxid" \
  | gzip -c \
  > "${outdir}/features.tsv.gz"
gzip -dc "${infile}" \
  | tail -n+2 \
  | cut -f1-4 \
  | gzip -c \
  >> "${outdir}/features.tsv.gz"

# Write sequences
echo -e "id\tseq" \
  | gzip -c \
  > "${outdir}/sequences.tsv.gz"
gzip -dc "${infile}" \
  | tail -n+2 \
  | cut -f1,9 \
  | gzip -c \
  >> "${outdir}/sequences.tsv.gz"

# Write names
echo -e "id\tname" \
  | gzip -c \
  > "${outdir}/names.tsv.gz"
gzip -dc "${infile}" \
  | tail -n+2 \
  | awk \
      -vFS='\t' \
      -vOFS='\t' \
      '{split($5,names, "|"); for(i in names) {print $1,names[i]}}' \
  | gzip -c \
  >> "${outdir}/names.tsv.gz"

# Write terms in long format
echo -e "id\term\tsource" \
  | gzip -c \
  > "${outdir}/terms.tsv.gz"
gzip -dc "${infile}" \
  | tail -n+2 \
  | awk \
      -vFS='\t' \
      -vOFS='\t' \
     '{split($6,ipr, "|");
     split($7,pfam, "|");
     split($8,go, "|");
     for(i in ipr) {print $1,ipr[i],"IPR"};
     for(i in pfam) {print $1,pfam[i],"PFAM"};
     for(i in go) {print $1,go[i],"GO"}}' \
  | gzip -c \
  >> "${outdir}/terms.tsv.gz"
