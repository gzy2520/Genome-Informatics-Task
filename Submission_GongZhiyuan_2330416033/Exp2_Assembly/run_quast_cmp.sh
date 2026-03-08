# ~/Geneinfo/run_quast_cmp.sh
set -euo pipefail

# ===== 路径：按需修改 =====
BASE=~/Geneinfo
REF_FA=$BASE/fna.fna      # <- 改成你的参考fasta
REF_GFF=$BASE/gff.gff    # <- 改成你的GFF/GTF（可去掉 -g 选项）
SOAP=$BASE/soap_out
OUT=$BASE/quast_cmp
THREADS=8

mkdir -p "$OUT"

# ===== 样本：按你的文件名就好 =====
C3="$SOAP/3_soap_out.contig"
C10="$SOAP/10_soap_out.contig"
C33="$SOAP/33_soap_out.contig"
S3="$SOAP/3_soap_out.scafSeq"
S10="$SOAP/10_soap_out.scafSeq"
S33="$SOAP/33_soap_out.scafSeq"

# 简单存在性检查（缺哪个就报错）
for f in "$C3" "$C10" "$C33" "$S3" "$S10" "$S33"; do
  [[ -s "$f" ]] || { echo "[ERR] missing: $f"; exit 1; }
done

# ===== 跑 QUAST（一次性多输入 + 自定义标签）=====
quast.py -o "$OUT" -t "$THREADS" -r "$REF_FA" -g "$REF_GFF" \
  -l "3x_contig,10x_contig,33x_contig,3x_scaffold,10x_scaffold,33x_scaffold" \
  "$C3" "$C10" "$C33" "$S3" "$S10" "$S33"

echo "[OK] QUAST done -> $OUT"

