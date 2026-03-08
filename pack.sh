#!/bin/bash

# 1. 定义提交目录名称 (格式：Submission_姓名_学号)
SUB_DIR="Submission_GongZhiyuan_2330416033"

# 2. 创建目录结构 (实验1-6)
mkdir -p $SUB_DIR/Exp1_Simulation
mkdir -p $SUB_DIR/Exp2_Assembly
mkdir -p $SUB_DIR/Exp3_Homology
mkdir -p $SUB_DIR/Exp4_DeNovo
mkdir -p $SUB_DIR/Exp5_DNA_Element
mkdir -p $SUB_DIR/Exp6_Visualization

echo "正在整理文件..."

# === 实验1: 测序模拟 ===
# 复制统计表和PDF图 (假设在当前目录，如果没有会跳过)
cp data.txt $SUB_DIR/Exp1_Simulation/ 2>/dev/null
cp *.pdf $SUB_DIR/Exp1_Simulation/ 2>/dev/null
# 将命令行记录写入文本文件 (因为是命令行操作的)
echo "art_illumina commands used in report." > $SUB_DIR/Exp1_Simulation/commands_log.txt

# === 实验2: 序列组装 ===
# 复制脚本和配置文件
cp 02_Assembly/lib.cfg $SUB_DIR/Exp2_Assembly/
cp 02_Assembly/*.sh $SUB_DIR/Exp2_Assembly/
# 复制 FastQC HTML 报告 (排除 zip)
cp 02_Assembly/fastqc_out/*.html $SUB_DIR/Exp2_Assembly/
# 复制 QUAST 最终报告
cp 02_Assembly/quast_cmp/report.html $SUB_DIR/Exp2_Assembly/quast_report.html
cp 02_Assembly/quast_cmp/report.pdf $SUB_DIR/Exp2_Assembly/quast_report.pdf
# (可选) 复制最终的 Scaffold 序列文件
cp 02_Assembly/soap_out/*.scafSeq $SUB_DIR/Exp2_Assembly/

# === 实验3: 同源搜索 ===
# 复制 GFF3 结果
cp 03_Homology_Search/Sc_perl_modified.gff3 $SUB_DIR/Exp3_Homology/
cp 03_Homology_Search/Sc_python_modified.gff3 $SUB_DIR/Exp3_Homology/
# 复制评估 Stats
cp 03_Homology_Search/*.stats $SUB_DIR/Exp3_Homology/
# 复制 bamstats 图片
cp 03_Homology_Search/plot-bamstats_out/*.png $SUB_DIR/Exp3_Homology/

# === 实验4: 从头预测 ===
# 复制 GFF3 和 提取的序列
cp 04_DeNovo_Prediction/augustus_out.gff3 $SUB_DIR/Exp4_DeNovo/
cp 04_DeNovo_Prediction/*.fa $SUB_DIR/Exp4_DeNovo/
# 复制 BLASTP 验证表和评估结果
cp 04_DeNovo_Prediction/blastp.outfmt6 $SUB_DIR/Exp4_DeNovo/
cp 04_DeNovo_Prediction/*.stats $SUB_DIR/Exp4_DeNovo/

# === 实验5: DNA元件预测 (重点是代码和图) ===
# 复制 Python 和 R 脚本
cp 05_DNA_Element/*.py $SUB_DIR/Exp5_DNA_Element/
cp 05_DNA_Element/*.R $SUB_DIR/Exp5_DNA_Element/
# 复制 图片
cp 05_DNA_Element/*.png $SUB_DIR/Exp5_DNA_Element/
# 复制 结果 GFF3 和 txt
cp 05_DNA_Element/*.gff3 $SUB_DIR/Exp5_DNA_Element/
cp 05_DNA_Element/score_distance.txt $SUB_DIR/Exp5_DNA_Element/

# === 实验6: 可视化 ===
# 复制 JBrowse 的配置文件 (如果存在)
if [ -f "06_Visualization/jbrowse/data/test/tracks.conf" ]; then
    cp 06_Visualization/jbrowse/data/test/tracks.conf $SUB_DIR/Exp6_Visualization/
fi
# 复制打包好的 jbrowse 数据包 (如果有)
cp 06_Visualization/jbrowse_data.zip $SUB_DIR/Exp6_Visualization/ 2>/dev/null

echo "文件整理完成，正在压缩..."

# 3. 打包为 ZIP
zip -r ${SUB_DIR}.zip $SUB_DIR

echo "======================================================"
echo "打包完成！文件名: ${SUB_DIR}.zip"
echo "请将此文件下载到本地，连同你的 Word/PDF 报告一起提交。"
echo "======================================================"
