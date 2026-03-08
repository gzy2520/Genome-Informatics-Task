import sys

# 1. 读取基因注释 (gff.gff)
# 我们只关心 "gene" 特征，且需要区分正负链
genes = {} # 结构: {'chrom': {'+': [start1, start2...], '-': [end1, end2...]}}

print("Loading gene annotations...")
with open("../gff.gff", 'r') as f:
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split('\t')
        if len(parts) < 9: continue
        
        feat_type = parts[2]
        if feat_type != 'gene': continue # 只提取 gene [cite: 265]
        
        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        strand = parts[6]
        
        # 提取 GeneID (从属性列)
        attr = parts[8]
        gene_id = "unknown"
        if "ID=" in attr:
            gene_id = attr.split("ID=")[1].split(";")[0]
        elif "Name=" in attr:
            gene_id = attr.split("Name=")[1].split(";")[0]

        if chrom not in genes:
            genes[chrom] = {'+': [], '-': []}
            
        # 存储格式: (位置, gene_id)
        # 对于正链，转录起始位点(TSS)约等于 start
        # 对于负链，转录起始位点(TSS)约等于 end
        if strand == '+':
            genes[chrom]['+'].append((start, gene_id))
        elif strand == '-':
            genes[chrom]['-'].append((end, gene_id))

# 对基因位置进行排序，以便二分查找或顺序遍历 [cite: 264]
for chrom in genes:
    genes[chrom]['+'].sort(key=lambda x: x[0])
    genes[chrom]['-'].sort(key=lambda x: x[0])

# 2. 读取 TATA-box 预测结果并计算距离
tata_file = "TATA_box_prediction.gff3"
out_gff = "TATA_annotated.gff3"
out_data = "score_distance.txt"

print("Annotating TATA boxes...")
with open(tata_file, 'r') as f, open(out_gff, 'w') as og, open(out_data, 'w') as od:
    od.write("Score\tDistance\n") # 写入表头
    og.write("##gff-version 3\n")
    
    for line in f:
        if line.startswith("#"): continue
        parts = line.strip().split('\t')
        
        chrom = parts[0]
        start = int(parts[3])
        end = int(parts[4])
        score = float(parts[5])
        strand = parts[6]
        attributes = parts[8]
        
        min_dist = 999999
        nearest_gene = "None"
        
        # 查找逻辑 [cite: 268-280]
        # 正链 TATA: 找 start > TATA_start 的最小 gene_start
        # 距离 = Gene_Start - TATA_Start
        if strand == '+' and chrom in genes:
            gene_list = genes[chrom]['+']
            # 简单的顺序遍历查找 (可以优化，但对于小基因组足够)
            for g_pos, g_id in gene_list:
                dist = g_pos - start
                # 我们寻找下游基因 (dist > -500, 允许在基因内部一点点或上游)
                # 通常 TATA 在 TSS 上游，所以 dist 应该是正数且较小
                if dist > -200: 
                    min_dist = dist
                    nearest_gene = g_id
                    break # 因为排过序，找到第一个就是最近的
                    
        # 负链 TATA: 找 end < TATA_end 的最大 gene_end
        # 距离 = TATA_End - Gene_End
        elif strand == '-' and chrom in genes:
            gene_list = genes[chrom]['-']
            # 倒序遍历，找第一个小于 TATA end 的
            for g_pos, g_id in reversed(gene_list):
                dist = end - g_pos
                if dist > -200:
                    min_dist = dist
                    nearest_gene = g_id
                    break

        # 过滤掉太远的 (例如距离超过 2000bp 可能就不是启动子了)
        if nearest_gene != "None" and min_dist < 5000:
            # 写入新 GFF
            new_attr = f"{attributes};Adjacent_GeneID={nearest_gene};distance={min_dist}"
            og.write(f"{chrom}\tPWM_Scan\tTATA_box\t{start}\t{end}\t{score}\t{strand}\t.\t{new_attr}\n")
            
            # 写入统计文件供 R 绘图
            od.write(f"{score}\t{min_dist}\n")

print("Analysis finished.")
