import math
import random
import sys

# 1. 定义 TATA-box 的频率矩阵 (Frequency Matrix) [cite: 79]
# 数据来源：教程表 1-4
# 为了防止 log(0)，我们添加一个极小的伪计数 0.001
freq_matrix = [
    {'A': 31.6, 'C': 24.5, 'G': 15.3, 'T': 28.6}, # Pos 1
    {'A': 16.3, 'C': 60.2, 'G': 10.2, 'T': 13.3}, # Pos 2
    {'A':  2.0, 'C':  3.0, 'G':  0.0, 'T': 94.9}, # Pos 3
    {'A': 90.8, 'C':  2.1, 'G':  2.0, 'T':  5.1}, # Pos 4
    {'A':  0.0, 'C':  0.0, 'G':  1.0, 'T': 99.0}, # Pos 5
    {'A': 94.9, 'C':  0.0, 'G':  0.0, 'T':  5.1}, # Pos 6
    {'A': 57.1, 'C':  0.0, 'G':  0.0, 'T': 42.9}, # Pos 7
    {'A': 100.0,'C':  0.0, 'G':  0.0, 'T':  0.0}, # Pos 8
    {'A': 27.6, 'C':  0.0, 'G':  2.0, 'T': 70.4}, # Pos 9
    {'A': 69.4, 'C':  3.1, 'G': 13.3, 'T': 14.3}, # Pos 10
    {'A': 11.2, 'C': 39.8, 'G': 37.8, 'T': 11.2}, # Pos 11
    {'A': 24.5, 'C': 52.0, 'G': 21.4, 'T':  2.1}  # Pos 12
]

# 背景概率 (假设基因组 GC含量均等)
bg = 0.25
pseudo = 0.001

# 2. 构建 PWM (Position Weight Matrix) - Log-Odds Score
pwm = []
for pos in freq_matrix:
    pwm_row = {}
    total = sum(pos.values())
    for base, percent in pos.items():
        # 计算概率 p
        p = (percent / 100.0)
        if p == 0: p = pseudo
        # 计算 Log-Odds: log2(p / background)
        score = math.log(p / bg, 2)
        pwm_row[base] = score
    pwm.append(pwm_row)

motif_len = len(pwm)

# 辅助函数：读取 FASTA
def read_fasta(file_path):
    genome = {}
    name = ""
    seq = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name: genome[name] = "".join(seq).upper()
                name = line[1:].split()[0] # 获取序列ID
                seq = []
            else:
                seq.append(line)
        if name: genome[name] = "".join(seq).upper()
    return genome

# 辅助函数：计算序列得分
def score_seq(sequence):
    s = 0
    if len(sequence) != motif_len: return -999
    for i in range(motif_len):
        base = sequence[i]
        if base not in ['A','C','G','T']: return -999 # 忽略 N
        s += pwm[i][base]
    return s

# 辅助函数：获取反向互补链
def rev_comp(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return "".join([complement.get(base, base) for base in seq[::-1]])

# 3. 主程序
genome_file = "../fna.fna"
output_file = "TATA_box_prediction.gff3"

print(f"Reading genome from {genome_file}...")
genome = read_fasta(genome_file)

# 设定一个初始硬阈值以减少计算量 (经验值，例如最大得分的 80%)
# 理论最大得分
max_possible_score = sum([max(row.values()) for row in pwm])
threshold = max_possible_score * 0.75 
print(f"Max score: {max_possible_score:.2f}, Pre-filter threshold: {threshold:.2f}")

results = []

print("Scanning genome (this may take a while)...")
with open(output_file, 'w') as out:
    out.write("##gff-version 3\n")
    
    for chrom, sequence in genome.items():
        seq_len = len(sequence)
        # 正链扫描
        for i in range(seq_len - motif_len + 1):
            segment = sequence[i : i+motif_len]
            if 'N' in segment: continue
            
            raw_score = score_seq(segment)
            
            # 只有得分够高才进行耗时的 P-value 计算 (Shuffle 测试) 
            if raw_score > threshold:
                # 简化的 P-value 计算：随机打乱当前片段 100 次
                # 如果打乱后的得分比原始得分高的次数越多，说明原始得分越不显著
                shuffle_count = 100
                better_score_count = 0
                seg_list = list(segment)
                for _ in range(shuffle_count):
                    random.shuffle(seg_list)
                    shuffled_seq = "".join(seg_list)
                    if score_seq(shuffled_seq) >= raw_score:
                        better_score_count += 1
                
                p_value = (better_score_count + 1) / (shuffle_count + 1)
                
                # 筛选条件：P-value < 0.05
                if p_value < 0.05:
                     # GFF3 格式: seqid source type start end score strand phase attributes
                    attr = f"pvalue={p_value:.4f};raw_score={raw_score:.2f}"
                    out.write(f"{chrom}\tPWM_Scan\tTATA_box\t{i+1}\t{i+motif_len}\t{raw_score:.2f}\t+\t.\t{attr}\n")

        # 负链扫描 (在反向互补链上找)
        # 策略：取反向互补链，位置对应关系要注意
        rc_sequence = rev_comp(sequence)
        for i in range(seq_len - motif_len + 1):
            segment = rc_sequence[i : i+motif_len]
            if 'N' in segment: continue
            raw_score = score_seq(segment)
            
            if raw_score > threshold:
                shuffle_count = 100
                better_score_count = 0
                seg_list = list(segment)
                for _ in range(shuffle_count):
                    random.shuffle(seg_list)
                    shuffled_seq = "".join(seg_list)
                    if score_seq(shuffled_seq) >= raw_score:
                        better_score_count += 1
                p_value = (better_score_count + 1) / (shuffle_count + 1)
                
                if p_value < 0.05:
                    # 负链坐标转换：
                    # 反向链上的索引 i 对应正链的：seq_len - (i + motif_len) + 1 到 seq_len - i
                    start = seq_len - (i + motif_len) + 1
                    end = seq_len - i
                    attr = f"pvalue={p_value:.4f};raw_score={raw_score:.2f}"
                    out.write(f"{chrom}\tPWM_Scan\tTATA_box\t{start}\t{end}\t{raw_score:.2f}\t-\t.\t{attr}\n")

print(f"Done! Results saved to {output_file}")
