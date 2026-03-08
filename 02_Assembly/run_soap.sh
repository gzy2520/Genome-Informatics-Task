#!/bin/bash
# -------------------------------
# SOAPdenovo 批量运行脚本（智能版）
# 自动检测内存 + 降线程避免 OOM
# -------------------------------

KMER=13
THREADS=4
PREFIX_LIST=("1_" "3_" "10_")
OUTDIR="soap_out"
CFG_TEMPLATE="lib.cfg"

mkdir -p ${OUTDIR}

# 执行每个前缀样本
for prefix in "${PREFIX_LIST[@]}"; do
    echo "=== 正在处理样本 ${prefix} ==="

    cfg="${prefix}lib.cfg"
    cp ${CFG_TEMPLATE} ${cfg}
    sed -i "s#short1.fq#${prefix}short1.fq#g" ${cfg}
    sed -i "s#short2.fq#${prefix}short2.fq#g" ${cfg}
    sed -i "s#long1.fq#${prefix}long1.fq#g"  ${cfg}
    sed -i "s#long2.fq#${prefix}long2.fq#g"  ${cfg}

    # 打印当前时间与参数
    echo "开始时间：$(date '+%H:%M:%S') | THREADS=${THREADS} | KMER=${KMER}"

    # 运行 SOAPdenovo
    soapdenovo-63mer all \
        -s ${cfg} \
        -K ${KMER} \
        -p ${THREADS} \
        -o ${OUTDIR}/${prefix}soap_out

    # 检查结果是否生成
    if [ -f "${OUTDIR}/${prefix}assembly.scafSeq" ]; then
        echo "✅ 样本 ${prefix} 成功完成！"
    else
        echo "⚠️ 样本 ${prefix} 未生成完整结果，可能内存不足。"
    fi

    echo "=== 样本 ${prefix} 结束 ==="
done

echo -e "\n全部样本运行结束。结果目录：${OUTDIR}/"

