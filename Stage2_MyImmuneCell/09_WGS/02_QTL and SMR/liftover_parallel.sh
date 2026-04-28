#!/bin/bash

set -euo pipefail

input_dir="/media/AnalysisDisk2/Yangshichen/2_My-Onek/WGS/Data/GWAS_ma/bbj-2/"
rename="/media/AnalysisDisk2/Yangshichen/2_My-Onek/WGS/Data/chr_rename.txt"
chain="/media/AnalysisDisk2/Yangshichen/2_My-Onek/WGS/Data/hg19ToHg38.over.chain.gz"
ref="/media/AnalysisDisk2/Yangshichen/2_My-Onek/WGS/Data/Homo_sapiens_assembly38.fasta.gz"

threads=8

for vcf in "$input_dir"/*.vcf.gz; do
    [[ "$vcf" == *.hg38.vcf.gz ]] && continue

    (
        base=$(basename "$vcf" .vcf.gz)

        chr_vcf="${input_dir}/${base}.chr.vcf.gz"
        hg38_vcf="${input_dir}/${base}.hg38.vcf.gz"
        reject="${input_dir}/${base}.rejected.vcf"

        echo "Processing $base"

        # rename chr
        bcftools annotate \
            --rename-chrs "$rename" \
            "$vcf" \
            -Oz -o "$chr_vcf"

        bcftools index "$chr_vcf"

        # liftover
        gatk LiftoverVcf \
            -I "$chr_vcf" \
            -O "$hg38_vcf" \
            -CHAIN "$chain" \
            -R "$ref" \
            --REJECT "$reject" \
            --RECOVER_SWAPPED_REF_ALT true

        # clean
        rm -f "$chr_vcf" "${chr_vcf}.csi" "$reject"

        echo "$base done"
    ) &

    # 控制并发数
    while (( $(jobs -r | wc -l) >= threads )); do
        sleep 1
    done

done

wait

echo "All finished."