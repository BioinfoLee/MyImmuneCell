### author by yangshichen
### 注意：脚本仅供参考，使用前请仔细阅读

import argparse
import os 
import pandas as pd

smr = '/media/scPBMC1_AnalysisDisk1/huangzhuoli/Script_HPC/software_gaoyue/SMR/smr-1.3.1-linux-x86_64/smr-1.3.1'

SMR_output_dir = '/media/AnalysisDisk2/Yangshichen/2_My-Onek/WES/Germline/Result/pQTL_GWAS_SMR/'

GWAS_dir = '/media/AnalysisDisk2/Yangshichen/2_My-Onek/WES/Germline/Data/BBJ_GWAS_for_SMR_OPERA/ma_file/'
pQTL_output_dir = '/media/AnalysisDisk2/Yangshichen/2_My-Onek/WES/Germline/Data/SMR_pQTL_besd/besd/'
bim = '/media/AnalysisDisk2/Yangshichen/2_My-Onek/WES/Germline/Data/Genetics/10.maf01'

trait_all = pd.read_csv('/media/AnalysisDisk2/Yangshichen/2_My-Onek/WES/Germline/Data/BBJ_GWAS_for_SMR_OPERA/BBJ_disease.txt',sep='\t',header=None)
trait_all = trait_all[0].to_list()
trait_all = trait_all[11:]

for trait in trait_all:
    print(trait)
    
    # GWAS 文件只处理一次，输出到临时文件
    gwas_file = f"{GWAS_dir}hum0197.v3.BBJ.{trait}.v1.ma"
    gwas_tmp = f"{SMR_output_dir}/{trait}_processed.ma"
    os.system(f"awk 'NR==1{{print $0}} NR>1{{split($1,a,\"_\"); chrom=a[1]; sub(/^chr/,\"\",chrom); print chrom\":\"a[2]\":\"$2\":\"$3, $2, $3, $4, $5, $6, $7, $8}}' {gwas_file} > {gwas_tmp}")
    
    pQTL = f"{pQTL_output_dir}bulk"
    gi = f'"{smr}" --bfile "{bim}" --gwas-summary {gwas_tmp} --beqtl-summary "{pQTL}" --thread-num 1 --out "{SMR_output_dir}/pQTLsmr_{trait}"'
    os.system(gi)