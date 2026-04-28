### author by yangshichen
### 注意：脚本仅供参考，使用前请仔细阅读

import argparse
import os 
import pandas as pd

CT_list = ['Adaptive NK cells','ALPL- MARCKS- NDNs','ASDC','Atypical naïve B cells','Basophils','CCR4- CD8+ Tcm','CCR4+ CD8+ Tcm','CD14+ cDC2','CD177int iLDNs','CD1C+ cDC2','CD27- IgD- atypical memory B cells','CD27- MAIT','CD27- Th1','CD27- Th17','CD27+ IgD- atypical memory B cells','CD27+ IgD+ atypical memory B cells','CD27+ MAIT','CD27+ Th1','CD27+ Th17','CD279+ SOX4+ Vδ1+ T cells','CD4+ Temra','CD56+ MAIT','CD56bright NK cells','CD56dim NK cells','CD62Lhi GZMK+ Vδ2+ T cells','CD8+ Temra','CD8+ Treg','CD95 memory B cells','cDC1','CILCP','CLP','Core classical monocytes','Core NDNs','CXCL8- PTGS2+ NDNs','DN T cells','Early memory B cells','FOS- NDNs','GBP1+ classical monocytes','GZMB+ CD4+ terminal effector T cells','GZMB+ CD8+ Tem','GZMB+ Vδ2+ T cells','GZMK+ CD8+ Tem','GZMK+ effector Vδ1+ T cells','GZMK+ Vδ2+ T cells','HLA-DRhi CD4+ terminal effector T cells','HLA-DRhi CD4+ Treg','HLA-DRhi CD8+ Tem','HSC_MPP','IFIT2- RNF213- NDNs','IGKChi Plasma cells','IGLL5hi Plasma cells','ILC2','ILCP','iNKT','Intermediate monocytes','IRF1- GBP2- NDNs','ISG+ classical monocytes','KLRB1+ CD4+ Treg','KLRC2+ effector Vδ1+ T cells','LAMP3+ DC','MC_MCP','Memory CD4+ Treg','MEP','MkP','MMP8+ CD177+ iLDNs','MMP8+ CD177+ mLDNs','MMP9+ CD177+ iLDNs','MPO- CD177- iLDNs','MPO+ CD177- iLDNs','MPO+ mLDNs','MT-ATP6- MT-CO2- NDNs','Naïve B cells','Naïve CD4+ T cells','Naïve CD4+ Treg','Naïve CD8+ T cells','Naïve Vδ1+ T cells','Non-classical monocytes','Non-switched memory B cells','pDCs','Plamsablasts','Platelets','Proliferative CD4+ Treg','Proliferative CD8+ memory T cells','Proliferative cytotoxic CD4+ T cells','Proliferative DN T cells','Proliferative help memory T cells','Proliferative iLDNs','Proliferative mLDNs','Proliferative NK cells','Proliferative Vδ2+ T cells','RGS2- NDNs','SOX4+ Vδ1+ T cells','Switched Memory B cells','Tfh','Th1_Th17','Th2','Th22','VIM- FLNA- NDNs','vNKT']

smr = '/media/scPBMC1_AnalysisDisk1/huangzhuoli/Script_HPC/software_gaoyue/SMR/smr-1.3.1-linux-x86_64/smr-1.3.1'
SMR_output_dir = '/media/AnalysisDisk2/Yangshichen/2_My-Onek/WGS/Result/eQTL_pQTL_SMR/'
GWAS_dir = '/media/AnalysisDisk2/Yangshichen/2_My-Onek/WES/Germline/Data/pQTL_ma/'
eQTL_output_dir = '/media/AnalysisDisk2/Yangshichen/2_My-Onek/WGS/Data/SMR_eQTL_besd/besd/'
bim = '/media/AnalysisDisk2/Yangshichen/2_My-Onek/WGS/Data/Genetics/10.maf01'
trait_all = pd.read_csv('/media/scPBMC1_AnalysisDisk1/huangzhuoli/CIMA_r1/Data/pQTL_ma/pQTL_list.txt',sep='\t',header=None)
trait_all = trait_all[0].to_list()

for trait in trait_all:
    print(trait)
    # 创建 trait 输出目录
    trait_dir = f'{SMR_output_dir}{trait}'
    if not os.path.exists(trait_dir):
        os.makedirs(trait_dir)

    # 只处理一次 GWAS 文件，生成临时文件
    gwas_file = f"{GWAS_dir}{trait}.ma"
    processed_gwas_file = f"{trait_dir}/{trait}_processed.ma"
    awk_command = (
        f"awk 'NR==1{{print $0}} NR>1{{split($1,a,\"_\"); "
        f"chrom=a[1]; sub(/^chr/,\"\",chrom); "
        f"print chrom\":\"a[2]\":\"$3\":\"$2, $2, $3, $4, $5, $6, $7, $8}}' "
        f"{gwas_file} > {processed_gwas_file}"
    )
    os.system(awk_command)

     # 定义函数：执行单个 celltype 的 SMR 命令
    from concurrent.futures import ThreadPoolExecutor, as_completed
    def run_smr(celltype):
        eQTL = f"{eQTL_output_dir}{celltype}"
        out_prefix = f"{trait_dir}/pQTLsmr_{celltype}"
        smr_command = (
            f'"{smr}" --bfile "{bim}" '
            f'--gwas-summary "{processed_gwas_file}" '
            f'--beqtl-summary "{eQTL}" '
            f'--thread-num 20 --diff-freq-prop 0.5 '
            f'--out "{out_prefix}"'
        )
        return os.system(smr_command)

    # 并行执行多个 celltype
    max_threads = 25
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        futures = {executor.submit(run_smr, celltype): celltype for celltype in CT_list}
        for fut in as_completed(futures):
            cell = futures[fut]
            code = fut.result()
            print(f"{trait} - {cell} finished with code {code}")