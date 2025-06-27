
import os,sys
import glob
import matplotlib.pyplot as plt
import polars as pl
import pandas as pd

consider_noncoding = True;

base_folder = './'

sv_folder = 'sv/'
#sv_folder = 'sv_2025Apr/'



with open(base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_coding.txt', 'r') as mr:
    genes5_list =[l.strip() for l in mr.readlines()];
    print("hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_coding:", len(genes5_list) )

if consider_noncoding:
   with open(base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_noncoding.txt', 'r') as mr:
       ncgenes5_list =[l.strip() for l in mr.readlines()];
       print("hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_noncoding: ", len(ncgenes5_list) )

with open(base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_coding.txt', 'r') as mr:
    genes3_list =[l.strip() for l in mr.readlines()];
    print("hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_coding: ", len(genes3_list) )

if consider_noncoding:
   with open(base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_noncoding.txt', 'r') as mr:
      ncgenes3_list =[l.strip() for l in mr.readlines()];
      print("hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_noncoding: ", len(ncgenes3_list))

with open(base_folder+sv_folder+'/merge_sv_region_coding.txt', 'r') as mr:
    genesSV_list = [l.strip() for l in mr.readlines()];
    print(sv_folder+"/merge_sv_region_coding: ", len(genesSV_list) )

if consider_noncoding:
    with open(base_folder+sv_folder+'/merge_sv_region_noncoding.txt', 'r') as mr:
       ncgenesSV_list = [l.strip() for l in mr.readlines()];
       print(sv_folder+"/merge_sv_region_noncoding: ", len(ncgenesSV_list) )


gwas_snp_df = pl.read_csv(base_folder+'snp/gwas_Bellenguez.sign.csv', separator='\t').to_pandas();


#################

if consider_noncoding:
   g_list = [genes5_list, ncgenes5_list, genes3_list, ncgenes3_list, genesSV_list, ncgenesSV_list]
   g_sav_f = [base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_coding_gwas.txt', base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_noncoding_gwas.txt', base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_coding_gwas.txt', base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_noncoding_gwas.txt', base_folder+sv_folder+'/merge_sv_region_coding_gwas.txt', base_folder+sv_folder+'/merge_sv_region_noncoding_gwas.txt']
else:
   g_list = [genes5_list, genes3_list, genesSV_list]
   g_sav_f = [base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_coding_gwas.txt', base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_coding_gwas.txt', base_folder+sv_folder+'/merge_sv_region_coding_gwas.txt']


gwas_snp_df_gene = gwas_snp_df[['REGION', 'CHR_ID', 'CHR_POS','MAPPED_GENE', 'UPSTREAM_GENE_DISTANCE','DOWNSTREAM_GENE_DISTANCE', 'SNPS','INTERGENIC', 'P-VALUE']]

snp_gene_list = gwas_snp_df_gene['MAPPED_GENE']
for genes_list_ind in range(len(g_list)):
   print()
   print(g_sav_f[genes_list_ind])
   genes_list = g_list[ genes_list_ind ]
   for ir in range(gwas_snp_df_gene.shape[0]):
      snp_genes = snp_gene_list.iloc[ir].split(' - ')
      #if gwas_snp_df_gene['CHR_ID'].iloc[ir] in ['12',12]:
      #   print(gwas_snp_df_gene.iloc[ir])
      for _sg in snp_genes:
         if _sg in genes_list: print('\t', gwas_snp_df_gene.iloc[ir], _sg);

   




