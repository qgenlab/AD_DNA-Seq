
import os,sys
import glob
import matplotlib.pyplot as plt
import polars as pl
import pandas as pd

import copy


plot_v_rna = False
consider_noncoding = False;
#consider_noncoding = True;

base_folder = './'
consider_early = True;
consider_early = False;

# previously, this value is 0.2
logfc_thre = 0.3

sv_folder = 'sv/'
#sv_folder = 'sv_2025Apr/'


def plot_volcono(df, fig_name, fc_c, p_c, padj_t = 0.05, fc_t = 1):
   df['-log10(P-adj)'] = -np.log10(df[p_c])
   fc_threshold = fc_t      # absolute log2FC threshold
   padj_threshold = padj_t # p-adj significance threshold
   plt.figure(figsize=(10, 7))
   #
   # Scatter points: color by significance
   colors = df.apply(lambda row: 
                  'red' if (row[fc_c] > fc_threshold and row['P-adj'] < padj_threshold) else
                  'blue' if (row[fc_c] < -fc_threshold and row['P-adj'] < padj_threshold) else
                  'grey', axis=1)
   #
   plt.scatter(df[fc_c], df['-log10(P-adj)'], c=colors, alpha=0.6, edgecolors='k')
   # Add axis labels and title
   plt.xlabel('log2 Fold Change', fontsize=14)
   plt.ylabel('-log10 Adjusted P-value', fontsize=14)
   plt.title('Volcano Plot', fontsize=16)
   # Add reference lines
   plt.axhline(-np.log10(padj_threshold), color='black', linestyle='--', linewidth=1)
   plt.axvline(fc_threshold, color='black', linestyle='--', linewidth=1)
   plt.axvline(-fc_threshold, color='black', linestyle='--', linewidth=1)
   #
   plt.tight_layout()
   plt.savefig(fig_name+'.png', dpi=600)
   plt.close('all')


if consider_early:
   rna_seq_file = {'EC_all_adj.txt':'Entorhinal Cortex','HP_all_adj.txt':'Hippocampus'}
else:
   rna_seq_file = {'EC_all_adj.txt':'Entorhinal Cortex', 'FC_all_adj.txt':'Frontal Cortex', 'HP_all_adj.txt':'Hippocampus', 'TC_all_adj.txt':'Temporal Cortex'}

ovlp_genes = {}
RNA_seq_dict = {}
All_genes = set()
sig_genes = set()

for _fn in rna_seq_file:
   RNA_seq_dict[rna_seq_file[_fn] ]= pl.read_csv(base_folder+'/brain_rna/'+_fn, separator='\t').to_pandas();
   RNA_seq_dict[rna_seq_file[_fn] ]['P-adj'] = RNA_seq_dict[rna_seq_file[_fn] ]['P-adj'].astype(float)

   All_genes.update( RNA_seq_dict[rna_seq_file[_fn] ]['GENE'].unique() )
   sign_df = RNA_seq_dict[rna_seq_file[_fn] ][((RNA_seq_dict[rna_seq_file[_fn] ]['P-adj']<=0.05 ) & ((RNA_seq_dict[rna_seq_file[_fn] ]['log2FC']<=-logfc_thre) | (RNA_seq_dict[rna_seq_file[_fn] ]['log2FC']>=logfc_thre)))]
   sig_genes.update( sign_df['GENE'].unique() )

   print("RNA-seq data:", rna_seq_file[_fn], (RNA_seq_dict[rna_seq_file[_fn] ]['P-adj']<=0.05).sum(), sign_df.shape)
   for c_g in sign_df['GENE']:
      if c_g not in ovlp_genes: ovlp_genes[ c_g ] = set()
      ovlp_genes[ c_g ].add( rna_seq_file[_fn] );
print("RNA-seq data:", "All-genes", len(All_genes), "P-adj<=0.05 & abs('log2FC')>=logfc_thre", len(sig_genes))

ovlp_n = [0, 0, 0, 0]
for c_g in ovlp_genes:
   r_n = len(ovlp_genes[c_g]);
   ovlp_n[r_n-1] += 1
print("ovlp-sign-genes", ovlp_n, ovlp_n[3]+ovlp_n[1]+ovlp_n[2] )

def filter_for_coding(mlist):
   f_list = []
   for _c_g in mlist:
      if len(_c_g)>3 and _c_g[:3]=='MIR': continue;
      cg_sp = _c_g.split('-')
      if len(cg_sp)>1 and len(cg_sp[1])>=3 and cg_sp[1][:2]=='AS':
         continue;
      f_list.append( _c_g )
   return f_list

with open(base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_coding.txt', 'r') as mr:
    genes5_list = filter_for_coding([l.strip() for l in mr.readlines()]);
    print("hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_coding: #genes, #genes&all-genes, #genes&sign-genes", len(genes5_list), len(set(genes5_list) & All_genes), len(set(genes5_list) & sig_genes) )

if consider_noncoding:
   with open(base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_noncoding.txt', 'r') as mr:
       ncgenes5_list =[l.strip() for l in mr.readlines()];
       print("hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_noncoding: #genes, #genes&all-genes, #genes&sign-genes", len(ncgenes5_list), len(set(ncgenes5_list) & All_genes), len(set(ncgenes5_list) & sig_genes) )

with open(base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_coding.txt', 'r') as mr:
    genes3_list =filter_for_coding([l.strip() for l in mr.readlines()]);
    print("hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_coding: #genes, #genes&all-genes, #genes&sign-genes", len(genes3_list), len(set(genes3_list) & All_genes), len(set(genes3_list) & sig_genes))

if consider_noncoding:
   with open(base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_noncoding.txt', 'r') as mr:
      ncgenes3_list =[l.strip() for l in mr.readlines()];
      print("hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_noncoding: #genes, #genes&all-genes, #genes&sign-genes", len(ncgenes3_list), len(set(ncgenes3_list) & All_genes), len(set(ncgenes3_list) & sig_genes))

with open(base_folder+sv_folder+'/merge_sv_region_coding.txt', 'r') as mr:
    genesSV_list = filter_for_coding([l.strip() for l in mr.readlines()]);
    print(sv_folder+"/merge_sv_region_coding: #genes, #genes&all-genes, #genes&sign-genes", len(genesSV_list), len(set(genesSV_list) & All_genes), len(set(genesSV_list) & sig_genes) )

if consider_noncoding:
    with open(base_folder+sv_folder+'/merge_sv_region_noncoding.txt', 'r') as mr:
       ncgenesSV_list = [l.strip() for l in mr.readlines()];
       print(sv_folder+"/merge_sv_region_noncoding: #genes, #genes&all-genes, #genes&sign-genes", len(ncgenesSV_list), len(set(ncgenesSV_list) & All_genes), len(set(ncgenesSV_list) & sig_genes) )


#################

if consider_noncoding:
   g_list = [genes5_list, ncgenes5_list, genes3_list, ncgenes3_list, genesSV_list, ncgenesSV_list]
   g_sav_f = [base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_coding_inrna.txt', base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_noncoding_inrna.txt', base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_coding_inrna.txt', base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_noncoding_inrna.txt', base_folder+sv_folder+'/merge_sv_region_coding_inrna.txt', base_folder+sv_folder+'/merge_sv_region_noncoding_inrna.txt']
else:
   g_list = [genes5_list, genes3_list, genesSV_list]
   g_sav_f = [base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_5_region_coding_inrna.txt', base_folder+'/meth/hg38_destrand2025Mar_meth_merge_nocancerCov10_3_region_coding_inrna.txt', base_folder+sv_folder+'/merge_sv_region_coding_inrna.txt']

for c_gl_ind in range(len(g_list)): # for a annotation-gene list
   sig_ge_oc = [0, 0, 0, 0]
   c_gl = g_list[c_gl_ind]
   sign_list = [];
   print('For', g_sav_f[c_gl_ind] )
   for _g_c in c_gl: # for each gene
      c_sig = {}
      c_sig["Gene" ] = _g_c
      is_sign = 0;
      for t_k in RNA_seq_dict:
         sig_df = RNA_seq_dict[t_k][RNA_seq_dict[t_k]['GENE']==_g_c]
         if (not (sig_df.empty)):
            if sig_df.shape[0]>1: print("Error: multiplge gene", t_k, _g_c)
            c_sig[t_k+'-padj'] = sig_df['P-adj'].iloc[0]
            c_sig[t_k+'-logfc'] = sig_df['log2FC'].iloc[0]
            if sig_df['P-adj'].iloc[0]<=0.05 and sig_df['log2FC'].abs().iloc[0]>=logfc_thre:
      	      is_sign += 1;
         else: c_sig[t_k] = float('nan') 
      if is_sign>0: 
         sign_list.append( c_sig )
         sig_ge_oc[is_sign-1] += 1
         if is_sign>=2: #(3 if not consider_early else 2): #len(RNA_seq_dict):
            print('\t', is_sign, _g_c);
   if len( sign_list )>0:
      print('\tsign_list', pd.DataFrame(sign_list).shape)
      if not consider_early:
         pl.from_pandas( pd.DataFrame(sign_list)).write_csv( g_sav_f[c_gl_ind], float_precision=7 );
      print('\tcount-ovlp', sig_ge_oc, sig_ge_oc[1]+sig_ge_oc[2]+sig_ge_oc[3])

print()
##################################
if consider_early:
   for c_gl_ind in range(len(g_list)):  # for a annotation-gene list
      c_gl = g_list[c_gl_ind]  
      org_dfs = []; # orginal_dataframe
      sig_dfs = []; # dataframe with significance
      ovl_dfs = []; # dataframe with olvp significance
      for t_k in RNA_seq_dict:
         org_dfs.append( RNA_seq_dict[t_k] );
         sig_df = RNA_seq_dict[t_k][( (RNA_seq_dict[t_k]['P-adj']<=0.05) & ( RNA_seq_dict[t_k]['log2FC'].abs()>=logfc_thre) ) ]
         sig_dfs.append( sig_df );
         ovl_df = sig_df[sig_df['GENE'].isin(c_gl)]
         ovl_dfs.append( ovl_df );
         print( t_k, RNA_seq_dict[t_k].shape[0], sig_df.shape[0], ovl_df.shape[0]);
      org_ovlp = (set(org_dfs[0]['GENE'])).intersection(set(org_dfs[1]['GENE']))
      sig_ovlp = (set(sig_dfs[0]['GENE'])).intersection(set(sig_dfs[1]['GENE']))
      ovl_ovlp = (set(ovl_dfs[0]['GENE'])).intersection(set(ovl_dfs[1]['GENE']))
      print("original", len(org_ovlp), org_dfs[0].shape[0]+org_dfs[1].shape[0]-len(org_ovlp))
      print("signific", len(sig_ovlp), sig_dfs[0].shape[0]+sig_dfs[1].shape[0]-len(sig_ovlp))
      print("ovlerlap", len(ovl_ovlp), ovl_dfs[0].shape[0]+ovl_dfs[1].shape[0]-len(ovl_ovlp), ovl_dfs[0].shape[0]-len(ovl_ovlp), ovl_dfs[1].shape[0]-len(ovl_ovlp))
   
sv_genes = copy.deepcopy(genesSV_list);
if consider_noncoding: sv_genes.extend( ncgenesSV_list )
genes5 = copy.deepcopy(genes5_list)
if consider_noncoding: genes5.extend( ncgenes5_list )
genes3 = copy.deepcopy(genes3_list)
if consider_noncoding: genes3.extend( ncgenes3_list )

set_svg = set(sv_genes)
set_gs5 = set(genes5)
set_gs3 = set(genes3)

print('ovlp 5&3', len(set_gs5 & set_gs3), len(set_gs5), len(set_gs3), len(set(All_genes) & set_gs5 & set_gs3) )
print('ovlp5&sv', len(set_gs5 & set_svg), len(set_gs5), len(set_svg), len(set(All_genes) & set_gs5 & set_svg) )
print('ovlp3&sv', len(set_gs3 & set_svg), len(set_gs3), len(set_svg), len(set(All_genes) & set_svg & set_gs3) )

print( 'ovlp5&sv', ','.join(list((set_gs5 & set_svg))), ','.join(list((set(All_genes) & set_svg & set_gs5))) )

print( 'ovlp3&sv', ','.join(list((set_gs3 & set_svg))), ','.join(list((set(All_genes) & set_svg & set_gs3))) ) 


print( "_5", len(set_gs5), len(set(All_genes) & set_gs5), ','.join(set(All_genes) & set_gs5), '\n', ','.join(set_gs5))
print()
print( "_3", len(set_gs3), len(set(All_genes) & set_gs3), ','.join(set(All_genes) & set_gs3), '\n', ','.join(set_gs3))
print()
print( "sv", len(set_svg), len(set(All_genes) & set_svg), ','.join(set(All_genes) & set_svg), '\n', ','.join(set_svg))
print()


