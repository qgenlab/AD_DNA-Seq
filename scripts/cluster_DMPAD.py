import os,sys
import gzip
from scipy.stats import fisher_exact

import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import polars as pl

import time
import copy

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
import ast  # For safer evaluation

from com_fun import *

def plot_Manhattan(df, mfile, orderedchrnum, bidirection=True, fdrname='adj.P.Val'):
    plt.figure(); # 'Chrm', 'Pos', 'MethyPerc', 'Index', 'P-Value', 'MethDeviation'
    # -log_10(pvalue)
    df['minuslog10pvalue'] = -np.log10(df[fdrname])
    df.loc[df['minuslog10pvalue'] > 7, 'minuslog10pvalue'] = 7
    if bidirection:
       df.loc[ df['Meth_dif']<0 , 'minuslog10pvalue'] = -df['minuslog10pvalue'][ df['Meth_dif']<0 ]
    # How to plot gene vs. -log10(pvalue) and colour it by chromosome?
    #
    df = copy.deepcopy(df[~(df['Chr'].isin(['X','Y','M','chrX','chrY','chrM']))])
    df['Chr_order'] = df['Chr'].str.replace('chr', '').astype(int); #pd.Categorical(df['Chr'], categories=orderedchrnum, ordered=True)
    df = df.sort_values(['Chr_order', 'Pos'])
    df['clors'] = 'steelblue'
    #for num, (name, _) in enumerate(df_grouped):
    #for num in range(len(orderedchrnum)):
    #   if num%2==1:
    #      name=orderedchrnum[ num ]
    #      df.loc[df['Chr']==name, 'clors' ] = 'black'
    for i, chr in enumerate(orderedchrnum):
        df.loc[df['Chr'] == chr, 'clors'] = 'steelblue' if i % 2 == 0 else 'seagreen' 
    #df.loc[df[fdrname]<0.05, 'clors'] = 'red'
    #
    # manhattan plot
    df['ind'] = range(len(df))
    df_grouped = df.groupby(('Chr_order'))
    #
    fig = plt.figure(figsize=(16, 7)) # Set the figure size
    ax = fig.add_subplot(111)
    #colors = ['darkred','darkgreen','darkblue', 'gold']
    x_labels = []
    x_labels_pos = []
    for num, (name, group) in enumerate(df_grouped):
        group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=group['clors'], ax=ax, s=0.4)
        x_labels.append(str(name))
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels, rotation=90)
    #
    # set axis limits
    ax.set_xlim([0, len(df)])
    minlogp = df['minuslog10pvalue'].min();
    maxlogp = df['minuslog10pvalue'].max();
    if minlogp<-10: minlogp = -10;
    if maxlogp>10: maxlogp = 10;
    plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')
    plt.axhline(y=np.log10(0.05), color='r', linestyle='--')
    #ax.set_xlim([minlogp, maxlogp])
    #
    # x axis label
    ax.set_xlabel('Chromosome', fontsize=20)
    ax.set_ylabel('Log10(q-value)             -Log10(q-value)', fontsize=18)
    plt.tick_params(axis='both', labelsize=18)
    plt.tight_layout()
    #plt.show()
    plt.savefig(mfile, dpi=600)
    plt.close('all')

def cluster_differential_methylation(meth_df, dms_df, min_pos=3, neural_change_limit=7.5, neurl_perc=30, opposite_perc=10):
    meth_df = meth_df.sort_index()
    dms_df = dms_df.sort_index()
    clustered_regions = []
    unclustered_dms = []
    clustered_dms = []
    for chrom, dms_chr in dms_df.groupby(level=0):  
      positions = dms_chr.index.get_level_values(1).to_numpy()
      meth_diffs = dms_chr["Meth_dif"].to_numpy()
      cluster_start = None
      cluster_end = None
      cluster_sites = []
      cluster_diffs = []
      sign_2 = [0, 0]
      i = -1;
      last_add = i;
      cluster_add_sites = []
      while i<len(positions)-1:
          i = i + 1;
          if cluster_start is None:
              cluster_start = positions[i]
              cluster_end = positions[i]
              cluster_sites.append(positions[i])
              if last_add<=i: cluster_add_sites.append(positions[i])
              cluster_diffs.append(meth_diffs[i])
              if meth_diffs[i]>0: sign_2[1] = sign_2[1] + 1;
              else: sign_2[0] = sign_2[0] + 1;
              continue
          # Check if the site is within 1kb of the last DMS in cluster
          if (positions[i] - cluster_end <= 1000) and ((sign_2[1]>0 and sign_2[0]==0 and meth_diffs[i]>0) or (sign_2[0]>0 and sign_2[1]==0 and meth_diffs[i]<0) ):
              cluster_sites.append(positions[i])
              if last_add<=i: cluster_add_sites.append(positions[i])
              cluster_diffs.append(meth_diffs[i])
              if meth_diffs[i]>0: sign_2[1] = sign_2[1] + 1;
              else: sign_2[0] = sign_2[0] + 1;
              cluster_end = cluster_sites[0] if len(cluster_sites)<min_pos else cluster_sites[-min_pos]
          else:
              # Step 2: Process and validate cluster
              if len(cluster_sites) >= min_pos:
                  trend = np.sign(np.mean(cluster_diffs))
                  if not (all(np.sign(cluster_diffs) == trend)):  # Ensure same trend
                     print('Wrong same sign', chrom, cluster_sites, trend, np.sign(cluster_diffs), cluster_diffs)
                  else:
                      cpg_diff = meth_df.loc[(chrom, slice(cluster_start, cluster_sites[-1])), "Meth_dif"]
                      if cpg_diff.empty:
                          print('Wrong: empty region', chrom, cluster_start, cluster_sites[-1])
                      else:
                         pos_trend_pct = (cpg_diff>neural_change_limit).sum()*100/cpg_diff.shape[0];
                         neg_trend_pct = (cpg_diff<-neural_change_limit).sum()*100/cpg_diff.shape[0];
                         neutral_pct = 100 - pos_trend_pct - neg_trend_pct
                         if neutral_pct < neurl_perc and ( (sign_2[1]>0 and neg_trend_pct<opposite_perc) or (sign_2[0]>0 and pos_trend_pct<opposite_perc)):
                              clustered_regions.append({
                                  "chromosome": chrom,
                                  "start": cluster_start,
                                  "end": cluster_sites[-1],
                                  "num_dms": len(cluster_sites),
                                  "avg_sign_meth": np.mean(cluster_diffs),
                                  "avg_meth": np.mean(cpg_diff),
                                  "std_meth": np.std(cpg_diff),
                                  "pos_trend_pct": pos_trend_pct,
                                  "neg_trend_pct": neg_trend_pct,
                                  "neutral_pct": neutral_pct
                              })
                              clustered_dms.append(dms_chr.loc[(chrom, cluster_sites), :])
                              # 0, 1, 2, 3, 4; i=5;
                              # 
                              if ((sign_2[1]>0 and sign_2[0]==0 and meth_diffs[i]>0) or (sign_2[0]>0 and sign_2[1]==0 and meth_diffs[i]<0) ):
                                 last_add = i
                                 i = i - (min_pos-1)
                              #if chrom=='chr11':
                              #   print('ADD', chrom, cluster_sites, cluster_start, cluster_end, len(cluster_sites), neutral_pct, sign_2, neg_trend_pct, pos_trend_pct)
                         else:
                              print('cannot addd', neutral_pct, sign_2, neg_trend_pct, pos_trend_pct, chrom,  cluster_sites, cluster_start, cluster_end, len(cluster_sites))
                              if len(cluster_add_sites)>0:
                                 unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                              last_add = i;
              else:
                   if len(cluster_add_sites)>0:
                      unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                   last_add = i;
                   #if chrom=='chr11':
                   #  print(chrom, cluster_sites, cluster_start, cluster_end, len(cluster_sites), sign_2 )
              # Reset for new cluster
              cluster_start = positions[i]
              cluster_end = positions[i]
              cluster_sites = [positions[i]]
              if last_add<=i: cluster_add_sites.append(positions[i])
              cluster_diffs = [meth_diffs[i]]
              sign_2 = [0, 0]
              if meth_diffs[i]>0: sign_2[1] = sign_2[1] + 1;
              else: sign_2[0] = sign_2[0] + 1;
      #  
      if len(cluster_sites) >= min_pos:
          trend = np.sign(np.mean(cluster_diffs))
          if not (all(np.sign(cluster_diffs) == trend)):  # Ensure same trend
             print('Wrong same sign',chrom, cluster_sites, trend, np.sign(cluster_diffs), cluster_diffs)
          else:
              cpg_diff = meth_df.loc[(chrom, slice(cluster_start, cluster_sites[-1])), "Meth_dif"]
              if cpg_diff.empty:
                  print('Wrong: empty region', chrom, cluster_start, cluster_sites[-1])
              else:
                 pos_trend_pct = (cpg_diff>neural_change_limit).sum()*100/cpg_diff.shape[0];
                 neg_trend_pct = (cpg_diff<-neural_change_limit).sum()*100/cpg_diff.shape[0];
                 neutral_pct = 100 - pos_trend_pct - neg_trend_pct
                 if neutral_pct < neurl_perc and ( (sign_2[1]>0 and neg_trend_pct<opposite_perc) or (sign_2[0]>0 and pos_trend_pct<opposite_perc)):
                      clustered_regions.append({
                          "chromosome": chrom,
                          "start": cluster_start,
                          "end": cluster_sites[-1],
                          "num_dms": len(cluster_sites),
                          "avg_sign_meth": np.mean(cluster_diffs),
                          "avg_meth": np.mean(cpg_diff),
                          "std_meth": np.std(cpg_diff),
                          "pos_trend_pct": pos_trend_pct,
                          "neg_trend_pct": neg_trend_pct,
                          "neutral_pct": neutral_pct
                      })
                      clustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                      # 0, 1, 2, 3, 4; i=5;
                      if ((sign_2[1]>0 and sign_2[0]==0 and meth_diffs[i]>0) or (sign_2[0]>0 and sign_2[1]==0 and meth_diffs[i]<0) ):
                         last_add = i
                         i = i - (min_pos-1)
                 else:
                    print('cannot addd_', neutral_pct, sign_2, neg_trend_pct, pos_trend_pct, chrom,  cluster_sites, cluster_start, cluster_end, len(cluster_sites)) 
                    if len(cluster_add_sites)>0:
                         unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
                    last_add = i;
      else:
           if len(cluster_add_sites)>0:
              unclustered_dms.append(dms_chr.loc[(chrom, cluster_add_sites), :])
    #
    #print(clustered_dms)
    cluster_df = pd.DataFrame(clustered_regions)
    unclustered_dms_df = pd.concat(unclustered_dms) if unclustered_dms else pd.DataFrame()
    clustered_dms_df = pd.concat(clustered_dms) if clustered_dms else pd.DataFrame()
    return cluster_df, unclustered_dms_df, clustered_dms_df


if __name__=='__main__':
   input_prefx = sys.argv[1]
   orderedchr = []
   orderedchrnum = []
   for ic in range(1, 23):
      orderedchr.append('chr'+str(ic));
      orderedchrnum.append(str(ic));

   orderedchr.append('chrX'); orderedchrnum.append('X')
   orderedchr.append('chrY'); orderedchrnum.append('Y')
   orderedchr.append('chrM'); orderedchrnum.append('M')

   min_value = 10**(-3)
   filter_samples_ratio=0.6
   readThresh=sys.argv[2];

   filter2_meth_df, filter2_covg_df, df_meth_dmp, ad_df_list, ct_df_list, others = filter_input(input_prefx, readThresh, min_value, filter_samples_ratio)

   min_pos_input = int(sys.argv[3]);
   if True: # not os.path.isfile(input_prefx+'_Manhattan_'+str(min_pos_input)+'.png'):
      df_meth_dmp_noindex = df_meth_dmp.reset_index()
      plot_Manhattan(df_meth_dmp_noindex, input_prefx+'_Manhattan_'+str(min_pos_input)+'.png', orderedchr, bidirection=True, fdrname='adj.P.Val')
   os.exit()

   region_df, unclustered_dms_df, cluster_dms_df = cluster_differential_methylation(filter2_meth_df, df_meth_dmp[df_meth_dmp['adj.P.Val']<0.05], min_pos=min_pos_input)

   pl.from_pandas(region_df).write_csv( input_prefx+'_'+str(min_pos_input)+'_region.csv' );
   pl.from_pandas(unclustered_dms_df.reset_index()).write_csv( input_prefx+'_'+str(min_pos_input)+'_unclustered_dms.csv' );
   pl.from_pandas(cluster_dms_df.reset_index()).write_csv( input_prefx+'_'+str(min_pos_input)+'_cluster_dms.csv' );


