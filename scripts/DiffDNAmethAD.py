
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
#import warnings
#warnings.filterwarnings('ignore', category=DeprecationWarning)
pd.set_option('future.no_silent_downcasting', True)


from com_fun import *


#
# /mnt/labshare/share/project_data/AD_data/DNAseq/merg_meth_perc.py

def plot_volcano(df, x_col, y_col, x_label, y_label, fig_name):
   plt.figure(figsize=(8, 6))
   sns.scatterplot(data=df, x=x_col, y=y_col, s=14, alpha=0.7)

   pvalue_cutoff = 0.05
   fc_cutoff = 0.25
   # Add lines for significance thresholds
   plt.axhline(y=-np.log10(pvalue_cutoff), color='red', linestyle='--')
   #plt.axvline(x=fc_cutoff, color='red', linestyle='--')
   #plt.axvline(x=-fc_cutoff, color='red', linestyle='--')

   # Label points that meet significance thresholds
   #significant_genes = df[(df['log2FoldChange'] > fc_cutoff) & (df['p_value'] < pvalue_cutoff)] 
   #for i, gene in significant_genes.iterrows():
   #    plt.text(gene['log2FoldChange'], gene['-log10_pvalue'], gene['Gene'], 
   #             fontsize=8, ha='left')

   # Customize plot
   plt.xlabel(x_label, fontsize=16)
   plt.ylabel(y_label, fontsize=16)
   plt.tick_params(axis='both', labelsize=16)

   plt.savefig( fig_name, dpi=600, bbox_inches='tight')
   plt.close('all')



def meth_differential(df, sample_category, row_ad_mean, row_ctrl_mean, csvfile):
   pandas2ri.activate()
   limma = importr('limma')
   df_r = pandas2ri.py2rpy( df )

   if 'Age' in sample_category[0].columns:   
      design_matrix = pd.DataFrame({
             'Group': [sample_category[0].loc[_cn[:-len('_Perc')], 'AD'] for _cn in df.columns ] ,
             'Age': [sample_category[0].loc[_cn[:-len('_Perc')], 'Age'] for _cn in df.columns ] ,
             'Intercept': [1] * (len(df.columns)),
          })
   else:
      design_matrix = pd.DataFrame({
             'Group': [sample_category[0].loc[_cn[:-len('_Perc')], 'AD'] for _cn in df.columns ] ,
             'Intercept': [1] * (len(df.columns)),
          })

   r_design = pandas2ri.py2rpy( design_matrix )

   fit = limma.lmFit(df_r, r_design) 
   r.assign("fit", fit)
   r_summary = """
summary(apply(fit$coefficients, 1, sd))
"""
   summary = r(r_summary)
   print(summary)

   contrast_matrix = limma.makeContrasts(Group='Group', levels=r_design)
   fit = limma.contrasts_fit(fit, contrast_matrix)
   r.assign("fit", fit)
   summary = r(r_summary)
   print(summary)

   #print(fit)
   #print(list(py_fit.names))  # Get the names of elements in the ListVector
   #for AD 
   py_fit = pandas2ri.rpy2py(fit);
   stdev_unscaled = np.array(py_fit.rx2("stdev.unscaled"))
   stdev_unscaled[stdev_unscaled < 0.1] = 0.1
   py_fit.rx2["stdev.unscaled"] = stdev_unscaled
   fit = pandas2ri.py2rpy( py_fit )
   #print(py_fit)

   ### for AD
   fit2 = limma.eBayes(fit, trend=False, robust=False, proportion=1,  stdev_coef_lim= r['c'](0.1, np.inf) )
   #
   #fit2 = limma.treat(fit, lfc=0.2, trend=False, robust=False)
   # For MONO and Prasun's data
   #fit2 = limma.treat(fit, lfc=0.15, trend=False, robust=False)
   # For MONO and Prasun's data
   #fit2 = limma.treat(fit, trend=False, robust=False)
   #fit2 = limma.eBayes(fit, stdev_coef_lim= r['c'](0.1, 2) )
   ##fit2 = limma.eBayes(fit, stdev_coef_lim= r['c'](0.05, 2) )
   #r.assign("fit2", fit2)
   #r_summary = """
   #summary(apply(fit2$coefficients, 1, sd))
   #"""
   #summary = r(r_summary)
   #print(summary)
   #
   results = limma.topTable(fit2, adjust_method='BH', number=df.shape[0], sort_by="p")
   #results = limma.topTable(fit, adjust_method='BH', number=df.shape[0])
   results_df = pandas2ri.rpy2py(results)

   results_df.index = pd.MultiIndex.from_tuples(results_df.index.to_series().apply(ast.literal_eval), names=['Chr', 'Pos', 'Strand']);

   df['Meth_dif'] = row_ad_mean - row_ctrl_mean

   df = pd.merge( df, results_df, left_index=True, right_index=True) ;

   df['-log10_adj.P.Val'] = -np.log10( df['adj.P.Val'] )

   pl.from_pandas(df.reset_index()).write_csv( csvfile );

   return df;

def _cal_mean_(series):
    # Drop NaN values
    series = series.dropna()
    if len( series )==0:
       return min_value
    return series.mean()

def mean_exclude_extremes(series):
    # Drop NaN values
    series = series.dropna()

    # Check if there are at least 3 values to exclude extremes
    if len(series) > 2:
        sorted_series = series.sort_values()
        # Exclude the smallest and largest values
        trimmed_series = sorted_series.iloc[1:-1]
        return trimmed_series.mean()
    else:
        # If not enough values to exclude extremes, return NaN
        return np.nan

# Function to generate random values in [-v, -v/3] or [v/3, v]
def generate_random_values(v):
    if np.random.rand() < 0.5:  # Randomly choose between the two ranges
        return np.random.uniform(-v, -v/3)
    else:
        return np.random.uniform(v/3, v)


if __name__=='__main__':
   np.random.seed(1)

   sample_list = read_sample_list(sys.argv[1])

   is_sep = sys.argv[2];
   file_name = sys.argv[3]
   input_ext = file_name.split('.')
   fig_name = input_ext[0] if len(input_ext)==1 else '.'.join(input_ext[:-1]);

   #sample_name = int(sys.argv[3]);
   #fil_cov = 10;

   min_list_cov = sys.argv[4].split(',')
   min_cov = int( min_list_cov[0] )
   #min_cov = int(sys.argv[4]) 
   if is_sep=='0':
      start_time = time.time();
      df_meth = pl.read_csv( sys.argv[5], separator=',', infer_schema_length=100000).to_pandas();
      df_meth.set_index(['Chr', 'Pos', 'Strand'], inplace=True)
      
      df_covg = pl.read_csv( sys.argv[6], separator=',', infer_schema_length=100000).to_pandas();
      df_covg.set_index(['Chr', 'Pos', 'Strand'], inplace=True)
      print("{} {:.3f} {} {}".format("Read meth_perc & cov", time.time() - start_time, df_meth.shape, df_covg.shape)); sys.stdout.flush();

   else:
   #elif is_sep=='1':
      #min_cov = 7;
      #min_cov = 10
      sample_name = int(sys.argv[5])
      sample_s = 0
      for this_bedf in sys.argv[6:]:
         sample_s = sample_s + 1
         print("Process:", this_bedf); sys.stdout.flush();
         start_time = time.time();
         ad_sample_name = '_'.join(this_bedf.split('/')[-sample_name].split('_')[:2])
         t_df_meth = pl.read_csv(this_bedf, separator='\t', has_header=False, columns=[0,1,5,9,10], infer_schema_length=100000).to_pandas();
         t_df_meth.columns = ['Chr', 'Pos', 'Strand', ad_sample_name+'_Cov', ad_sample_name+'_Perc']
         t_df_meth.set_index(['Chr', 'Pos', 'Strand'], inplace=True)

         #t_df_meth = t_df_meth[ t_df_meth[ad_sample_name+'_Cov' ]>min_cov ] ; # for 7 only
         # t_df_meth = t_df_meth[ (t_df_meth[ad_sample_name+'_Cov' ]>=min_cov) & (t_df_meth[ad_sample_name+'_Perc' ]=='nan') ]
         t_df_meth = t_df_meth[ t_df_meth[ad_sample_name+'_Cov' ]>=min_cov ]

         if sample_s==1:
            df_meth = (t_df_meth[[ad_sample_name+'_Perc']]).copy(deep=True)
            df_covg = (t_df_meth[[ad_sample_name+'_Cov']]).copy(deep=True)
         else:
            df_meth = pd.merge(df_meth, t_df_meth[[ad_sample_name+'_Perc']], left_index=True, right_index=True, suffixes=(None, None), how='outer');
            df_covg = pd.merge(df_covg, t_df_meth[[ad_sample_name+'_Cov']],  left_index=True, right_index=True, suffixes=(None, None), how='outer');
         print("{} {:.3f}".format(ad_sample_name, time.time() - start_time)); sys.stdout.flush();

      pl.from_pandas(df_meth.reset_index()).write_csv( fig_name+'_meth_perc_'+str(len(sys.argv[6:]))+'_cov'+str(min_cov)+'.csv' )
      pl.from_pandas(df_covg.reset_index()).write_csv( fig_name+'_cov_'+str(len(sys.argv[6:]))+'_cov'+str(min_cov)+'.csv' )

   #print(sample_list)
   nc_all = df_meth.columns;
   ad_df_list = [[],[]]
   ct_df_list = [[],[]]
   for c in nc_all:
      c_sp = c.split('_')
      if len(c_sp)>1 and c_sp[-1] in ['Cov', 'Perc']:
         if '_'.join(c_sp[:(2 if len(c_sp)>2 else 1)]) in sample_list[1]:
            ad_df_list[0].append(c)
            ad_df_list[1].append( '_'.join(c_sp[:(2 if len(c_sp)>2 else 1)])+'_Cov' )
         elif '_'.join(c_sp[:(2 if len(c_sp)>2 else 1)]) in sample_list[2]:
            ct_df_list[0].append(c)
            ct_df_list[1].append( '_'.join(c_sp[:(2 if len(c_sp)>2 else 1)])+'_Cov' )
         else: print("Not support", c)
      else:
         pass;
         #ad_df_list[0].append( c) 
         #ct_df_list[0].append( c)
         #ad_df_list[1].append( c)
         #ct_df_list[1].append( c)
   print( len(ad_df_list[0]), ad_df_list )
   print( len(ct_df_list[0]), ct_df_list )

   df_covg.fillna(0)
   mprange_count = (df_meth > 1.5).sum().sum()
   if mprange_count > (df_meth.shape[0] * df_meth.shape[1])*0.01:
      df_meth = df_meth/100.0

   row_std = df_meth.std(axis=1, skipna=True)
   min_std = 0.075;
   min_std = 0.1
   #min_std = 0.05

   readThresh = min_cov
   if len(min_list_cov)>1:
      readThresh = int(min_list_cov[1] )
   filter_samples_ratio = 0.6;
   minSamp_num = [ len(ad_df_list[0])*filter_samples_ratio, len(ct_df_list[0])*filter_samples_ratio ]
   for _i in range(len(minSamp_num)):
      if minSamp_num[_i]<2: 
         minSamp_num[_i]=2;
  
   start_time = time.time(); 
   #filtered_meth_df = copy.deepcopy(df_meth[( ( (df_covg>=readThresh).sum(axis=1)>=(minSamp_num[0]+minSamp_num[1]) ) & (row_std>min_std ) ) ])
   #filtered_covg_df = copy.deepcopy(df_covg[( ( (df_covg>=readThresh).sum(axis=1)>=(minSamp_num[0]+minSamp_num[1]) ) & (row_std>min_std ) ) ])

   # Filter: total samples with enough coverage;
   filtered_meth_df = copy.deepcopy(df_meth[( ( (df_covg>=readThresh).sum(axis=1)>=(minSamp_num[0]+minSamp_num[1]) ) ) ] )
   filtered_covg_df = copy.deepcopy(df_covg[( ( (df_covg>=readThresh).sum(axis=1)>=(minSamp_num[0]+minSamp_num[1]) ) ) ] )

   print('Orignal:', df_meth.shape, df_covg.shape, 'filtered', filtered_meth_df.shape, filtered_covg_df.shape, minSamp_num, len(ad_df_list[0]), len(ct_df_list[0]))
   print("{} {:.3f}".format("Loading and filter1", time.time() - start_time)); sys.stdout.flush();

   # before 2025 Jan 24;
   #ctrl_cov_cond = filtered_covg_df.loc[:, ct_df_list[1]]>=readThresh
   #case_cov_cond = filtered_covg_df.loc[:, ad_df_list[1]]>=readThresh
   ##group_filer = ( (ctrl_cov_cond.sum(axis=1)>=minSamp_num[1]) | (case_cov_cond.sum(axis=1)>=minSamp_num[0]) )
   ## Filter: samples in each group with enough coverage;
   #group_filer = ( (ctrl_cov_cond.sum(axis=1)>=minSamp_num[1]) & (case_cov_cond.sum(axis=1)>=minSamp_num[0]) )


   ctrl_cov_cond = (filtered_covg_df.loc[:, ct_df_list[1]]>=readThresh).sum(axis=1)
   case_cov_cond = (filtered_covg_df.loc[:, ad_df_list[1]]>=readThresh).sum(axis=1)
   group_filer = ( ( ctrl_cov_cond>=minSamp_num[1]) & ( case_cov_cond>=minSamp_num[0]) )

   # use apply is too slow;
   #row_ctrl_mean = filtered_meth_df.loc[:, ct_df_list[0]].apply(_cal_mean_, axis=1); #
   #row_ad_mean = filtered_meth_df.loc[:, ad_df_list[0]].apply(_cal_mean_, axis=1); #
   row_ctrl_mean = filtered_meth_df.loc[:, ct_df_list[0]].mean(axis=1, skipna=True).fillna(min_value)
   row_ad_mean = filtered_meth_df.loc[:, ad_df_list[0]].mean(axis=1, skipna=True).fillna(min_value)

   # Filter                                  ctrl has more na, ctrl mean less than 0.2; or ad have more na with ad mean less than 0.2
   small_mean = 0.2

   if minSamp_num[0]+minSamp_num[1] < len(ad_df_list[0]) + len(ct_df_list[0]): 
      group_filer = group_filer | ( ( ( case_cov_cond>=minSamp_num[0] ) & ( ctrl_cov_cond<minSamp_num[1]) & ( row_ctrl_mean<small_mean ) ) |  ( ( ctrl_cov_cond>=minSamp_num[1]) & ( case_cov_cond<minSamp_num[0] ) & ( row_ad_mean<small_mean ) ) )   

   #ad_na_sum = filtered_meth_df.loc[:, ad_df_list[0]].isna().sum(axis=1)
   #ctrl_na_s = filtered_meth_df.loc[:, ct_df_list[0]].isna().sum(axis=1)
   #group_filer = group_filer | ( ( ( ad_na_sum<len(ad_df_list[0])-minSamp_num[0] ) & ( ctrl_na_s>=minSamp_num[1] ) & ( row_ctrl_mean<small_mean ) ) |  ( ( ctrl_na_s<len(ct_df_list[0])-minSamp_num[1] ) & ( ad_na_sum>=minSamp_num[0] ) & ( row_ad_mean<small_mean ) ) )

   # previous version: before 2025 Jan 24;
   # group_filer = group_filer | ( ( ( filtered_meth_df.loc[:, ad_df_list[0]].isna().sum(axis=1)<len(ad_df_list[0])-minSamp_num[0] ) & ( filtered_meth_df.loc[:, ct_df_list[0]].isna().sum(axis=1)>=minSamp_num[1] ) & ( row_ctrl_mean<small_mean ) ) |  ( ( filtered_meth_df.loc[:, ct_df_list[0]].isna().sum(axis=1)<len(ct_df_list[0])-minSamp_num[1] ) & ( filtered_meth_df.loc[:, ad_df_list[0]].isna().sum(axis=1)>=minSamp_num[0] ) & ( row_ad_mean<small_mean ) ) )
   # 
   #group_filer = group_filer | ( ( ( filtered_meth_df.loc[:, ct_df_list[0]].isna().sum(axis=1)>=len(ct_df_list[0])-minSamp_num[1] ) & ( row_ad_mean-row_ctrl_mean>0.3 ) ) |  ( ( filtered_meth_df.loc[:, ad_df_list[0]].isna().sum(axis=1)>=len(ad_df_list[0])-minSamp_num[0] ) & ( row_ctrl_mean - row_ad_mean>0.3 ) ) )

   filter2_meth_df = copy.deepcopy(filtered_meth_df[group_filer])
   print("{} {:.3f}".format("Filter2", time.time() - start_time)); sys.stdout.flush();

   # fill with row mean;
   #row_mean = filter2_meth_df.apply(_cal_mean_, axis=1);
   #filter2_meth_df.fillna(row_mean);
   
   # fill with group mean;
   # /mnt/analysis/qgenlab/projects/AD_data_analysis/AD_unmapped/LLM_differential.py
   ## This process is slow;
   #row_ctrl_mean = filter2_meth_df.loc[:, ct_df_list[0]].apply(_cal_mean_, axis=1); #
   #row_ad_mean = filter2_meth_df.loc[:, ad_df_list[0]].apply(_cal_mean_, axis=1); #
   row_ctrl_mean = copy.deepcopy(row_ctrl_mean[group_filer])
   row_ad_mean = copy.deepcopy(row_ad_mean[group_filer])

   ## calculate within-group std for randomness;
   add_randomness = False;
   if add_randomness:
      ctrl_row_std = filter2_meth_df[ ct_df_list[0] ].std(axis=1, skipna=True).fillna(0)
      adpt_row_std = filter2_meth_df[ ad_df_list[0] ].std(axis=1, skipna=True).fillna(0)
      #
      ctrl_row_std = min_std - ctrl_row_std
      adpt_row_std = min_std - adpt_row_std
      ctrl_row_std[ctrl_row_std<0] = 0
      adpt_row_std[adpt_row_std<0] = 0
   ###########
   
   #filter2_meth_df.loc[:, ct_df_list[0]].fillna( row_ctrl_mean)
   #filter2_meth_df.loc[:, ad_df_list[0]].fillna( row_ad_mean)
   for  _c_n in filter2_meth_df.columns:
      if _c_n in ct_df_list[0]: filter2_meth_df[ _c_n] = filter2_meth_df[ _c_n].fillna(row_ctrl_mean)
      # wrong; make correction
      #elif _c_n in ad_df_list[0]: filter2_meth_df[ _c_n] = filter2_meth_df[ _c_n].fillna(row_ctrl_mean)
      elif _c_n in ad_df_list[0]: filter2_meth_df[ _c_n] = filter2_meth_df[ _c_n].fillna(row_ad_mean)
      else: print("Warning!! unexpected row", _c_n)

   print("{} {:.3f}".format("Filter2", time.time() - start_time)); sys.stdout.flush();
   row_std = filter2_meth_df.std(axis=1, skipna=True)
   filter3_meth_df = copy.deepcopy( filter2_meth_df[ row_std>=min_std ] )

   print("{} {:.3f} {} {}".format("Filter3", time.time() - start_time, filter2_meth_df.shape, filter3_meth_df.shape)); sys.stdout.flush();

   ##### add randomness;
   if add_randomness:
      ctrl_row_std = ctrl_row_std[ row_std>=min_std ]
      adpt_row_std = adpt_row_std[ row_std>=min_std ]
      random_pdf = pd.DataFrame(index=filter2_meth_df[ row_std>=min_std ].index)
      for  _c_n in filter2_meth_df.columns:
         if _c_n in ct_df_list[0]: 
             #random_pdf[ _c_n ] = np.random.uniform(-ctrl_row_std, ctrl_row_std) 
             random_pdf[ _c_n ] = ctrl_row_std.apply(generate_random_values)
         elif _c_n in ad_df_list[0]: 
             #random_pdf[ _c_n ] = np.random.uniform(-adpt_row_std, adpt_row_std) ;
             random_pdf[ _c_n ] = adpt_row_std.apply(generate_random_values)
         else: print("Warning!! unexpected row", _c_n)
         
      for  _c_n in filter2_meth_df.columns:
         filter3_meth_df[ _c_n] = filter3_meth_df[ _c_n] + random_pdf[ _c_n ];
      print ('randomness', time.time() - start_time, random_pdf)  

   #row_ctrl_mean = filter3_meth_df.loc[:, ct_df_list[0]].apply(_cal_mean_, axis=1); #
   #row_ad_mean = filter3_meth_df.loc[:, ad_df_list[0]].apply(_cal_mean_, axis=1); #

   res_df = meth_differential(filter3_meth_df, sample_list, row_ad_mean, row_ctrl_mean, fig_name+'.csv')
   print (time.time() - start_time,  res_df.columns)
   
   #res_df['-log10_adj.P.Val'] = -np.log10( res_df['adj.P.Val'] )
   plot_volcano(res_df, 'logFC', '-log10_adj.P.Val', 'logFC(AD-Ctrl)', '-log10(padj) for AF', fig_name+'_percentage.png' );
   print (time.time() - start_time)

   #row_ctrl_mean = res_df.loc[:, ct_df_list[0]].apply(_cal_mean_, axis=1); #
   #row_ad_mean = res_df.loc[:, ad_df_list[0]].apply(_cal_mean_, axis=1); #
   #res_df['Meth_dif'] = row_ad_mean - row_ctrl_mean

   plot_volcano(res_df, 'Meth_dif', '-log10_adj.P.Val', 'Mean methylation difference (AD-Ctrl)', '-log10(padj)', fig_name+'_volcano.png' );

   print( res_df[ res_df['adj.P.Val']<=0.05 ].shape)

   #with pd.option_context('display.max_columns', None):
   #   print('AD', ad_df_list)
   #   print(df[df['adj.P.Val']<0.05][ad_df_list])
   #   print('ctrl', ct_df_list)
   #   print(df[df['adj.P.Val']<0.05][ct_df_list])
   

  
