
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
import pysam
import re

import time

import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
import ast  # For safer evaluation
#import warnings
#warnings.filterwarnings('ignore', category=DeprecationWarning)
pd.set_option('future.no_silent_downcasting', True)

import gzip
import copy

from com_fun import *


#                         rep_ext_bp=20, seg_ext_bp
def _segment_genome(bed_file, orderedchr, rep_ext_bp = 20, seg_ext_bp=20, seg_bin=100, unusedchr=None):
   # Load BED file using pandas
   bed_df = pl.read_csv(bed_file, separator='\t', has_header=True, columns=[0,1,2,4], infer_schema_length=10000).to_pandas();
   bed_df = bed_df[bed_df["chr"].isin(orderedchr)]
   #print( bed_df )
   # Adjust repeat regions by subtracting 20bp from start and adding 20bp to end
   bed_df["start"] = bed_df["start"] - rep_ext_bp
   bed_df["end"] = bed_df["end"] + rep_ext_bp

   # Ensure start is not negative
   bed_df["start"] = bed_df["start"].apply(lambda x: max(x, 0))

   # Sort repeat regions by chromosome and start position
   bed_df = bed_df.sort_values(by=["chr", "start"]).reset_index(drop=True)

   # Create list of final genome segments
   genome_segments = []
   # Iterate over repeat regions and create non-repeat segments
   previous_end = 0
   for i, row in bed_df.iterrows():
       if (not unusedchr==None) and row["chr"] in unusedchr: continue;
       if row["chr"] not in orderedchr:
          continue;
       chrom, start, end, repeat_cat = row["chr"], row["start"], row["end"], row["repeat_cat"]
       # Handle non-repeat region
       if start > previous_end:
           non_repeat_start = previous_end
           non_repeat_end = start + rep_ext_bp
           non_repeat_length = non_repeat_end - non_repeat_start
           rougth_seg_n = non_repeat_length/seg_bin
           if rougth_seg_n>0.5: 
              this_bin_size = int(non_repeat_length/int(rougth_seg_n + 0.5) )
           else: this_bin_size = seg_bin
 
           while non_repeat_start + this_bin_size < non_repeat_end:
              genome_segments.append([chrom, non_repeat_start, non_repeat_start + this_bin_size, "Non-repeat"])
              non_repeat_start += (this_bin_size - seg_ext_bp)
           non_repeat_length = non_repeat_end - non_repeat_start
           if non_repeat_length < this_bin_size/2:
               adjust_half = non_repeat_length // 2
               genome_segments[-1][2] = genome_segments[-1][2] + adjust_half
               genome_segments.append([chrom, start - (non_repeat_length-adjust_half), end, repeat_cat])
           else:
               genome_segments.append([chrom, non_repeat_start, non_repeat_end, "Non-repeat"])
               genome_segments.append([chrom, start, end, repeat_cat])
       else:
          # Add repeat region as a segment
          genome_segments.append([chrom, start, end, repeat_cat])
       previous_end = end - rep_ext_bp  # Update previous_end for next iteration

   segments_by_chr = {}
   cur_chr = ''
   for _c_seg in genome_segments:
      if (not _c_seg[0]==cur_chr):
         segments_by_chr[ _c_seg[0] ] = []
         cur_chr =  _c_seg[0]
      segments_by_chr[ _c_seg[0] ].append( copy.deepcopy(_c_seg[1:3] ))
   print("segments_by_chr.keys()" , segments_by_chr.keys() )
   segments_df = pd.DataFrame(genome_segments, columns=["chr", "start", "end", "repeat_cat"])
   segments_df = segments_df.sort_values(by=["chr", "start"]).reset_index(drop=True)
   return (segments_df, genome_segments, segments_by_chr)

def _read_sv_(genome_segments, vcf_files, output_file, sample_ind, orderedchr, unusedchr):
   # Process SVs from VCF files
   sample_sv_counts = {}
   sample_sv_bases = {}
  
   for _ad_sv_type in Def_sv_types:
      sample_sv_counts[_ad_sv_type] = {}
      sample_sv_bases[_ad_sv_type] = {} 
   #segments_df_by_chr = {chr_name: df for chr_name, df in genome_segments[0].groupby("chr")}
   segments_by_chr = genome_segments[2]
   start_time = time.time();

   sv_count_dict = {}
   sv_base_dict = {}

   sv_add_dict = {}

   testi = 0; 
   for vcf_file in vcf_files:
       #print("process: ", vcf_file); sys.stdout.flush()
       sample_name = vcf_file.split("/")[-sample_ind]  # Extract sample name
       for _ad_sv_type in Def_sv_types:
          sample_sv_counts[_ad_sv_type][sample_name] = {}
          sample_sv_bases[_ad_sv_type][sample_name] = {}

       sv_count_dict[sample_name] = {}
       sv_base_dict[sample_name] = {}

       sv_add_dict[sample_name] = {}

       cur_chr = ''
       chr_ind = 0;
       with open(vcf_file, 'r') as vcfread:
          lines = vcfread.readlines()
       for sv in lines:
           sv = sv.strip();
           if len(sv)==0: continue;
           if sv[0]=='#': continue;
           sv = sv.split('\t')

           chrom, pos, m_filter = sv[0], int(sv[1])-1, sv[6]
           if (not unusedchr==None) and chrom in unusedchr: continue;
           if chrom  not in orderedchr: continue;
           if not m_filter=='PASS': continue;          

           info_sp = sv[7].split(';')
           svtype = None;
           svlen= None;
           support = None;
           sv_end = None;
           STDEV_POS = None;
           find_n = 0;
           for _isp in info_sp:
             if find_n>=5: break;
             t_type = _isp.split('=')
             if t_type[0]=='SVTYPE': svtype = t_type[1]; find_n=find_n+1
             elif t_type[0]=='SVLEN': svlen = abs(int(t_type[1])); find_n=find_n+1
             elif t_type[0]=='SUPPORT': support = int(t_type[1]); find_n=find_n+1
             elif t_type[0]=='END': sv_end  = int(t_type[1])-1; find_n=find_n+1
             elif t_type[0]=='STDEV_POS': STDEV_POS= float(t_type[1]); find_n=find_n+1
           if (sv_end==None and (not svlen==None)) or ((not sv_end==None) and svlen==None) or (STDEV_POS==None): # sv_end==None or svlen==None:
             print("No end or svlen", sv, svtype, svlen, sv_end, STDEV_POS)
           if svlen==None:
             svlen = 10000

           if support<5: continue;
           if svlen<50:
              print("Warning: SV length less than 50", sv)

           if not (chrom==cur_chr): 
              cur_chr = chrom
              chr_ind = 0;
           if sv_end==None: 
              sv_end = pos
              breakpoints = [[pos,pos]]
           else:
              breakpoints = [[pos,pos], [sv_end, sv_end]]

           if svtype in Def_sv_types[:-1]: SVTYPE = svtype
           else: 
              SVTYPE = 'OTHER'
       
           if svtype not in sv_count_dict[sample_name]:
              sv_count_dict[sample_name][svtype] = 0
              sv_base_dict[sample_name][svtype] = 0
           sv_count_dict[sample_name][svtype] = sv_count_dict[sample_name][svtype] + 1
           sv_base_dict[sample_name][svtype] = sv_base_dict[sample_name][svtype] + svlen

           while chr_ind<len(segments_by_chr[cur_chr]) and pos > segments_by_chr[cur_chr][chr_ind][1]:
              chr_ind = chr_ind + 1
           
           add_segi = 0
           for bk1 in breakpoints:
              pbs = bk1[0]
              svbend = bk1[1]
              used_ind = chr_ind
              while used_ind<len(segments_by_chr[cur_chr]):
                 seg_start, seg_end = segments_by_chr[cur_chr][used_ind][0], segments_by_chr[cur_chr][used_ind][1]
                 if svbend>=seg_start and pbs<=seg_end:
                    segment_key = f"{chrom}:{seg_start}-{seg_end}"
                    if segment_key not in sample_sv_counts[SVTYPE][sample_name]:
                        sample_sv_counts[SVTYPE][sample_name][segment_key] = 0
                        sample_sv_bases[SVTYPE][sample_name][segment_key] = 0
                    sample_sv_counts[SVTYPE][sample_name][segment_key] += 1
                    sample_sv_bases[SVTYPE][sample_name][segment_key] += svlen
                    if SVTYPE not in sv_add_dict[sample_name]:
                       sv_add_dict[sample_name][SVTYPE] = 0
                    sv_add_dict[sample_name][SVTYPE] += 1
                    add_segi += 1
                 elif svbend<seg_start: break;
                 used_ind = used_ind + 1
           if svtype not in Def_sv_types[:-2]: 
              if add_segi>4 and used_ind<len(segments_by_chr[cur_chr]) and segments_by_chr[cur_chr][used_ind][1]-segments_by_chr[cur_chr][chr_ind][0]>5000:
                  print('\t', SVTYPE, 'add_segi='+str(add_segi), chr_ind, used_ind, sv[:2], sv_end); #sv[7:])
           testi += 1
       print("process: ", sample_name, testi, "{:.2f}".format( time.time()- start_time ), sv_count_dict[sample_name], sv_base_dict[sample_name], sv_add_dict[sample_name]); sys.stdout.flush()
       start_time = time.time();
   
   segments_list = []
   for _t_row in genome_segments[1]:
      if _t_row==genome_segments[1][0]: print("_t_row", _t_row )
      segments_list.append( "{}:{}-{}".format( _t_row[0], _t_row[1], _t_row[2]  ) )

   filtered_counts_df = {}
   filtered_bases_df = {}
   for _ad_sv_type in Def_sv_types:
      count_df = pd.DataFrame(index=segments_list, columns=[vcf.split("/")[-sample_ind] for vcf in vcf_files]).fillna(0)
      # Fill DataFrame with computed values
      for sample, segment_data in sample_sv_counts[_ad_sv_type].items():
          start_time = time.time()
          for segment, sv_bases in segment_data.items():
              count_df.loc[segment, sample] = sv_bases
      #print("process: ",sample, "{:.2f}".format( time.time()- start_time )); sys.stdout.flush()
      print('count_df:', count_df.shape, _ad_sv_type, "{:.2f}".format( time.time()- start_time ) ); sys.stdout.flush()
      filtered_count_df = count_df.drop(index=count_df.index[(count_df == 0).all(axis=1)])
      print('filtered_df:', filtered_count_df.shape)
      for sample in sample_sv_counts[_ad_sv_type]:
         filtered_count_df[sample] = filtered_count_df[sample].astype(int)
      pl.from_pandas(filtered_count_df.reset_index()).write_csv( output_file+_ad_sv_type+'_counts.csv' );
      filtered_counts_df[ _ad_sv_type ] = filtered_count_df

      base_df = pd.DataFrame(index=segments_list, columns=[vcf.split("/")[-sample_ind] for vcf in vcf_files]).fillna(0)
      # Fill DataFrame with computed values
      for sample, segment_data in sample_sv_bases[_ad_sv_type].items():
          start_time = time.time()
          for segment, sv_bases in segment_data.items():
              base_df.loc[segment, sample] = sv_bases
      #print("process: ",sample, "{:.2f}".format( time.time()- start_time )); sys.stdout.flush()
      print('base_df:', base_df.shape, _ad_sv_type, "{:.2f}".format( time.time()- start_time ) ); sys.stdout.flush()
      filtered_base_df = base_df.drop(index=base_df.index[(base_df == 0).all(axis=1)])
      print('filtered_base_df:', filtered_base_df.shape)
      for sample in sample_sv_bases[_ad_sv_type]:
         filtered_base_df[sample] = filtered_base_df[sample].astype(int)
      pl.from_pandas(filtered_base_df.reset_index()).write_csv( output_file+_ad_sv_type+'_bases.csv' );
      filtered_bases_df[ _ad_sv_type ] =  filtered_base_df
   
   print(f"SV matrix saved to {output_file}")
   return filtered_counts_df, filtered_bases_df


def sv_differential(df, af_list, sample_category, csvfile):
   pandas2ri.activate()
   limma = importr('limma')
   df_r = pandas2ri.py2rpy( df[af_list] )
   
   design_matrix = pd.DataFrame({
          'Group': [sample_category[0].loc['_'.join(_cn.split('_')[:(2 if len(_cn.split('_'))>=2 else 1)]), 'AD'] for _cn in af_list] ,
          'Intercept': [1] * (len(af_list)),
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
   #print(summary)

   fit2 = limma.eBayes(fit, stdev_coef_lim= r['c'](0.1, 10) )

   results = limma.topTable(fit2, adjust_method='BH', number=df.shape[0])
   results_df = pandas2ri.rpy2py(results)

   df = pd.merge( df, results_df, left_index=True, right_index=True) ;
   
   df.sort_values(by='adj.P.Val', inplace=True);

   pl.from_pandas(df.reset_index()).write_csv( csvfile );

   return df;


def compute_effect_size(df, group1, group2):
    g1_mean = df[group1].mean(axis=1)
    g2_mean = df[group2].mean(axis=1)
    g1_std = df[group1].std(axis=1)
    g2_std = df[group2].std(axis=1)
    
    # Pooled std (average of two stds)
    pooled_std = np.sqrt((g1_std**2 + g2_std**2) / 2)
    effect_size = (g1_mean - g2_mean) / pooled_std.replace(0, np.nan)
    
    return effect_size.abs()  # Use absolute effect size


if __name__=='__main__':
   t_op = sys.argv[1]
   group_file = sys.argv[2]; #3]; 
   sample_list = read_sample_list(group_file)

   sup_str = sys.argv[3].split(',')
   min_sup_thr = float(sup_str[0]); #5]);
   min_uniqsup_thr = min_sup_thr
   if len(sup_str)>1:
     min_uniqsup_thr = float(sup_str[1]);

   if t_op in [1, '1']:
      segment_file = sys.argv[4]; #2];
      segment_bin = int(sys.argv[5]); #4])

      extend_bp = 20
      orderedchr, orderedchrnum = get_chr_human(contain_xy=True);
      unusedchr=['chrY', 'chrM']
      if extend_bp < segment_bin*0.5:
         extend_bp = int(segment_bin*0.5 + 0.5)
  
      start_time = time.time();
      genome_segments = _segment_genome(segment_file, orderedchr, rep_ext_bp=20, seg_ext_bp=extend_bp, seg_bin=segment_bin, unusedchr=unusedchr)
      print("segment: time: {:.2f}".format( time.time()- start_time ))

   if t_op in [1, '1']:
      outp_f_format = ( sys.argv[6]+('ovlpBK2025_p{}_{}_s{}'.format(extend_bp,segment_bin, len(sys.argv[8:]))) )
      sample_id = int(sys.argv[7])
      vcf_file_list = sys.argv[8:]
      print("Operation: from separate vcf files")
      print("segment_file=", segment_file);
      print("group_file=", group_file);
      print("segment_bin={}; min_sup_thr={}; extend_bp={} sample_id={}".format(segment_bin, min_sup_thr, extend_bp, sample_id))
      print("outp_f_format=", outp_f_format)
      print("first VCF=", vcf_file_list[:2])
      svs_cout_df, svs_base_df =  _read_sv_(genome_segments, vcf_file_list, outp_f_format, sample_id, orderedchr, unusedchr)
   else:
      outp_f_format = sys.argv[4]
      print("Operation: from merged files")
      print("group_file=", group_file);
      print("min_sup_thr={} ".format(min_sup_thr ))
      print("outp_f_format=", outp_f_format)
      print("from files (count, then base)", sys.argv[5], sys.argv[6], sys.argv[7]);

      sv_ts = sys.argv[5].split(',');
      cout_files = sys.argv[6].split(',');
      base_files = sys.argv[7].split(',');
      segment_bin = int(sys.argv[8])
      extend_bp = 20
      if extend_bp < segment_bin*0.5: extend_bp = int(segment_bin*0.5 + 0.5)
      if len(sys.argv)>10:
         sample_id = int(sys.argv[9])
         vcf_file_list = sys.argv[10:]
      else:
         sample_id = None;
         vcf_file_list = None;
      svs_cout_df = {}
      svs_base_df = {}
      for _c_f_ind in range(len(cout_files)):
         print("read", sv_ts[_c_f_ind], cout_files[_c_f_ind], base_files[_c_f_ind])
         svs_cout_df[sv_ts[_c_f_ind]] = pl.read_csv( cout_files[_c_f_ind], separator=',', infer_schema_length=100000).to_pandas();
         svs_cout_df[sv_ts[_c_f_ind]].set_index(['index'], inplace=True)

         svs_base_df[sv_ts[_c_f_ind]] = pl.read_csv(base_files[_c_f_ind], separator=',', infer_schema_length=100000).to_pandas();
         svs_base_df[sv_ts[_c_f_ind]].set_index(['index'], inplace=True)

   adsv_df_list, ctsv_df_list, other_c = get_groupd_columns(svs_cout_df[ list(svs_cout_df.keys())[0] ].columns, sample_list)
   print("AD samples", adsv_df_list)
   print("Ctrl samples", ctsv_df_list)

   if not sample_id==None:
      stdev_data = {}
      for c_vcf in vcf_file_list:
         stdev_data[c_vcf.split('/')[-sample_id]] = read_vcf_to_dataframe_dynamic_columns(c_vcf, orderedchr, unusedchr).groupby(('CHROM'))
   else: stdev_data = {}
 
   cd_thr = 0.5
   final_sv_res = {}
   for _ak in ["chr", "start", "end", "sample", 'svtype', 'svlen', 'group']:
      final_sv_res[ _ak ] = []
   for _sv_type in svs_cout_df:
      outp_ft_format = outp_f_format+_sv_type
      sv_cout_df = svs_cout_df[_sv_type]
      sv_base_df = svs_base_df[_sv_type]
      print("SV type", _sv_type)
      #print("Check:", sv_cout_df[adsv_df_list].sum(), sv_cout_df[ctsv_df_list].sum() )
      #print("bf", sv_cout_df.shape, sv_base_df.shape)

      sv_cout_dfadsv_NonZero = (sv_cout_df[adsv_df_list] != 0).sum(axis=1)
      sv_cout_dfctsv_NonZero = (sv_cout_df[ctsv_df_list] != 0).sum(axis=1)
      sv_cout_dfadsv_Percent = sv_cout_dfadsv_NonZero / len(adsv_df_list)
      sv_cout_dfctsv_Percent = sv_cout_dfctsv_NonZero / len(ctsv_df_list)
      group_filter = ((sv_cout_dfadsv_Percent >= min_sup_thr) | (sv_cout_dfctsv_Percent >= min_sup_thr))
      #print("af_group_filter", sv_cout_df[group_filter].shape, sv_base_df[group_filter].shape)  

      if t_op not in [1, '1', 0, '0']:
         sv_cout_dif = sv_cout_df
         sv_base_dif = sv_base_df
      else:
         #print(sv_cout_df.shape); print(sv_cout_df.dtypes)
         dif_count_filter = compute_effect_size(sv_cout_df, adsv_df_list, ctsv_df_list);
         #print(sv_base_df.shape); print(sv_base_df.dtypes)
         dif_bases_filter = compute_effect_size(sv_base_df, adsv_df_list, ctsv_df_list);

         #print( (group_filter & (dif_count_filter>= cd_thr)).sum() );
         #print( (group_filter & (dif_bases_filter>= cd_thr)).sum() );

         #sv_cout_dif = sv_differential(copy.deepcopy(sv_cout_df[group_filter]), adsv_df_list+ctsv_df_list, sample_list, outp_ft_format+'_counts_dif.csv')
         sv_cout_dif = sv_differential(copy.deepcopy(sv_cout_df[ (group_filter & (dif_count_filter>= cd_thr)) ]), adsv_df_list+ctsv_df_list, sample_list, outp_ft_format+'_counts_dif.csv')
         #print(sv_cout_dif)

         #sv_base_dif = sv_differential(copy.deepcopy(sv_base_df[group_filter]), adsv_df_list+ctsv_df_list, sample_list, outp_ft_format+'_base_dif.csv')
         sv_base_dif = sv_differential(copy.deepcopy(sv_base_df[ (group_filter & (dif_bases_filter>= cd_thr)) ]), adsv_df_list+ctsv_df_list, sample_list, outp_ft_format+'_base_dif.csv')
         #print(sv_base_dif)

         '''if 'chr12:37240944-37241899' in sv_cout_dif.index:
            m_ind = sv_cout_dif.index.get_loc('chr12:37240944-37241899')
            print( sv_cout_dif.iloc[ m_ind-1 ] )
            print( sv_cout_dif.iloc[ m_ind ] )
            print( sv_cout_dif.iloc[ m_ind+1 ] )
         if 'chr12:37240944-37241899' in sv_base_dif.index:
            #print( sv_base_dif.loc['chr12:37240944-37241899'] )
            m_ind = sv_base_dif.index.get_loc('chr12:37240944-37241899')
            print( sv_base_dif.iloc[ m_ind-1 ] )
            print( sv_base_dif.iloc[ m_ind ] )
            print( sv_base_dif.iloc[ m_ind+1 ] ) '''

         sv_cout_sum = [sv_cout_df[ ((sv_cout_dfadsv_Percent >= min_sup_thr) & (sv_cout_dfctsv_Percent <0.01)) | ((sv_cout_dfadsv_Percent <0.01) & (sv_cout_dfctsv_Percent >= min_sup_thr)) ][adsv_df_list].sum(), sv_cout_df[ ((sv_cout_dfadsv_Percent >= min_sup_thr) & (sv_cout_dfctsv_Percent <0.01)) | ((sv_cout_dfadsv_Percent <0.01) & (sv_cout_dfctsv_Percent >= min_sup_thr)) ][ctsv_df_list].sum()]
         sv_base_sum = [sv_base_df[ ((sv_cout_dfadsv_Percent >= min_sup_thr) & (sv_cout_dfctsv_Percent <0.01)) | ((sv_cout_dfadsv_Percent <0.01) & (sv_cout_dfctsv_Percent >= min_sup_thr)) ][adsv_df_list].sum(), sv_base_df[ ((sv_cout_dfadsv_Percent >= min_sup_thr) & (sv_cout_dfctsv_Percent <0.01)) | ((sv_cout_dfadsv_Percent <0.01) & (sv_cout_dfctsv_Percent >= min_sup_thr)) ][ctsv_df_list].sum()]
  
         sv_base_sum[0] = np.where(sv_base_sum[0] > 10000000, 10000000, sv_base_sum[0])
         sv_base_sum[1] = np.where(sv_base_sum[1] > 10000000, 10000000, sv_base_sum[1])

         fig, axes = plt.subplots(1, 2, figsize=(10, 5))
         axes[0].violinplot([sv_cout_sum[0].values.tolist(), sv_cout_sum[1].values.tolist()]); #, labels=["AD", "Ctrl"])
         axes[0].set_xticks([1, 2])  # Set x-axis positions
         axes[0].set_xticklabels(["AD", "Ctrl"])  # Set x-axis labels
         axes[0].set_title("SV Count Sum")
         axes[1].violinplot([sv_base_sum[0], sv_base_sum[1]]); #labels=["AD", "Ctrl"])
         axes[1].set_xticks([1, 2])
         axes[1].set_xticklabels(["AD", "Ctrl"])
         axes[1].set_title("SV Base Sum")
         # Show plots
         plt.tight_layout()
         plt.savefig(outp_ft_format+"_unique.png", dpi=600, bbox_inches='tight')
         plt.close('all')
  
      sv_cout_dfadsv_NonZero = (sv_cout_dif[adsv_df_list] != 0).sum(axis=1)
      sv_cout_dfctsv_NonZero = (sv_cout_dif[ctsv_df_list] != 0).sum(axis=1)
      sv_cout_dfadsv_Percent = sv_cout_dfadsv_NonZero / len(adsv_df_list)
      sv_cout_dfctsv_Percent = sv_cout_dfctsv_NonZero / len(ctsv_df_list)
      uniq_index =( ((sv_cout_dfadsv_Percent >= min_uniqsup_thr) & (sv_cout_dfctsv_Percent <0.01)) | ((sv_cout_dfadsv_Percent <0.01) & (sv_cout_dfctsv_Percent >= min_uniqsup_thr)) )
      res_cout_df = sv_cout_dif[ uniq_index | ( (sv_cout_dif['adj.P.Val']<=0.05) & ( (sv_cout_dfctsv_Percent<0.15) | (sv_cout_dfadsv_Percent<0.15) ) ) ]

      sv_base_dfadsv_NonZero = (sv_base_dif[adsv_df_list] != 0).sum(axis=1)
      sv_base_dfctsv_NonZero = (sv_base_dif[ctsv_df_list] != 0).sum(axis=1)
      sv_base_dfadsv_Percent = sv_base_dfadsv_NonZero / len(adsv_df_list)
      sv_base_dfctsv_Percent = sv_base_dfctsv_NonZero / len(ctsv_df_list)
      res_base_df = sv_base_dif[ ( sv_base_dif.index.isin( sv_cout_dif[uniq_index].index )  ) | ( (sv_base_dif['adj.P.Val']<=0.05)  & ( (sv_base_dfadsv_Percent<0.15) | (sv_base_dfctsv_Percent<0.15) ) ) ]

      merg_indexes = [ res_cout_df.index, res_base_df[ ~(res_base_df.index.isin( res_cout_df.index ) ) ].index ]
      merged_loci = {}
      for _t_index in merg_indexes:
         for c_ind_loc in _t_index:
            c_loc = c_ind_loc.split(':')
            c_loc[1] = c_loc[1].split('-')
            if c_loc[0] not in merged_loci: merged_loci[ c_loc[0] ] = []
            merged_loci[ c_loc[0] ].append( [int(c_loc[1][0])-extend_bp, int(c_loc[1][1])+extend_bp] )

      save_sv_keys = {}
      g_2 = [adsv_df_list, ctsv_df_list]
      for c_k in merged_loci: # for each chromoseom
         merged_loci[ c_k ] = sorted(merged_loci[ c_k ])
         for pos_inf in merged_loci[ c_k ]: # for each region;
            for _gi in range(len(g_2)): #for two groups;
               for _s_k in g_2[_gi]: # for each sample
                  chr_df = stdev_data[_s_k].get_group(c_k);
                  _in_region = chr_df[((chr_df['POS'] <= pos_inf[1]) & (chr_df['POS'] >= pos_inf[0])) | ((chr_df['SVEND'] <= pos_inf[1]) & (chr_df['SVEND'] >= pos_inf[0]))  ]
                  if _sv_type=='OTHER':
                     _in_region = _in_region[_in_region['SVTYPE'].isin(Def_sv_types[4:])]
                  else:  _in_region = _in_region[_in_region['SVTYPE']==_sv_type]
                  if not (_in_region.empty):
                     for _, t_region in _in_region.iterrows():
                        if (_s_k, c_k, t_region['POS'], t_region['SVEND']) not in save_sv_keys:
                           #if c_k=='chr12' and ((t_region['POS']>=37240942 and t_region['SVEND']<=37245721)):# or ( t_region['POS']>=37245699 and t_region['SVEND']<=37245716 )):
                           #   print ("need check", c_k+':'+str(pos_inf[0])+'-'+str(pos_inf[1]) )
                           #final_sv_res.append([ c_k, t_region['POS'], t_region['SVEND'], _s_k, _in_region['SVTYPE'], _in_region['SVLEN' ] ])
                           final_sv_res["chr"].append( c_k )
                           final_sv_res["start"].append( t_region['POS'] )
                           final_sv_res["end"].append( t_region['SVEND'] )
                           final_sv_res["sample"].append( _s_k )
                           final_sv_res["svtype"].append( t_region['SVTYPE'] )
                           final_sv_res["svlen"].append( t_region['SVLEN' ] )
                           final_sv_res["group"].append("AD" if _gi==0 else "Ctrl")
                           save_sv_keys[ (_s_k, c_k, t_region['POS'], t_region['SVEND']) ] = 0

   #final_sv_df = pd.DataFrame(final_sv_res, columns=["chr", "start", "end", "sample", 'svtype', 'svlen']).sort_values(by=['chr', 'start', 'end'])
   final_sv_df = pd.DataFrame(final_sv_res ).sort_values(by=['chr', 'start', 'end'])
   final_sv_df["start"] = final_sv_df["start"].astype(int)
   final_sv_df["end"] = final_sv_df["end"].astype(int)
   final_sv_df["chr"] = final_sv_df["chr"].astype(str)
   final_sv_df["sample"] = final_sv_df["sample"].astype(str)
   final_sv_df['svtype'] = final_sv_df["svtype"].astype(str)
   final_sv_df['svlen'] = final_sv_df["svlen"].astype(int)

   rem_ind_list = []
   res_region = []
   for mchr, t_svs in final_sv_df.groupby(('chr')):
      prev_i = 0;
      start_i = -1;
      for _t_i in range(len(t_svs)):
          if t_svs.iloc[_t_i]['start'] - t_svs.iloc[prev_i]['start']<extend_bp:
             if start_i==-1: start_i = prev_i
             prev_i = _t_i
          else:
             group_size = t_svs.iloc[start_i:_t_i].groupby(("group"))['sample'].nunique()
             if 'AD' in group_size:
                t_ad_p = group_size['AD'] / len(adsv_df_list);
             else: t_ad_p = 0
             if 'Ctrl' in group_size:
                t_ct_p = group_size['Ctrl'] / len(ctsv_df_list);
             else: t_ct_p = 0
             if (t_ad_p<0.15 and t_ct_p>=min_sup_thr) or (t_ad_p>=min_sup_thr and t_ct_p<0.15):
                print(mchr, t_svs.iloc[start_i]['start'], t_svs.iloc[prev_i]['start'] );
                res_region.append([ mchr, t_svs.iloc[start_i:_t_i]['start'].median(), t_svs.iloc[start_i:_t_i]['end'].median(),t_svs.iloc[start_i:_t_i]['start'].min(), t_svs.iloc[start_i:_t_i]['start'].max(), t_svs.iloc[start_i:_t_i]['end'].min(), t_svs.iloc[start_i:_t_i]['end'].max(), t_ad_p, t_ct_p, '/'.join(f"{k}:{v}" for k, v in (t_svs.iloc[start_i:_t_i]['svtype'].value_counts().to_dict().items()) ) ])
                print(t_svs.iloc[start_i:_t_i] );
             else:
                rem_ind_list.extend( t_svs.iloc[start_i:_t_i].index.tolist())
             prev_i = _t_i
             start_i = prev_i
      group_size = t_svs.iloc[start_i:_t_i].groupby(("group"))['sample'].nunique()
      if 'AD' in group_size:
         t_ad_p = group_size['AD'] / len(adsv_df_list);
      else: t_ad_p = 0
      if 'Ctrl' in group_size:
        t_ct_p = group_size['Ctrl'] / len(ctsv_df_list);
      else: t_ct_p = 0
      if (t_ad_p<0.15 and t_ct_p>=min_sup_thr) or (t_ad_p>=min_sup_thr and t_ct_p<0.15):
         print(mchr, t_svs.iloc[start_i]['start'], t_svs.iloc[prev_i]['start'] );
         #res_region.append([ mchr, t_svs.iloc[start_i:_t_i]['start'].min(), t_svs.iloc[start_i:_t_i]['start'].max(), t_svs.iloc[start_i:_t_i]['end'].min(), t_svs.iloc[start_i:_t_i]['end'].max(), t_ad_p, t_ct_p, '/'.join(f"{k}:{v}" for k, v in (t_svs.iloc[start_i:_t_i]['svtype'].value_counts().to_dict().items()) ) ])
         res_region.append([ mchr, t_svs.iloc[start_i:_t_i]['start'].median(), t_svs.iloc[start_i:_t_i]['end'].median(),t_svs.iloc[start_i:_t_i]['start'].min(), t_svs.iloc[start_i:_t_i]['start'].max(), t_svs.iloc[start_i:_t_i]['end'].min(), t_svs.iloc[start_i:_t_i]['end'].max(), t_ad_p, t_ct_p, '/'.join(f"{k}:{v}" for k, v in (t_svs.iloc[start_i:_t_i]['svtype'].value_counts().to_dict().items()) ) ])
         print(t_svs.iloc[start_i:_t_i] );
      else:
         rem_ind_list.extend( t_svs.iloc[start_i:_t_i].index.tolist())

   outp_ft_format = outp_f_format+'_res'
   pl.from_pandas( pd.DataFrame(res_region, columns=["chr", "start", "end", "start_min", "start_max", "end_min", "end_max", "AD_freq", 'Ctrl_freq', 'svinfo']).sort_values(by=['chr', 'start'])  ).write_csv( outp_ft_format+"_region.csv" );
   print(final_sv_df.shape)
   final_sv_df = final_sv_df.drop(rem_ind_list)
   print(final_sv_df.shape)
   pl.from_pandas( final_sv_df ).write_csv( outp_ft_format+".csv" );
   generate_sv_count_violinplots(final_sv_df, adsv_df_list, ctsv_df_list, outp_ft_format)
   generate_svlen_violinplots(final_sv_df, adsv_df_list, ctsv_df_list, outp_ft_format)

