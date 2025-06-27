
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
           if sv_end==None: sv_end = pos

           if svtype in Def_sv_types[:-1]: SVTYPE = svtype
           else: 
              #print("OTHER: ", sv)
              SVTYPE = 'OTHER'
       
           if svtype not in sv_count_dict[sample_name]:
              sv_count_dict[sample_name][svtype] = 0
              sv_base_dict[sample_name][svtype] = 0
           sv_count_dict[sample_name][svtype] = sv_count_dict[sample_name][svtype] + 1
           sv_base_dict[sample_name][svtype] = sv_base_dict[sample_name][svtype] + svlen

           while chr_ind<len(segments_by_chr[cur_chr]) and pos > segments_by_chr[cur_chr][chr_ind][1]:
              chr_ind = chr_ind + 1
           used_ind = chr_ind
           add_segi = 0
           while used_ind<len(segments_by_chr[cur_chr]):
              seg_start, seg_end = segments_by_chr[cur_chr][used_ind][0], segments_by_chr[cur_chr][used_ind][1]
              if sv_end>=seg_start and pos<=seg_end:
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
              elif sv_end<seg_start: break;
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

   min_sup_thr = float(sys.argv[3]); #5]);

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
      outp_f_format = ( sys.argv[6]+('ovlpSep_p{}_{}_s{}'.format(extend_bp,segment_bin, len(sys.argv[8:]))) )
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
      if len(sys.argv)>9:
         sample_id = int(sys.argv[8])
         vcf_file_list = sys.argv[9:]
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
   for _sv_type in svs_cout_df:
      outp_ft_format = outp_f_format+_sv_type
      sv_cout_df = svs_cout_df[_sv_type]
      sv_base_df = svs_base_df[_sv_type]
      print("SV type", _sv_type)
      print("bf", sv_cout_df.shape, sv_base_df.shape)

      sv_cout_dfadsv_NonZero = (sv_cout_df[adsv_df_list] != 0).sum(axis=1)
      sv_cout_dfctsv_NonZero = (sv_cout_df[ctsv_df_list] != 0).sum(axis=1)
      sv_cout_dfadsv_Percent = sv_cout_dfadsv_NonZero / len(adsv_df_list)
      sv_cout_dfctsv_Percent = sv_cout_dfctsv_NonZero / len(ctsv_df_list)
      group_filter = ((sv_cout_dfadsv_Percent >= min_sup_thr) | (sv_cout_dfctsv_Percent >= min_sup_thr))
      print("af_group_filter", sv_cout_df[group_filter].shape, sv_base_df[group_filter].shape)  

      if t_op not in [1, '1', 0, '0']:
         sv_cout_dif = sv_cout_df
         sv_base_dif = sv_base_df
      else:
         dif_count_filter = compute_effect_size(sv_cout_df, adsv_df_list, ctsv_df_list);
         dif_bases_filter = compute_effect_size(sv_base_df, adsv_df_list, ctsv_df_list);

         print( (group_filter & (dif_count_filter>= cd_thr)).sum() );
         print( (group_filter & (dif_bases_filter>= cd_thr)).sum() );

         sv_cout_dif = sv_differential(copy.deepcopy(sv_cout_df[ (group_filter & (dif_count_filter>= cd_thr)) ]), adsv_df_list+ctsv_df_list, sample_list, outp_ft_format+'_counts_dif.csv')
         print(sv_cout_dif)

         sv_base_dif = sv_differential(copy.deepcopy(sv_base_df[ (group_filter & (dif_bases_filter>= cd_thr)) ]), adsv_df_list+ctsv_df_list, sample_list, outp_ft_format+'_base_dif.csv')
         print(sv_base_dif)

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
      uniq_index =( ((sv_cout_dfadsv_Percent >= min_sup_thr) & (sv_cout_dfctsv_Percent <0.01)) | ((sv_cout_dfadsv_Percent <0.01) & (sv_cout_dfctsv_Percent >= min_sup_thr)) )
      res_cout_df = sv_cout_dif[ uniq_index | (sv_cout_dif['adj.P.Val']<=0.05) ]
      res_base_df = sv_base_dif[ ( sv_base_dif.index.isin( sv_cout_dif[uniq_index].index ) ) | (sv_base_dif['adj.P.Val']<=0.05) ]



