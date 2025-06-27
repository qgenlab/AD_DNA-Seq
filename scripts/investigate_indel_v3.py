
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

def _read_indel_(genome_segments, vcf_files, output_file, sample_ind, orderedchr, unusedchr):
   # Process INDELs from VCF files
   sample_indel_counts = {}
   sample_indel_bases = {}
  
   for _ad_indel_type in Def_indel_types:
      sample_indel_counts[_ad_indel_type] = {}
      sample_indel_bases[_ad_indel_type] = {} 
   segments_by_chr = genome_segments[2]
   start_time = time.time();

   indel_count_dict = {}
   indel_base_dict = {}

   indel_add_dict = {}

   testi = 0; 
   for vcf_file in vcf_files:
       #print("process: ", vcf_file); sys.stdout.flush()
       sample_name = vcf_file.split("/")[-sample_ind]  # Extract sample name
       for _ad_indel_type in Def_indel_types:
          sample_indel_counts[_ad_indel_type][sample_name] = {}
          sample_indel_bases[_ad_indel_type][sample_name] = {}

       indel_count_dict[sample_name] = {}
       indel_base_dict[sample_name] = {}

       indel_add_dict[sample_name] = {}

       cur_chr = ''
       chr_ind = 0;
       with open(vcf_file, 'r') as vcfread:
          lines = vcfread.readlines()
       for indel in lines:
           indel = indel.strip();
           if len(indel)==0: continue;
           if indel[0]=='#': continue;
           indel = indel.split('\t')

           chrom, pos, m_filter = indel[0], int(indel[1])-1, indel[6]
           if (not unusedchr==None) and chrom in unusedchr: continue;
           if chrom  not in orderedchr: continue;
           if not m_filter=='PASS': continue;          

           ref_seq, alt_seq = indel[3], indel[4]
           alt_seq = alt_seq.split(',')

           if not chrom==cur_chr: 
              cur_chr = chrom
              chr_ind = 0;

           if ',' in ref_seq:
              print('Warning mutiple sequences in Ref', ref_seq);
           for _alt_s in alt_seq:
              if len(_alt_s)<len(ref_seq):
                 indeltype = 'del'
                 indellen = len(ref_seq)-len(_alt_s)
                 indel_end = pos + indellen
              elif len(_alt_s)>len(ref_seq):
                 indeltype = 'ins'
                 indellen = len(_alt_s)-len(ref_seq)
                 indel_end = pos
              else:
                 if ',' not in indel[4]:
                    print("Alt is equal REF", indel)
                 continue;

              INDELTYPE = indeltype
              if indeltype not in indel_count_dict[sample_name]:
                 indel_count_dict[sample_name][indeltype] = 0
                 indel_base_dict[sample_name][indeltype] = 0
              indel_count_dict[sample_name][indeltype] = indel_count_dict[sample_name][indeltype] + 1
              indel_base_dict[sample_name][indeltype] = indel_base_dict[sample_name][indeltype] + indellen

              while chr_ind<len(segments_by_chr[cur_chr]) and pos > segments_by_chr[cur_chr][chr_ind][1]:
                 chr_ind = chr_ind + 1
              used_ind = chr_ind
              while used_ind<len(segments_by_chr[cur_chr]):
                 seg_start, seg_end = segments_by_chr[cur_chr][used_ind][0], segments_by_chr[cur_chr][used_ind][1]
                 if indel_end>=seg_start and pos<=seg_end:
                    segment_key = f"{chrom}:{seg_start}-{seg_end}"
                    if segment_key not in sample_indel_counts[INDELTYPE][sample_name]:
                        sample_indel_counts[INDELTYPE][sample_name][segment_key] = 0
                        sample_indel_bases[INDELTYPE][sample_name][segment_key] = 0
                    sample_indel_counts[INDELTYPE][sample_name][segment_key] += 1
                    sample_indel_bases[INDELTYPE][sample_name][segment_key] += indellen
                    if INDELTYPE not in indel_add_dict[sample_name]:
                       indel_add_dict[sample_name][INDELTYPE] = 0
                    indel_add_dict[sample_name][INDELTYPE] += 1
                 elif indel_end<seg_start: break;
                 used_ind = used_ind + 1
           testi += 1
       print("process: ", sample_name, testi, "{:.2f}".format( time.time()- start_time ), indel_count_dict[sample_name], indel_base_dict[sample_name], indel_add_dict[sample_name]); sys.stdout.flush()
       start_time = time.time();
   
   segments_list = []
   for _t_row in genome_segments[1]:
      if _t_row==genome_segments[1][0]: print("_t_row", _t_row )
      segments_list.append( "{}:{}-{}".format( _t_row[0], _t_row[1], _t_row[2]  ) )

   filtered_counts_df = {}
   filtered_bases_df = {}
   for _ad_indel_type in Def_indel_types:
      count_df = pd.DataFrame(index=segments_list, columns=[vcf.split("/")[-sample_ind] for vcf in vcf_files]).fillna(0)
      # Fill DataFrame with computed values
      for sample, segment_data in sample_indel_counts[_ad_indel_type].items():
          start_time = time.time()
          for segment, indel_bases in segment_data.items():
              count_df.loc[segment, sample] = indel_bases
      #print("process: ",sample, "{:.2f}".format( time.time()- start_time )); sys.stdout.flush()
      print('count_df:', count_df.shape, _ad_indel_type, "{:.2f}".format( time.time()- start_time ) ); sys.stdout.flush()
      filtered_count_df = count_df.drop(index=count_df.index[(count_df == 0).all(axis=1)])
      print('filtered_df:', filtered_count_df.shape)
      for sample in sample_indel_counts[_ad_indel_type]:
         filtered_count_df[sample] = filtered_count_df[sample].astype(int)
      pl.from_pandas(filtered_count_df.reset_index()).write_csv( output_file+_ad_indel_type+'_counts.csv' );
      filtered_counts_df[ _ad_indel_type ] = filtered_count_df

      base_df = pd.DataFrame(index=segments_list, columns=[vcf.split("/")[-sample_ind] for vcf in vcf_files]).fillna(0)
      # Fill DataFrame with computed values
      for sample, segment_data in sample_indel_bases[_ad_indel_type].items():
          start_time = time.time()
          for segment, indel_bases in segment_data.items():
              base_df.loc[segment, sample] = indel_bases
      #print("process: ",sample, "{:.2f}".format( time.time()- start_time )); sys.stdout.flush()
      print('base_df:', base_df.shape, _ad_indel_type, "{:.2f}".format( time.time()- start_time ) ); sys.stdout.flush()
      filtered_base_df = base_df.drop(index=base_df.index[(base_df == 0).all(axis=1)])
      print('filtered_base_df:', filtered_base_df.shape)
      for sample in sample_indel_bases[_ad_indel_type]:
         filtered_base_df[sample] = filtered_base_df[sample].astype(int)
      pl.from_pandas(filtered_base_df.reset_index()).write_csv( output_file+_ad_indel_type+'_bases.csv' );
      filtered_bases_df[ _ad_indel_type ] =  filtered_base_df
   
   print(f"INDEL matrix saved to {output_file}")
   return filtered_counts_df, filtered_bases_df


def indel_differential(df, af_list, sample_category, csvfile):
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
      genome_segments = _segment_genome(segment_file, orderedchr, rep_ext_bp=10, seg_ext_bp=extend_bp, seg_bin=segment_bin, unusedchr=unusedchr)
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
      indels_cout_df, indels_base_df =  _read_indel_(genome_segments, vcf_file_list, outp_f_format, sample_id, orderedchr, unusedchr)
   else:
      outp_f_format = sys.argv[4]
      print("Operation: from merged files")
      print("group_file=", group_file);
      print("min_sup_thr={} ".format(min_sup_thr ))
      print("outp_f_format=", outp_f_format)
      print("from files (count, then base)", sys.argv[5], sys.argv[6], sys.argv[7]);

      indel_ts = sys.argv[5].split(',');
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
      indels_cout_df = {}
      indels_base_df = {}
      for _c_f_ind in range(len(cout_files)):
         print("read", indel_ts[_c_f_ind], cout_files[_c_f_ind], base_files[_c_f_ind])
         indels_cout_df[indel_ts[_c_f_ind]] = pl.read_csv( cout_files[_c_f_ind], separator=',', infer_schema_length=100000).to_pandas();
         indels_cout_df[indel_ts[_c_f_ind]].set_index(['index'], inplace=True)

         indels_base_df[indel_ts[_c_f_ind]] = pl.read_csv(base_files[_c_f_ind], separator=',', infer_schema_length=100000).to_pandas();
         indels_base_df[indel_ts[_c_f_ind]].set_index(['index'], inplace=True)

   adindel_df_list, ctindel_df_list, other_c = get_groupd_columns(indels_cout_df[ list(indels_cout_df.keys())[0] ].columns, sample_list)
   print("AD samples", adindel_df_list)
   print("Ctrl samples", ctindel_df_list)

   if not sample_id==None:
      stdev_data = {}
      for c_vcf in vcf_file_list:
         stdev_data[c_vcf.split('/')[-sample_id]] = read_clair_vcf_to_dataframe(c_vcf, orderedchr, unusedchr).groupby(('CHROM'))
   else: stdev_data = {}
 
   cd_thr = 0.5
   final_indel_res = {}
   for _ak in ["chr", "start", "end", "sample", 'indeltype', 'indellen', 'group']:
      final_indel_res[ _ak ] = []
   for _indel_type in indels_cout_df:
      outp_ft_format = outp_f_format+_indel_type
      indel_cout_df = indels_cout_df[_indel_type]
      indel_base_df = indels_base_df[_indel_type]
      if indel_cout_df.shape[0]<1 or indel_base_df.shape[0]<1:
         print("INDEL type", _indel_type, 'Zero', indel_cout_df.shape, indel_base_df.shape)
         continue
      print("INDEL type", _indel_type)

      indel_cout_dfadindel_NonZero = (indel_cout_df[adindel_df_list] != 0).sum(axis=1)
      indel_cout_dfctindel_NonZero = (indel_cout_df[ctindel_df_list] != 0).sum(axis=1)
      indel_cout_dfadindel_Percent = indel_cout_dfadindel_NonZero / len(adindel_df_list)
      indel_cout_dfctindel_Percent = indel_cout_dfctindel_NonZero / len(ctindel_df_list)
      group_filter = ((indel_cout_dfadindel_Percent >= min_sup_thr) | (indel_cout_dfctindel_Percent >= min_sup_thr))

      if t_op not in [1, '1', 0, '0']:
         indel_cout_dif = indel_cout_df
         indel_base_dif = indel_base_df
      else:
         dif_count_filter = compute_effect_size(indel_cout_df, adindel_df_list, ctindel_df_list);
         dif_bases_filter = compute_effect_size(indel_base_df, adindel_df_list, ctindel_df_list);

         indel_cout_dif = indel_differential(copy.deepcopy(indel_cout_df[ (group_filter & (dif_count_filter>= cd_thr)) ]), adindel_df_list+ctindel_df_list, sample_list, outp_ft_format+'_counts_dif.csv')
         indel_base_dif = indel_differential(copy.deepcopy(indel_base_df[ (group_filter & (dif_bases_filter>= cd_thr)) ]), adindel_df_list+ctindel_df_list, sample_list, outp_ft_format+'_base_dif.csv')
         #print(indel_base_dif)

         indel_cout_sum = [indel_cout_df[ ((indel_cout_dfadindel_Percent >= min_sup_thr) & (indel_cout_dfctindel_Percent <0.01)) | ((indel_cout_dfadindel_Percent <0.01) & (indel_cout_dfctindel_Percent >= min_sup_thr)) ][adindel_df_list].sum(), indel_cout_df[ ((indel_cout_dfadindel_Percent >= min_sup_thr) & (indel_cout_dfctindel_Percent <0.01)) | ((indel_cout_dfadindel_Percent <0.01) & (indel_cout_dfctindel_Percent >= min_sup_thr)) ][ctindel_df_list].sum()]
         indel_base_sum = [indel_base_df[ ((indel_cout_dfadindel_Percent >= min_sup_thr) & (indel_cout_dfctindel_Percent <0.01)) | ((indel_cout_dfadindel_Percent <0.01) & (indel_cout_dfctindel_Percent >= min_sup_thr)) ][adindel_df_list].sum(), indel_base_df[ ((indel_cout_dfadindel_Percent >= min_sup_thr) & (indel_cout_dfctindel_Percent <0.01)) | ((indel_cout_dfadindel_Percent <0.01) & (indel_cout_dfctindel_Percent >= min_sup_thr)) ][ctindel_df_list].sum()]
  
         indel_base_sum[0] = np.where(indel_base_sum[0] > 10000000, 10000000, indel_base_sum[0])
         indel_base_sum[1] = np.where(indel_base_sum[1] > 10000000, 10000000, indel_base_sum[1])

         fig, axes = plt.subplots(1, 2, figsize=(10, 5))
         axes[0].violinplot([indel_cout_sum[0].values.tolist(), indel_cout_sum[1].values.tolist()]); #, labels=["AD", "Ctrl"])
         axes[0].set_xticks([1, 2])  # Set x-axis positions
         axes[0].set_xticklabels(["AD", "Ctrl"])  # Set x-axis labels
         axes[0].set_title("INDEL Count Sum")
         axes[1].violinplot([indel_base_sum[0], indel_base_sum[1]]); #labels=["AD", "Ctrl"])
         axes[1].set_xticks([1, 2])
         axes[1].set_xticklabels(["AD", "Ctrl"])
         axes[1].set_title("INDEL Base Sum")
         # Show plots
         plt.tight_layout()
         plt.savefig(outp_ft_format+"_unique.png", dpi=600, bbox_inches='tight')
         plt.close('all')
  
      indel_cout_dfadindel_NonZero = (indel_cout_dif[adindel_df_list] != 0).sum(axis=1)
      indel_cout_dfctindel_NonZero = (indel_cout_dif[ctindel_df_list] != 0).sum(axis=1)
      indel_cout_dfadindel_Percent = indel_cout_dfadindel_NonZero / len(adindel_df_list)
      indel_cout_dfctindel_Percent = indel_cout_dfctindel_NonZero / len(ctindel_df_list)
      uniq_index =( ((indel_cout_dfadindel_Percent >= min_uniqsup_thr) & (indel_cout_dfctindel_Percent <0.01)) | ((indel_cout_dfadindel_Percent <0.01) & (indel_cout_dfctindel_Percent >= min_uniqsup_thr)) )
      res_cout_df = indel_cout_dif[ uniq_index | ( (indel_cout_dif['adj.P.Val']<=0.05) & ( (indel_cout_dfctindel_Percent<0.15) | (indel_cout_dfadindel_Percent<0.15) ) ) ]

      indel_base_dfadindel_NonZero = (indel_base_dif[adindel_df_list] != 0).sum(axis=1)
      indel_base_dfctindel_NonZero = (indel_base_dif[ctindel_df_list] != 0).sum(axis=1)
      indel_base_dfadindel_Percent = indel_base_dfadindel_NonZero / len(adindel_df_list)
      indel_base_dfctindel_Percent = indel_base_dfctindel_NonZero / len(ctindel_df_list)
      res_base_df = indel_base_dif[ ( indel_base_dif.index.isin( indel_cout_dif[uniq_index].index )  ) | ( (indel_base_dif['adj.P.Val']<=0.05)  & ( (indel_base_dfadindel_Percent<0.15) | (indel_base_dfctindel_Percent<0.15) ) ) ]

      merg_indexes = [ res_cout_df.index, res_base_df[ ~(res_base_df.index.isin( res_cout_df.index ) ) ].index ]
      merged_loci = {}
      for _t_index in merg_indexes:
         for c_ind_loc in _t_index:
            c_loc = c_ind_loc.split(':')
            c_loc[1] = c_loc[1].split('-')
            if c_loc[0] not in merged_loci: merged_loci[ c_loc[0] ] = []
            merged_loci[ c_loc[0] ].append( [int(c_loc[1][0])-extend_bp, int(c_loc[1][1])+extend_bp] )

      save_indel_keys = {}
      g_2 = [adindel_df_list, ctindel_df_list]
      for c_k in merged_loci: # for each chromoseom
         merged_loci[ c_k ] = sorted(merged_loci[ c_k ])
         for pos_inf in merged_loci[ c_k ]: # for each region;
            for _gi in range(len(g_2)): #for two groups;
               for _s_k in g_2[_gi]: # for each sample
                  chr_df = stdev_data[_s_k].get_group(c_k);
                  _in_region = chr_df[((chr_df['POS'] <= pos_inf[1]) & (chr_df['POS'] >= pos_inf[0])) | ((chr_df['INDELEND'] <= pos_inf[1]) & (chr_df['INDELEND'] >= pos_inf[0]))  ]
                  if _indel_type=='OTHER':
                     _in_region = _in_region[_in_region['INDELTYPE'].isin(Def_indel_types[2:])]
                  else:  _in_region = _in_region[_in_region['INDELTYPE']==_indel_type]
                  if not (_in_region.empty):
                     for _, t_region in _in_region.iterrows():
                        if (_s_k, c_k, t_region['POS'], t_region['INDELEND']) not in save_indel_keys:
                           final_indel_res["chr"].append( c_k )
                           final_indel_res["start"].append( t_region['POS'] )
                           final_indel_res["end"].append( t_region['INDELEND'] )
                           final_indel_res["sample"].append( _s_k )
                           final_indel_res["indeltype"].append( t_region['INDELTYPE'] )
                           final_indel_res["indellen"].append( t_region['INDELLEN' ] )
                           final_indel_res["group"].append("AD" if _gi==0 else "Ctrl")
                           save_indel_keys[ (_s_k, c_k, t_region['POS'], t_region['INDELEND']) ] = 0

   final_indel_df = pd.DataFrame(final_indel_res ).sort_values(by=['chr', 'start', 'end'])
   final_indel_df["start"] = final_indel_df["start"].astype(int)
   final_indel_df["end"] = final_indel_df["end"].astype(int)
   final_indel_df["chr"] = final_indel_df["chr"].astype(str)
   final_indel_df["sample"] = final_indel_df["sample"].astype(str)
   final_indel_df['indeltype'] = final_indel_df["indeltype"].astype(str)
   final_indel_df['indellen'] = final_indel_df["indellen"].astype(int)

   pl.from_pandas( final_indel_df ).write_csv( outp_ft_format+"_bffilter.csv" );

   rem_ind_list = []
   res_region = []
   for mchr, t_indels in final_indel_df.groupby(('chr')):
      prev_i = 0;
      start_i = -1;
      for _t_i in range(len(t_indels)):
          if t_indels.iloc[_t_i]['start'] - t_indels.iloc[prev_i]['start']<extend_bp:
             if start_i==-1: start_i = prev_i
             prev_i = _t_i
          else:
             group_size = t_indels.iloc[start_i:_t_i].groupby(("group"))['sample'].nunique()
             if 'AD' in group_size:
                t_ad_p = group_size['AD'] / len(adindel_df_list);
             else: t_ad_p = 0
             if 'Ctrl' in group_size:
                t_ct_p = group_size['Ctrl'] / len(ctindel_df_list);
             else: t_ct_p = 0
             if (t_ad_p<0.15 and t_ct_p>=min_sup_thr) or (t_ad_p>=min_sup_thr and t_ct_p<0.15):
                print(mchr, t_indels.iloc[start_i]['start'], t_indels.iloc[prev_i]['start'] );
                res_region.append([ mchr, t_indels.iloc[start_i:_t_i]['start'].median(), t_indels.iloc[start_i:_t_i]['end'].median(),t_indels.iloc[start_i:_t_i]['start'].min(), t_indels.iloc[start_i:_t_i]['start'].max(), t_indels.iloc[start_i:_t_i]['end'].min(), t_indels.iloc[start_i:_t_i]['end'].max(), t_ad_p, t_ct_p, '/'.join(f"{k}:{v}" for k, v in (t_indels.iloc[start_i:_t_i]['indeltype'].value_counts().to_dict().items()) ) ])
                print(t_indels.iloc[start_i:_t_i] );
             else:
                rem_ind_list.extend( t_indels.iloc[start_i:_t_i].index.tolist())
             prev_i = _t_i
             start_i = prev_i
      group_size = t_indels.iloc[start_i:_t_i].groupby(("group"))['sample'].nunique()
      if 'AD' in group_size:
         t_ad_p = group_size['AD'] / len(adindel_df_list);
      else: t_ad_p = 0
      if 'Ctrl' in group_size:
        t_ct_p = group_size['Ctrl'] / len(ctindel_df_list);
      else: t_ct_p = 0
      if (t_ad_p<0.15 and t_ct_p>=min_sup_thr) or (t_ad_p>=min_sup_thr and t_ct_p<0.15):
         print(mchr, t_indels.iloc[start_i]['start'], t_indels.iloc[prev_i]['start'] );
         res_region.append([ mchr, t_indels.iloc[start_i:_t_i]['start'].median(), t_indels.iloc[start_i:_t_i]['end'].median(),t_indels.iloc[start_i:_t_i]['start'].min(), t_indels.iloc[start_i:_t_i]['start'].max(), t_indels.iloc[start_i:_t_i]['end'].min(), t_indels.iloc[start_i:_t_i]['end'].max(), t_ad_p, t_ct_p, '/'.join(f"{k}:{v}" for k, v in (t_indels.iloc[start_i:_t_i]['indeltype'].value_counts().to_dict().items()) ) ])
         print(t_indels.iloc[start_i:_t_i] );
      else:
         rem_ind_list.extend( t_indels.iloc[start_i:_t_i].index.tolist())

   outp_ft_format = outp_f_format+'_res'
   pl.from_pandas( pd.DataFrame(res_region, columns=["chr", "start", "end", "start_min", "start_max", "end_min", "end_max", "AD_freq", 'Ctrl_freq', 'indelinfo']).sort_values(by=['chr', 'start'])  ).write_csv( outp_ft_format+"_region.csv" );
   print(final_indel_df.shape)
   final_indel_df = final_indel_df.drop(rem_ind_list)
   print(final_indel_df.shape)
   if final_indel_df.shape[0]>0:
      pl.from_pandas( final_indel_df ).write_csv( outp_ft_format+".csv" );
      generate_sv_count_violinplots(final_indel_df, adindel_df_list, ctindel_df_list, outp_ft_format, col_t='indeltype')
      generate_svlen_violinplots(final_indel_df, adindel_df_list, ctindel_df_list, outp_ft_format, col_t='indetype', col_l='indellen')

