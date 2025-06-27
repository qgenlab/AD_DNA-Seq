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
pd.set_option('future.no_silent_downcasting', True)

#

min_value = 10**(-3)

default_plot_ext = ['png', 'jpeg', 'jpg', 'tiff', 'pdf']

Def_sv_types = ['DEL', 'INS', 'DUP', 'INV', 'OTHER']
Def_indel_types = ['del', 'ins', 'other']

def read_sample_list(sample_file='sample_category.txt'):
   df = pd.read_csv(sample_file, sep='\t');
   if 'Cancer' in df.columns:
      ad_list = (df[(df['AD']==1)&(df['Cancer']==0)])['Sample'].to_list()
      print( ad_list )
      ctrlist = (df[(df['AD']==0)&(df['Cancer']==0)])['Sample'].to_list()
      print( ctrlist )
   else:
      ad_list = (df[(df['AD']==1)])['Sample'].to_list()
      print( ad_list )
      ctrlist = (df[(df['AD']==0)])['Sample'].to_list()
      print( ctrlist )
   #print(df.head())
   df = df.set_index(['Sample'])
   #print(df.head())
   return (df, ad_list, ctrlist);

def filter_input(input_prefx, readThresh='7,10', min_value = 10**(-3), filter_samples_ratio=0.6, sample_file='sample_category.txt'):
   readThresh_sp = readThresh.split(',')
   readThresh = int(readThresh_sp[0])
   df_methp = pl.read_csv(input_prefx+'_meth_perc_20_cov'+str(readThresh)+'.csv', separator=',', infer_schema_length=100000).to_pandas();
   df_methp.set_index(['Chr', 'Pos', 'Strand'], inplace=True)
   #
   df_covag = pl.read_csv(input_prefx+'_cov_20_cov'+str(readThresh)+'.csv', separator=',', infer_schema_length=100000).to_pandas();
   df_covag.set_index(['Chr', 'Pos', 'Strand'], inplace=True)
   #
   df_meth_dmp = pl.read_csv(input_prefx+'.csv', separator=',', infer_schema_length=100000).to_pandas();
   df_meth_dmp.set_index(['Chr', 'Pos', 'Strand'], inplace=True)
   #
   sample_list = read_sample_list(sample_file)
   #
   nc_all = df_methp.columns;
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
         others.append( c)
   #
   others = []
   for c in df_meth_dmp.columns:
      if c not in nc_all:
          others.append(c)

   readThresh = int(readThresh_sp[1])
   df_covag.fillna(0)
   minSamp_num = [ len(ad_df_list[0])*filter_samples_ratio, len(ct_df_list[0])*filter_samples_ratio ]
   for _i in range(len(minSamp_num)):
      if minSamp_num[_i]<2:
         minSamp_num[_i]=2;
   #
   filtered_meth_df = copy.deepcopy(df_methp[( ( (df_covag>=readThresh).sum(axis=1)>=(minSamp_num[0]+minSamp_num[1]) ) ) ] )
   filtered_covg_df = copy.deepcopy(df_covag[( ( (df_covag>=readThresh).sum(axis=1)>=(minSamp_num[0]+minSamp_num[1]) ) ) ] )
   print('Orignal:', df_methp.shape, df_covag.shape, 'filtered', filtered_meth_df.shape, filtered_covg_df.shape, minSamp_num, len(ad_df_list[0]), len(ct_df_list[0]))
   #
   ctrl_cov_cond = (filtered_covg_df.loc[:, ct_df_list[1]]>=readThresh).sum(axis=1)
   case_cov_cond = (filtered_covg_df.loc[:, ad_df_list[1]]>=readThresh).sum(axis=1)
   group_filer = ( ( ctrl_cov_cond>=minSamp_num[1]) & ( case_cov_cond>=minSamp_num[0]) )
   row_ctrl_mean = filtered_meth_df.loc[:, ct_df_list[0]].mean(axis=1, skipna=True).fillna(min_value)
   row_ad_mean = filtered_meth_df.loc[:, ad_df_list[0]].mean(axis=1, skipna=True).fillna(min_value)
   small_mean = 0.2
   #
   if minSamp_num[0]+minSamp_num[1] < len(ad_df_list[0]) + len(ct_df_list[0]):
      group_filer = group_filer | ( ( ( case_cov_cond>=minSamp_num[0] ) & ( ctrl_cov_cond<minSamp_num[1]) & ( row_ctrl_mean<small_mean ) ) |  ( ( ctrl_cov_cond>=minSamp_num[1]) & ( case_cov_cond<minSamp_num[0] ) & ( row_ad_mean<small_mean ) ) )
   #
   filter2_meth_df = copy.deepcopy(filtered_meth_df[group_filer])
   #
   row_ctrl_mean = filter2_meth_df.loc[:, ct_df_list[0]].mean(axis=1, skipna=True).fillna(min_value)
   row_ad_mean = filter2_meth_df.loc[:, ad_df_list[0]].mean(axis=1, skipna=True).fillna(min_value)
   #
   filter2_meth_df['Meth_dif'] = row_ad_mean - row_ctrl_mean
   print(filter2_meth_df.shape)
   #              0                        1                                   2           3           4         5
   return (filter2_meth_df, copy.deepcopy(filtered_covg_df[group_filer]), df_meth_dmp, ad_df_list, ct_df_list, others)


def parse_annotation(annotation):
    """
    Parse the annotation field to categorize into gene intron, gene exon, IG, and ENCODE types.
    """
    categories = set()
    encode_types = set()
    repeat = set()
    cpg_epic = set()

    enhance_list = dict();
    gene_list = dict();

    nb_suf = 'nb100'
    #EPIC not consider
    annotations = annotation.split('/')
    for entry in annotations:
        if ':' in entry:
            key, value = entry.split(':', 1)

            if key == 'IG':
                categories.add('IG')
            elif key == 'ENCODE' or key == 'CCRE':
                encode_type = value.split('@')[1] if '@' in value else value
                if encode_type[-len(nb_suf):]==nb_suf:
                   encode_type = encode_type.split('_')[0]+'_nb'
                encode_types.add(encode_type)
                en_id = value.split('@')[0] if '@' in value else value
                if en_id not in gene_list:
                   gene_list[ en_id ] = set()
                gene_list[ en_id ].add(encode_type)
            elif key == 'REPEAT':
                repeat.add( value);
            elif key=='EPIC' or key=='EPICv2' or key=='HM450':
                cpg_epic.add( value)
            else:
                if 'intron' in value.lower():
                    categories.add('Gene_Intron')
                    if 'intron' not in gene_list: gene_list[ 'intron'] = set()
                    gene_list[ 'intron'].add(key);
                elif 'exon' in value.lower():
                    categories.add('Gene_Exon')
                    if 'exon' not in gene_list: gene_list[ 'exon'] = set()
                    gene_list[ 'exon'].add(key);

    # Ensure IG is not counted if others exist
    if ('Gene_Intron' in categories or 'Gene_Exon' in categories or encode_types):
        categories.discard('IG')
    if 'enhD' in encode_types or 'enhP' in encode_types or 'prom' in encode_types or 'K4m3' in encode_types:
       categories.discard('Gene_Intron')
    if 'enhD' in encode_types:
       encode_types.discard('enhD');
       encode_types.add('enh')
    if 'enhP' in encode_types:
       encode_types.discard('enhP');
       encode_types.add('enh')
    merg_enh = []
    for _et in encode_types:
       if '_' in _et and _et.split('_')[0] in ['enhD', 'enhP']:
          merg_enh.append( _et );
    for _et in merg_enh:
       encode_types.discard( _et );
       encode_types.add('enh'+'_'+_et.split('_')[1] )

    return categories, encode_types, repeat, cpg_epic, gene_list
    


def convert_dict_to_str(m_dict):
   str_list = []
   for dk in sorted(list(m_dict.keys())):
      str_list.append("{}:{}".format( str(dk), str(m_dict[dk]) ))
   return '/'.join(str_list)


def get_chr_human(contain_xy=False, contain_m=False):
   orderedchr = []
   orderedchrnum = []
   for ic in range(1, 23):
      orderedchr.append('chr'+str(ic));
      orderedchrnum.append(str(ic));

   if contain_xy:
      orderedchr.append('chrX'); orderedchrnum.append('X')
      orderedchr.append('chrY'); orderedchrnum.append('Y')
   if contain_m:
      orderedchr.append('chrM'); orderedchrnum.append('M')

   return (orderedchr, orderedchrnum)

def read_vcf_to_dataframe_dynamic_columns(filename, orderedchr, unusedchr, target_region=None):
    data= []
    columns = None  # Initialize columns as None

    with open(filename, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#') and not line.startswith('##'):  # Header line
                columns = line.strip()[1:].split('\t')  # Extract column names (skip the initial '#')
                break  # Stop after reading the header
        if columns is None:
            raise ValueError("No column header line found in the VCF file.")

        vcf_file.seek(0) #reset the file pointer to the beginning.

        for line in vcf_file:
            if not line.startswith('#'):  # Skip header lines
                row = line.strip().split('\t')
                record = dict(zip(columns, row))  # Create dictionary with column names
                if (not unusedchr==None) and record['CHROM'] in unusedchr: continue;
                if record['CHROM'] not in orderedchr: continue;
                if not record['FILTER']=='PASS': continue;
                record['POS'] = int(record['POS'])-1
                record['QUAL'] = int(record['QUAL'])

                support = None;
                sv_pos_std = None;
                sv_end = None;
                sv_length = None;
                SVTYPE = None;
                SVLEN = None;
                for _isp in record['INFO'].split(';'):
                   t_type = _isp.split('=')
                   if t_type[0]=='STDEV_POS': sv_pos_std = float(t_type[1])
                   elif t_type[0]=='SUPPORT': support = int(t_type[1])
                   elif t_type[0]=='END': sv_end  = int(t_type[1])-1
                   elif t_type[0]=='SVTYPE': SVTYPE = t_type[1]
                   elif t_type[0]=='SVLEN': SVLEN = int(t_type[1])
                if support<5: continue;
                if sv_end==None: sv_end = int(record['POS'])

                if not target_region==None:
                    if not (record['CHROM'] == target_region[0] and record['POS']>=target_region[1] and record['POS']<=target_region[2]):
                       continue;
                record['STDEV_POS'] = sv_pos_std
                record['SVEND'] = sv_end
                record['SVTYPE'] = SVTYPE
                record['SVLEN'] = SVLEN

                data.append({a_k:record[a_k] for a_k in ['CHROM', 'POS', 'QUAL','STDEV_POS', 'SVEND', 'SVTYPE', 'SVLEN']})

    df = pd.DataFrame(data)
    return df

def read_clair_vcf_to_dataframe(c_vcf, orderedchr, unusedchr):
    data= []
    columns = None  # Initialize columns as None

    with open(c_vcf, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#') and not line.startswith('##'):  # Header line
                columns = line.strip()[1:].split('\t')  # Extract column names (skip the initial '#')
                break  # Stop after reading the header
        if columns is None:
            raise ValueError("No column header line found in the VCF file.")

        vcf_file.seek(0) #reset the file pointer to the beginning.

        for line in vcf_file:
            if not line.startswith('#'):  # Skip header lines
                row = line.strip().split('\t')
                record = dict(zip(columns, row))  # Create dictionary with column names
                if (not unusedchr==None) and record['CHROM'] in unusedchr: continue;
                if record['CHROM'] not in orderedchr: continue;
                if not record['FILTER']=='PASS': continue;
                record['POS'] = int(record['POS'])-1
                record['QUAL'] = float(record['QUAL'])

                ref_seq, alt_seq = record['REF'], record['ALT']
                alt_seq = alt_seq.split(',')
                if ',' in ref_seq:
                    print('Warning mutiple sequences in Ref', ref_seq);
                for _alt_s in alt_seq:
                    if len(_alt_s)<len(ref_seq):
                       indeltype = 'Del'
                       indellen = len(ref_seq)-len(_alt_s)
                       indelend = record['POS'] + indellen
                    elif len(_alt_s)>len(ref_seq):
                       indeltype = 'Ins'
                       indellen = len(_alt_s)-len(ref_seq)
                       indelend = record['POS']
                    else:
                       if ',' not in record['ALT']:
                          print("Alt is equal REF", indel)
                       continue;
                    add_r = {}
                    for a_k in ['CHROM', 'POS', 'QUAL']:
                       add_r[ a_k ] = record[a_k]
                    add_r['INDELTYPE'] = indeltype
                    add_r['INDELLEN'] = indellen
                    add_r['INDELEND'] = indelend
                    data.append( add_r )

    df = pd.DataFrame(data)
    return df



def get_groupd_columns(nc_all, sample_list):
   adsv_df_list = []
   ctsv_df_list = []
   other_c = []
   for c in nc_all:
         c_sp = c.split('_')
         if '_'.join(c_sp[:(2 if len(c_sp)>=2 else 1)]) in sample_list[1]:
            adsv_df_list.append(c)
         elif '_'.join(c_sp[:(2 if len(c_sp)>=2 else 1)]) in sample_list[2]:
            ctsv_df_list.append(c)
         else: other_c.append(c); #print("Not support", c)
   #print("AD samples", adsv_df_list)
   #print("Ctrl samples", ctsv_df_list)
   return adsv_df_list, ctsv_df_list, other_c


def generate_sv_count_violinplots(final_sv_df, adsv_df_list, ctsv_df_list, outp_ft_format, col_t='svtype'):
    svtypes = final_sv_df[col_t].unique()
    fig, axes = plt.subplots(len(svtypes), 1, figsize=(8, 5 * len(svtypes)))

    for i, svtype in enumerate(svtypes):
        ad_counts = [len(final_sv_df[(final_sv_df['sample'] == sample) & (final_sv_df[col_t] == svtype)]) for sample in adsv_df_list]
        ct_counts = [len(final_sv_df[(final_sv_df['sample'] == sample) & (final_sv_df[col_t] == svtype)]) for sample in ctsv_df_list]

        sns.violinplot(data=[ad_counts, ct_counts], ax=axes[i])
        axes[i].set_title(f"SV-type: ({svtype})", fontsize=18)
        axes[i].set_xticks([0, 1])
        axes[i].set_xticklabels(['AD', 'Control'], fontsize=16)
        axes[i].set_ylabel(f"{svtype} Count", fontsize=16)
        axes[i].tick_params(axis='both', labelsize=16)

    plt.tight_layout()
    plt.savefig(outp_ft_format+"_num.png", dpi=600, bbox_inches='tight')
    plt.close('all')

def generate_svlen_violinplots(final_sv_df, adsv_df_list, ctsv_df_list, outp_ft_format, min_len=-2000, max_len=1000, col_t='svtype', col_l='svlen'):
    svtypes = final_sv_df[col_t].unique()
    fig, axes = plt.subplots(len(svtypes), 1, figsize=(8, 5 * len(svtypes)))
    for i, svtype in enumerate(svtypes):
        ad_svlens = pd.concat([
            final_sv_df[(final_sv_df['sample'] == sample) & (final_sv_df[col_t] == svtype)][col_l]
            for sample in adsv_df_list
        ])
        ct_svlens = pd.concat([
            final_sv_df[(final_sv_df['sample'] == sample) & (final_sv_df[col_t] == svtype)][col_l]
            for sample in ctsv_df_list
        ])

        # Construct labeled DataFrame
        plot_df = pd.concat([
            pd.DataFrame({"svlen": ad_svlens, "group": "AD"}),
            pd.DataFrame({"svlen": ct_svlens, "group": "Control"})
        ], ignore_index=True)
        plot_df["svlen"] = plot_df["svlen"].clip(lower=min_len, upper=max_len)

        sns.violinplot(x="group", y="svlen", data=plot_df, ax=axes[i]); #, palette={"AD": "blue", "Control": "red"})
        axes[i].set_title(f"SV Length ({svtype})", fontsize=18)
        axes[i].set_ylabel(f"{svtype} Length", fontsize=16)
        axes[i].tick_params(axis='both', labelsize=16)
    plt.tight_layout()
    #plt.show()
    plt.savefig(outp_ft_format+"_len.png", dpi=600, bbox_inches='tight')
    plt.close('all')

def read_gtf_gene_exon(gtf_file, header_rown=5): # ./reference_genome/human/hg38/gencode/gencode.v42.chr_patch_hapl_scaff.annotation.gtf
   gtf_df = pd.read_csv(gtf_file, sep='\t', header=None, skiprows=header_rown);
   gene_dict = []
   exon_dict = []
   t_exon = [];
   for _, t_row in gtf_df.iterrows():
      if t_row[2]=='gene':
         if len(t_exon)>0:
            t_exon = sorted(t_exon)
            m_exon = [t_exon[0] ];
            for _n_ex in t_exon[1:]:
               if _n_ex[1]<=m_exon[-1][2]:
                  m_exon[-1][1]=min( m_exon[-1][1], _n_ex[1] )
                  m_exon[-1][2]=max( m_exon[-1][2], _n_ex[2] )
               else:
                  m_exon.append( _n_ex )
            for _n_ex in m_exon:
               exon_dict.append( _n_ex )
         t_exon = [];
         gene_name = ""
         for fieldi in t_row[8].split(';'):
            f_lsp = fieldi.split()
            if f_lsp[0]=='gene_name':
               gene_name = f_lsp[1][1:-1]
               break;
         if gene_name=="":
            print("Warning!!! Cannot find gene name for --"+t_row)
         #                     0             1              2           
         gene_dict.append( [t_row[0], int(t_row[3])-1, int(t_row[4]), gene_name] )
      elif t_row[2]=='transcript':
         pass;
      else: #                                0            1             2
         t_exon.append( [t_row[2], int(t_row[3])-1, int(t_row[4]), gene_name] )
    
   if len(t_exon)>0:
      t_exon = sorted(t_exon)
      m_exon = [t_exon[0] ];
      for _n_ex in t_exon[1:]:
         if _n_ex[1]<=m_exon[-1][2]:
            m_exon[-1][1]=min( m_exon[-1][1], _n_ex[1] )
            m_exon[-1][2]=max( m_exon[-1][2], _n_ex[2] )
         else:
            m_exon.append( _n_ex )
      for _n_ex in m_exon:
         exon_dict.append( _n_ex )

   return pd.DataFrame(gene_dict, columns=['chr', 'start', 'end', 'gene_name']), pd.DataFrame(exon_dict, columns=['chr', 'start', 'end', 'gene_name'])
 
