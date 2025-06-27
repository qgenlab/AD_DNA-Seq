
import os,sys

import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.cm as cm
import numpy as np

from com_fun import *


def calculate_region_distance(start1, end1, start2, end2):
    if end1 >= start2 and end2 >= start1:
        return 0  # Regions overlap

    # Calculate distance if no overlap
    if end1 < start2:
        return start2 - end1  # Region 1 is before Region 2
    else:
        return start1 - end2  # Region 2 is before Region 1


def process_regions(region_file, gtf_file, gene_bed_file, ccre_file):
    # Read region file
    #regions_df = pd.read_csv(region_file, sep=',', header=0, names=['chr', 'start', 'end'], usecols=[0,1,2])
    regions_df = pd.read_csv(region_file, sep=',', header=0)
    if ('chromosome' in regions_df.columns) and ('chr' not in regions_df.columns):
       regions_df = regions_df.rename(columns={'chromosome':'chr'})
    
    gene_bed = pd.read_csv(gene_bed_file, sep='\t', header=None, names=['chr', 'start', 'end', 'gene', 'strand'], usecols=[0,1,2,3,5])
    ccre_bed = pd.read_csv(ccre_file, sep='\t', header=None, names=['chr', 'start', 'end', 'ccre_type', 'ccre'], usecols=[0,1,2,12, 13])

    gene_df, exon_df = read_gtf_gene_exon(gtf_file)

    exon_genes = defaultdict(int)
    intro_genes  = defaultdict(int)
    prom_genes = defaultdict(int)
    enhP_genes = defaultdict(int)
    enhD_genes = defaultdict(int)

    for _, region in regions_df.iterrows():
        matched_ccre = ccre_bed[
            (ccre_bed['chr'] == region['chr']) & 
            (ccre_bed['start'] <= region['end']) & 
            (ccre_bed['end'] >= region['start'])
        ]

        for _, annotation in matched_ccre.iterrows():
            _de = annotation['ccre_type']
            if _de=='enhD': t_dict = enhD_genes
            elif _de=='enhP': t_dict = enhP_genes
            elif _de=='prom': t_dict = prom_genes
            else:
               print("does not support", _de);
               continue;
            if _de in ['enhD', 'enhP', 'prom']:
                t_dict[ annotation['ccre'] ] += 1;

    for _, region in regions_df.iterrows():
        matched_g = gene_df[
            (gene_df['chr'] == region['chr']) &
            (gene_df['start'] <= region['end']) &
            (gene_df['end'] >= region['start'])
        ]
        if matched_g.empty: continue;
        matched_exons = exon_df[
            (exon_df['chr'] == region['chr']) &
            (exon_df['start'] <= region['end']) &
            (exon_df['end'] >= region['start'])
        ]
        if matched_exons.empty:
           for _, annotation in matched_g.iterrows():
              _ge = annotation['gene_name']
              intro_genes[ _ge ] += 1
        for _, annotation in matched_exons.iterrows():
           _ge = annotation['gene_name']
           exon_genes[ _ge ] += 1

    promF_genes = defaultdict(int)
    enhPF_genes = defaultdict(int)
    enhDF_genes = defaultdict(int)
   
    enhd_thr = 500000
    enhp_thr = 50000
    prom_thr = 2000
    
    thresholds = [enhd_thr, enhp_thr, prom_thr]
    gene_list = [enhD_genes, enhP_genes, prom_genes]
    res_list = [enhDF_genes, enhPF_genes, promF_genes]
    for _gi in range(len(thresholds)):
       t_thr = thresholds[ _gi ]
       t_ccre = gene_list[ _gi ]
       t_genes = res_list[ _gi ]
       for t_e in t_ccre:
          ccre_region = ccre_bed[ ccre_bed['ccre']==t_e ] 
          if not ccre_region.shape[0]==1:
             print("wrong id", t_e, ccre_region.shape, ccre_region)
             continue;
          res_gene_def = gene_bed[ (gene_bed['chr']==ccre_region['chr'].iloc[0]) & ( gene_bed['start']-t_thr <=ccre_region['end'].iloc[0] ) & ( gene_bed['end']+t_thr>=ccre_region['start'].iloc[0]  )   ]
          if _gi==0: this_t = "enhD"
          elif _gi==1: this_t = "enhP"
          else: this_t = "Prom"
          outpt_list = [t_e+'\t'+this_t]
          if _gi in [0, 1]:
             for _, g_row in res_gene_def.iterrows():
                t_genes[ g_row['gene'] ] += 1;
                outpt_list.append("\n\tANNOg:"+g_row['gene']+':'+str(calculate_region_distance(ccre_region['start'].iloc[0], ccre_region['end'].iloc[0], g_row['start'],  g_row['end']) ))
             print("ANNOe:"+''.join( outpt_list ) );
          else:
             for _, g_row in res_gene_def.iterrows():
                if g_row['strand']=='+':
                   if g_row['start']>ccre_region['start'].iloc[0]: 
                      t_genes[ g_row['gene'] ] += 1;
                      outpt_list.append("\n\tANNOg:"+g_row['gene']+':'+str( g_row['start']-ccre_region['start'].iloc[0] ) )
                else:
                   if g_row['end'] < ccre_region['start'].iloc[0]: 
                      t_genes[ g_row['gene'] ] += 1;
                      outpt_list.append("\n\tANNOg:"+g_row['gene']+':'+str( g_row['end'] - ccre_region['start'].iloc[0] ) )
             print("ANNOe:"+''.join( outpt_list ) );

    
    return exon_genes, intro_genes, promF_genes, enhPF_genes, enhDF_genes




if __name__=='__main__':
   # Input files (modify paths if needed)
   region_file = sys.argv[1]
   gtf_file = sys.argv[2]

   # Process data
   exon_genes, intro_genes, prom_genes, enhP_genes, enhD_genes = process_regions(region_file, gtf_file, sys.argv[3], sys.argv[4])

   print('exon', convert_dict_to_str(exon_genes ))
   print('intron', convert_dict_to_str(intro_genes ))
   ccre_genes = defaultdict(int)
   for c_g_list in [prom_genes, enhP_genes, enhD_genes]:
      for _c_k in c_g_list:
         ccre_genes[ _c_k] += c_g_list[ _c_k ]
   print('ccre', convert_dict_to_str(ccre_genes ))

   final_list = set();
   rem_list = defaultdict(int)
   for c_g_list in [prom_genes, enhP_genes, enhD_genes, exon_genes, intro_genes]:
      for _c_k in c_g_list:
         if len(_c_k)>4 and _c_k[:4] in ['ENSG', 'LINC']:
            rem_list[ _c_k] += 1;
         elif len(_c_k.split('-AS'))>1:
            rem_list[ _c_k] += 1;
            final_list.add( _c_k.split('-AS')[0] );
         else:
            final_list.add( _c_k )

   #print( rem_list )
   #print( final_list )
   with open(region_file[:-4]+'_noncoding.txt', 'w') as file:
    for item in rem_list:
        file.write(str(item) + '\n')  # Convert to string and add newline

   with open(region_file[:-4]+'_coding.txt', 'w') as file:
    for item in final_list:
        file.write(str(item) + '\n')  # Convert to string and add newline

   for _fg in final_list:
     pass; #print(_fg)


