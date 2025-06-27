
import os,sys

import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.cm as cm
import numpy as np

from com_fun import *

def process_regions(region_file, ccre_file, gtf_file, rep_file):
    """
    Reads region file and annotation file, finds matches, and counts occurrences.
    """
    # Read region file
    regions_df = pd.read_csv(region_file, sep=',', header=0, names=['chr', 'start', 'end'], usecols=[0,1,2])

    ccre_bed = pd.read_csv(ccre_file, sep='\t', header=None, names=['chr', 'start', 'end', 'ccre_type', 'ccre'], usecols=[0,1,2,12, 13])
    gene_df, exon_df = read_gtf_gene_exon(gtf_file)

    default_col = [5,6,7,11]
    default_new = ['chr', 'start', 'end', 'labels']
    repeat_regions_df = pl.read_csv(rep_file, separator='\t', has_header=False, columns=default_col, new_columns=default_new, infer_schema_length=10000).to_pandas();
    if repeat_regions_df[default_new[-1]].str.contains('#').any() or repeat_regions_df[default_new[-1]].str.contains('/').any():
       repeat_regions_df['labels'] = repeat_regions_df[default_new[-1]].apply(split_and_get_second)

    total_counts_1 = defaultdict(int)
    total_counts_2 = defaultdict(int)
    total_counts_3 = defaultdict(int)

    for _, region in regions_df.iterrows():
        matched_rep = repeat_regions_df[
            (repeat_regions_df['chr'] == region['chr']) &
            (repeat_regions_df['start'] <= region['end']) &
            (repeat_regions_df['end'] >= region['start'])
        ]
        reg_size = region['end'] - region['start']
        max_rep_type = [0, '']
        for _, annotation in matched_rep.iterrows():
            t_size = min(region['end'], annotation['end']) - max(region['start'], annotation['start'])
            if t_size >= reg_size - t_size:
               if t_size>max_rep_type[0]:
                  max_rep_type[0] = t_size
                  max_rep_type[1] = annotation['labels']
        if max_rep_type[0] == 0:
           total_counts_3['No_repeat'] += 1
        else:
           total_counts_3[max_rep_type[1]] += 1

    for _, region in regions_df.iterrows():
        matched_ccre = ccre_bed[
            (ccre_bed['chr'] == region['chr']) &
            (ccre_bed['start'] <= region['end']) &
            (ccre_bed['end'] >= region['start'])
        ]
        ovlp_bp_dict = defaultdict(int)
        has_ccre = False;
        for _, annotation in matched_ccre.iterrows():
            t_size = min(region['end'], annotation['end']) - max(region['start'], annotation['start'])
            ovlp_bp_dict[ annotation['ccre_type'] ] += t_size
            ovlp_bp_dict[ 'Total_ccre' ] += t_size
            has_ccre = True;


        matched_g = gene_df[
            (gene_df['chr'] == region['chr']) &
            (gene_df['start'] <= region['end']) &
            (gene_df['end'] >= region['start'])
        ]
        if matched_g.empty:
            ovlp_bp_dict['Total_IG'] += region['end'] - region['start']
        else:
           matched_exons = exon_df[
               (exon_df['chr'] == region['chr']) &
               (exon_df['start'] <= region['end']) &
               (exon_df['end'] >= region['start'])
           ]
           if matched_exons.empty:
              for _, annotation in matched_g.iterrows():
                 ovlp_bp_dict['Intron'] += min(region['end'], annotation['end']) - max(region['start'], annotation['start'])
           else:
              for _, annotation in matched_exons.iterrows():
                 ovlp_bp_dict['Exon'] += min(region['end'], annotation['end']) - max(region['start'], annotation['start'])

        if has_ccre:
           for ent in ['enhD', 'enhP', 'prom', 'K4m3', 'enh', 'CTCF']:
               if ent in ovlp_bp_dict:
                  if ent in  ['enhD', 'enhP', 'enh']:
                     total_counts_1[ 'Enh' ] += 1;
                     total_counts_2[ 'Enh' ] += 1;
                  #elif ent in  ['prom']:
                  #   total_counts_1[ ent ] += 1;
                  #   total_counts_2[ ent ] += 1;
                  else:
                     total_counts_1[ 'Other cCRE' ] += 1;
                     total_counts_2[ 'Other cCRE' ] += 1;
           if 'Exon' in ovlp_bp_dict:
               total_counts_1[ 'Exon' ] += 1
               total_counts_2['Gene_Body'] += 1
        else:
           if 'Exon' in ovlp_bp_dict:
              total_counts_1[ 'Exon' ] += 1
              total_counts_2['Gene_Body'] += 1
           elif 'Intron' in ovlp_bp_dict:
              total_counts_1[ 'Intron' ] += 1
              total_counts_2['Gene_Body'] += 1
           else:
              total_counts_1['IG'] += 1
              total_counts_2['IG'] += 1


    return total_counts_1, total_counts_2, total_counts_3

# Function to format autopct to show both raw count and percentage
def autopct_format(pct, all_values):
    absolute = int(round(pct/100. * sum(all_values)))  # Convert percentage to raw count
    return f"{absolute} ({pct:.1f}%)"  # Format as 'count (percent%)'

def plot_pie_chart(data, title, fig_name):
    """
    Plot and save a pie chart.
    """
    labels = list(data.keys())
    sizes = list(data.values())
    labels = sorted(list(data.keys()) )
    sizes = []
    for _lb_ in labels:
       sizes.append( data[_lb_] )

    #cmap = cm.get_cmap('viridis', len(sizes)) # Get the viridis colormap with the number of slices
    #colors = [cmap(i) for i in range(len(sizes))]  # Get colors from the colormap
    colors = plt.get_cmap('Blues')(np.linspace(0.2, 0.7, len(sizes)))

    plt.figure(figsize=(7, 7))
    #plt.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=140, colors=colors)
    plt.pie(sizes, labels=labels, autopct='%1.0f%%', startangle=140, colors=colors, textprops={'fontsize': 16})
    #plt.pie(sizes, labels=labels, autopct=lambda pct: autopct_format(pct, sizes), startangle=140, colors=colors)
    plt.title(title)
    #plt.show()
    plt.savefig( fig_name, dpi=600, bbox_inches='tight')
    plt.close('all')



if __name__=='__main__':
   # Input files (modify paths if needed)
   region_file = sys.argv[1]
   # Process data
   total_counts_1, total_counts_2, total_counts_3 = process_regions(region_file, sys.argv[2], sys.argv[3], sys.argv[4])

   # Plot pie charts
   plot_pie_chart(total_counts_1, "Fraction of Occurrences: Gene Intron, Gene Exon, IG, and ENCODE Types", sys.argv[1][:-4]+'_SVinex2.png')
   plot_pie_chart(total_counts_2, "Fraction of Occurrences: Gene Body, IG, and ENCODE Types", sys.argv[1][:-4]+'_SVbody2.png')
   plot_pie_chart(total_counts_3, "", sys.argv[1][:-4]+'_SVrepeat2.png')



