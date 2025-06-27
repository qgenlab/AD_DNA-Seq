
import os,sys

import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import matplotlib.cm as cm
import numpy as np

from com_fun import *

def process_regions(region_file, annotation_file):
    """
    Reads region file and annotation file, finds matches, and counts occurrences.
    """
    # Read region file
    regions_df = pd.read_csv(region_file, sep=',', header=0, names=['chr', 'start', 'end'], usecols=[0,1,2])

    # Read annotation file
    annotation_df = pd.read_csv(annotation_file, sep='\t', header=None, names=['chr', 'start', 'end', 'id', 'annotation'])

    total_counts_1 = defaultdict(int)
    total_counts_2 = defaultdict(int)
    total_counts_3 = defaultdict(int)
    total_counts_4 = defaultdict(int)
    total_counts_p4 = defaultdict(int)

    # Iterate through regions and match annotations
    for _, region in regions_df.iterrows():
        matched_annotations = annotation_df[
            (annotation_df['chr'] == region['chr']) & 
            (annotation_df['start'] <= region['end']) & 
            (annotation_df['end'] >= region['start'])
        ]
        this_counts_1 = defaultdict(int)
        this_counts_2 = defaultdict(int)
        en_list = defaultdict(int)
        this_counts_3 = defaultdict(int)
        this_counts_4 = defaultdict(int)
        gend_counts = defaultdict(int)

        for _, annotation in matched_annotations.iterrows():
            categories, encode_types, repeat, cpg_epic, gene_list = parse_annotation(annotation['annotation'])

            for _ge in gene_list:
               for _de in gene_list[ _ge ]:
                  if (_ge + ':'+_de) not in gend_counts: gend_counts[  _ge + ':'+_de ] = 1;
                  else: gend_counts[  _ge + ':'+_de ] += 1
            # Count each category only once per region
            for cat in categories:
                this_counts_1[cat] += 1
                this_counts_2['Gene_Body' if 'Gene_' in cat else cat] += 1
            
            if len(repeat)==0: this_counts_3['No_repeat'] += 1;
            else:
                for rep_c in repeat:
                    this_counts_3[rep_c] += 1

            if len(cpg_epic)==0: 
                #total_counts_4['non-EPIC'] += 1;
                this_counts_4['non-EPIC'] += 1;
            else: 
                #total_counts_4['EPIC'] += 1;
                this_counts_4['EPIC'] += 1;

            for encode_type in encode_types:
                this_counts_1[encode_type] += 1
                this_counts_2[encode_type] += 1
                en_list[encode_type] += 1

        de_gede = []
        for _gede in gend_counts:
           if _gede[-3:]=='_nb' and _gede[:-3] in gend_counts:
              gend_counts[ _gede[:-3] ] += gend_counts[ _gede ]
              de_gede.append( _gede)
        for _gede in de_gede:
           del gend_counts[ _gede ]
        if len(this_counts_1)>0:
           print(region['chr']+':'+str(region['start'])+'-'+str(region['end'])+'/'+str(region['end']-region['start']), convert_dict_to_str(this_counts_1), convert_dict_to_str(gend_counts) )
        if 'Gene_Intron' in this_counts_1 or 'Gene_Exon' in this_counts_1 or len(en_list)>0:
           if 'IG' in this_counts_1: del this_counts_1['IG']
           if 'IG' in this_counts_2: del this_counts_2['IG']
        if 'enhD' in en_list or 'enhP' in en_list or 'prom' in en_list or 'K4m3' in en_list or 'enh' in en_list:
           if 'Gene_Intron' in this_counts_1: del this_counts_1['Gene_Intron']

        if 'No_repeat' in this_counts_3 and this_counts_3['No_repeat']< matched_annotations.shape[0] - this_counts_3['No_repeat']:
           del this_counts_3['No_repeat']

        if 'Gene_Exon' in this_counts_1:
           print("Exon", region)

        en_dis = []
        en_cod_type = list(this_counts_1.keys())
        for _c in en_cod_type:
           if '_' in _c and _c.split('_')[0] in ['enhD', 'enhP', 'prom', 'K4m3', 'enh'] and _c.split('_')[0] in en_cod_type:
              en_dis.append( _c )
        for _c in en_dis:
           if _c in this_counts_1: del this_counts_1[_c]
           if _c in this_counts_2: del this_counts_2[_c]

        for _c in this_counts_1:
           total_counts_1[ _c ] += 1
        for _c in this_counts_2:
           total_counts_2[ _c ] += 1
        for _c in this_counts_3:
           total_counts_3[ _c ] += 1
        for _c in this_counts_4:
           total_counts_4[ _c ] += 1
        for _c in this_counts_4:
           total_counts_p4[ _c ] += this_counts_4[ _c ]

    return total_counts_1, total_counts_2, total_counts_3, total_counts_4, total_counts_p4

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

    colors = plt.get_cmap('Blues')(np.linspace(0.2, 0.7, len(sizes)))

    plt.figure(figsize=(7, 7))
    plt.pie(sizes, labels=labels, autopct='%1.0f%%', startangle=140, colors=colors, textprops={'fontsize': 16})
    plt.savefig( fig_name, dpi=600, bbox_inches='tight')
    plt.close('all')



if __name__=='__main__':
   # Input files (modify paths if needed)
   region_file = sys.argv[1]
   annotation_file = sys.argv[2]

   # Process data
   total_counts_1, total_counts_2, total_counts_3, total_counts_4, total_counts_p4 = process_regions(region_file, annotation_file)

   save_fsp = sys.argv[3].strip().split('/')
   if save_fsp[-1]=='':
      save_f  = '/'.join(save_fsp);
   else:
      save_f  = '/'.join(save_fsp[:-1])
   if not os.path.isdir( save_f ):
      os.system('mkdir '+save_f)

   # Plot pie charts
   plot_pie_chart(total_counts_1, "Fraction of Occurrences: Gene Intron, Gene Exon, IG, and ENCODE Types", sys.argv[3]+'_inex2.png')
   plot_pie_chart(total_counts_2, "Fraction of Occurrences: Gene Body, IG, and ENCODE Types", sys.argv[3]+'_body2.png')
   plot_pie_chart(total_counts_3, "", sys.argv[3]+'_repeat2.png')
   plot_pie_chart(total_counts_4, "", sys.argv[3]+'_epic2.png')
   plot_pie_chart(total_counts_p4, "", sys.argv[3]+'_epicP2.png')



