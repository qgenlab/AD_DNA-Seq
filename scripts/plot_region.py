
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

repeat_colors = {
    "SINE": "purple",
    "LINE": "blue",
    "LTR": "green",
    "Simple_repeat": "cyan",
    "Satellite": "brown",
    "Low_complexity": "pink",
    "DNA": "gray"
}

enhancer_colors = {
    "enhP": "red",
    "enhD": "orange",
    "prom": "gold"
}

#def sliding_window_avg(df, start_col, window_size=100, step_size=50, sample_start_ind=3):
def sliding_window_avg(df, start_col, window_size=50, step_size=25, sample_start_ind=3):
    """
    Compute sliding window average methylation for each sample.
    df: DataFrame with ['chromosome', 'start'] + sample columns.
    start_col: Column name for start positions.
    window_size: Size of the window (default=50bp).
    step_size: Step size between windows (default=25bp).
    Returns a new DataFrame with windowed average methylation.
    """
    new_rows = []
    min_pos, max_pos = df[start_col].min(), df[start_col].max()

    for start in range(min_pos, max_pos - window_size + 1, step_size):
        end = start + window_size
        window_data = df[(df[start_col] >= start) & (df[start_col] < end)]
        if not window_data.empty:
            avg_values = window_data.iloc[:, sample_start_ind:].mean().to_dict()  # Compute mean for each sample
            avg_values[start_col] = start + window_size // 2  # Use center position
            new_rows.append(avg_values)

    return pd.DataFrame(new_rows)

def plot_methylation_curve(region_df, methylation_df, group1_samples, group2_samples, base_plt_folder, repeat_regions_df, enhancer_promoter_df, window_size=50, step_size=25):
    """
    Generate smoothed methylation plots using sliding window averaging.
    - Left subplot: Individual sample curves.
    - Right subplot: Mean ± STD curves for two groups.
    
    region_df: BED format DataFrame ['chromosome', 'start', 'end'].
    methylation_df: DataFrame with ['chromosome', 'start'] + sample columns.
    group1_samples, group2_samples: Lists of sample names for each group.
    """
    annotation_height = 2
    add_y = 5
    min_std_area = 7.5;
    sep_area_dict = []
    
    for _, region in region_df.iterrows():
        chrom, start, end = region['chromosome'], region['start'], region['end']
        base_fn = "{}-{}-{}".format(chrom, start, end)
        reg_len = end - start;
        if reg_len<400:
           reg_len = int((400-reg_len)/2)
           start -= reg_len
           if start<0: start = 0
           end += reg_len

        # Extract methylation data for the region
        region_data = methylation_df[
            (methylation_df['Chr'] == chrom) & 
            (methylation_df['Pos'] >= start) & 
            (methylation_df['Pos'] <= end)
        ].copy()

        if region_data.empty:
            print(f"No data for region: {chrom}:{start}-{end}")
            continue

        # Compute sliding window averages
        windowed_data = sliding_window_avg(region_data, 'Pos', window_size, step_size, sample_start_ind=3)
        #print(windowed_data)
        if windowed_data.empty:
            print(f"No data for window_data: {chrom}:{start}-{end}")
            print( region_data )
            continue

        # Prepare figure with two subplots
        fig, axes = plt.subplots(1, 2, figsize=(15, 5))

        # === Left Subplot: Individual Sample Curves ===
        for sample in group1_samples:
            if sample in windowed_data.columns:
                axes[0].plot(windowed_data['Pos'], windowed_data[sample], label=sample, color='blue', alpha=0.5)

        for sample in group2_samples:
            if sample in windowed_data.columns:
                axes[0].plot(windowed_data['Pos'], windowed_data[sample], label=sample, color='red', alpha=0.5)

        axes[0].set_xlabel("Genomic Coordinate", fontsize=13)
        axes[0].set_ylabel("Methylation Level", fontsize=13)
        axes[0].tick_params(axis='both', labelsize=13)

        # === Right Subplot: Mean ± STD Curves for Groups ===
        windowed_data['Group1_Mean'] = windowed_data[group1_samples].mean(axis=1)
        windowed_data['Group1_STD'] = windowed_data[group1_samples].std(axis=1)
        windowed_data['Group2_Mean'] = windowed_data[group2_samples].mean(axis=1)
        windowed_data['Group2_STD'] = windowed_data[group2_samples].std(axis=1)

        x_range = [windowed_data['Pos'].iloc[0], windowed_data['Pos'].iloc[-1] ]

        axes[1].plot(windowed_data['Pos'], windowed_data['Group1_Mean'], label="AD", color='blue', linestyle='-')
        axes[1].fill_between(windowed_data['Pos'], windowed_data['Group1_Mean'] - windowed_data['Group1_STD'], 
                             windowed_data['Group1_Mean'] + windowed_data['Group1_STD'], color='blue', alpha=0.2)

        axes[1].plot(windowed_data['Pos'], windowed_data['Group2_Mean'], label="Control", color='red', linestyle='-')
        axes[1].fill_between(windowed_data['Pos'], windowed_data['Group2_Mean'] - windowed_data['Group2_STD'], 
                             windowed_data['Group2_Mean'] + windowed_data['Group2_STD'], color='red', alpha=0.2)

        axes[1].set_xlabel("Genomic Coordinate", fontsize=13)
        axes[1].set_ylabel("Average Methylation ± STD", fontsize=13)
        # Optionally, increase tick label font size for better visibility
        axes[1].tick_params(axis='both', labelsize=13)
        #axes[1].set_title(f"Smoothed Methylation: {chrom}:{start}-{end}")
        def_loc = 'center right'
        low_part = 20;
        is_upper = False
        if windowed_data['Group1_Mean'].iloc[-1]-windowed_data['Group1_STD'].iloc[-1]>low_part and  windowed_data['Group2_Mean'].iloc[-1]-windowed_data['Group2_STD'].iloc[-1]>low_part:
           def_loc = 'lower right'
        elif windowed_data['Group1_Mean'].iloc[0]-windowed_data['Group1_STD'].iloc[0]>low_part and  windowed_data['Group2_Mean'].iloc[0]-windowed_data['Group2_STD'].iloc[0]>low_part:
           def_loc = 'lower left'
        elif windowed_data['Group1_Mean'].iloc[0]+windowed_data['Group1_STD'].iloc[0]<100-low_part and  windowed_data['Group2_Mean'].iloc[0]+windowed_data['Group2_STD'].iloc[0]<100-low_part:
           def_loc = 'upper right'
           is_upper = True
        elif windowed_data['Group1_Mean'].iloc[-1]+windowed_data['Group1_STD'].iloc[-1]<100-low_part and  windowed_data['Group2_Mean'].iloc[-1]+windowed_data['Group2_STD'].iloc[-1]<100-low_part:
           def_loc = 'upper left'
           is_upper = True
        if is_upper:
           bbox = plt.gca().get_window_extent().transformed(plt.gca().transData.inverted())
           bbox_height_inches = bbox.height
           bbox_height_points = bbox_height_inches * 72  # 1 inch = 72 points
           # Calculate the offset in axes coordinates
           offset_points = 10
           offset_axes_y = -offset_points / bbox_height_points
           axes[1].legend(fontsize=11, loc=def_loc, ncol=1, bbox_to_anchor=(0, 0 + offset_axes_y), bbox_transform=plt.gca().transAxes)
        else:
           axes[1].legend(fontsize=11, loc=def_loc, ncol=1)

        xmin, xmax = max(axes[0].get_xlim()[0], axes[1].get_xlim()[0]), min(axes[0].get_xlim()[1], axes[1].get_xlim()[1])

        if not (repeat_regions_df is None):
           #print( axes[0].get_xlim(), axes[1].get_xlim() )
           repeat_in_region = repeat_regions_df[
                  (repeat_regions_df['chromosome'] == chrom) &
                  (repeat_regions_df['start'] <= end) &
                  (repeat_regions_df['end'] >= start)
           ]
           print(base_fn, repeat_in_region.shape)
           for ax_i in range(2):
              ymin, ymax = axes[ax_i].get_ylim()
              annotation_y = ymax + add_y  # Move above the highest methylation signal
              for _, repeat in repeat_in_region.iterrows():
                  repeat_type = repeat.get('labels', "Other")  # Default if type missing
                  color = repeat_colors.get(repeat_type, "gray")  # Default to gray
                  plot_start = xmin if repeat['start']<xmin else repeat['start']
                  plot_end = xmax if repeat['end']>xmax else repeat['end']
                  axes[ax_i].add_patch(plt.Rectangle((plot_start, annotation_y), plot_end - plot_start, annotation_height, 
                                                  color=color, alpha=0.7, label=repeat['labels']))
                  axes[ax_i].text( (plot_end+plot_start) / 2,  # X-center
                       annotation_y + annotation_height,  # Y-center
                       repeat['labels'],
                       ha='center', va='center', fontsize=10, color='black', fontweight='bold'
                   )
              if repeat_in_region.shape[0]>0:
                  axes[ax_i].set_ylim(0, annotation_y + add_y)
        if not (enhancer_promoter_df is None):
           enhancer_in_region = enhancer_promoter_df[
                  (enhancer_promoter_df['chromosome'] == chrom) &
                  (enhancer_promoter_df['start'] <= end) &
                  (enhancer_promoter_df['end'] >= start)
           ]
           print( base_fn, enhancer_in_region.shape)
           for ax_i in range(2):
              ymin, ymax = axes[ax_i].get_ylim()
              annotation_y = ymax + add_y  # Move above the highest methylation signal
              for _, enhancer in enhancer_in_region.iterrows():
                  enhancer_type = enhancer.get('cat', "Other")  # Default if missing
                  color = enhancer_colors.get(enhancer_type, "black")  # Default to black
                  plot_start = xmin if enhancer['start']<xmin else enhancer['start']
                  plot_end = xmax if enhancer['end']>xmax else enhancer['end']
                  axes[ax_i].add_patch(plt.Rectangle((plot_start, annotation_y), plot_end - plot_start, annotation_height, 
                                                  color=color, alpha=0.7, label=enhancer['labels']))
                  axes[ax_i].text((plot_end+plot_start) / 2,  # X-center
                       annotation_y + annotation_height,  # Y-center
                       enhancer['labels'],
                       ha='center', va='center', fontsize=10, color='black', fontweight='bold'
                   )
              if enhancer_in_region.shape[0]>0:
                 axes[ax_i].set_ylim(0, annotation_y + add_y)

        plt.tight_layout()
        plt.savefig(base_plt_folder+'/long_region_plot_'+base_fn+'.png', dpi=600)
        plt.close('all')

        t_n_ear = 0;
        for _i_sa in range(len(windowed_data['Group1_Mean'])):
           t_n_ear += (windowed_data['Group1_Mean'].iloc[_i_sa] - windowed_data['Group2_Mean'].iloc[_i_sa])/( (windowed_data['Group1_STD'].iloc[_i_sa] if windowed_data['Group1_STD'].iloc[_i_sa]>min_std_area else min_std_area)+ (windowed_data['Group2_STD'].iloc[_i_sa] if windowed_data['Group2_STD'].iloc[_i_sa]>min_std_area else min_std_area) )
        t_n_ear = t_n_ear / (len(windowed_data['Group1_Mean']) if len(windowed_data['Group1_Mean']) >=10 else 10)
        sep_area_dict.append([abs(t_n_ear), t_n_ear, base_fn])

    return sorted(sep_area_dict)

def split_and_get_second(value):
    """
    Splits a string by '#' or '/' and returns the second part, or the original string if no split.
    """
    if '#' in value and '/' in value:
       hash_index = value.find('#')
       slash_index = value.find('/')
       return value.split('#')[1].split('/')[0]
    elif '#' in value:
       return value.split('#')[1]
    elif '/' in value:
       return value.split('/')[0]
    else: return value;

if __name__=='__main__':
   input_prefx = sys.argv[1]
   min_value = 10**(-3)
   filter_samples_ratio=0.6
   readThresh=sys.argv[2];

   filter2_meth_df, filter2_covg_df, df_meth_dmp, ad_df_list, ct_df_list, others = filter_input(input_prefx, readThresh, min_value, filter_samples_ratio)

   # Read region file
   regions_df = pd.read_csv(sys.argv[3], sep=',', header=0, names=['chromosome', 'start', 'end'], usecols=[0,1,2])
   base_plt_folder = sys.argv[3][:-4]
   if not os.path.isdir( base_plt_folder ):
       os.system('mkdir '+base_plt_folder)

   if len(sys.argv)>4:
      rep_info = sys.argv[4].split(',')
      default_col = [5,6,7,11]
      default_new = ['chromosome', 'start', 'end', 'labels']
      if len(rep_info)>2:
         default_col = []
         default_new = []
         read_numcol = len(rep_info)//2
         for _i in range(read_numcol):
            default_col.append( int(rep_info[1+_i]) )
            default_new.append(     rep_info[1+_i+read_numcol] )
      print( default_col, default_new)
      repeat_regions_df = pl.read_csv(rep_info[0], separator='\t', has_header=False, columns=default_col, new_columns=default_new, infer_schema_length=10000).to_pandas();
      if repeat_regions_df[default_new[-1]].str.contains('#').any() or repeat_regions_df[default_new[-1]].str.contains('/').any():
         repeat_regions_df['labels'] = repeat_regions_df[default_new[-1]].apply(split_and_get_second)
   else:
      repeat_regions_df = None;

   if len(sys.argv)>5:
      enp_info = sys.argv[5].split(',')
      default_col = [0,1,2,12,13]
      default_new = ['chromosome', 'start', 'end', 'cat', 'id']
      if len(enp_info)>2:
         default_col = []
         default_new = []
         read_numcol = len(enp_info)//2
         for _i in range(read_numcol):
            default_col.append( int(enp_info[1+_i]) )
            default_new.append(     enp_info[1+_i+read_numcol] )
      print( default_col, default_new)
      enhancer_promoter_df = pl.read_csv(enp_info[0], separator='\t', has_header=False, columns=default_col, new_columns=default_new, infer_schema_length=10000).to_pandas(); 
      enhancer_promoter_df['labels'] =  enhancer_promoter_df['cat'] + '_' + enhancer_promoter_df['id']
   else:
      enhancer_promoter_df = None;

   sep_area_dict = plot_methylation_curve(regions_df, filter2_meth_df.reset_index(), ad_df_list[0], ct_df_list[0], base_plt_folder, repeat_regions_df, enhancer_promoter_df)
   #print( sep_area_dict )
   sprint = []
   ri = 0;
   for _sa in sep_area_dict[::-1]:
      ri = ri + 1;
      sprint.append("[{:.3f},{:.3f},{}]".format( _sa[0], _sa[1], _sa[2] ))
      # plt.savefig(base_plt_folder+'/long_region_plot_'+base_fn+'.png', dpi=600)
      os.system('cp {} {}'.format( base_plt_folder+'/long_region_plot_'+_sa[2]+'.png', base_plt_folder+'/R'+('{:02d}'.format(ri))+'long_region_plot_'+_sa[2]+'.png' ))
      if len(sprint)>0 and len(sprint)%4==0:
         print(' '.join( sprint ));
         sprint = []
   print(' '.join( sprint ));
