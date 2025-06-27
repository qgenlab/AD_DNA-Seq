
import os,sys
import pandas as pd
import polars as pl
import numpy as np


def find_closest_region_endpoints(row_a, df_b):
    """
    Finds the closest region in df_b to a given row (region) in df_a.

    Args:
        row_a (pd.Series): A row from df_a with 'chr', 'start', 'end' columns.
        df_b (pd.DataFrame): The other dataframe with 'chr', 'start', 'end' columns.

    Returns:
        pd.Series or None: The closest row from df_b, or None if no match on 'chr'.
    """
    chr_a = row_a['chr']
    start_a = row_a['start']
    end_a = row_a['end']

    df_b_chr = df_b[df_b['chr'] == chr_a].copy()

    if df_b_chr.empty:
        return None

    def calculate_distance(row_b):
        start_b = row_b['start']
        end_b = row_b['end']
        distances = [
            abs(start_a - start_b),
            abs(start_a - end_b),
            abs(end_a - start_b),
            abs(end_a - end_b)
        ]
        return min(distances)


    df_b_chr['distance'] = df_b_chr.apply(calculate_distance, axis=1)
    closest_row_b = df_b_chr.loc[df_b_chr['distance'].idxmin()]
    return closest_row_b, df_b_chr['distance'].min()

def _find_close(df1, df2):
   close_regin_list = []
   for df1i in range(df1.shape[0]):
      closerow_ = find_closest_region_endpoints(df1.iloc[df1i], df2);
      if closerow_==None:
          continue;
      t_dis = closerow_[1]
      closerow_ = closerow_[0]
      if t_dis<500000:
         close_regin_list.append([ "{}:{}-{}".format(df1.iloc[df1i]['chr'],df1.iloc[df1i]['start'], df1.iloc[df1i]['end'] ), "{}".format(t_dis), "{}:{}-{}".format(closerow_['chr'],closerow_['start'], closerow_['end'] ) ] )
   print(len(close_regin_list))
   for _r in close_regin_list:
      print(' '.join(_r))

   print()

if __name__=='__main__':
   dmr_region_file = sys.argv[2]
   svr_region_file = sys.argv[1]

   dmr_regions_df = pd.read_csv(dmr_region_file, sep=',', header=0, names=['chr', 'start', 'end'], usecols=[0,1,2])
   svr_regions_df = pd.read_csv(svr_region_file, sep=',', header=0, names=['chr', 'start', 'end'], usecols=[0,1,2])

   _find_close(dmr_regions_df, svr_regions_df)
   _find_close(svr_regions_df, dmr_regions_df)


