
import os,sys
import pandas as pd
import polars as pl
from intervaltree import Interval, IntervalTree
import copy


def check_ovlp(this_set, check_df):
   ovlp_set = set()
   for t_r in this_set:
      t_r_sp = t_r.split(':');
      t_r_sp[1] = t_r_sp[1].split('-')
      t_r_sp[1][0] = int(t_r_sp[1][0])
      t_r_sp[1][1] = int(t_r_sp[1][1]);
      check_i = 0;
      for _t_df in check_df:
         t_c_df = _t_df[ _t_df['chr']==t_r_sp[0] ]
         if t_c_df.empty: continue;
         else:
             for _i in range(t_c_df.shape[0]):
                max_start = max(t_r_sp[1][0], t_c_df.iloc[_i]['start'] );
                min_end = min( t_r_sp[1][1], t_c_df.iloc[_i]['end'] );
                if min_end - max_start >0: 
                   check_i += 1;
                   break;
      if check_i==len(check_df):
          ovlp_set.add( t_r )
   return ovlp_set

def merge_dataframes_and_overlap(df1, df2, df3, shared_all, identical_share):
    """
    Merges rows from three DataFrames based on a shared key, then merges overlapping rows.

    Args:
        df1, df2, df3 (pd.DataFrame): The three DataFrames to merge.
        shared_all (set): A set of keys to filter rows from the DataFrames.

    Returns:
        pd.DataFrame: A merged DataFrame with overlapping rows merged.
    """

    # Filter rows based on the shared key
    filtered_df1 = df1[df1['key'].isin(shared_all)]
    filtered_df2 = df2[df2['key'].isin(shared_all)]
    filtered_df3 = df3[df3['key'].isin(shared_all)]

    # Merge the filtered DataFrames
    merged_df = pd.concat([filtered_df1, filtered_df2, filtered_df3], ignore_index=True)

    # Function to merge overlapping rows within a chromosome
    def merge_overlapping_rows(group):
        group = group.sort_values(by=['start'])  # Assuming 'start' and 'end' columns
        merged_rows = []
        i = 0
        while i < len(group):
            current_row = group.iloc[i].copy()
            j = i + 1
            while j < len(group):
                next_row = group.iloc[j]
                if current_row['end'] >= next_row['start']:
                    if (group.iloc[j]['key'] in identical_share ) and (group.iloc[j-1]['key'] in identical_share ):
                       j += 1
                       continue;
                    current_row['start'] = min(current_row['start'], next_row['start'])
                    current_row['end'] = max(current_row['end'], next_row['end'])
            
                    current_row['start_min'] = min(current_row['start_min'], next_row['start_min'])
                    current_row['end_min'] = min(current_row['end_min'], next_row['end_min'])
                    current_row['start_max'] = max(current_row['start_max'], next_row['start_max'])
                    current_row['end_max'] = max(current_row['end_max'], next_row['end_max'])
                   
                    current_row['AD_freq'] = (current_row['AD_freq'] + next_row['AD_freq'] )/2
                    current_row['Ctrl_freq'] = (current_row['Ctrl_freq'] + next_row['Ctrl_freq'])/2
                    current_row['svinfo'] = current_row['svinfo']  + '/' + next_row['svinfo']
 
                    j += 1
                else:
                    break
            merged_rows.append(current_row)
            i = j
        return pd.DataFrame(merged_rows)

    merged_df = merged_df.groupby('chr').apply(merge_overlapping_rows).reset_index(drop=True)

    return merged_df


def add_key(df):
    return df.assign(key=df['chr'] + ':' + df['start'].astype(str) + '-' + df['end'].astype(str))

def get_ovlp(df1, df2, df3, file1, file2, file3):
   df1 = add_key(df1)
   df2 = add_key(df2)
   df3 = add_key(df3)
   
   # Sets for comparison
   s1, s2, s3 = set(df1['key']), set(df2['key']), set(df3['key'])
   
   # Find overlaps
   shared_all = s1 & s2 & s3
   shared_12 = (s1 & s2) - shared_all
   shared_13 = (s1 & s3) - shared_all
   shared_23 = (s2 & s3) - shared_all
   unique_1 = s1 - s2 - s3
   unique_2 = s2 - s1 - s3
   unique_3 = s3 - s1 - s2

   identical_share = copy.deepcopy( shared_all )

   print(f"Shared by all 3: {len(shared_all)}")
   print(f"Shared by {file1} and {file2} only: {len(shared_12)}")
   print(f"Shared by {file1} and {file3} only: {len(shared_13)}")
   print(f"Shared by {file2} and {file3} only: {len(shared_23)}")
   print(f"Unique to {file1}: {len(unique_1)}")
   print(f"Unique to {file2}: {len(unique_2)}")
   print(f"Unique to {file3}: {len(unique_3)}")
   print()

   check_list = [ [shared_12, [df3] ], \
                  [shared_13, [df2] ], \
                  [shared_23, [df1] ], \
                  [unique_1,  [df2, df3] ], \
                  [unique_2,  [df1, df3] ], \
                  [unique_3,  [df1, df2] ] ]

   for _c_l_i in range(len(check_list)):
       part_ovlp = check_ovlp(check_list[_c_l_i][0], check_list[_c_l_i][1]  )
       if len( part_ovlp )>0:
           shared_all.update( part_ovlp )
           #check_list[_c_l_i][0] = check_list[_c_l_i][0] - part_ovlp
           check_list[_c_l_i][0].difference_update( part_ovlp )
   
   print(f"Shared by all 3: {len(shared_all)}")
   print(f"Shared by {file1} and {file2} only: {len(shared_12)}")
   print(f"Shared by {file1} and {file3} only: {len(shared_13)}")
   print(f"Shared by {file2} and {file3} only: {len(shared_23)}")
   print(f"Unique to {file1}: {len(unique_1)}")
   print(f"Unique to {file2}: {len(unique_2)}")
   print(f"Unique to {file3}: {len(unique_3)}")

   merge_df = merge_dataframes_and_overlap(df1, df2, df3, shared_all, identical_share)
   pl.from_pandas( merge_df ).write_csv( "sv_2025Apr/merge_sv_region.csv" );

   return df1, df2, df3, shared_all, shared_12, shared_13, shared_23, unique_1, unique_2, unique_3

def build_interval_tree(df):
    tree = {}
    for _, row in df.iterrows():
        tree.setdefault(row['chr'], IntervalTree()).addi(row['start'], row['end'])
    return tree

def get_partial_ovlp(shared_12, shared_13, shared_23, unique_1, unique_2, unique_3, df1, df2, df3):
   tree1 = build_interval_tree(df1)
   tree2 = build_interval_tree(df2)
   tree3 = build_interval_tree(df3)
   #print(tree1)
   
   # Check overlaps for regions unique to file1
   def check_overlap(key, df_lookup, tree_other):
       row = df_lookup[df_lookup['key'] == key].iloc[0]
       #print('In check_overlap', key )
       #print('In check_overlap', tree_other[ row['chr'] ]  )
       overlaps = tree_other[row['chr']].overlap(row['start'], row['end']) if row['chr'] in tree_other else []
       #print ('In check_overlap', overlaps)
       return overlaps
   
   check_list = [['shared_12', shared_12, [ tree3 ], df2], \
                 ['shared_13', shared_13, [ tree2 ], df1], \
                 ['shared_23', shared_23, [ tree1 ], df3], \
                 ['unique_1', unique_1, [tree2, tree3], df1], \
                 ['unique_2', unique_2, [tree1, tree3], df2], \
                 ['unique_3', unique_3, [tree1, tree2], df3] ]

   for _t_cl in check_list:
      print( _t_cl[0], len(_t_cl[1]) );
      for key in list(_t_cl[1]):  # just show a few for brevity
         ovlpi = 0
         for _t_tree in _t_cl[2]:
            ovlpi = ovlpi + 1
            t_ovlp = check_overlap(key, _t_cl[3], _t_tree)
            if len( t_ovlp )>0:
               print(f"\t{key}: overlaps with {ovlpi}: {t_ovlp}")
            else: pass; #print("\t!!!No ovlp", key)

if __name__=='__main__':

   df1 = pd.read_csv(sys.argv[1], sep=',')
   df2 = pd.read_csv(sys.argv[2], sep=',')
   df3 = pd.read_csv(sys.argv[3], sep=',')

   print('rows:', df1.shape[0], df2.shape[0], df3.shape[0] )

   df1['start'] = df1['start'].astype(int)
   df2['start'] = df2['start'].astype(int)
   df3['start'] = df3['start'].astype(int)
   df1['end'] = df1['end'].astype(int)
   df2['end'] = df2['end'].astype(int)
   df3['end'] = df3['end'].astype(int)

   df1['end'] = df1['end'] + 1
   df2['end'] = df2['end'] + 1
   df3['end'] = df3['end'] + 1
  
   df1, df2, df3, shared_all, shared_12, shared_13, shared_23, unique_1, unique_2, unique_3 = get_ovlp(df1, df2, df3, sys.argv[1], sys.argv[2], sys.argv[3])

   get_partial_ovlp(shared_12, shared_13, shared_23, unique_1, unique_2, unique_3, df1, df2, df3)     



