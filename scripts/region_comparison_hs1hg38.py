
import sys,os

import pandas as pd
from intervaltree import Interval, IntervalTree

import polars as pl


def find_overlaps(file1_path, file2_path, file3_path):
    """
    Finds overlaps between genomic regions defined in two files.

    Args:
        file1_path (str): Path to the first CSV file (chr,start,end).
        file2_path (str): Path to the second file (may have #Split lines).

    Returns:
        pandas.DataFrame: DataFrame containing the overlap information.
    """

    # 1. Read the first file
    df1 = pd.read_csv(file1_path, header=0, names=["chr", "start", "end"], usecols=[0,1,2])
    df1['start'] = df1['start'].astype(int)
    df1['end'] = df1['end'].astype(int)
    print(file1_path, df1.shape)
    df2 = pd.read_csv(file2_path, sep='\t', header=None, names=["chr", "start", "end"], usecols=[0,1,2])
    df2['start'] = df2['start'].astype(int)
    df2['end'] = df2['end'].astype(int)
    print(file2_path, df2.shape)

    # 2. Read and process the second file
    regions3 = []
    with open(file3_path, 'r') as f:
        current_chr = None
        current_start = None
        for line in f:
            line = line.strip()
            if line.startswith("#"):
                continue  # Skip #Split lines
            parts = line.split()
            if len(parts) == 3:
                chr2 = parts[0]
                start2 = int(parts[1])
                end2 = int(parts[2])
                regions3.append({"chr": chr2, "start": start2, "end": end2})
    df3 = pd.DataFrame(regions3)
    print(file3_path, df3.shape)

    # 3. Find overlaps
    not_in_df = []
    overlaps12 = []
    df_2_dict = {}
    not1_count = 0;
    for row1 in df1.itertuples():
        chr1, start1, end1 = row1.chr, row1.start, row1.end
        ovlp_t = 0;
        for row2 in df2.itertuples():
            chr2, start2, end2 = row2.chr, row2.start, row2.end
            t_r_str = "{}:{}-{}".format(chr2, start2, end2)
            if t_r_str not in df_2_dict: df_2_dict[ t_r_str ] = 0;
            if chr1 == chr2:  # Check chromosome match
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                if overlap_start < overlap_end:  # Check for actual overlap
                    ovlp_t += 1;
                    df_2_dict[ t_r_str ] += 1;
                    overlaps12.append({
                        "chr": chr1,
                        "start1": start1,
                        "end1": end1,
                        "start2": start2,
                        "end2": end2,
                        "overlap_start": overlap_start,
                        "overlap_end": overlap_end,
                        "overlap_length": overlap_end - overlap_start
                    })
        if ovlp_t==1: pass
        else:
           not1_count += 1 
           print('df1', ovlp_t, not1_count, chr1+":{}-{}".format( start1, end1))
        if ovlp_t==0:
           not_in_df.append( {"chr": chr1, "start": start1, "end": end1 }) 
    print('overlaps12', len(overlaps12))
    not1_count = 0
    for df2_key in sorted(list(df_2_dict.keys())):
       if df_2_dict[ df2_key ] ==1: pass;
       else:
           not1_count += 1 
           print('df2', not1_count, df_2_dict[ df2_key ], df2_key)
           not1_count += 1
    pl.from_pandas( pd.DataFrame(not_in_df) ).write_csv( file1_path[:-4]+'_notinhg38.csv', separator='\t' );

    overlaps13 = []
    for row1 in df1.itertuples():
        chr1, start1, end1 = row1.chr, row1.start, row1.end
        for row2 in df3.itertuples():
            chr2, start2, end2 = row2.chr, row2.start, row2.end
            if chr1 == chr2:  # Check chromosome match
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                if overlap_start < overlap_end:  # Check for actual overlap
                    overlaps13.append({
                        "chr": chr1,
                        "start1": start1,
                        "end1": end1,
                        "start2": start2,
                        "end2": end2,
                        "overlap_start": overlap_start,
                        "overlap_end": overlap_end,
                        "overlap_length": overlap_end - overlap_start
                    })
    print('overlaps13', len(overlaps13))

    overlaps23 = []
    for row1 in df2.itertuples():
        chr1, start1, end1 = row1.chr, row1.start, row1.end
        for row2 in df3.itertuples():
            chr2, start2, end2 = row2.chr, row2.start, row2.end
            if chr1 == chr2:  # Check chromosome match
                overlap_start = max(start1, start2)
                overlap_end = min(end1, end2)
                if overlap_start < overlap_end:  # Check for actual overlap
                    overlaps23.append({
                        "chr": chr1,
                        "start1": start1,
                        "end1": end1,
                        "start2": start2,
                        "end2": end2,
                        "overlap_start": overlap_start,
                        "overlap_end": overlap_end,
                        "overlap_length": overlap_end - overlap_start
                    })
    print( 'overlaps23', len(overlaps23))

    return pd.DataFrame(overlaps12), pd.DataFrame(overlaps13)

# Example usage:
overlap_df = find_overlaps(sys.argv[1], sys.argv[2], sys.argv[3])
print()
print(overlap_df)

