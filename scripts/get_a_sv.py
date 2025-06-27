
from com_fun import *



def read_vcf_to_dataframe_dynamic_columns(filename, target_region=None):
    """Reads a VCF file and returns a pandas DataFrame with dynamic column names."""
    data = []
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
                if not target_region==None:
                    if not (record['CHROM'] == target_region[0] and int(record['POS'])>=target_region[1] and int(record['POS'])<=target_region[2]):
                       continue;
                data.append(record)

    df = pd.DataFrame(data)
    if df.shape[0]>0:
       df['POS'] = df['POS'].astype(int)
       df['QUAL'] = df['QUAL'].astype(int)
    return df

if __name__=='__main__':
   mchr = sys.argv[1]
   start = int(sys.argv[2])
   end = int(sys.argv[3])

   pos_list = []
   for c_vcf in sys.argv[4:]:
      # Example usage (replace 'your_vcf_file.vcf' with your file)
      vcf_df = read_vcf_to_dataframe_dynamic_columns(c_vcf, [mchr, start, end])

      if vcf_df.shape[0]>0:      
         # Print the first few rows of the DataFrame
         print(c_vcf, vcf_df)
         if 'SVTYPE=DEL' in vcf_df.iloc[0]['INFO']:
            pos_list.append(( vcf_df.iloc[0]['POS'], '-'+vcf_df.iloc[0]['REF'] ) )
         else:
            pos_list.append(( vcf_df.iloc[0]['POS'], '+'+vcf_df.iloc[0]['ALT'] ) )

   for _pk in sorted(pos_list):
      print(_pk[0], _pk[1])


