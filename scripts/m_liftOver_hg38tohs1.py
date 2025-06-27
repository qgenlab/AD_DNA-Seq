
import os,sys

from com_fun import *

liftover_program = sys.argv[1]
chan_file = sys.argv[2]
hs_1_inputs = sys.argv[3].split(',')

for _csv_ in hs_1_inputs:
   if '_cluster_dms.csv' in _csv_:
      df_meth_dmp_forlift = pl.read_csv(_csv_, separator=',', infer_schema_length=100000).to_pandas();
      sign_df_meth_dmp_forlift = df_meth_dmp_forlift[df_meth_dmp_forlift['adj.P.Val']<=0.05]
      bed_df = pd.DataFrame({'chr':sign_df_meth_dmp_forlift['Chr'], 'start':sign_df_meth_dmp_forlift['Pos'], 'end':sign_df_meth_dmp_forlift['Pos']+1 });
   else:
      sign_df_meth_dmp_forlift = pl.read_csv(_csv_, separator=',', infer_schema_length=100000).to_pandas();
      bed_df = pd.DataFrame({'chr':sign_df_meth_dmp_forlift['chromosome'], 'start':sign_df_meth_dmp_forlift['start'], 'end':sign_df_meth_dmp_forlift['end'] });
      
   bed_df.to_csv(_csv_[:-4]+'.bed', sep='\t', index=False, header=False)
   #                                   0                   1          2
   cmd = " {0} {1} {2} {3} {4}".format(liftover_program, _csv_[:-4]+'.bed', chan_file, _csv_[:-4]+('_tohg38.bed' if 'hs1_' in _csv_ else '_tohs1.bed'), _csv_[:-4]+'_unmapped.bed')
   print(cmd); sys.stdout.flush()
   os.system(cmd);


