import os
import pandas as pd
import openpyxl
import sys


file_name = sys.argv[1]
antigen_name = file_name.split("_")[1]

cd= os.getcwd()

develop_dir = '%s/%s/%s/develop' %(cd,antigen_name,file_name)
pyig_dir = '%s/pyig_result' %(develop_dir)
tap_dir = '%s/tap_result' %(develop_dir)


pyig_files= [os.path.join(pyig_dir,x) for x in os.listdir(pyig_dir) if x.endswith("_pyig.tsv")]

tap_files= [os.path.join(tap_dir,x) for x in os.listdir(tap_dir) if x.endswith("_tap.tsv")]


pyig_col = ['PDB', 'pyig_H1','pyig_H2','pyig_H3','pyig_L1','pyig_L2','pyig_L3','H1_flag','H2_flag','H3_flag','L1_flag','L2_flag','L3_flag','pyig_flag']
pyig_total_df = pd.DataFrame()

for pyig_file in pyig_files:
	pyig = pd.read_csv(pyig_file, sep='\t', header= None)
	pyig_total_df = pd.concat([pyig_total_df, pyig], ignore_index=True)

pyig_total_df.columns = pyig_col


tap_col = ['PDB','tap_H1','tap_H2','tap_H3','tap_L1','tap_L2','tap_L3','CDR_len','Hydro_score','PPC','PNC','SFvCSP','CDR_len_DBscore','Hydro_DBscore','PPC_DBscore','PNC_DBscore','SFvCSP_DBscore']
tap_total_df = pd.DataFrame()


for tap_file in tap_files:
	tap = pd.read_csv(tap_file, sep='\t',header=None)
	tap_total_df = pd.concat([tap_total_df, tap], ignore_index=True)

tap_total_df.columns = tap_col


pyig_total_result='%s/%s_pyig_total.tsv' %(pyig_dir,file_name)
pyig_total_df.to_csv(pyig_total_result,sep ='\t',index=False)


tap_total_result='%s/%s_tap_total.tsv' %(tap_dir,file_name) 
tap_total_df.to_csv(tap_total_result,sep='\t',index=False)

df= pd.merge(pyig_total_df,tap_total_df, on = 'PDB' , how = 'outer') 

pyig_tap_total = '%s/%s_pyig_tap_total.tsv' %(develop_dir,file_name)
df.to_csv(pyig_tap_total, sep= '\t',index=False)



#filtered_df = df[(df['H3_flag'] == 'X') & 
#				 (df['H2_flag'] == 'O') &
#				 (df['H1_flag'] == 'X') &
#				 (df['L3_flag'] == 'O') &
#				 (df['L2_flag'] == 'O') &
#				 (df['L1_flag'] == 'O') &
#                 (df['CDR_len_DBscore'] >= 1) & 
#                 (df['Hydro_DBscore'] >= 0) & 
#                 (df['PPC_DBscore'] >= 1) & 
#                 (df['PNC_DBscore'] >= 1) & 
#                 (df['SFvCSP_DBscore'] >= 1)]


####for test#####

filtered_df = df[(df['CDR_len_DBscore'] >= 0) &
                 (df['Hydro_DBscore'] >= 0) &
                 (df['PPC_DBscore'] >= 0) &
                 (df['PNC_DBscore'] >= 0) &
                 (df['SFvCSP_DBscore'] >= 0)]





#print(df)
print(filtered_df)

pyig_tap_filtered = '%s/%s_pyig_tap_filtered.tsv' %(develop_dir,file_name)
filtered_df.to_csv(pyig_tap_filtered, sep= '\t',index=False)
