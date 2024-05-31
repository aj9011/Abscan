import os
import pandas as pd
import sys
import subprocess as sp
import shutil

def rmsd_cal(ori, deep):
	rmsd_cmd = '%s/rmsd_total_v2 %s %s bb H > %s_hchain_rmsd.txt' %(GEAR,ori, deep, deep.replace('.deepab_new_numbering.pdb', ''))
	sp.call(rmsd_cmd, shell = True)
	rmsd_cmdL = '%s/rmsd_total_v2 %s %s bb L > %s_lchain_rmsd.txt' %(GEAR,ori, deep, deep.replace('.deepab_new_numbering.pdb', ''))
	sp.call(rmsd_cmdL, shell = True)


file_name= sys.argv[1]
antigen_name= file_name.split("_")[1]

cd = os.getcwd()

GEAR = '/lwork01/neoscan_gear'
input_dir = '%s/%s/%s/input'%(cd,antigen_name,file_name)
pdb_dir='%s/%s/%s/PDB' %(cd,antigen_name,file_name)
develop_dir = '%s/%s/%s/develop' %(cd,antigen_name,file_name)
deepab_dir = '%s/%s/%s/DeepAb/result' %(cd,antigen_name,file_name)
ori_dir= [os.path.join(pdb_dir,x) for x in os.listdir(pdb_dir) if x.endswith('_HL.pdb')][0]
ori_pdb = [os.path.join(input_dir,x) for x in os.listdir(input_dir) if x.endswith('_antigen_filled.pdb')][0]

tool_path= '/tools2/abscan_tools'


filtered_df =  pd.read_csv('%s/%s_pyig_tap_filtered.tsv'%(develop_dir,file_name),sep='\t')

filtered_df['antigen_cdr'] = filtered_df['PDB'].str.replace('.deepab_new_numbering.pdb', '', regex=False)
#print(filtered_df)

filtered_df['cdr_index'] = filtered_df['antigen_cdr'].str.split("_").str[1]

print(filtered_df)


rmsd_dir= '%s/%s_sec_final_rmsd_stats.txt' %(deepab_dir,antigen_name)

rmsd_df =pd.read_csv(rmsd_dir ,sep='\t',names=['epitope','cdr_index','rmsd'])

merged_temp= pd.merge(filtered_df,rmsd_df, how ='left', on = ['cdr_index'])
#print(merged_temp)

merged = merged_temp.drop_duplicates(subset='PDB', keep='first')

top50_rmsd= merged.sort_values('rmsd').head(50)

top50_pdb= top50_rmsd['PDB'].tolist()


top50_dir= [os.path.join(develop_dir,x) for x in top50_pdb]

for top50 in top50_dir:
	rmsd_cal(ori_dir,top50)


h_rmsd_list= []
l_rmsd_list= []

for top50 in top50_dir:
	h_rmsd_file = '%s_hchain_rmsd.txt' %(top50.replace('.deepab_new_numbering.pdb', ''))
	l_rmsd_file = '%s_lchain_rmsd.txt' %(top50.replace('.deepab_new_numbering.pdb', ''))
	with open(h_rmsd_file,'r') as hr:
		h_line = hr.readlines()[0]
		h_rmsd_list.append(h_line.strip())
	with open(l_rmsd_file,'r') as lr:
		l_line = lr.readlines()[0]
		l_rmsd_list.append(l_line.strip())

h_result_dir= '%s/%s_hchain_rmsd.txt' %(develop_dir,antigen_name)
l_result_dir= '%s/%s_lchain_rmsd.txt' %(develop_dir,antigen_name)

with open(h_result_dir, 'w' ) as hw:
	hw.write('\n'.join(h_rmsd_list))

with open(l_result_dir, 'w') as lw:
	lw.write('\n'.join(l_rmsd_list))





Lpass_pdb=[]
Hpass_pdb=[]

with open(l_result_dir, 'r') as lr:
	l_lines = lr.readlines()
	for l_line in l_lines:
		rmsd= float(l_line.strip().split('\t')[-1])
		if rmsd < 8:  #for test
		#if rmsd < 2:
			Lpass_pdb.append(l_line.strip().split('\t')[0])

with open(h_result_dir, 'r') as hr:
	h_lines = hr.readlines()
	for h_line in h_lines:
		rmsd= float(h_line.strip().split('\t')[-1])
		if rmsd < 8:
			Hpass_pdb.append(h_line.strip().split('\t')[0])

print(Hpass_pdb)
print(Lpass_pdb)
s1=set(Hpass_pdb)
s2=set(Lpass_pdb)
pass_pdb = s1.intersection(s2)



snugdock_input_list_temp= list(pass_pdb)
snugdock_input_list = [os.path.join(deepab_dir,x.split("/")[-1].replace(".deepab_new_numbering.pdb",".deepab_tr.pdb")) for x in pass_pdb]

print ('XXXXX')
print(snugdock_input_list)

if not os.path.exists('%s/snugdock_input' %(develop_dir)):
	os.makedirs('%s/snugdock_input'%(develop_dir))
	

for snugdock_input in snugdock_input_list:
	shutil.copy(snugdock_input, '%s/snugdock_input'%(develop_dir))
	deep = os.path.basename(snugdock_input)
	L_lines=[]
	H_lines=[]
	A_lines=[]
	with open('%s/snugdock_input/%s'%(develop_dir,deep),'r') as r:
		lines= r.readlines()
		for line in lines:
			if line.startswith('ATOM') and line[21]=='L':
				L_lines.append(line.strip())
			elif line.startswith('ATOM') and line[21]=='H':
				H_lines.append(line.strip())

	with open(ori_pdb,'r') as a_r:
		antigen_lines= a_r.readlines()
		for antigen_A in antigen_lines:
			if antigen_A.startswith('ATOM') and antigen_A[21]=='A':
				A_lines.append(antigen_A.strip())
	L_lines.append('TER\n')
	H_lines.append('TER\n')
	with open('%s/snugdock_input/%s' %(develop_dir,deep.replace(".deepab_tr.pdb", ".LHA.pdb")),'w') as w:
		L= '\n'.join(L_lines)
		H= '\n'.join(H_lines)
		A= '\n'.join(A_lines)
		final_lines = L + H + A
		w.write(final_lines)
	
	LHA_pdb= '%s' %(deep.replace(".deepab_tr.pdb", ".LHA.pdb"))
	with open('%s/snugdock_input/%s'%(develop_dir,LHA_pdb.replace("LHA.pdb","LHA_sorted.pdb")),'w') as ff:
		with open('%s/snugdock_input/%s'%(develop_dir,LHA_pdb),'r') as ff1:
			lines = ff1.readlines()
			count=1
			for line in lines:
				if line.startswith('ATOM'):
					new_line= line[0:6]+ str(count).rjust(5,' ') + line[11:]
					ff.write(new_line)
					count += 1
				else:
					ff.write(line)


uniq_list=filtered_df['PDB'].tolist()
filtered_df.to_csv('%s/develop_pyig_tap_info.tsv' %(develop_dir),sep='\t' ,index =False)
print(filtered_df['PDB'])
pass_index = [x.split("/")[-1] for x in snugdock_input_list ]
top50_df = filtered_df[filtered_df['PDB'].isin(pass_index)]
#print("dddddddddddddddddddddddddddddddddddddd",top50_df)
top50_df.to_csv('%s/develop_top50.tsv' %(develop_dir),sep='\t' ,index =False)
