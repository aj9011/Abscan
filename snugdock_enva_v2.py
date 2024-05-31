import os
import sys
import glob
import shutil
import subprocess as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr


def change_cdrh3_num(parsed_pdb):
	antigens_list= []
	antigen_num = 0
	cdr_num=[]
	with open(parsed_pdb, 'r') as an:
		antigen_nums= an.readlines()
		for i, aa in enumerate(antigen_nums):
			if aa.startswith('ATOM') > 0 :
				if aa[21]=='A':
					if i == len(antigen_nums)-3:
						antigen_num= int(aa[6:11].strip())
						antigens_list.append(aa.rstrip())
					else:
						antigens_list.append(aa.rstrip())
				elif aa[21]=='B':
					cdr_num.append(aa[22:27].strip())
		cdr_region = sorted(list(set(cdr_num)))
		cdr_regions=list(set(cdr_region))
		cdr_regions.sort()
		cdr_regions.sort(key=lambda x : int(x) if len(x) <=2 else (int(x[:-1]) if x[0] == '9'  else (int(x) if len(x) ==3 else int(x[:-1]) )))
		dic2={}
		for nn, n in enumerate(cdr_regions):
			dic2[n] = nn+1
	new_lines_list=[]
	with open(parsed_pdb,'r' ) as cdr_r:
		cdr_lines=cdr_r.readlines()
		for i, cdr_line in enumerate(cdr_lines):
			if cdr_line.startswith('ATOM') > 0 :
				if cdr_line[21]=='B':
					new_numbering = dic2[cdr_line[22:27].strip()]
					new_lines= cdr_line[:6] + str(i+antigen_num+1).rjust(5,' ') + cdr_line[11:22] + str(new_numbering).rjust(4,' ') + " " + cdr_line[27:].rstrip()
					new_lines_list.append(new_lines)
				else:
					new_lines_list.append(cdr_line.strip())
	cdrh3_result = '%s/%s' %(parsed_dock_res1,os.path.basename(parsed_pdb))
	with open(cdrh3_result,'w')  as w_ag_cdr_r:
		w_ag_cdr_r.write('\n'.join(new_lines_list))


tpdb = sys.argv[1]
file_name= sys.argv[2]
antigen_name= file_name.split("_")[1]
pipeline = sys.argv[3]

cd= os.getcwd()

snugdock_dir = '%s/%s/%s/%s/snugdock_input'%(cd,antigen_name,file_name,pipeline)
pdb_dir = '%s/%s/%s/PDB'%(cd,antigen_name,file_name)

GEAR = '/lwork01/neoscan_gear'
tmp = '/lwork01/%s_abs_tmp_data'%(file_name)
snugdock_output = '%s/%s/%s/%s/snugdock_output' %(cd,antigen_name,file_name,pipeline)
os.makedirs(snugdock_output,exist_ok=True)
parsed_dock_res = '%s/%s/parsed_dock_res'%(tmp,tpdb)
parsed_dock_res1 = '%s/%s/parsed_dock_res1'%(tmp,tpdb)
os.makedirs('%s/%s/parsed_dock_res'%(tmp,tpdb),exist_ok=True)
os.makedirs('%s/%s/parsed_dock_res1'%(tmp,tpdb),exist_ok=True)
pdbs = glob.glob('%s/%s/dock_res/*.pdb'%(tmp,tpdb))
bes = []
confs = []
for pdb in pdbs:
	pdb1 = os.path.basename(pdb).replace('.pdb','')
	confs.append(pdb1)
	with open('%s/%s.pdb'%(parsed_dock_res,pdb1),'w') as f:
		with open(pdb,'r') as f1:
			lines = f1.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 :
					if line[21]=='A':
						f.write(line)
					elif line[21]=='H' and (int(line[22:26].strip()) >= 95 and int (line[22:26].strip()) <=102):
						line = line[:21] + 'B' + line[22:]
						f.write(line)
	change_cdrh3_num('%s/%s.pdb'%(parsed_dock_res,pdb1))
	os.system('%s/enva5.0 -k %s/%s.pdb > %s/%s_energy.out'%(GEAR,parsed_dock_res1,pdb1,parsed_dock_res1,pdb1))
	if os.path.exists('%s/%s_energy.out'%(parsed_dock_res1,pdb1)) and os.path.getsize('%s/%s_energy.out'%(parsed_dock_res1,pdb1)) > 0 :
		with open('%s/%s_energy.out'%(parsed_dock_res1,pdb1),'r') as f2:
			lines = f2.readlines()
			for line in lines:
				score= float(line.strip().split(" ")[2])
				bes.append(score)
	else:
		bes.append(0)

df = pd.DataFrame()
df['description'] = confs
df['BE'] = bes
df.to_csv('%s/%s/total_env_score.txt'%(tmp,tpdb),sep='\t',index=False)

pdbs = glob.glob('%s/%s/dock_res/%s.LHA_sorted_0001_*.pdb'%(tmp,tpdb,tpdb))

with open('%s/%s/isc_irms_score.txt'%(tmp,tpdb),'w') as f1:
	f1.write('description\tI_sc\tIrms\tCAPRI_rank\n')
	for pdb in pdbs:
		pdb1 = os.path.basename(pdb).replace('.pdb','')
		with open(pdb,'r') as f:
			lines = f.readlines()
			i_sc_count=0
			rmsd_count =0
			capri_count = 0
			for line in lines:
				if line.startswith('I_sc') > 0 :
					try:
						i_sc =line.strip().split(" ")[1]
						i_sc_count += 1
					except ValueError:
						pass
				elif line.startswith('Irms ') > 0 :
					try:
						rmsd = line.strip().split(" ")[1]
						rmsd_count += 1
					except ValueError:
						pass
				elif line.startswith('CAPRI_rank') > 0:
					try:
						capri = line.strip().split(" ")[1]
						capri_count += 1
					except ValueError:
						pass	
			if i_sc_count == 0:
				i_sc = 'NaN'
			if rmsd_count == 0:
				rmsd = 'NaN'
			if capri_count == 0:
				capri = 'NaN'
		f1.write('%s\t%s\t%s\t%s\n'%(pdb1,i_sc,rmsd,capri))

df1 = pd.read_csv('%s/%s/isc_irms_score.txt'%(tmp,tpdb),sep='\t')
df2 = pd.merge(df,df1)
df2.to_csv('%s/%s/isc_irms_score_BE.txt'%(tmp,tpdb),sep='\t',index=False)

ori_dic={}
isc_dic={}

with open('%s/ori_env_isc_score.txt'%(pdb_dir),'r') as f2:
	ori_lines= f2.readlines()
	for ori_line in ori_lines:
		anti=ori_line.strip().split('\t')[0]
		score=float(ori_line.strip().split('\t')[1])
		i_score= float(ori_line.strip().split('\t')[2])
		ori_dic[anti] = score
		isc_dic[anti] = i_score

pd.set_option('display.max_rows', None)

df2['I_sc'] = pd.to_numeric(df2['I_sc'], errors='coerce')
df2['Irms'] = pd.to_numeric(df2['Irms'], errors='coerce')

enva_capri= df2[(df2['BE'] < 0) & (df2['Irms'] <20)]
enva_capri= enva_capri.dropna(subset=['CAPRI_rank'])
enva_capri = enva_capri.reset_index(drop=True)
df3 = df2[(df2['I_sc'] < 0) & (df2['Irms'] <20)]
df3 = df3.reset_index(drop=True)

colors = {0: 'black', 1: 'yellow', 2: 'red', 3: 'green'}
plt.scatter(df3['Irms'], df3['I_sc'], c=df3['CAPRI_rank'].apply(lambda x: colors[x]), alpha=0.5)
top10 = df3.sort_values(['I_sc'], ascending=True).head(10)
top10.to_csv('%s/%s/isc_top10.txt' %(tmp,tpdb),sep='\t',index=False)
for i, row in top10.iterrows():
	plt.annotate(row['description'].split('_')[-1], (row['Irms'], row['I_sc']), fontsize=5)
plt.xlabel('Irms')
plt.ylabel('I_sc')
handles = [plt.plot([], [], marker='o', ls="", color=c, alpha=0.5)[0] for c in colors.values()]
labels = list(colors.keys())
plt.legend(handles, labels, title='CAPRI_rank', loc='upper left', fontsize=8)

plt.axhline(isc_dic[anti], color='r', linewidth=1)
plt.savefig('%s/%s/isc.png' %(tmp,tpdb))

plt.clf()
plt.scatter(enva_capri['Irms'], enva_capri['BE'], c=enva_capri['CAPRI_rank'].apply(lambda x: colors[x]), alpha=0.5)
top10_enva = enva_capri.sort_values('BE', ascending=True).head(10)
top10_enva.to_csv('%s/%s/enva_top10.txt' %(tmp,tpdb),sep='\t',index=False)
for i, row in top10.iterrows():
	if row['description'] in enva_capri['description'].values:
		index = enva_capri[enva_capri['description'] == row['description']].index[0]
		plt.annotate(row['description'].split('_')[-1], (enva_capri.loc[index, 'Irms'], enva_capri.loc[index, 'BE']), fontsize=5)

plt.xlabel('Irms')
plt.ylabel('BE')

handles = [plt.plot([], [], marker='o', ls="", color=c, alpha=0.5)[0] for c in colors.values()]
labels = list(colors.keys())
plt.legend(handles, labels, title='CAPRI_rank', loc='upper left', fontsize=8)

plt.axhline(ori_dic[anti], color='r', linewidth=1)
plt.savefig('%s/%s/enva.png' %(tmp,tpdb))

shutil.copytree('%s/%s'%(tmp,tpdb),'%s/%s'%(snugdock_output,tpdb))
