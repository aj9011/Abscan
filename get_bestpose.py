import os, glob, multiprocessing, sys
import pandas as pd
import subprocess as sp
import shutil

file_name = sys.argv[1]
antigen_name= file_name.split("_")[-1]
md_file = sys.argv[2]
pipeline = sys.argv[3]


wdir = os.getcwd()

md_input = '%s/%s/%s/%s/MD_input' %(wdir,antigen_name,file_name,pipeline)
tmp_work = '/lwork01/abs_tmp_work'
md_res = '%s/%s/%s/%s/MD_result'%(wdir,antigen_name,file_name,pipeline)
GEAR = '/lwork01/neoscan_gear'
#try:
#	if not os.path.exists(md_res):
#		os.mkdir(md_res,0777)
#except OSError:
#	pass

#if not os.path.exists('%s/bestpose'%(md_res)):
#	os.mkdir('%s/bestpose'%(md_res),0777)

if not os.path.exists(md_res):
	os.makedirs(md_res)

if not os.path.exists('%s/bestpose'%(md_res)):
	os.makedirs('%s/bestpose'%(md_res))







os.chdir('%s/%s/dock_res'%(tmp_work,md_file))

pdbs = glob.glob('*.pdb')

tpdbs = []
bes = []
for i in range(5000):
	idx = 0 
	tpdb = md_file + '-' + str(i+1)
	with open('%s_r.pdb'%(tpdb),'w') as f:
		with open('%s.pdb'%(tpdb),'r') as f1:
			lines =f1.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 and line[21]=='H':
					line = line[:21] + 'B' + line[22:]
					f.write(line)
				else:
					f.write(line)
	os.system('%s/enva5.0 -k %s_r.pdb > %s_energy.out'%(GEAR,tpdb,tpdb))
	os.remove('%s_r.pdb'%(tpdb))
	if i >= 3999:
		tpdbs.append(tpdb)		
		if os.path.exists('%s_energy.out'%(tpdb)) and os.path.getsize('%s_energy.out'%(tpdb)) > 0 :
			with open('%s_energy.out'%(tpdb),'r') as ff4:
				lines = ff4.readlines()
				for line in lines:
					if line.startswith('ENERGY') > 0 :
						envs = ' '.join(line.split()).split(' ')
						bes.append(float(envs[1]))
						idx += 1
			if idx == 0 :
				bes.append(0)

df = pd.DataFrame()
df['PDB'] = tpdbs
df['BE'] = bes
be_ave = df['BE'].mean()
df['del'] = abs(float(be_ave) - df['BE'])
df = df.sort_values(['del'],ascending=True)
df = df.reset_index(drop=True)
df.to_csv('%s/%s_8-10ns_confs_BE.txt'%(md_res,md_file),sep='\t',index=False)
shutil.copy('%s.pdb'%(df.at[0,'PDB']),'%s/bestpose'%(md_res))
shutil.copytree('%s/%s'%(tmp_work,md_file),'%s/%s'%(md_res,md_file))
os.chdir(wdir)
