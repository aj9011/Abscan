import os,sys
import pandas as pd
import matplotlib.pyplot as plt
import glob
import shutil

def mk_renumbered_pdb(pdb):
	h_nums = []
	a_nums = []
	with open(pdb,'r') as p_r:
		p_lines=p_r.readlines()
		for p in p_lines:
			if p.startswith('ATOM') > 0 and p[21] =='H':
				if p[22:26].strip() not in h_nums:
					h_nums.append(p[22:26].strip())
			if p.startswith('ATOM') > 0 and p[21] =='A':
				if p[22:26].strip() not in a_nums:
					a_nums.append(p[22:26].strip())

	h_range= range(1,len(h_nums) + 1)
	a_range= range(1,len(a_nums) + 1)
	dic_H= {name:value for name , value in zip(h_nums,h_range) }
	dic_A= {name:value for name , value in zip(a_nums,a_range) }
	old_ch = ''
	with open(pdb.replace('.pdb', '_new_numbering.pdb'),'w') as f2:
		with open(pdb,'r') as p_r:
			p_lines=p_r.readlines()
			ter_index = 0
			for p_line in p_lines:
				if p_line.startswith('ATOM') > 0 and p_line[21]=='H':
					new_line = p_line[:22] + str(dic_H[p_line[22:26].strip()]).rjust(4,' ') + " " + p_line[27:]
					if old_ch != p_line[21] and old_ch != '':
						f2.write('TER\n')
					f2.write(new_line)
					old_ch = p_line[21]
				elif p_line.startswith('ATOM') > 0 and p_line[21]=='A':
					new_line = p_line[:22] + str(dic_A[p_line[22:27].strip()]).rjust(4,' ') + " " + p_line[27:]
					if old_ch != p_line[21] and old_ch != '':
						f2.write('TER\n')
					f2.write(new_line)
					old_ch = p_line[21]
					


#os.chdir('./230407_2bdn_data_without_native')

cd= os.getcwd()
pdb = sys.argv[1]
file_name =sys.argv[2]
antigen_name = file_name.split("_")[1]

ori_isc= 0
ori_enva_isc_result='%s/%s/%s/PDB/ori_env_isc_score.txt' %(cd,antigen_name,file_name)
with open(ori_enva_isc_result,'r') as ori_r:
	ori_lines= ori_r.readlines()
	for ori_line in ori_lines:
		if antigen_name in ori_line:
			ori_isc = float(ori_line.strip().split('\t')[-1])
					
pd.set_option('display.max_rows', None)


snugdock_output = '%s/%s/%s/develop/snugdock_output' %(cd,antigen_name,file_name)
md_dir= '%s/%s/%s/develop/MD_input' %(cd,antigen_name,file_name)

if not os.path.exists(md_dir):
	os.makedirs(md_dir)

isc_pd = pd.read_csv('%s/%s/isc_irms_score_BE.txt'%(snugdock_output,pdb),sep='\t')
top10 = isc_pd.sort_values('I_sc', ascending=True).head(10)
below_ori = isc_pd[isc_pd['I_sc'] < ori_isc]
condition_top10 = top10['CAPRI_rank'] <= 1
condition_below_ori = below_ori['CAPRI_rank'] <= 1
if condition_top10.all() and condition_below_ori.all():
	all_zero_top10 = (top10['CAPRI_rank'] == 0).all()
	all_zero_below_ori = (below_ori['CAPRI_rank'] == 0).all()
	if all_zero_top10 and all_zero_below_ori:
		print ('fail:%s'%(pdb))
	else:
		isc_rms=isc_pd[(isc_pd['Irms'] <6) & (isc_pd['CAPRI_rank'] ==1)]
		isc_sam_temp=isc_rms.loc[isc_rms['I_sc'].idxmin()]
		isc_sam=isc_sam_temp['description']
		mk_renumbered_pdb('%s/%s/dock_res/%s.pdb'%(snugdock_output,pdb,isc_sam))
		with open('%s/%s.pdb'%(md_dir,isc_sam),'w') as ff:
			with open('%s/%s/dock_res/%s_new_numbering.pdb'%(snugdock_output,pdb,isc_sam),'r') as r:
				lines = r.readlines()
				count = 0 
				an_count = 0
				for line in lines:
					if line.startswith('ATOM') > 0 and line[21]=='H':
						if line[77]=='H':
							continue
						ff.write(line)
						count += 1
					elif line.startswith('ATOM') > 0 and line[21]=='A':
						if line[77]=='H':
							continue
						new_line = line[0:4] + str(count+ an_count + 2 ).rjust(7,' ') + line[11:]
						ff.write(new_line)
						an_count += 1
					elif line.startswith('TER') > 0 :
						ff.write(line)
else:
	isc_rms2=isc_pd[(isc_pd['Irms'] <7) & (isc_pd['CAPRI_rank'] >= 2) ]
	isc_sam_temp2=isc_rms2.loc[isc_rms2['I_sc'].idxmin()]
	isc_sam2=isc_sam_temp2['description']
	print('capri:',isc_sam2)
	mk_renumbered_pdb('%s/%s/dock_res/%s.pdb'%(snugdock_output,pdb,isc_sam2))
	with open('%s/%s.pdb'%(md_dir,isc_sam2),'w') as ff:
		with open('%s/%s/dock_res/%s_new_numbering.pdb'%(snugdock_output,pdb,isc_sam2),'r') as r:
			lines = r.readlines()
			count = 0
			an_count = 0
			for line in lines:
				if line.startswith('ATOM') > 0 and line[21]=='H':
					if line[77]=='H':
						continue
					ff.write(line)
					count += 1
				elif line.startswith('ATOM') > 0 and line[21]=='A':
					if line[77]=='H':
						continue
					new_line = line[0:4] + str(count+ an_count + 2 ).rjust(7,' ') + line[11:]
					ff.write(new_line)
					an_count += 1
				elif line.startswith('TER') > 0 :
					ff.write(line)
