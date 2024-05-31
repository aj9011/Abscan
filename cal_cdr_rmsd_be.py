import os,sys
import subprocess as sp
import glob
import shutil

def sheba_run(ref,org) :
	cmd3 = GEAR + '/sheba_01 -x ' + ref + '.pdb ' + org + '.pdb'
	sp.call(cmd3,shell=True)
	cmd4 = GEAR + '/sheba_01 -t ' + org + '.trf ' + org + '.pdb'
	sp.call(cmd4,shell=True)
	shutil.move(org + '.pdb.pdb',org +'_tr.pdb')

def get_cdr(x):

	cdr_lines=[]
	with open('%s' %(x),'r') as pdb:
		pd=pdb.readlines()
		start_line=0
		end_line= 0
		cdr_region=[]
		for i ,p_c in enumerate(pd):
			if p_c[0:4] == 'ATOM' :
				if p_c[21] == 'H':
					if int(p_c[22:26].strip()) <= 94:
						start_line = i+1
					elif int(p_c[22:26].strip()) == 103:
						if end_line == 0:
							end_line = i
		cdr_dic={}
		cdr_seq=''
		for v2 ,p_c2 in enumerate(pd):
			if v2 >=start_line and v2 <end_line:
				cdr_region.append(p_c2[22:27].strip())
				cdr_dic[p_c2[22:27].strip()] = one_letter[p_c2[16:20].strip()]
		print(cdr_dic)
		for key, value in cdr_dic.items():
			cdr_seq += value
			cdr_regions=list(set(cdr_region))
			cdr_regions.sort()
			cdr_regions.sort(key=lambda x : int(x) if len(x) <=2 else (int(x[:-1]) if x[0] == '9'  else (int(x) if len(x) ==3 else int(x[:-1]) )))

		
		for v3, p_c3 in enumerate(pd):
			if p_c3[0:4] =='ATOM':
				if p_c3[21] =='H':
					if p_c3[22:26].strip() in cdr_regions:
						cdr_lines.append(p_c3.strip())


		cdrh3_lines = '\n'.join(cdr_lines)

		if x.endswith('_HL.pdb'):
			with open('%s' %(x.replace('_HL.pdb','_HL_H3.pdb')),'w') as w:
				w.write(cdrh3_lines)
		elif x.endswith('deepab_tr.pdb'):
			with open('%s' %(x.replace('.deepab_tr.pdb','_H3.pdb')),'w') as w:
				w.write(cdrh3_lines)

							

def rmsd_cal(ori, deep):
	#rmsd_sh = '/STG24-3/abscan/abscan_github/rmsd_calculator'

	
	rmsd_cmd = '%s/rmsd_total_v2 %s %s bb H > %s/%s_%s_rmsd.txt' %(GEAR, ori, deep, rmsd_result_dir, deep.split("/")[-1].split("_")[0], deep.split("/")[-1].split("_")[1])
#	print (rmsd_cmd)
	sp.call(rmsd_cmd, shell = True)




def prepare_mini(deep_h3):
	line_clus=[]
	with open(deep_h3,'r') as deep_h3r:
		deep_h3_lines= deep_h3r.readlines()
		for deep_h3_l in deep_h3_lines:
			new_deep_h3 = deep_h3_l[:21] + 'B' + deep_h3_l[22:].rstrip()
			line_clus.append(new_deep_h3)
	deep_h3_new_lines = '\n'.join(line_clus)
	
	with open('%s' %(deep_h3.replace("_H3.pdb","_CDRH3.pdb")),'w') as  h3_w:
		h3_w.write(deep_h3_new_lines)

def merge_ag_cdrh3(ag,cdrh3):
	antigens_list= []
	with open(ag, 'r') as an:
		antigen_nums= an.readlines()
		antigen_num=0
		for i, aa in enumerate(antigen_nums):
			if aa[0:4] == 'ATOM':
				if i == len(antigen_nums)-3:
					antigen_num= int(aa[6:11].strip())
					antigens_list.append(aa.rstrip())
				else:
					antigens_list.append(aa.rstrip())

		cdr_num=[]
		with open(cdrh3, 'r' ) as sam_cdr_r:
			sam_cdr= sam_cdr_r.readlines()
			for sam in sam_cdr:
				cdr_num.append(sam[22:27].strip())
		cdr_region = sorted(list(set(cdr_num)))
		cdr_regions=list(set(cdr_region))
		cdr_regions.sort()
		cdr_regions.sort(key=lambda x : int(x) if len(x) <=2 else (int(x[:-1]) if x[0] == '9'  else (int(x) if len(x) ==3 else int(x[:-1]) )))
		dic2={}
		for nn, n in enumerate(cdr_regions):
			dic2[n] = nn+1
		new_lines_list=[]
		with open(cdrh3,'r' ) as cdr_r:
			cdr_lines=cdr_r.readlines()
			for i, cdr_line in enumerate(cdr_lines):
				new_numbering = dic2[cdr_line[22:27].strip()]
				new_lines= cdr_line[:6] + str(i+antigen_num+1).rjust(5,' ') + cdr_line[11:22] + str(new_numbering).rjust(4,' ') + " " + cdr_line[27:].rstrip()
				new_lines_list.append(new_lines)
		new_lines_list.insert(0,'TER')
		final_list= antigens_list + new_lines_list
		ddir = '/'.join(cdrh3.split("/")[0:-1])
		with open('%s/final_%s_%s.pdb' %(ddir,  ag.split("/")[-1].replace('_antigen_filled.pdb','') , cdrh3.split("/")[-1].replace('_CDRH3.pdb','')),'w') as w_ag_cdr_r:
			w_ag_cdr_r.write('\n'.join(final_list))





one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',\
			'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',\
			'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',\
			'GLY':'G', 'PRO':'P', 'CYS':'C'}


if __name__ == "__main__":

	##rmsd & enva##

	###prepare rmsd : have to get original heavy light chain of antibody structure (ex:5xxy_HL.pdb)
	##we calculate only cdrh3 after superimpose candidates with original antibody structure
	file_name = sys.argv[1]
	antigen_name= file_name.split("_")[1]
	c_dir=os.getcwd()
	input_dir = '%s/%s/%s/input/' %(c_dir,antigen_name,file_name)
	pdb  =  '%s/%s/%s/PDB' %(c_dir,antigen_name,file_name)
	ori_dir = [os.path.join(pdb,x) for x in os.listdir(pdb) if x.endswith('_HL.pdb')][0]
	ori_pdb = ori_dir.split("/")[-1]
	deepab_result_dir = '%s/%s/%s/DeepAb/result'  %(c_dir, antigen_name, file_name)
	env_dir= '%s/env/'%(deepab_result_dir)
	os.makedirs(env_dir , exist_ok=True)

	GEAR = '/lwork01/neoscan_gear'
	script = '/KHIT1/ABS_model/script'
#	tool_dir='/tools2/abscan_tools' 
#	sheba_dir = '%s/sheba_01_20200602' %(tool_dir)
#	rmsd_dir = '%s/rmsd_calculator' %(tool_dir)
#	enva_dir ='%s/enva4.5w' %(tool_dir)


	rmsd_result_dir = '%s/rmsd' %(deepab_result_dir)

	if not os.path.exists(rmsd_result_dir):
		os.makedirs(rmsd_result_dir)
	
#	deepab_list= [os.path.join(deepab_result_dir,x) for x in os.listdir(deepab_result_dir) if x.endswith('deepab.pdb')]
	deepab_list = [os.path.join(root, x) for root, _, files in os.walk(deepab_result_dir) for x in files if x.endswith('deepab.pdb')]
	print(deepab_list)


	ori_copy = '%s/%s' %(c_dir,ori_pdb)
	if not os.path.exists(ori_copy):
		shutil.copy(ori_dir, ori_copy)

	for deepab in deepab_list:
		deep_pdb = deepab.split("/")[-1]
		deep_copy = '%s/%s' %(c_dir, deep_pdb)
		if not os.path.exists(deep_copy):
			shutil.copy(deepab, deep_copy)
			
		sheba_run(ori_pdb.replace('.pdb','') , deep_pdb.replace('.pdb',''))

		trf_file = deep_pdb.replace(".pdb",".trf")
		trans_deepab= deep_pdb.replace(".pdb","_tr.pdb")
		shutil.copy(trans_deepab, deepab_result_dir)

		os.remove(trans_deepab)
		os.remove(trf_file)
		os.remove(deep_copy)
		
	os.remove(ori_copy)
	

	super_deepabs= [os.path.join(deepab_result_dir,x) for x in os.listdir(deepab_result_dir) if x.endswith('deepab_tr.pdb')]
	print(super_deepabs)
	get_cdr(ori_dir)

	for super_deep in super_deepabs:
		print("ssssssssss",super_deep)
		get_cdr(super_deep)
	

	ori_cdrh3_dir= ori_dir.replace('_HL.pdb','_HL_H3.pdb')

	deepab_h3s= [os.path.join(deepab_result_dir,x) for x in os.listdir(deepab_result_dir) if x.endswith('_H3.pdb')]

	print(ori_cdrh3_dir,deepab_h3s)

	for deepab_h3 in deepab_h3s:
		rmsd_cal(ori_cdrh3_dir, deepab_h3)
		prepare_mini(deepab_h3)

	

	rmsd_list= [os.path.join(rmsd_result_dir,x) for x in os.listdir(rmsd_result_dir) if  x.endswith('_rmsd.txt') and not x.endswith('_total_rmsd.txt')]

	rmsd_lines= []
	
	for rmsd in rmsd_list:
		with open(rmsd,'r') as rmsd_r:
			rmsd_line = rmsd_r.readlines()[0].strip()
			rmsd_lines.append(rmsd_line)

	total_rmsd_file= '%s/%s_total_rmsd.txt' %(rmsd_result_dir, antigen_name)

	total_rmsd_lines= '\n'.join(rmsd_lines)

	with open(total_rmsd_file,'w') as rmsd_w:
		rmsd_w.write(total_rmsd_lines)


			

	chagned_cdr_chs= [os.path.join(deepab_result_dir,x) for x in os.listdir(deepab_result_dir) if x.endswith('_CDRH3.pdb')]

	enva_antigen= '%s/%s_antigen_filled.pdb' %(input_dir, ori_pdb.replace('_HL.pdb',''))

	for chagned_cdr_ch in chagned_cdr_chs:
		merge_ag_cdrh3(enva_antigen,chagned_cdr_ch)


	run_minimize_cmd = 'python %s/minimize_pep_io.py %s' %(script, file_name)
	sp.call(run_minimize_cmd , shell =True)

	
	mini_patt= '**/final_*_0001.pdb'
	mini_pdbs= [os.path.abspath(x) for x in glob.glob(os.path.join(deepab_result_dir, mini_patt))]
	enva_result_file = '%s/%s_enva_result.txt' %(env_dir,antigen_name)
#	print ('*******')
#	print (mini_pdbs)

	energy_line= []
	for mini_pdb in mini_pdbs:
		mini_pdb = mini_pdb.replace('.pdb','')
		only_pdb = os.path.basename(mini_pdb)
		os.system('%s/enva5.0 -K %s.pdb > %s_energy.out'%(GEAR,mini_pdb,mini_pdb))
		if os.path.getsize('%s_energy.out'%(mini_pdb))!=0:
			with open('%s_energy.out'%(mini_pdb),'r') as r:
				lines = r.readlines()
				for line in lines:
					if line.startswith('ENERGY'):
						line_list=line.strip().split(" ")
						line_l= [i for i in line_list if i]
						energy_line.append(only_pdb +'\t' + line_l[1])
		else:
			energy_line.append(only_pdb +'\t0.000')
	
	with open('%s/%s_enva_result.txt' %(env_dir,antigen_name),'w') as w:
		w.write('\n'.join(energy_line))

	dic={}
	with open(total_rmsd_file, 'r') as r:
		lines = r.readlines()
		for line in lines:
			h3_pdb=line.strip().split('\t')[0].split('/')[-1]
			rmsd= line.strip().split('\t')[-1]
			dic[h3_pdb] = rmsd
	enva_dic={}
	with open(enva_result_file,'r' ) as er:
		e_lines = er.readlines()
		for e_l in e_lines:
			enva_h3_pdb = e_l.strip().split('\t')[0].split('/')[-1]
			enva_energy= e_l.strip().split("\t")[-1]
			enva_dic[enva_h3_pdb] = enva_energy
	
	print (enva_dic)

	for_epi_dir='%s/%s/%s/cdrgen_result' %(c_dir,antigen_name,file_name)

	set_list=[os.path.join(for_epi_dir,x) for x in os.listdir(for_epi_dir) if x.endswith('_epi_cdr_set.txt')]

	final_list=[]
	final_enva_list=[]
	for  se in set_list:
		with open(se,'r') as s_r:
			see=s_r.readlines()
			for s in see:
				epi= s.strip().split("\t")[0]
				cdr=s.strip().split("\t")[1]
				ind=  '%s_%s_H3.pdb' %(antigen_name,cdr)
			#	print (ind)
				enva_ind=  'final_%s_%s_%s_0001' %(ori_pdb.replace('_HL.pdb',''), antigen_name, cdr)
			#	print (enva_ind)
				if ind in dic :
					new_line= '%s\t%s\t%s' %(epi,cdr,dic[ind])
					final_list.append(new_line)
				if enva_ind in enva_dic:
					enva_new_line = '%s\t%s\t%s' %(epi,cdr,enva_dic[enva_ind])
					final_enva_list.append(enva_new_line)

	ww= "\n".join(final_list)
	e_w= "\n".join(final_enva_list)
	final_rmsd='%s/%s_final_rmsd_stats.txt' %(deepab_result_dir,antigen_name)
	final_enva='%s/%s_final_enva_stats.txt' %(deepab_result_dir,antigen_name)
	with open(final_rmsd,'w') as w:
		w.write(ww)
	with open(final_enva,'w') as ew:
		ew.write(e_w)

