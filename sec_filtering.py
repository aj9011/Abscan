import os,sys
import glob
import shutil


def count_deep_sec_num(name):
	sec_count=0
	with open(name,'r') as r_envs:
		envs=r_envs.readlines()[1:]
		env_length=len(envs)
		for e in envs:
			if e[21] ==' ':
				pass
			else:
				sec_count +=1
	return env_length, sec_count

def count_ori_sec_num(names,length):
	ori_sec_count=0
	new_lines=[]
	with open(names,'r') as r_envs:
		envs=r_envs.readlines()[1:]
		env_length=len(envs)
		for i, e in enumerate(envs):
			if i < length:
				new_lines.append(e.strip())
				if e[21] ==' ':
					pass
				else:
					ori_sec_count +=1
	with open(names+"2",'w') as w:
		w.write("\n".join(new_lines))
	return ori_sec_count



def get_sec_filtered_file(deepab_info):
	rmsd_sec_result= []
	enva_sec_result= []
	with open(deepab_info,'r') as dp_r:
		dp_lines = dp_r.readlines()
		for dp in dp_lines:
			sp = dp.strip().split("\t")
			deep_dir = sp[0]
			flag= sp[-1]
			if flag =='pass':
				rmsd_result=  '%s/%s_final_rmsd_stats.txt' %(deepab_dir,antigen_name)
				with open(rmsd_result,'r') as fr:
					fr_lines = fr.readlines()
					for fr_line in fr_lines:
						cdr= fr_line.strip().split("\t")[1]
						if cdr in deep_dir:
							rmsd_sec_result.append(fr_line.strip())

				enva_result= '%s/%s_final_enva_stats.txt' %(deepab_dir,antigen_name)
				with open(enva_result,'r') as fe:
					fe_lines = fe.readlines()
					for fe_line in fe_lines:
						cdr = fe_line.strip().split("\t")[1]
						if cdr in deep_dir:
							enva_sec_result.append(fe_line.strip())
	sec_rmsd_info = '%s/%s_sec_final_rmsd_stats.txt' %(deepab_dir,antigen_name)
	sec_env_info = '%s/%s_sec_final_enva_stats.txt' %(deepab_dir,antigen_name)			
	with open(sec_rmsd_info,'w') as rw:
		rw.write('\n'.join(rmsd_sec_result))
	with open(sec_env_info,'w') as ew:
		ew.write('\n'.join(enva_sec_result))

def mk_renumbered_pdb(pdb,develop_dir):
	h_nums = []
	l_nums = []
	with open(pdb,'r') as p_r:
		p_lines=p_r.readlines()
		for p in p_lines:
			if p.startswith('ATOM') > 0 and p[12:16].strip()=='CA' and p[21]=='H':
				h_nums.append(p[22:27].strip())
			if p.startswith('ATOM') > 0 and p[12:16].strip()=='CA' and p[21]=='L':
				l_nums.append(p[22:27].strip())

	h_range= range(1,len(h_nums) + 1)
	l_range= range(1,len(l_nums) + 1)
	dic_H= {name:value for name , value in zip(h_nums,h_range) }
	dic_L= {name:value for name , value in zip(l_nums,l_range) }
	new_lines=[]
	with open(pdb,'r') as p_r:
		p_lines=p_r.readlines()
		for p in p_lines:
			if p.startswith('ATOM') > 0 and p[21]=='H':
				new_line = p[:22] + str(dic_H[p[22:27].strip()]).rjust(4,' ') + " " + p[27:].rstrip()
				new_lines.append(new_line)
			elif p.startswith('ATOM') > 0 and p[21]=='L':
				new_line = p[:22] + str(dic_L[p[22:27].strip()]).rjust(4,' ') + " " + p[27:].rstrip()
				new_lines.append(new_line)	

	pdb = os.path.basename(pdb)
	with open('%s/%s' %(develop_dir,pdb.replace('_tr.pdb', '_new_numbering.pdb')),'w') as wn:
		wn.write('\n'.join(new_lines))
					
	 
file_name = sys.argv[1]
antigen_name = file_name.split("_")[1]


cd = os.getcwd()
GEAR = '/lwork01/neoscan_gear'
deepab_dir= '%s/%s/%s/DeepAb/result'%(cd,antigen_name,file_name)
pdb_dir = '%s/%s/%s/PDB'%(cd,antigen_name,file_name)
develop_dir = '%s/%s/%s/develop'%(cd,antigen_name,file_name)
pyig_result = '%s/pyig_result' %(develop_dir)
tap_result = '%s/tap_result' %(develop_dir)

sec_dir= '%s/secondary' %(deepab_dir)
os.makedirs(sec_dir , exist_ok=True)
os.makedirs(develop_dir,exist_ok=True)
os.makedirs(pyig_result,exist_ok=True)
os.makedirs(tap_result,exist_ok=True)
pdb_info= '%s/%s/%s/input/%s_pdb_information.txt' %(cd,antigen_name,file_name,antigen_name)
ben_dic={}
with open(pdb_info, 'r' ) as ben_r:
	bens= ben_r.readlines()[1:]
	for ben in bens:
		ben_sp=ben.strip().split("\t")
		ben_dic[ben_sp[0]]= ben_sp[1]+"_"+ben_sp[3]+"_"+ben_sp[4]+"_"+ben_sp[5]


#print(ben_dic)


for k, v in ben_dic.items():
	#print(v)
	ori_pdb_name = '%s/%s/%s/PDB/%s_HL.pdb' %(cd,antigen_name,file_name,k)
	env_chain= v.split("_")[1]
	env_chain_L= v.split("_")[2]
	os.system('%s/enva5.0 -s %s %s > %s/%s.sec_env'%(GEAR,ori_pdb_name,env_chain,pdb_dir,k+'_HL'))
	os.system('%s/enva5.0 -s %s %s > %s/%s.sec_envL'%(GEAR,ori_pdb_name,env_chain_L,pdb_dir,k+'_HL'))
	antigen_file_name= v.split("_")[0]

	deepabs=glob.glob(os.path.join(deepab_dir, '*deepab_tr.pdb'))
	sec_num = 0
	sec_numL= 0
	passed_num=0
	fail_num=0
	deep_summary=[]
	deep_num= 0
	for deepab in deepabs:
		deep_num=len(deepabs)
		deepab = os.path.basename(deepab)
#		only_deep = deep.split("/")[-1]
#		if not os.path.exists('%s/%s' %(ori_dir, only_deep)): 
#			shutil.copy(deep,ori_dir)
		os.system('%s/enva5.0 -s %s/%s H > %s/%s.sec_env'%(GEAR,deepab_dir,deepab,deepab_dir,deepab.replace('_deepab_tr.pdb','')))
		os.system('%s/enva5.0 -s %s/%s L > %s/%s.sec_envL'%(GEAR,deepab_dir,deepab,deepab_dir,deepab.replace('_deepab_tr.pdb','')))
		deep_enva_name='%s/%s.sec_env' %(deepab_dir,deepab.replace('_deepab_tr.pdb',''))
		deep_enva_length, deep_enva_count = count_deep_sec_num(deep_enva_name)

		#deep_enva_nameL= deep.split("/")[-1].replace(".deepab.pdb.pdb",".sec_envL")
		deep_enva_nameL='%s/%s.sec_envL' %(deepab_dir,deepab.replace('_deepab_tr.pdb',''))
		deep_enva_lengthL, deep_enva_countL = count_deep_sec_num(deep_enva_nameL)

		#print(deep_enva_count,sec_num)
		if sec_num == 0:
			sec_num = count_ori_sec_num('%s/%s.sec_env' %(pdb_dir,k+'_HL'),deep_enva_length)
		if sec_numL ==0:
			sec_numL= count_ori_sec_num('%s/%s.sec_envL' %(pdb_dir,k+'_HL'),deep_enva_lengthL)


		#if deep_enva_count/sec_num >= 0.7 and deep_enva_countL/sec_numL >= 0.5  : #for_production
		if deep_enva_count/sec_num >= 0.7 and deep_enva_countL/sec_numL >= 0.5  :  #for_test
			passed_num +=1
			shutil.copy('%s/%s'%(deepab_dir,deepab), sec_dir)
			mk_renumbered_pdb('%s/%s'%(deepab_dir,deepab),develop_dir)
			deep_line= '%s\t%d\t%d\t%d\tpass' %(deepab, deep_enva_count, sec_num, deep_enva_length)
			deep_summary.append(deep_line)
		else:
			fail_num +=1
			deep_line= '%s\t%d\t%d\t%d\tfail' %(deepab, deep_enva_count, sec_num, deep_enva_length)
			deep_summary.append(deep_line)

	#	os.remove(deep_enva_name)
	#	os.remove(deep_enva_nameL)

	deep_info = '%s/%s_deepab_info.txt' %(deepab_dir,antigen_file_name)
	total_line= "%d\t%d\t%d" %(deep_num, passed_num,fail_num)
	with open('%s/%s_total_sec_result.txt' %(deepab_dir,antigen_file_name),'w' ) as t:
		t.write(total_line)
	with open(deep_info,'w') as y:
		y.write('\n'.join(deep_summary))

	

	get_sec_filtered_file(deep_info)
