import os,sys
import glob
import subprocess as sp
import shutil


def change_cdrh3_num(parsed_pdb):
	antigens_list= []
	antigen_num = 0
	cdr_num=[]
	with open(parsed_pdb, 'r') as an:
		antigen_nums= an.readlines()
		for i, aa in enumerate(antigen_nums):
			if aa.startswith('ATOM') > 0 :
				if 'A' == aa[21]: 
					if i == len(antigen_nums)-3:
						antigen_num= int(aa[6:11].strip())
						antigens_list.append(aa.rstrip())
					else:
						antigens_list.append(aa.rstrip())
				elif 'B' ==aa[21]:
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
				if 'B' == cdr_line[21]:
					new_numbering = dic2[cdr_line[22:27].strip()]
					new_lines= cdr_line[:6] + str(i+antigen_num+1).rjust(5,' ') + cdr_line[11:22] + str(new_numbering).rjust(4,' ') + " " + cdr_line[27:].rstrip()
					new_lines_list.append(new_lines)
				else: 
					new_lines_list.append(cdr_line.strip())
	cdrh3_result = '%s' %(parsed_pdb.replace('_parsed.pdb','_cdrh3_numbering.pdb'))
	with open(cdrh3_result,'w')  as w_ag_cdr_r:
		w_ag_cdr_r.write('\n'.join(new_lines_list))		










#cd= os.getcwd()
#antigen=sys.argv[1]
file_name= sys.argv[1]
antigen_name= file_name.split("_")[1]
GEAR = '/lwork01/neoscan_gear'

cd= os.getcwd()


ori_dir = '%s/%s/%s/PDB' %(cd,antigen_name,file_name)

isc_score_pdb = [os.path.join(ori_dir,x) for x in os.listdir(ori_dir) if x.endswith('_0001.pdb')][0]
isc_score=0
with open(isc_score_pdb,'r') as isc_r:
	isc_pdbs = isc_r.readlines()
	for isc in isc_pdbs:
		if 'dG_separated ' in isc:
			isc_score = float(isc.split(" ")[1])





ori_pdb = [os.path.join(ori_dir,x) for x in os.listdir(ori_dir) if x.endswith('HLA.pdb')][0]

ss=''
antigen_list=[]
cdr_list=[]
with open(ori_pdb,'r') as rc:
	cc= rc.readlines()
	for co in cc:
		if co.startswith('ATOM') > 0:
			if 'A' == co[21]:
				antigen_list.append(co.strip())
			elif 'H' == co[21] and (int(co[23:26].strip()) >=95 and int(co[23:26].strip()) <= 102):
				cdr_new= co[:21] + 'B' + co[22:]
				cdr_list.append(cdr_new.strip())
final_line= '\n'.join(antigen_list) + '\nTER\n' + '\n'.join(cdr_list)
#print(final_line)
with open(ori_pdb.replace(".pdb","_parsed.pdb"),'w') as w:
	w.write(final_line)
parsed= ori_pdb.replace(".pdb","_parsed.pdb")
change_cdrh3_num(parsed)
ch_cdr=  ori_pdb.replace(".pdb","_cdrh3_numbering.pdb")
only_parsed= ch_cdr.split("/")[-1]
env= ori_pdb.replace(".pdb","_cdrh3_numbering.env")
		
if not os.path.exists(only_parsed):
	shutil.copy(ch_cdr,only_parsed)
os.system('%s/enva5.0 -k %s > %s' %(GEAR,only_parsed,env))
with open(env,'r') as envr:
	envs=envr.readlines()
	for en in envs:
		#print(en)
		score= float(en.strip().split(" ")[2])
		#print(score)
		env_isc_lin= '%s\t%f\t%f' %(antigen_name ,score,isc_score )
		ss = env_isc_lin
		#print(score)
os.remove(only_parsed)


with open('%s/ori_env_isc_score.txt' %(ori_dir),'w') as w2:
	w2.write(ss)


