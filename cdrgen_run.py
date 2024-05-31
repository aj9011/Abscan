import os,sys
import subprocess as sp
import shutil

file_name=sys.argv[1]

antigen_name = sys.argv[1].split("_")[1]



cd = os.getcwd()
antigen_dir= '%s/%s' %(cd,antigen_name)
antigen_file_dir= '%s/%s/%s' %(cd,antigen_name,file_name)
input_dir= '%s/%s/%s/input' %(cd,antigen_name,file_name)
script = '/KHIT1/ABS_model/script'


info = '%s/input_prepare/%s_pdb_information.txt' %(cd,antigen_name)
epi_file = '%s/input_prepare/%s_epitope' %(cd, antigen_name)



if not  os.path.exists(antigen_dir):
	os.makedirs(antigen_dir)


if not  os.path.exists(antigen_file_dir):
	os.makedirs(antigen_file_dir)

if not  os.path.exists(input_dir):
    os.makedirs(input_dir)



shutil.copy(info, input_dir)
shutil.copy(epi_file, input_dir)


epi_dir= '%s/%s_epitope' %(input_dir,antigen_name)


epi_list=[]


with open(epi_dir,'r') as r:
	epi_lines=r.readlines()
	for epi_line in epi_lines:
		antigen_seq = epi_line.strip().split(":")[0]
		epi_list.append(antigen_seq)


for s in epi_list:
	cdr_gen_line='time python %s/CDR_generator.py  221107_non_H_clean_data.json  %s  %s' %(script,s, file_name)
	sp.call(cdr_gen_line,shell=True)

