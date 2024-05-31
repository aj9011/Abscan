import os
import sys
import shutil
import glob
import pandas as pd
import argparse

GEAR = '/lwork01/neoscan_gear'
ABS_MODEL = '/KHIT1/ABS_model/model'
SCRIPT = '/KHIT1/ABS_model/script'

wdir = os.getcwd()

parser = argparse.ArgumentParser(description='MD run in AB-ARS')
parser.add_argument('--pdb',dest='pdb',help='Ag-Ab complex')
parser.add_argument('--sample',dest='sample',help='sample_id')
parser.add_argument('--pipeline',dest='pipeline',help='sample_id')
parser.add_argument('--gpu-id', dest='fixed_gpu_id', type=str, help='Fix gpu device id (device id: 0,1,..), default: None', default=0)
args=parser.parse_args()

sam = args.pdb
file_name = args.sample
sam1 = args.pipeline
sam2 = args.fixed_gpu_id
antigen_name = file_name.split("_")[1]

#print type(sam1)

md_input = os.path.join(wdir,antigen_name,file_name,sam1,'MD_input')
#if sam1 == 'develop':
#	tpdb = glob.glob('%s.LHA_sorted_0001_*.pdb'%(sam))[0]
#elif sam1 == 'binding':
#	tpdb = sam		


if sam1 == 'develop':
	tpdb = sam
elif sam1 == 'binding':
	tpdb = glob.glob('%s.LHA_sorted_0001_*.pdb'%(sam))[0]
	


tmp_work = '/lwork01/abs_tmp_work'
tmp_sam_dir = '%s/%s' %(tmp_work,sam)
#try:
#	if not os.path.exists(tmp_work):
#		os.mkdir(tmp_work,0777)
#except OSError:
#	pass

#try:
#	if not os.path.exists('%s/%s'%(tmp_work,sam)):
#		os.mkdir('%s/%s'%(tmp_work,sam),0777)
#except OSError:
#	pass

if not os.path.exists(tmp_work):
	os.makedirs(tmp_work)
if not os.path.exists(tmp_sam_dir):
	os.makedirs(tmp_sam_dir)


shutil.copy('%s/%s.pdb'%(md_input,tpdb),'%s/%s/%s.pdb'%(tmp_work,sam,sam))
os.chdir('%s/%s'%(tmp_work,sam))

with open(sam + '_H.pdb','w') as ff:
	with open(sam + '.pdb','r') as rr:
		lines = rr.readlines()
		for line in lines:
			if line.startswith('ATOM') > 0 and line[21]=='H' and line[77]!='H':
				ff.write(line)
	ff.write('TER\n')

with open(sam + '_A.pdb','w') as ff1:
	with open(sam + '.pdb','r') as rr:
		lines = rr.readlines()
		for line in lines:
			 if line.startswith('ATOM') > 0 and line[21]=='A' and line[77]!='H':
				ff1.write(line)
	ff1.write('END\n')


#os.system('python %s/run_MD_peptide.py -ir %s.pdb -il %s.pdb -t %s -n 1 -fix-gpu-id %s -heat-time 0.02 -equil-time 0.05 -prod-total-time 1 -prod-freq-time 0.1 -wrt-pd-crd-freq 0.002 -wrt-pd-out-freq 0.002'%(SCRIPT,tpdb+'_H',tpdb+'_A',tpdb,sam2)) # for test
os.system('python %s/run_MD_peptide.py -ir %s.pdb -il %s.pdb -t %s -n 1 -fix-gpu-id %s -heat-time 0.02 -equil-time 0.5 -prod-total-time 10 -prod-freq-time 1 -wrt-pd-crd-freq 0.002 -wrt-pd-out-freq 0.002'%(SCRIPT,sam+'_H',sam+'_A',sam,sam2)) # for production

os.chdir(wdir)
