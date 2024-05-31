import os
import sys
import shutil
import glob
import pandas as pd
import argparse

RES = [ 'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL','HIE','HID']

def conv(sam,A_res_conv,A_nres_conv,H_res_conv,H_nres_conv,rch,pch,folder):
	os.chdir(folder)
	pdbs = glob.glob('*.pdb.*')
	orn = ''
	rn = ''
	idx = 0
	idx1 = 0
	pdbs.sort()
	for pdb in pdbs:
		ter = 0
		ser = pdb.split('.')[3]
		with open(sam + '-' + ser + '_temp.pdb','w') as ff1:
			with open(pdb,'r') as ff:
				lines = ff.readlines()
				for line in lines:
					if line.startswith('ATOM') > 0 and line[17:20].strip() in RES :
						if orn == '' :
							if ter == 0 :
								TEXT = line[:22] + H_nres_conv[line[22:26]] + line[26:]
							else:
								TEXT = line[:22] + A_nres_conv[line[22:26]] + line[26:]
						elif orn != line[22:26] :
							idx += 1
							if ter == 0 :
								TEXT = line[:22] + H_nres_conv[line[22:26]] + line[26:]
							else:
								TEXT = line[:22] + A_nres_conv[line[22:26]] + line[26:]
						else:
							if ter == 0 :
								TEXT = line[:22] + H_nres_conv[line[22:26]] + line[26:]
							else:
								TEXT = line[:22] + A_nres_conv[line[22:26]] + line[26:]
						if TEXT.find('HIE') > 0 :
							TEXT = TEXT.replace('HIE','HIS')
						elif TEXT.find('HID') > 0 :
							TEXT = TEXT.replace('HID','HIS')
						if ter == 0 :
							TEXT1 = TEXT[:21] + rch + TEXT[22:]
						else:
							TEXT1 = TEXT[:21] + pch + TEXT[22:]
						ff1.write(TEXT1)
						orn = line[22:26]
					elif line.startswith('TER') > 0 :
						ff1.write('TER\n')
						idx = 0
						ter += 1
					elif line.startswith('END') > 0 :
						ff1.write(line)
		os.system('reduce -Trim %s-%s_temp.pdb > ../../dock_res/%s-%s.pdb'% (sam,ser,sam,ser))
		os.remove(sam + '-' + ser + '_temp.pdb')	

GEAR = '/lwork01/neoscan_gear'
MODEL = '/KHIT1/NEO_model/model'
WEIGHT = '/KHIT1/NEO_model/weight'
LIB = '/KHIT1/NEO_model/lib'
SCRIPT = '/KHIT1/NEO_model/script'

wdir = os.getcwd()

parser = argparse.ArgumentParser(description='snapshot generation from MD in AB-ARS')
parser.add_argument('--pdb',dest='pdb',help='Ag-Ab complex')
parser.add_argument('--sample',dest='sample',help='sample_name')
parser.add_argument('--pipeline',dest='pipeline',help='sample_id')
args=parser.parse_args()

pdb = args.pdb
file_name = args.sample
antigen_name = file_name.split("_")[1]
sam1 = args.pipeline


tmp_work = '/lwork01/abs_tmp_work'
md_res = '%s/%s/%s/%s/MD_result'%(wdir,antigen_name,file_name,sam1)
#try:
#	if not os.path.exists(md_res):
#		os.mkdir(md_res,0777)
#except OSError:
#	pass


if not os.path.exists(md_res):
	os.makedirs(md_res)


os.chdir('%s/%s'%(tmp_work,pdb))
#try:
#	if not os.path.exists('dock_res'):
#		os.mkdir('dock_res',0777)
#except OSError:
#	pass



if not os.path.exists('dock_res'):
	os.makedirs('dock_res')


H_res = []
H_nres = []
A_res = []
A_nres = []
with open(pdb + '_H.pdb','r') as rr:
	lines = rr.readlines()
	for line in lines:
		if line.startswith('ATOM') > 0 and line[12:16].strip() == 'CA':
			H_nres.append(line[22:26])
			H_res.append(line[17:20])

with open(pdb + '_A.pdb','r') as ff1:
	lines = ff1.readlines()
	for line in lines:
		if line.startswith('ATOM') > 0 and line[12:16].strip() == 'CA':
			A_nres.append(line[22:26])
			A_res.append(line[17:20])

ipdb = glob.glob('prep/*_initial_solv.pdb')[0]
H_res_conv = {}
H_nres_conv = {}
A_res_conv = {}
A_nres_conv = {}
idx = 0
ter = 0
with open(ipdb,'r') as ff2 :
	lines = ff2.readlines()
	for line in lines:
		if line[17:20] in RES and line[12:16].strip() == 'CA':
			if ter == 0:
				H_nres_conv[line[22:26]] = H_nres[idx]
				H_res_conv[line[17:20]] = H_res[idx]
				idx += 1
			else:
				A_nres_conv[line[22:26]] = A_nres[idx]
				A_res_conv[line[17:20]] = A_res[idx]
				idx += 1
		elif line.startswith('TER') > 0 :
			ter += 1
			idx = 0

os.system('python %s/anal_MD_v1.py -t %s -feat none -nofeature -for-cnn'%(SCRIPT,pdb))
conv(pdb,A_res_conv,A_nres_conv,H_res_conv,H_nres_conv,'H','A','traj_1/pdb_from_prod')
os.chdir(tmp_work)
#print '%s/%s'%(md_res,pdb)
#shutil.copytree(pdb,'%s/%s'%(md_res,pdb))

os.chdir(wdir)
