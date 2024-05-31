import os
import sys
import shutil
import pandas as pd

def refine(sam,rpdb):
	RES1 = 'loop_refine_res'
#	try:
#		if not os.path.exists(RES1):
#			os.mkdir(RES1,0777)
#	except OSError:
#		pass
	if not os.path.exists(RES1):
		os.makedirs(RES1)
	loop_st = []
	loop_end = []
	seqs = []
	with open('%s.remodel'%(sam),'r') as f1:
		lines = f1.readlines()
		for line in lines:
			cols = line[:-1].split(' ')
			seqs.append(cols[1])
	for i in range(len(seqs)):
		if i > 0 and seqs[i]=='X' and seqs[i-1]!='X':
			loop_st.append(i+1)
		elif i > 0 and seqs[i]!='X' and seqs[i-1]=='X':
			loop_end.append(i)
	with open('%s.loops'%(sam),'w') as f:
		for i in range(len(loop_st)):
			f.write('LOOP %d %d 0 0 1\n'%(int(loop_st[i]),int(loop_end[i])))
	with open('%s_refine_flags'%(sam),'w') as f1:
		f1.write('-in:file:s %s/%s.pdb\n'%(RES,rpdb))
		f1.write('-in:file:fullatom\n')
		f1.write('-loops:loop_file %s.loops\n'%(sam))
		f1.write('-nstruct 10\n')
		f1.write('-loops:fast\n')
		f1.write('-loops:max_kic_build_attempts 250\n')
		f1.write('-loops:remodel perturb_kic\n')
		f1.write('-loops:refine refine_kic\n')
		f1.write('-ex1\n')
		f1.write('-ex2\n')
		f1.write('-packing:use_input_sc\n')
		f1.write('-out:file:fullatom\n')
		f1.write('-out:path:all %s\n'%(RES1))
		f1.write('-out:file:scorefile %s_refine.sc'%(sam))
	os.system('mpirun -np 20 %s/loopmodel.mpi.linuxgccrelease @%s_refine_flags > %s_refine.log'%(ROSETTA_BIN,sam,sam))
	#os.system('%s/loopmodel.default.linuxgccrelease @%s_refine_flags > %s_refine.log' %(ROSETTA_BIN,sam,sam)) 

	tscore = []
	rpdbs = []
	with open('%s/%s_refine.sc'%(RES1,sam),'r') as f3:
		lines = f3.readlines()
		for line in lines:
			if line.startswith('SCORE') > 0:
				envs = ' '.join(line.split()).split(' ')
				if envs[1]!='total_score':
					tscore.append(float(envs[1]))
					rpdbs.append(envs[len(envs)-1])

	df = pd.DataFrame()
	df['RPDB'] = rpdbs
	df['total_score'] = tscore
	df = df.sort_values(['total_score'],ascending=True)
	df = df[df['total_score']!=0]
	df = df.reset_index(drop=True)
	if df.shape[0]!=0:
		if float(df.at[0,'total_score']) < 0 :
			os.system('%s/reduce -Trim %s/%s.pdb > %s/%s_antigen_filled.pdb'%(GEAR,RES1,df.at[0,'RPDB'],input_dir,sam))		
		else:
			print('Warning:refinemnet not properly done check the result')
			os.system('%s/reduce -Trim %s/%s.pdb > %s/%s_antigen_filled.pdb'%(GEAR,RES1,df.at[0,'RPDB'],input_dir,sam))
	else:
		print('refinement failed')

REMODEL_TOOLS = '/lwork01/rosetta_src_2019.40.60963_bundle/tools/remodel'
ROSETTA_BIN = '/lwork01/rosetta_src_2019.40.60963_bundle/main/source/bin'
GEAR = '/lwork01/neoscan_gear'


file_name=sys.argv[1]
antigen_name = sys.argv[1].split("_")[1]
cd = os.getcwd()

input_dir= '%s/%s/%s/input' %(cd,antigen_name,file_name)


if not  os.path.exists(input_dir):
    os.makedirs(input_dir)


pdb_info = '%s/%s_pdb_information.txt' %(cd,antigen_name)


with open(pdb_info,'r') as pd_r:
	pd_lines= pd_r.readlines()[1:][0]
	sam = pd_lines.split("\t")[0]
		

sam_dir = '%s/input_prepare' %(cd)

os.chdir(sam_dir)


mis_loop = 0
nres_old = 0
ch = []

with open('%s.pdb' %(sam),'r') as f:
	lines = f.readlines()
	for line in lines:
		if line.startswith('ATOM') > 0 and line[12:16].strip()=='CA':
			nres = int(line[22:26].strip())
			ch.append(line[21])
			diff = nres-nres_old
			if diff > 1 and nres_old != 0 :
				mis_loop += 1
			nres_old = nres

ch = list(set(ch))
#print ch
#pdb = sam.split('_')[1]
pdb = sam

if mis_loop > 0 :
	print('%s has missing loop so we will model missing loop'%(sam))
	#try:
	#	if not os.path.exists('loop_remodel_res'):
	#		os.mkdir('loop_remodel_res',0777)
	#except OSError:
	#	pass
	if not os.path.exists('loop_remodel_res'):
		os.makedirs('loop_remodel_res')
	RES = 'loop_remodel_res'
	if not os.path.exists('%s.tsv'%(pdb[0:4])):
		print('curl --silent ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/%s.xml.gz | gunzip | python %s/parse_sifts.py > %s.tsv'%(pdb[0:4],GEAR,pdb[0:4]))
		os.system('curl --silent ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/%s.xml.gz | gunzip | python %s/parse_sifts.py > %s.tsv'%(pdb[0:4],GEAR,pdb[0:4]))
	pdb_seq = ''
	um_seq = ''
	pdb_st = 0
	pdb_end = 0 
	um = ''
	um_st = 0
	um_end = 0 
	um_ser = []
	with open('%s.tsv'%(pdb[0:4]),'r') as f1:
		lines = f1.readlines()
		for line in lines:
			envs = ' '.join(line.split()).split(' ')
			if envs[1]==ch[0]:
				if len(um) < 1:
					um = envs[4]
				if envs[3]!='null':
					if pdb_st == 0 :
						pdb_st = int(envs[3])
						um_st = int(envs[6])
					else:
						pdb_end = int(envs[3])
						um_end = int(envs[6])
	shutil.move('%s.pdb'%(sam),'%s_org.pdb'%(sam))
	with open('%s.pdb'%(sam),'w') as f3:
		with open('%s_org.pdb'%(sam),'r') as f4:
			lines = f4.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 :
					resn = int(line[22:26].strip())
					if resn >= pdb_st and resn <= pdb_end:
						f3.write(line)
				else:
					f3.write(line)
	os.system('%s/getBluePrintFromCoords.pl -pdbfile %s.pdb -chain %s > %s_temp.remodel'%(REMODEL_TOOLS,sam,ch[0],sam))
	with open('%s.tsv'%(pdb[0:4]),'r') as f2:
		lines = f2.readlines()
		for line in lines:
			envs = ' '.join(line.split()).split(' ')
			if envs[1]== ch[0] and int(envs[6]) >= um_st and int(envs[6]) <= um_end:
				if envs[3]=='null':
					envs[2]='-'	
				pdb_seq = pdb_seq + envs[2]
				um_seq = um_seq + envs[5]
				um_ser.append(int(envs[6]))
#	print pdb_st
#	print pdb_end
#	print len(pdb_seq)
#	print len(um_seq)
	mis_loop_st = []
	mis_loop_st_ser = []
	mis_loop_end = []
	mis_loop_end_ser = []
	mis_loop_seq = []
	for i in range(len(pdb_seq)-1):
		if pdb_seq[i+1]=='-' and pdb_seq[i]!='-' and (i+1) < len(pdb_seq):
			mis_loop_st.append(um_ser[i+1])
			mis_loop_st_ser.append(i+1)
		if pdb_seq[i+1]!='-' and pdb_seq[i]=='-' and (i+1) < len(pdb_seq):
			mis_loop_end.append(um_ser[i])
			mis_loop_end_ser.append(i+2)

	os.system('wget https://rest.uniprot.org/uniprotkb/%s.fasta'%(um))
	seq = ''
	with open('%s.fasta'%(um),'r') as f:
		lines = f.readlines()
		for line in lines:
			if line[0]=='>':
				continue
			seq = seq + line[:-1]
	for i in range(len(mis_loop_st)):
		mis_loop_seq.append(seq[mis_loop_st[i]-1:mis_loop_end[i]])
	total_loop = 0
	for i in range(len(mis_loop_st_ser)):
		print('length of missing loop:%d'%(len(mis_loop_seq[i])))
		if len(mis_loop_seq[i]) > 15:
			print('WARNING:Too long missing loop')
			mis_loop_seq[i] = mis_loop_seq[i][0:14]
		if i > 0 :
			total_loop = total_loop + len(mis_loop_seq[i-1])
		#	mis_loop_st_ser[i] = mis_loop_st_ser[i]-len(mis_loop_seq[i-1])
			mis_loop_st_ser[i] = mis_loop_st_ser[i]-total_loop
#	print mis_loop_st
#	print mis_loop_st_ser
#	print mis_loop_seq 
#	print len(mis_loop_seq)
	idx1 = 0	
	with open('%s.remodel'%(sam),'w') as f1:
		with open('%s_temp.remodel'%(sam),'r') as f:
			lines = f.readlines()
			for line in lines:
				cols = line[:-1].split(' ')
				if idx1 < len(mis_loop_seq):
					if int(cols[0]) == mis_loop_st_ser[idx1]:
						f1.write('%d %s L PIKAA %s\n'%(mis_loop_st_ser[idx1],cols[1],cols[1]))
						for j in range(len(mis_loop_seq[idx1])):
							f1.write('0 X L PIKAA %s\n'%(mis_loop_seq[idx1][j]))
		#				idx1 += 1
					elif int(cols[0]) == mis_loop_st_ser[idx1] + 1 :
						f1.write('%s %s L PIKAA %s\n'%(cols[0],cols[1],cols[1]))
						idx1 += 1
					else:
						f1.write('%s %s %s\n'%(cols[0],cols[1],cols[2]))	
				else:
					f1.write('%s %s %s\n'%(cols[0],cols[1],cols[2]))

	with open('%s_flags'%(sam),'w') as f2:
		f2.write('-in:file:s %s_0001.pdb\n'%(sam))
		f2.write('-remodel:blueprint %s.remodel\n'%(sam))
		f2.write('-run:chain %s\n'%(ch[0]))
		f2.write('-remodel:num_trajectory 10\n')
		f2.write('-nstruct 10\n')
	#	f2.write('-remodel:RemodelLoopMover:allowed_closure_attempts 5\n')
		f2.write('-out:path:all %s\n'%(RES))
		f2.write('-out:file:scorefile %s.sc\n'%(sam))
#		f2.write('-chain A\n')

#	os.system('%s/score.opencl.linuxgccrelease -in:file:s %s.pdb -out:output -no_optH false -ignore_zero_occupancy false'%(ROSETTA_BIN,sam))
#	if not os.path.exists('outdir_1'):
#		os.makedirs('outdir_1')
	os.system('%s/score.default.linuxgccrelease -in:file:s %s.pdb -out:output  -no_optH false -ignore_zero_occupancy false'%(ROSETTA_BIN,sam))
#	os.system('mpirun -np 2 %s/score.mpi.linuxgccrelease -in:file:s %s.pdb -out:output -no_optH false -ignore_zero_occupancy false'%(ROSETTA_BIN,sam))
#	shutil.move('%s_0001.pdb' %(sam),'./outdir_1')
	os.system('mpirun -np 20 %s/remodel.mpi.linuxgccrelease @%s_flags > %s.log'%(ROSETTA_BIN,sam,sam))
	
	#os.system('%s/remodel.default.linuxgccrelease @%s_flags > %s.log'%(ROSETTA_BIN,sam,sam))
	if os.path.exists('ROSETTA_CRASH.log') :
		err_idx = 0
		with open('ROSETTA_CRASH.log','r') as ff:
			lines = ff.readlines()
			for line in lines:
				if line.find('the pose does not have residue with chain') > 0:
					err_idx += 1
		if err_idx > 0 :
			with open('%s_flags'%(sam),'a') as f2:	
				f2.write('-chain A\n')
			os.system('mpirun -np 20 %s/remodel.mpi.linuxgccrelease @%s_flags > %s.log'%(ROSETTA_BIN,sam,sam))
			#os.system('%s/remodel.default.linuxgccrelease @%s_flags > %s.log'%(ROSETTA_BIN,sam,sam))
	err = 0 
	tscore = []
	rpdbs = []
	with open('%s/%s.sc'%(RES,sam),'r') as f3:
		lines = f3.readlines()
		for line in lines:
			if line.startswith('SCORE') > 0:
				envs = ' '.join(line.split()).split(' ')
				if envs[1]!='total_score':
					tscore.append(float(envs[1]))
					rpdbs.append(envs[len(envs)-1])
	df = pd.DataFrame()
	df['RPDB'] = rpdbs
	df['total_score'] = tscore
	df = df.sort_values(['total_score'],ascending=True)
	df = df[df['total_score']!=0]
	df = df.reset_index(drop=True)
	if df.shape[0]==0:
		print('Error in loop modeling run again')
		shutil.move(RES,RES+'_error')
		#os.mkdir(RES,0777)
		if not os.path.exists(RES):
			os.makedirs(RES)
		os.system('mpirun -np 20 %s/remodel.mpi.linuxgccrelease @%s_flags > %s.log_v1'%(ROSETTA_BIN,sam,sam))
		#os.system('%s/remodel.default.linuxgccrelease @%s_flags > %s.log_v1'%(ROSETTA_BIN,sam,sam))
		err  =  0
		tscore = []
		rpdbs = []
		with open('%s/%s.sc'%(RES,sam),'r') as f3:
			lines = f3.readlines()
			for line in lines:
				if line.startswith('SCORE') > 0:
					envs = ' '.join(line.split()).split(' ')
					if envs[1]!='total_score':
						tscore.append(float(envs[1]))
						rpdbs.append(envs[len(envs)-1])
		df = pd.DataFrame()
		df['RPDB'] = rpdbs
		df['total_score'] = tscore
		df = df.sort_values(['total_score'],ascending=True)
		df = df[df['total_score']!=0]
		df = df.reset_index(drop=True)
		if df.shape[0]==0:
			print('Still error in loop modeling and check everything')
		else:
			if float(df.at[0,'total_score']) < 0 :
				os.system('%s/reduce -Trim %s/%s.pdb > %s/%s_antigen_fiiled.pdb'%(GEAR,RES,df.at[0,'RPDB'],input_dir,sam))
			elif float(df.at[0,'total_score']) > 0 :
				print('refinement required')
				refine(sam,df.at[0,'RPDB'])
	else:
		if float(df.at[0,'total_score']) < 0 :
			os.system('%s/reduce -Trim %s/%s.pdb > %s/%s_antigen_filled.pdb'%(GEAR,RES,df.at[0,'RPDB'],input_dir,sam))
		elif float(df.at[0,'total_score']) > 0 :
			print('refinement required')
			refine(sam,df.at[0,'RPDB'])
	os.remove('%s.tsv'%(pdb[0:4]))
	os.remove('%s.fasta'%(um))
	os.remove('%s_temp.remodel'%(sam))
