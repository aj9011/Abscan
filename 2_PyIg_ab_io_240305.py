import sys, os, json
import shutil
import subprocess as sp
import math


def clean (pdb):
	clean_pdb ='%s/clean_pdb.py' %(rosetta_dir)
	cmd_clean_H = 'python2 %s %s H' %(clean_pdb, pdb)
	sp.call(cmd_clean_H, shell = True)
	cmd_clean_L = 'python2 %s %s L' %(clean_pdb, pdb)
	sp.call(cmd_clean_L, shell = True)

def ANARCI (fasta, numbering, chain) :
	cmd1= 'ANARCI --scheme %s -i %s -o %s.anrc' %(numbering,fasta,fasta)
	sp.call(cmd1,shell=True)
	res_dic_H={'H1':'','H2':'','H3':''}
	res_dic_L={'L1':'','L2':'','L3':''}
	for line in open(fasta+'.anrc').readlines():
		if line[0] == 'H':
			pos=int(line[2:8].strip())
			AA=line[10]
			if pos in H1_pos:
				res_dic_H['H1']+=AA
			elif pos in H2_pos:
				res_dic_H['H2']+=AA
			elif pos in H3_pos:
				res_dic_H['H3']+=AA
		elif line[0] == 'L':
			pos=int(line[2:8].strip())
			AA=line[10]
			if pos in L1_pos:
				res_dic_L['L1']+=AA
			elif pos in L2_pos:
				res_dic_L['L2']+=AA
			elif pos in L3_pos:
				res_dic_L['L3']+=AA

	if chain == 'H':
		keys=sorted(res_dic_H.keys())
		res=[]
		for key in keys:
			res.append(res_dic_H[key].replace('-','').strip())
		rm = 'rm ' + fasta+'.anrc'
		sp.call(rm, shell = True)
		return( '\t'.join(res))
	elif chain == 'L':
		keys=sorted(res_dic_L.keys())
#		keys.sort()
		res=[]
		for key in keys:
			res.append(res_dic_L[key].replace('-','').strip())
		rm = 'rm ' + fasta+'.anrc'
		sp.call(rm, shell = True)
		return( '\t'.join(res))

def distance (q_list, s_list, angle_kind) :
	sum_d = 0
	for n, angle in enumerate(q_list[1:]) :
		d = (2 * (1 - math.cos(float(q_list[n]) - float(s_list[n]))))
		sum_d += d
                #print(d)
#	print(sum_d)
	return sum_d


def sheba_run(ref,org,ref_path,org_path) :
	if org_path != './':
		shutil.copy('%s/%s.pdb'%(org_path,org),'./%s.pdb'%(org))
	if ref_path != './':
		shutil.copy('%s/%s.pdb'%(ref_path,ref),'./%s.pdb'%(ref))
	cmd3 = GEAR + '/sheba_01 -x ' + ref + '.pdb ' + org + '.pdb'
	sp.call(cmd3,shell=True)
	cmd4 = GEAR + '/sheba_01 -t ' + org + '.trf ' + org + '.pdb'
	sp.call(cmd4,shell=True)
	shutil.move(org + '.pdb.pdb',ref+'_'+ org+'.pdb')
	os.remove(org+'.trf')
#	os.remove(org +'.pdb')

def rmsd (CDR, DB_CDR, kind) :
	chain = kind[0]
	clean_cmd1 = 'python2 %s/clean_pdb.py %s %s' %(rosetta_dir,CDR,chain)
	clean_cmd2 = 'python2 %s/clean_pdb.py %s %s' %(rosetta_dir,DB_CDR,chain)
	sp.call(clean_cmd1, shell = True)
	sp.call(clean_cmd2, shell = True)

	deep_pdb = '%s_%s.pdb' %(CDR.replace('.pdb', ''), chain)
	db_pdb = '%s_%s.pdb' %(DB_CDR.replace('.pdb', ''), chain)
	rmsd_cmd = '%s/rmsd_total_v2 %s %s bb %s > %s_%s_rmsd.txt' %(GEAR, deep_pdb, db_pdb, chain,deep_pdb.replace('.pdb',''),db_pdb.replace('.pdb',''))
	rmsd_result = '%s_%s_rmsd.txt' %(deep_pdb.replace('.pdb', ''), db_pdb.replace('.pdb', ''))
	sp.call(rmsd_cmd, shell = True)
	#print(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;",CDR.replace('.pdb', '') + '__' + DB_CDR.replace('.pdb', '') + '.pdb')
	rmsd_res = (open(rmsd_result).read()).strip().split('\t')[-1]
	#print("rmsd_res........................",rmsd_res)
	os.remove(rmsd_result)
	os.remove(deep_pdb)
	os.remove(db_pdb)
	os.remove(deep_pdb.replace('.pdb','')+'.fasta')
	os.remove(db_pdb.replace('.pdb','')+'.fasta')
	return rmsd_res

def CDR_str (pdb_str, kind, CDR_pos) :
	CDR_pdb = ''
	#print("pdb_str",pdb_str)
	CDR_str_name = pdb_str.replace('.pdb', '_') + kind + '.pdb'
	chain = kind[0]
	with open(CDR_str_name,'w') as ff:
		with open(pdb_str,'r') as ff1:
			lines = ff1.readlines()
			for line in lines:
				if line.startswith('ATOM') > 0 :
					pos = int(line[22:26].strip())
					if pos in CDR_pos:
						line = line[:21] + chain + line[22:]
						ff.write(line)

def comparison2 (kind, length, phi_list, psi_list, CDR, pdb, num) :
	CDR_str_res_name = pdb.replace('.pdb', '_') + kind + '.pdb'
	DB_pdb_key = list(DB_dic[kind][length][num].keys())  ##io
	DB_pdb_key = DB_pdb_key[0]
	DB_phi = DB_dic[kind][length][num][DB_pdb_key]['phi']
	DB_psi = DB_dic[kind][length][num][DB_pdb_key]['psi']
	DB_CDR_pos = DB_dic[kind][length][num][DB_pdb_key]['CDR_pos']
	sum_d_phi = distance (phi_list, DB_phi, 'phi')
	sum_d_psi = distance (psi_list, DB_psi, 'psi')
	sum_d_all = (sum_d_phi + sum_d_psi)
	#print(phi_list,DB_phi,psi_list,DB_psi)
	#sum_d_list.append(sum_d_all)
	sheba_run(pdb.replace('.pdb',''), DB_pdb_key,'./',pyig_model)
	#print("DB_pdb_key",DB_pdb_key,pyig_model)
	sheba1_res_name = pdb.replace('.pdb', '') + '_' + DB_pdb_key + '.pdb'
	CDR_str(sheba1_res_name, kind, DB_CDR_pos)
	DB_CDR_str_res_name = sheba1_res_name.replace('.pdb','') + '_' + kind + '.pdb'
	sheba_run(CDR_str_res_name.replace('.pdb',''), DB_CDR_str_res_name.replace('.pdb',''), './','./')
	sheba1_res_name2 = CDR_str_res_name.replace('.pdb', '') + '_' + DB_CDR_str_res_name
	rmsd_score = rmsd(CDR_str_res_name, sheba1_res_name2, kind)
	#print("rmsd_score.................................",rmsd_score)
	###print(CDR_str_res_name, sheba1_res_name2)
	if float(rmsd_score) <= 2.0 :
		sum_d_list.append(sum_d_all)
		if sum_d_all not in sum_d_dic :
			sum_d_dic[sum_d_all] = [kind + '-' + length + '-' + num]
			rmsd_dic[sum_d_all] = [str(rmsd_score)]
		else :
			sum_d_dic[sum_d_all].append(kind + '-' + length + '-' + num)
			rmsd_dic[sum_d_all].append(str(rmsd_score))
	os.remove(DB_CDR_str_res_name)
	os.remove(sheba1_res_name)
	os.remove(sheba1_res_name2)
	os.remove(DB_pdb_key+'.pdb')

def comparison (kind, length, phi_list, psi_list, CDR, pdb, CDR_pos, omega_r) :
	global sum_d_dic, sum_d_list, rmsd_dic
	CDR_str_res_name = pdb.replace('.pdb', '_') + kind + '.pdb'
	num_list = DB_dic[kind][length].keys()
	sum_d_dic = {}
	sum_d_list = []
	rmsd_dic = {}
	for num in num_list :
		if 'cis' in num :
			cis_pos = int(num[num.find('cis') + 3]) -1
			if CDR[cis_pos] != 'P' :
				pass
			else :
				cis_candi_AA_num = int(CDR_pos[cis_pos]) - 1
			#### for n, omega_line in enumerate(omega_r) :
			#### changed by io above line
				for k, omega_line in enumerate([v for v in omega_r.split("\n") if v]):
					if k <= len([v for v in omega_r.split("\n") if v ]) -2:
	#					print("omega_line:",omega_line)
						omega_num = int(omega_line[3:omega_line.find(':')])
						post_omega_AA = omega_r[k+1][:3]
						omega_angle = abs(float((omega_line[(omega_line.find(':') + 2) :]).strip()))
						if post_omega_AA == 'PRO' and cis_candi_AA_num == omega_num and omega_angle < 90:
							comparison2 (kind, length, phi_list, psi_list, CDR, pdb, num)
							break
				
		else :		
			comparison2 (kind, length, phi_list, psi_list, CDR, pdb, num)
#	rm = 'rm ' + CDR_str_res_name + '.omega'
#	sp.call(rm, shell = True)
	if sum_d_list != [] :
		best = min(sum_d_list)
		if best <= 40 :
			sum_d_list.sort()
			###for d in sum_d_list:
				###print (kind, CDR, sum_d_dic[d], d, rmsd_dic[d])
			return (sum_d_dic[best][0], best, rmsd_dic[best][0])
		else :
			return ('NA', 'NA', 'NA')
	else :
		return ('NA', 'NA', 'NA')


	
def enva (pdb, H, L, fasta_H, fasta_L):
	res_dic_H={'H1':[],'H2':[],'H3':[]}
	res_dic_L={'L1':[],'L2':[],'L3':[]}
	pdb = pdb.replace('.pdb', '')
	pdb_H = pdb + '_H.pdb'
	pdb_L = pdb + '_L.pdb'
	cmd_enva_H = '%s/enva5.0 -e %s H > %sH%s.env'%(GEAR,pdb_H,pdb[0:4],pdb_H[5:])
	cmd_enva_L = '%s/enva5.0 -e %s L > %sL%s.env'%(GEAR,pdb_L,pdb[0:4],pdb_L[5:])

	sp.call(cmd_enva_H, shell = True)
	sp.call(cmd_enva_L, shell = True)

	omega_H = 'dihed.pl -list omega %s > %s.omega'%(pdb_H,pdb_H)
	omega_L = 'dihed.pl -list omega %s > %s.omega'%(pdb_L,pdb_L)
	sp.call(omega_H, shell = True)
	sp.call(omega_L, shell = True)
	H1_omega, H2_omega, H3_omega = '', '', ''
	L1_omega, L2_omega, L3_omega = '', '', ''

#io`s command
	env_H = pdb[:4] + 'H' + pdb_H[5:] + '.env'
	env_L = pdb[:4] + 'L' + pdb_L[5:] + '.env'
	env_H_read = open(env_H).readlines()[1:]
	env_L_read = open(env_L).readlines()[1:]

	fasta_H_read = (''.join(open(fasta_H).readlines()[1:])).replace('\n','')
	fasta_L_read = (''.join(open(fasta_L).readlines()[1:])).replace('\n','')
	H = H.split('\t')
	L = L.split('\t')
	CDR_H1, CDR_H2, CDR_H3, CDR_L1, CDR_L2, CDR_L3 = H[0], H[1], H[2], L[0], L[1], L[2]
	CDR_H1_pos, CDR_H2_pos, CDR_H3_pos = range(fasta_H_read.find(CDR_H1) + 1, fasta_H_read.find(CDR_H1) + len(CDR_H1) + 1), range(fasta_H_read.find(CDR_H2) + 1, fasta_H_read.find(CDR_H2) + len(CDR_H2) + 1), range(fasta_H_read.find(CDR_H3) + 1, fasta_H_read.find(CDR_H3) + len(CDR_H3) + 1)
	CDR_L1_pos, CDR_L2_pos, CDR_L3_pos = range(fasta_L_read.find(CDR_L1) + 1, fasta_L_read.find(CDR_L1) + len(CDR_L1) + 1), range(fasta_L_read.find(CDR_L2) + 1, fasta_L_read.find(CDR_L2) + len(CDR_L2) + 1), range(fasta_L_read.find(CDR_L3) + 1, fasta_L_read.find(CDR_L3) + len(CDR_L3) + 1)
	C1, C2, C3, C4, C5, C6 = '', '', '', '', '', ''
	phi_H1, phi_H2, phi_H3, phi_L1, phi_L2, phi_L3 = [], [], [], [], [], []
	psi_H1, psi_H2, psi_H3, psi_L1, psi_L2, psi_L3 = [], [], [], [], [], []
	HP_H1, HP_H2, HP_H3, HP_L1, HP_L2, HP_L3= 0, 0, 0, 0, 0, 0
	pk_H1, pk_H2, pk_H3, pk_L1, pk_L2, pk_L3= 0, 0, 0, 0, 0, 0
	ac_H1, ac_H2, ac_H3, ac_L1, ac_L2, ac_L3= 0, 0, 0, 0, 0, 0
	po_H1, po_H2, po_H3, po_L1, po_L2, po_L3= 0, 0, 0, 0, 0, 0

	CDR_info= "%s\t%s\t%s\t%s\t%s\t%s\t" %(CDR_H1, CDR_H2, CDR_H3, CDR_L1, CDR_L2, CDR_L3)
	for omega_line in open(pdb_H + '.omega').readlines() :
		omega_num = int(omega_line[3:omega_line.find(':')])
#		print("omega_num:", omega_num)
		if omega_num in CDR_H1_pos :
			H1_omega += omega_line
		elif omega_num in CDR_H2_pos :
			H2_omega += omega_line
		elif omega_num in CDR_H3_pos :
			H3_omega += omega_line
	os.remove('%s.omega'%(pdb_H))
	for omega_line in open(pdb_L + '.omega').readlines() :
		omega_num = int(omega_line[3:omega_line.find(':')])
		if omega_num in CDR_L1_pos :
			L1_omega += omega_line
		elif omega_num in CDR_L2_pos :
			L2_omega += omega_line
		elif omega_num in CDR_L3_pos :
			L3_omega += omega_line
	os.remove('%s.omega'%(pdb_L))
	###fasta_H_read.find(CDR_H3) =
#	print("OMEGA_H&L:",H1_omega,H2_omega,H3_omega,L1_omega,L2_omega,L3_omega) 
	for line in env_H_read :
		if line.startswith('ATOM') > 0:
			pos = int(line[22:26].strip())
			AA = line[17:20]
			#HP = int(line[-45:-43].strip())
			#pk = int(line[-36:-34].strip())
			#ac = int(line[-32:-30].strip())
			#po = int(line[-29:-27].strip())
			if pos in CDR_H1_pos:
				phi = float(line[67:76].strip())
				psi = float(line[76:85].strip())
				phi_H1.append(phi)
				psi_H1.append(psi)
				C1 += AA_dic[AA]
				#print("confirm",AA_dic[AA], CDR_H1,phi,psi)
			elif pos in CDR_H2_pos:
				phi = float(line[67:76].strip())
				psi = float(line[76:85].strip())
				phi_H2.append(phi)
				psi_H2.append(psi)
				C2 += AA_dic[AA]
				#print("confirm",AA_dic[AA], CDR_H2,phi,psi)
			elif pos in CDR_H3_pos:
				phi = float(line[67:76].strip())
				psi = float(line[76:85].strip())
				phi_H3.append(phi)
				psi_H3.append(psi)
				C3 += AA_dic[AA]
				#print("confirm",AA_dic[AA], CDR_H3,phi,psi)
	for line in env_L_read :
		if line[:4] == 'ATOM' :
			pos = int(line[22:26].strip())
			AA = line[17:20]
			if pos in CDR_L1_pos:
				phi = float(line[67:76].strip())
				psi = float(line[76:85].strip())
				phi_L1.append(phi)
				psi_L1.append(psi)
				C4 += AA_dic[AA]
				#print("confirm",AA_dic[AA], CDR_L1,phi,psi)
			elif pos in CDR_L2_pos:
				phi = float(line[67:76].strip())
				psi = float(line[76:85].strip())
				phi_L2.append(phi)
				psi_L2.append(psi)
				C5 += AA_dic[AA]
				#print("confirm",AA_dic[AA], CDR_L2,phi,psi)
			elif pos in CDR_L3_pos:
				phi = float(line[67:75].strip())
				psi = float(line[76:85].strip())
				phi_L3.append(phi)
				psi_L3.append(psi)
				C6 += AA_dic[AA]
				#print("confirm",AA_dic[AA], CDR_L3,phi,psi)
	for kind in ['H1', 'H2', 'H3', 'L1', 'L2', 'L3']:
		if kind == 'H1' :
			CDR = CDR_H1
			phi_list = phi_H1
			psi_list = psi_H1
			pdb_in = pdb_H
			CDR_pos = CDR_H1_pos

			omega_r = H1_omega

		elif kind == 'H2' :
			CDR = CDR_H2
			phi_list = phi_H2
			psi_list = psi_H2
			pdb_in = pdb_H
			CDR_pos = CDR_H2_pos

			omega_r = H2_omega

		elif kind == 'H3' :
			CDR = CDR_H3
			phi_list = phi_H3
			psi_list = psi_H3
			pdb_in = pdb_H
			CDR_pos = CDR_H3_pos

			omega_r = H3_omega

		elif kind == 'L1' :
			CDR = CDR_L1
			phi_list = phi_L1
			psi_list = psi_L1
			pdb_in = pdb_L
			CDR_pos = CDR_L1_pos
			
			omega_r = L1_omega

		elif kind == 'L2' :
			CDR = CDR_L2
			phi_list = phi_L2
			psi_list = psi_L2
			pdb_in = pdb_L
			CDR_pos = CDR_L2_pos

			omega_r = L2_omega

		elif kind == 'L3' :
			CDR = CDR_L3
			phi_list = phi_L3
			psi_list = psi_L3
			pdb_in = pdb_L
			CDR_pos = CDR_L3_pos
			omega_r = L3_omega

		length = str(len(CDR))
		CDR_str (pdb_in, kind, CDR_pos)
		result_loop, result_distance, result_rmsd = comparison (kind, length, phi_list, psi_list, CDR, pdb_in, CDR_pos, omega_r)
		final_output.append(kind + '\t' + CDR + '\t' + result_loop + '\t' + str(result_distance) + '\t' + str(result_rmsd))
		os.remove(pdb_in.replace('.pdb', '_') + kind + '.pdb')
	os.remove(env_H)
	os.remove(env_L)

	return CDR_info

def TAP (r_pdb, q_pdb):
	enva (r_pdb, q_pdb)


if __name__ == '__main__':

	pdb = sys.argv[1]
	numbering = sys.argv[2]
#	numbering = 'aho'
	file_name=sys.argv[3]
	antigen_name = file_name.split("_")[1]
	cd= os.getcwd()
	GEAR = '/lwork01/neoscan_gear'
	develop_dir = '%s/%s/%s/develop'%(cd,antigen_name,file_name)
	pyig_result = '%s/pyig_result' %(develop_dir)
	script = '/KHIT1/ABS_model/script'
	os.chdir(pyig_result)

	pyig_model = '/KHIT1/ABS_model/PyIg_DB_enva4.5w'
	rosetta_dir='/lwork01/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts'
	DB_name = '%s/PyIgDB_angle_%s_enva4.5w.json' %(pyig_model,numbering)
	tool_path= '/tools2/abscan_tools'
	
	
	shutil.copy('%s/%s'%(develop_dir,pdb),pyig_result)
	
	with open(DB_name,'r') as indb:
		DB_dic=json.load(indb)

	final_output = []
	AA_dic = { "CYS":"C","SER":"S","THR":"T","PRO":"P","ALA":"A", "GLY":"G","ASN":"N","ASP":"D","GLU":"E","GLN":"Q", "HIS":"H","ARG":"R","LYS":"K","MET":"M","ILE":"I", "LEU":"L","VAL":"V","PHE":"F","TYR":"Y","TRP":"W"}
	
	if numbering == 'kabat' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = range(31, 36), range(50, 66), range(95, 103), range(24, 35), range(50, 57), range(89,98)
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = range(1, 31), range(36, 50), range(66, 95), range(103, 114), range(1, 24), range(35, 50), range(57, 89), range(98, 108)
	elif numbering == 'imgt' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = range(27, 39), range(56, 66), range(105, 118), range(27, 39), range(56, 66), range(105,118)
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = range(1, 27), range(39, 56), range(66, 105), range(118, 128), range(1, 27), range(39, 56), range(66, 105), range(118, 129)
	elif numbering == 'aho' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = range(24, 43), range(57, 70), range(107, 139), range(24, 43), range(57, 73), range(107,139)
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = range(1, 24), range(43, 57), range(70, 107), range(139, 149), range(1, 24), range(43, 57), range(70, 107), range(139, 150)
	elif numbering == 'chothia' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = range(26, 33), range(52, 57), range(95, 103), range(24, 35), range(50, 57), range(89,98)
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = range(1, 26), range(33, 52), range(57, 95), range(103, 114), range(1, 24), range(35, 50), range(57, 89), range(98, 108)
	elif numbering == 'chothia_S' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = range(27, 33), range(52, 57), range(95, 103), range(24, 35), range(50, 57), range(89,98)
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = range(1, 27), range(33, 52), range(57, 95), range(103, 114), range(1, 24), range(35, 50), range(57, 89), range(98, 108)
	
	clean(pdb)
	fasta_H = pdb.replace('.pdb', '_H.fasta')
	fasta_L = pdb.replace('.pdb', '_L.fasta')

	H = ANARCI (fasta_H, numbering, 'H')
	L = ANARCI (fasta_L, numbering, 'L')

	enva_res = enva(pdb, H, L, fasta_H, fasta_L)
	os.remove(fasta_H)
	os.remove(fasta_L)
	os.remove(fasta_H.replace('fasta', 'pdb'))
	os.remove(fasta_L.replace('fasta', 'pdb'))

	final_output2 = []
	no_of_pass = 0
	for res in final_output :
		res2 = res.split('\t')[3]
		if res2 != 'NA' :
			no_of_pass += 1
			final_output2.append('O')
		elif res2 == 'NA' :
			final_output2.append('X')

	os.remove(pdb)


#	print('\n'.join(final_output2))
	#print(output)
	pyig_out_dir = '%s/%s_pyig.tsv' %(pyig_result,pdb.replace('.deepab_new_numbering.pdb', ''))
	open(pyig_out_dir, 'w').write(pdb + '\t' + enva_res + '\t'.join(final_output2) + '\t' + str(no_of_pass) + '\n')
