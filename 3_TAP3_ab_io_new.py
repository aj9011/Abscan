import sys, os, math
import subprocess as sp
import shutil
import pandas as pd


def mk_list(col) :
    global per
    col_list = df[col].tolist()
    col_list.sort()
    col_min = min(col_list)
    col_max = max(col_list)
    col_amber_front = (col_list[:per])[-1]
    col_amber_rear = (col_list[-(per):])[0]
    return col_amber_front, col_amber_rear, col_min, col_max



def clean (pdb):
	cmd_clean_H = 'python2 %s/clean_pdb.py %s H' %(rosetta_dir,pdb)
	sp.call(cmd_clean_H, shell = True)
	cmd_clean_L = 'python2 %s/clean_pdb.py %s L' %(rosetta_dir,pdb)
	sp.call(cmd_clean_L, shell = True)
	h_pdb = pdb.replace('.pdb', '_H.pdb')
	l_pdb = pdb.replace('.pdb', '_L.pdb')
	for line in open(h_pdb) :
		line = line.strip()
		if line.startswith('ATOM') > 0 :
			num = int(line[22:26].strip())
			if num not in clean_pdb_dic['H'] :
				clean_pdb_dic['H'][num] = [line]
			else :
				clean_pdb_dic['H'][num].append(line)

	for line in open(l_pdb) :
		line = line.strip()
		if line.startswith('ATOM') > 0 :
			num = int(line[22:26].strip())
			if num not in clean_pdb_dic['L'] :
				clean_pdb_dic['L'][num] = [line]
			else :
				clean_pdb_dic['L'][num].append(line)

def ANARCI (fasta, numbering, chain) :
	#cmd1= 'ANARCI --scheme ' + numbering + ' -i '+fasta +' -o '+fasta+'.anrc'
	cmd1= 'ANARCI --scheme %s -i %s -o  %s.anrc' %(numbering,fasta,fasta)
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
		os.remove(fasta+'.anrc')
		return( '\t'.join(res))
	elif chain == 'L':
		keys=sorted(res_dic_L.keys())
		res=[]
		for key in keys:
			res.append(res_dic_L[key].replace('-','').strip())
		os.remove(fasta+'.anrc')
		return( '\t'.join(res))

def distance (x1, y1, z1, x2, y2, z2) :
	return math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)

def Hydrophobicity () :
	result = 0
	chain_pdb_nums = []
	check_list = []
	for chain in ['H', 'L'] :
		pdb_nums = clean_pdb_dic[chain].keys()
		for num in pdb_nums :
			chain_pdb_nums.append(chain + '_' + str(num))
  
	for n, r_chain_num in enumerate(chain_pdb_nums) :
		r_chain_num = r_chain_num.split('_')
		r_chain = r_chain_num[0]
		r_num = int(r_chain_num[1])
		if r_num in CDR_pos[r_chain] or r_num in anchor_pos[r_chain] or r_num in in4_pos[r_chain]:
			if r_num in SER_dic[r_chain] :
				for q_chain_num in chain_pdb_nums[:n] + chain_pdb_nums[n+1 :] :
					q_chain_num = q_chain_num.split('_')
					if set(['_'.join(r_chain_num), '_'.join(q_chain_num)]) not in check_list :
						q_chain = q_chain_num[0]
						q_num = int(q_chain_num[1])
						if q_num in CDR_pos[q_chain] or q_num in anchor_pos[q_chain] or q_num in in4_pos[q_chain]:
							if q_num in SER_dic[q_chain] :
								r_list, q_list = clean_pdb_dic[r_chain][r_num], clean_pdb_dic[q_chain][q_num]
								dis_list = []
								for r_line in r_list :
									if r_line.strip()[-1] in heavy_atom :
										x1, y1, z1 = float(r_line[27:38].strip()), float(r_line[38:46].strip()), float(r_line[46:54].strip())
										for q_line in q_list :
											if q_line.strip()[-1] in heavy_atom :
												x2, y2, z2 = float(q_line[27:38].strip()), float(q_line[38:46].strip()), float(q_line[46:54].strip())
												dis = distance (x1, y1, z1, x2, y2, z2)
												dis_list.append(dis)
								if min(dis_list) < 7.5 :
									r12 = min(dis_list)
									r_resi = r_line[17:20]
									if r_num in salt_bridge_dic[r_chain] :
										r_resi = 'GLY'
									q_resi = q_line[17:20]
									if q_num in salt_bridge_dic[q_chain] :
										q_resi = 'GLY'
									r_Hydro = nor_hydro_dic[r_resi]
									q_Hydro = nor_hydro_dic[q_resi]
									Hydro = (r_Hydro * q_Hydro)/(r12**2)
									#print(r_num,r_chain,q_num,q_chain,Hydro)
									result += Hydro
									check_list.append(set(['_'.join(r_chain_num), '_'.join(q_chain_num)]))
	return (result)




def SER (pdb, chain) :
	cmd = 'freesasa -f seq --shrake-rupley ' + pdb
	enva_res = (sp.getoutput(cmd))
	for line in enva_res.split('\n') :
		if line[:3] == 'SEQ' and line[4] == chain :
			ASA = float(line.strip()[-8:].strip())
			AA = line[12:15]
			pos = int(line[6:10].strip())
			if (ASA/Ala_X_Ala[AA]) * 100.0 >= 7.5 :
				SER_dic[chain].append(pos)

def salt_bridge (pdb) :
	for line in open(pdb).readlines() :
		if line.startswith('ATOM') > 0 :
			AA = line[17:20].strip()
			point = line[13:17].strip()
			if (AA == 'LYS' and point == 'NZ') or (AA == 'ARG' and point == 'NH1') :
				AB_dic['A'].append(line.strip())
			elif (AA == 'ASP' and point == 'OD2') or (AA == 'GLU' and point == 'OE2'):
				AB_dic['B'].append(line.strip())
	for A_line in AB_dic['A'] :
		x1, y1, z1 = float(A_line[27:38].strip()), float(A_line[38:46].strip()), float(A_line[46:54].strip())
		for B_line in AB_dic['B'] :
			x2, y2, z2 = float(B_line[27:38].strip()), float(B_line[38:46].strip()), float(B_line[46:54].strip())
			distance_res = distance (x1, y1, z1, x2, y2, z2)
			if distance_res < 3.2 :
				A_data = A_line[17:20].strip() + ' ' + A_line[22:26].strip()
				A_chain = A_line[21]
				B_data = B_line[17:20].strip() + ' ' + B_line[22:26].strip()
				B_chain = B_line[21]
				salt_bridge_dic[A_chain].append(int(A_line[22:26].strip()))
				salt_bridge_dic[B_chain].append(int(B_line[22:26].strip()))


				
def in4 (pdb) :
	pos_list = []
	dic = {}
	for line in open(pdb).readlines() :
		if line.startswith('ATOM') > 0 :
			chain = line[21]
			pos = int(line[22:26].strip())
			atom = line.strip()[-1]
			if pos in SER_dic[chain] :#and atom in heavy_atom:
				pos = str(pos)
				pos_list.append(chain + pos)
				if chain + pos not in dic :
					dic[chain + pos] = {atom : [line.strip()]}
				else :
					if atom not in dic[chain + pos] :
						dic[chain + pos][atom] = [line.strip()]
					else :
						dic[chain + pos][atom].append(line.strip())
	pos_list = list(set(pos_list))
	for chain_pos in pos_list :
		flag = 'x'
		atom_list = dic[chain_pos].keys()
		pos = int(chain_pos[1:].strip())
		chain = chain_pos[0]
		for atom in atom_list :
			r_lines = dic[chain_pos][atom]
			for r_line in r_lines :
				if r_line[13:17].strip() == 'CA' :
					x1, y1, z1 = float(r_line[27:38].strip()), float(r_line[38:46].strip()), float(r_line[46:54].strip())
					for chain_pos2 in pos_list :
						if chain_pos != chain_pos2 :
							for atom2 in list(dic[chain_pos2].keys()) :
								if atom2 in heavy_atom :
									q_lines = dic[chain_pos2][atom2]
									for q_line in q_lines :
										x2, y2, z2 = float(q_line[27:38].strip()), float(q_line[38:46].strip()), float(q_line[46:54].strip())
										dis = distance (x1, y1, z1, x2, y2, z2)
										if dis < 2.455:
											flag = 'o'
											break
                   
		if flag == 'o' :
			chain = chain_pos[0]
			pos = int(chain_pos[1:].strip())
			in4_pos[chain].append(pos)

def charge_P () :
	check_list = []
	PPC = 0
	PNC = 0
	chain_pdb_nums = []
	for chain in ['H', 'L'] :
		pdb_nums = clean_pdb_dic[chain].keys()
		for num in pdb_nums :
			chain_pdb_nums.append(chain + '_' + str(num))
	for n, r_chain_num in enumerate(chain_pdb_nums) :
		r_chain_num = r_chain_num.split('_')
		r_chain = r_chain_num[0]
		r_num = int(r_chain_num[1])
		if r_num in SER_dic[r_chain] :
			if r_num in CDR_pos[r_chain] or r_num in anchor_pos[r_chain] or r_num in in4_pos[r_chain] :
				#for q_chain_num in chain_pdb_nums[n+1 :] :
				for q_chain_num in chain_pdb_nums[:n] + chain_pdb_nums[n+1 :] :
					q_chain_num = q_chain_num.split('_')
					if set(['_'.join(r_chain_num), '_'.join(q_chain_num)]) not in check_list :
						q_chain = q_chain_num[0]
						q_num = int(q_chain_num[1])
						if q_num in SER_dic[q_chain] :
							if q_num in CDR_pos[q_chain] or q_num in anchor_pos[q_chain] or q_num in in4_pos[q_chain] :
								r_list, q_list = clean_pdb_dic[r_chain][r_num], clean_pdb_dic[q_chain][q_num]
								dis_list = []
								for r_line in r_list :
									x1, y1, z1 = float(r_line[27:38].strip()), float(r_line[38:46].strip()), float(r_line[46:54].strip())
									for q_line in q_list :
										x2, y2, z2 = float(q_line[27:38].strip()), float(q_line[38:46].strip()), float(q_line[46:54].strip())
										dis = distance (x1, y1, z1, x2, y2, z2)
										dis_list.append(dis)
								r_resi = r_line[17:20]
								q_resi = q_line[17:20] 
								r_Hydro = 0
								q_Hydro = 0
								#r_Hydro = nor_hydro_dic[r_resi]
								#q_Hydro = nor_hydro_dic[q_resi]
								if r_resi in charge_dic :
									r_Hydro = charge_dic[r_resi]
								if q_resi in charge_dic :
									q_Hydro = charge_dic[q_resi]
								if r_num in salt_bridge_dic[r_chain] :
									r_Hydro = 0
								if q_num in salt_bridge_dic[q_chain] :
									q_Hydro = 0
								if r_Hydro != 0 and q_Hydro != 0 and min(dis_list) < 7.5 and ((r_Hydro < 0 and q_Hydro < 0) or (r_Hydro > 0 and q_Hydro > 0)) :
									r12 = min(dis_list)
									if (r_Hydro > 0 and q_Hydro > 0) :
										Hydro = (r_Hydro * q_Hydro)/(r12**2)
										PPC += Hydro
									elif (r_Hydro < 0 and q_Hydro < 0) :
										Hydro = (r_Hydro * q_Hydro)/(r12**2)
										PNC += Hydro
									check_list.append(set(['_'.join(r_chain_num), '_'.join(q_chain_num)]))
	return (PPC, PNC) 
        #if r_num in SER_dic[r_chain] :
        #    if r_num in CDR_pos[r_chain] or r_num in anchor_pos[r_chain] or r_num in in4_pos[r_chain] :
 
def charge_SFvCSP () :
	H_charge = 0
	L_charge = 0
	chain_pdb_nums = []
	for chain in ['H', 'L'] :
		pdb_nums = clean_pdb_dic[chain].keys()
		for num in pdb_nums :
			chain_pdb_nums.append(chain + '_' + str(num))
	for n, r_chain_num in enumerate(chain_pdb_nums) :
		r_chain_num = r_chain_num.split('_')
		r_chain = r_chain_num[0]
		r_num = int(r_chain_num[1])
		r_resi = clean_pdb_dic[r_chain][r_num][0][17:20]
		if r_num in SER_dic[r_chain] : 
			if r_resi in charge_dic :
				r_charge = charge_dic[r_resi]
				if r_num in salt_bridge_dic[r_chain] :
					r_charge = 0
				if r_chain == 'H' :
					H_charge += r_charge
				elif r_chain == 'L' :
					L_charge += r_charge
	return(H_charge * L_charge)
        
#def cleaning(pdb) :
#    cmd = 'rm ' + pdb.replace('.pdb', '_H.pdb') + ' ' + pdb.replace('.pdb', '_L.pdb') + ' ' + pdb.replace('.pdb', '_H.fasta') + ' ' + pdb.replace('.pdb', '_L.fasta')
#    sp.call(cmd, shell = True) 



def check_DB (value, col_min, col_amber_front, col_amber_rear, col_max) :
	if col_amber_front < value < col_amber_rear :
		score = '1'
	elif (col_min <= value <= col_amber_front) or (col_amber_rear <= value <= col_max) :
		score = '0.5'
	elif (value < col_min) or (value > col_max) :
		score = '0'
	return score


if __name__ == '__main__':

	pdb = sys.argv[1]
	numbering = sys.argv[2]

	file_name = sys.argv[3]
	antigen_name = file_name.split("_")[1]


	cd =os.getcwd()



	ABS_model = '/KHIT1/ABS_model'
	DB = '%s/Therapeutic_Models_526_db_new_numbering.tsv' %(ABS_model)
	tools_dir= '/tools2/abscan_tools'
	clean_py = '%s/clean_pdb_JH.py' %(tools_dir)
	develop_dir = '%s/%s/%s/develop' %(cd,antigen_name,file_name)
	tap_result_dir= '%s/tap_result' %(develop_dir)
	rosetta_dir='/lwork01/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts'

	os.chdir(tap_result_dir)

	shutil.copy('%s/%s'%(develop_dir,pdb),tap_result_dir)

	df = pd.read_csv(DB, sep = '\t')
	per = int(round(df.shape[0] * 0.05))

	len_amber_front, len_amber_rear, len_min, len_max = mk_list('CDR_len')
	hydro_amber_front, hydro_amber_rear, hydro_min, hydro_max = mk_list('PSH')
	PPC_amber_front, PPC_amber_rear, PPC_min, PPC_max = mk_list('PPC')
	PNC_amber_front, PNC_amber_rear, PNC_min, PNC_max = mk_list('PNC')
	SFvCSP_amber_front, SFvCSP_amber_rear, SFvCSP_min, SFvCSP_max = mk_list('SFvCSP')

	hydro_dic = {'ALA' : 1.80, 'CYS' : 2.50, 'ASP' : -3.50, 'GLU' : -3.50, 'PHE' : 2.80, 'GLY' : -0.40, 'HIS' : -3.20, 'ILE' : 4.50, 'LYS' : -3.90, 'LEU' : 3.80, 'MET' : 1.90, 'ASN' : -3.50, 'PRO' : -1.60, 'GLN' : -3.50, 'ARG' : -4.50, 'SER' : -0.80, 'THR' : -0.70, 'VAL' : 4.20, 'TRP' : -0.90, 'TYR' : -1.30}
	Ala_X_Ala = {'ALA' : 110.2, 'CYS' : 140.4, 'ASP' : 144.1, 'GLU' : 174.7, 'PHE' : 200.7, 'GLY' : 78.7, 'HIS' : 181.9, 'ILE' : 185.0, 'LYS' : 205.7, 'LEU' : 183.1, 'MET' : 200.1, 'ASN' : 146.4, 'PRO' : 141.9, 'GLN' : 178.6, 'ARG' : 229.0, 'SER' : 117.2, 'THR' : 138.7, 'VAL' : 153.7, 'TRP' : 240.5, 'TYR' : 213.7}
	charge_dic = {'ASP' : -1, 'GLU' : -1, 'LYS' : +1, 'ARG' : +1, 'HIS' : +0.1}
	heavy_atom = {'C' : '', 'O' : '', 'N' : '', 'S' : ''}
	nor_hydro_dic = {}

	if numbering == 'kabat' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = range(31, 36), range(50, 66), range(95, 103), range(24, 35), range(50, 57), range(89,98)
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = range(1, 31), range(36, 50), range(66, 95), range(103, 114), range(1, 24), range(35, 50), range(57, 89), range(98, 108)
	elif numbering == 'imgt' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = list(range(27, 39)), list(range(56, 66)), list(range(105, 118)), list(range(27, 39)), list(range(56, 66)), list(range(105,118))
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = list(range(1, 27)), list(range(39, 56)), list(range(66, 105)), list(range(118, 128)), list(range(1, 27)), list(range(39, 56)), list(range(66, 105)), list(range(118, 129))
	elif numbering == 'aho' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = range(24, 43), range(57, 70), range(107, 139), range(24, 43), range(57, 73), range(107,139)
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = range(1, 24), range(43, 57), range(70, 107), range(139, 149), range(1, 24), range(43, 57), range(70, 107), range(139, 150)
	elif numbering == 'chothia' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = range(26, 33), range(52, 57), range(95, 103), range(24, 35), range(50, 57), range(89,98)
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = range(1, 26), range(33, 52), range(57, 95), range(103, 114), range(1, 24), range(35, 50), range(57, 89), range(98, 108)
	elif numbering == 'chothia_S' :
		H1_pos, H2_pos, H3_pos, L1_pos, L2_pos, L3_pos = range(27, 33), range(52, 57), range(95, 103), range(24, 35), range(50, 57), range(89,98)
		FRH1_pos, FRH2_pos, FRH3_pos, FRH4_pos, FRL1_pos, FRL2_pos, FRL3_pos, FRL4_pos = range(1, 27), range(33, 52), range(57, 95), range(103, 114), range(1, 24), range(35, 50), range(57, 89), range(98, 108)



	clean_pdb_dic = {'H' : {}, 'L' : {}}

	
	max_H = (max(hydro_dic.values()))
	min_H = (min(hydro_dic.values()))

	for key in hydro_dic :
		value = hydro_dic[key]
		new_value = ((value - min_H)/(max_H - min_H)) + 1
		nor_hydro_dic[key] = new_value


	CDR_pos = {'H' : [], 'L' : []}
	anchor_pos = {'H' : [], 'L' : []}
	in4_pos = {'H' : [], 'L' : []}
	SER_dic = {'H' : [], 'L' : []}
	AB_dic = {'A' : [], 'B' : []}
	salt_bridge_dic = {'H' : [], 'L' : []}


	clean(pdb)
	fasta_H = pdb.replace('.pdb', '_H.fasta')
	fasta_L = pdb.replace('.pdb', '_L.fasta')
	H = ANARCI (fasta_H, numbering, 'H')
	L = ANARCI (fasta_L, numbering, 'L')
	fasta_H_read = (''.join(open(fasta_H).readlines()[1:])).replace('\n','')
	fasta_L_read = (''.join(open(fasta_L).readlines()[1:])).replace('\n','')
	H = H.split('\t')
	L = L.split('\t')
	CDR_H1, CDR_H2, CDR_H3, CDR_L1, CDR_L2, CDR_L3 = H[0], H[1], H[2], L[0], L[1], L[2]
	CDR_H1_pos, CDR_H2_pos, CDR_H3_pos = list(range(fasta_H_read.find(CDR_H1) + 1, fasta_H_read.find(CDR_H1) + len(CDR_H1) + 1)), list(range(fasta_H_read.find(CDR_H2) + 1, fasta_H_read.find(CDR_H2) + len(CDR_H2) + 1)), list(range(fasta_H_read.find(CDR_H3) + 1, fasta_H_read.find(CDR_H3) + len(CDR_H3) + 1))
	CDR_L1_pos, CDR_L2_pos, CDR_L3_pos = list(range(fasta_L_read.find(CDR_L1) + 1, fasta_L_read.find(CDR_L1) + len(CDR_L1) + 1)), list(range(fasta_L_read.find(CDR_L2) + 1, fasta_L_read.find(CDR_L2) + len(CDR_L2) + 1)), list(range(fasta_L_read.find(CDR_L3) + 1, fasta_L_read.find(CDR_L3) + len(CDR_L3) + 1))
	CDR_pos['H'] = CDR_H1_pos + CDR_H2_pos + CDR_H3_pos
	CDR_pos['L'] = CDR_L1_pos + CDR_L2_pos + CDR_L3_pos
	anchor_pos['H'] = [(min(CDR_H1_pos) -1), (max(CDR_H1_pos) +1), (min(CDR_H2_pos) -1), (max(CDR_H2_pos) +1), (min(CDR_H3_pos) -1), (max(CDR_H3_pos) +1)]
	anchor_pos['L'] = [(min(CDR_L1_pos) -1), (max(CDR_L1_pos) +1), (min(CDR_L2_pos) -1), (max(CDR_L2_pos) +1), (min(CDR_L3_pos) -1), (max(CDR_L3_pos) +1)]




	SER (pdb, 'H')
	SER (pdb, 'L')
	in4(pdb)
	salt_bridge(pdb)
	PPC, PNC = (charge_P())
#	cleaning(pdb)
	#CDR_len = (len(CDR_pos['H']) + len( CDR_pos['L']))
	CDR_len = len(CDR_H1) + len(CDR_H2) + len(CDR_H3) + len(CDR_L1) + len(CDR_L2) + len(CDR_L3)
	Hydro_score = Hydrophobicity ()
	SFvCSP = charge_SFvCSP ()

	#print(len_min, len_amber_front, len_amber_rear, len_max)
	#print(hydro_min, hydro_amber_front, hydro_amber_rear, hydro_max)
	#print(PPC_amber_rear, PPC_max)
	#print(PNC_amber_rear, PNC_max)
	#print(SFvCSP_min, SFvCSP_amber_front)



	CDR_len_DBscore = check_DB (CDR_len, len_min, len_amber_front, len_amber_rear, len_max)
	Hydro_DBscore = check_DB (Hydro_score, hydro_min, hydro_amber_front, hydro_amber_rear, hydro_max)
	PPC_DBscore = check_DB (PPC, -10000, -10000, PPC_amber_rear, PPC_max)
	PNC_DBscore = check_DB (PNC, -10000, -10000, PNC_amber_rear, PNC_max)
	SFvCSP_DBscore = check_DB (SFvCSP, SFvCSP_min, SFvCSP_amber_front, 10000, 10000)

#	rm_pdb = 'rm ' + pdb
#	sp.call(rm_pdb, shell = True)
	os.remove(pdb)
	os.remove(pdb.replace('.pdb', '_H.pdb'))
	os.remove(pdb.replace('.pdb', '_L.pdb'))
	os.remove(pdb.replace('.pdb', '_H.fasta'))
	os.remove(pdb.replace('.pdb', '_L.fasta'))

	#print(pdb.replace('cp_', '') + '\t' + str(CDR_len) + '\t' + str(Hydro_score) + '\t' + str(PPC) + '\t' + str(PNC) + '\t' + str(SFvCSP) + '\t' + CDR_len_DBscore + '\t' + Hydro_DBscore + '\t' + PPC_DBscore + '\t' + PNC_DBscore + '\t' + SFvCSP_DBscore)

	CDR_info='%s\t%s\t%s\t%s\t%s\t%s\t' %(CDR_H1,CDR_H2,CDR_H3,CDR_L1,CDR_L2,CDR_L3)

	out_dir= '%s/%s_tap.tsv' %(tap_result_dir,pdb.replace('.deepab_new_numbering.pdb',''))
	#open(output, 'a').write(output_line + '\t' + CDR_H1 +'\t' + CDR_H2 + '\t' + CDR_H3 +'\n')
	open(out_dir, 'a').write(pdb + '\t' + CDR_info + str(CDR_len) + '\t' + str(Hydro_score) + '\t' + str(PPC) + '\t' + str(PNC) + '\t' + str(SFvCSP) + '\t' + CDR_len_DBscore + '\t' + Hydro_DBscore + '\t' + PPC_DBscore + '\t' + PNC_DBscore + '\t' + SFvCSP_DBscore + '\n')
