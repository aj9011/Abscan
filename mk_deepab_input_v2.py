import os,sys
import requests
import urllib.request
import logging
import re
import shutil
import subprocess as sp




def download_file(url, output_path):
    logging.basicConfig(format='%(asctime)s %(message)s', filename='pdb_log.txt', level=logging.ERROR)
    req = requests.get(url)
    if not req.ok:
        logging.error(f"Error code = {req.status_code}")
    else:
        with open(output_path, 'w') as f:
            f.write(req.content.decode('utf-8'))



def download_chothia_pdb_files(pdb_ids, antibody_database_path,
                               max_workers=10):
    """
    :param pdb_ids: A set of PDB IDs to download
    :type pdb_ids: set(str)
    :param antibody_database_path: Path to the directory to save the PDB files to.
    :type antibody_database_path: str
    :param max_workers: Max number of workers in the thread pool while downloading.
    :type max_workers: int
    """
    pdb_file_paths = [os.path.join(antibody_database_path, pdb + '.pdb') for pdb in pdb_ids]
    # Download PDBs using multiple threads
    download_url = 'http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/pdb/{}/?scheme=chothia'
    urls = [download_url.format(pdb) for pdb in pdb_ids]

    for args in zip(urls, pdb_file_paths):
        download_file(args[0], args[1])




def get_cdrh3(ids):
	pattern = "</a><br>H3: <a href.*?</a><br>"
	pattern2 = '<b>Sequence</b></td><td>.*'
	cdrh3_structures = {}
	long_seqs = list()
	for n, id_ in enumerate(ids):
		cdrh3_structures[id_] = list()
		url = f"http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/cdrsearch/?pdb={id_}&CDRdef_pdb=Chothia"
		with urllib.request.FancyURLopener({}).open(url) as conn:
			content = conn.read()
		matches = re.findall(pattern, content.decode())
		for m in matches:
			m = m[:-8]
			i = -1
			while m[i - 1] != '>':
				i -= 1
			complete = m[i - 12:i - 2] != 'incomplete'
			chain_pattern = 'chain=.'
			chain = re.findall(chain_pattern, m)[0][-1]
			if complete:
				m = m[i:]
				if '...' in m:
					break
					long_seqs.append(id_)
					url = f"http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/cdrviewer/?CDRdef=chothia&pdb={id_}&chain={chain}&loop=CDRH3"
					with urllib.request.FancyURLopener({}).open(url) as f:
						content = f.read()
					m = re.search(pattern2, content.decode())[0]
					m = m[:-10]
					i = -1
					while m[i - 1] != '>':
						i -= 1
				m = m[i:]
				cdrh3_structures[id_].append((chain, m))
			else:
				print(n, id_, (chain, m))
	        # print(n, id_, cdrh3_structures[id_])
	return cdrh3_structures






if __name__ == "__main__":

	cd= os.getcwd()
	file_name = sys.argv[1]
	antigen = file_name.split("_")[1]
	cdrgen_result_dir = '%s/%s/%s/cdrgen_result' %(cd,antigen, file_name)
	pdb_info_dir = '%s/%s/%s/input/%s_pdb_information.txt' %(cd,antigen,file_name,antigen)

	rosetta_bin='/lwork01/rosetta_src_2019.40.60963_bundle/main/source/bin'
	benchmark_dic={}
	#os.chdir('./epitope')
#	deepab_dir= '%s/DeepAb' %(cd)
	deepab_dir= '%s/%s/%s/DeepAb'%(cd,antigen,file_name)
	deepab_inp_dir = '%s/%s/%s/DeepAb/input'%(cd,antigen,file_name)
	with open(pdb_info_dir,'r') as b_r:
		b_lines=b_r.readlines()[1:]
		for b in b_lines:
			sp_line=b.strip().split("\t")
			pdb_name=sp_line[0]
			antigen_name=sp_line[1]
			heavy_chain=sp_line[3]
			light_chain=sp_line[4]
			antigen_chain =sp_line[5]
			cdr_cd_start=int(sp_line[8])
			cdr_cd_end=int(sp_line[9])
			bench_info='%s_%s_%s_%s_%s_%s' %(pdb_name, heavy_chain,light_chain,str(cdr_cd_start),str(cdr_cd_end),antigen_chain)
			benchmark_dic[antigen_name] = bench_info

#	deepab_data_file= '%s/data/%s/%s' %(deepab_dir,antigen_name, file_name)


	if not os.path.exists(deepab_dir):
		os.makedirs(deepab_dir)
		os.makedirs(deepab_inp_dir)

	
	pdb_dir= '%s/%s/%s/PDB' %(cd,antigen,file_name)
	if not os.path.exists(pdb_dir):
		os.makedirs(pdb_dir)
	
 
	cdr_result_files=[x for x in os.listdir(cdrgen_result_dir) if x.endswith('_E_total_result.txt')]

	for cdr_result_file in cdr_result_files:
		fantigen=cdr_result_file.strip().split("_")[0]
		fpdb_name= benchmark_dic[fantigen].split("_")[0]
		fheavy_chain= benchmark_dic[fantigen].split("_")[1]
		flight_chain=benchmark_dic[fantigen].split("_")[2]
		fcdr_cd_start= int(benchmark_dic[fantigen].split("_")[3])
		fcdr_cd_end=  int(benchmark_dic[fantigen].split("_")[4])
		fantigen_chain = benchmark_dic[fantigen].split("_")[5]
		fpdb_list=[]
		fpdb_list.append(fpdb_name)

		pdb_file_dir= '%s/%s.pdb' %(pdb_dir, fpdb_name)
		clean_pdb = '%s/%s_clean.pdb' %(pdb_dir,fpdb_name)
		HL_pdb = '%s/%s_HL.pdb' %(pdb_dir,fpdb_name)
		HLA_pdb = '%s/%s_HLA.pdb' %(pdb_dir,fpdb_name)
		if not os.path.isfile(pdb_file_dir):
			download_chothia_pdb_files(list(fpdb_list),pdb_dir,max_workers=10)

			output_pdb=[]
			output_H_pdb=[]
			output_L_pdb=[]
			output_A_pdb=[]
			with open(pdb_file_dir,'r') as r_p:
				r_lines=r_p.readlines()
				for r_line in r_lines:
					if r_line[:4] == 'ATOM':
						output_pdb.append(r_line.strip())
						if r_line[21] == fheavy_chain:
							output_H_pdb.append(r_line.strip())
						elif r_line[21] == flight_chain:
							output_L_pdb.append(r_line.strip())
						elif r_line[21] == fantigen_chain:
							output_A_pdb.append(r_line.strip())
			
			output_H_pdb.append('TER')
			output_HL_pdb = output_H_pdb + output_L_pdb
			output_HL_lines = '\n'.join(output_HL_pdb)

			output_HL_pdb.append('TER')
			output_HLA_pdb = output_HL_pdb + output_A_pdb
			output_HLA_lines = '\n'.join(output_HLA_pdb)

			print(output_HL_lines)			
			output_pdb = '\n'.join(output_pdb)
			with open(clean_pdb,'w') as w_c:
				w_c.write(output_pdb)

			with open(HL_pdb,'w') as w_hl:
				w_hl.write(output_HL_lines)


			with open(HLA_pdb,'w') as w_hll:
				w_hll.write(output_HLA_lines)


			

			interface_score_cmd =  '%s/InterfaceAnalyzer.linuxgccrelease -interface LH_A  -s %s -out:path:pdb  %s' %(rosetta_bin,HLA_pdb,pdb_dir)
			sp.call(interface_score_cmd,shell=True)
			
			rosetta_py='/lwork01/miniconda3/envs/python2/bin/python2'  ##directory adjustments are required
			get_H_fasta='%s /lwork01/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts/get_fasta_from_pdb.py %s %s %s/%s_%s_H.pdb.fasta' %(rosetta_py,clean_pdb,fheavy_chain,pdb_dir,fpdb_name,fheavy_chain)

			get_L_fasta='%s /lwork01/rosetta_src_2019.40.60963_bundle/tools/protein_tools/scripts/get_fasta_from_pdb.py %s %s %s/%s_%s_L.pdb.fasta' %(rosetta_py,clean_pdb,flight_chain,pdb_dir,fpdb_name,flight_chain)

			sp.call(get_H_fasta,shell=True)
			sp.call(get_L_fasta,shell=True)
			print(get_H_fasta,get_L_fasta)
		
		h_fasta =  '%s/%s_%s_H.pdb.fasta' %(pdb_dir,fpdb_name,fheavy_chain)
		l_fasta = '%s/%s_%s_L.pdb.fasta' %(pdb_dir,fpdb_name,flight_chain)

		with open(h_fasta,'r') as h_f:
			h_seq=h_f.readlines()[1]

		with open(l_fasta,'r') as l_f:
			l_seq=l_f.readlines()[1]


		cc= get_cdrh3(list(fpdb_list))
		fpdb_list=[]
		print(cc)
		cdrh3r= cc[fpdb_name]
		print (cdrh3r)

		for h3 in cdrh3r:
			if h3[0] == fheavy_chain:
				cdrh3=h3[1].strip()


		cdr_result_file_dir = '%s/%s' %(cdrgen_result_dir,cdr_result_file)

		epitope=cdr_result_file.split("_")[1]

		total_set_temp = []
		total_set_dir = '%s/total_%s_epi_cdr_set.txt' %(cdrgen_result_dir,epitope)

		with open(cdr_result_file_dir,'r') as cdr_r:
			cdrs=cdr_r.readlines()
			for cdr_candi in cdrs:
				#print(cdrh3[0:fcdr_cd_start-1], cdr_candi,cdrh3[fcdr_cd_end:])
				new_cdr=cdrh3[0:fcdr_cd_start-1] + cdr_candi.strip()+ cdrh3[fcdr_cd_end:].strip()
				new_H3= h_seq.replace(cdrh3,new_cdr).strip()
				fasta_line= ">:H\n%s\n>:L\n%s" %(new_H3,l_seq)
				epi_fasta_dir= '%s/%s_%s.fasta' %(deepab_inp_dir,fantigen,new_cdr)
				with open(epi_fasta_dir, 'w')  as fa_w:
					fa_w.write(fasta_line)

				nn= '%s\t%s' %(epitope,new_cdr)
				total_set_temp.append(nn)
			#	deepab_dir_new= '%s/data/%s/%s' %(deepab_dir,antigen_name,file_name)
			#	shutil.copy(epi_fasta_dir,deepab_dir_new)

		total_set = list(set(total_set_temp))

		total_set_w = '\n'.join(total_set)
		with open(total_set_dir, 'w' ) as tw:
			tw.write(total_set_w)
				
						
						



						
	
				






