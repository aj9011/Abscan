import sys, os, json
from functools import lru_cache  as cache
import pandas as pd
from tqdm import tqdm




@cache(maxsize=None)
def merge (outpep, direction, tail, pos, connected,outpep_hbond) :
	global outpep_dic_ant, ant_dic, ant_dic_pro, outpep_dic_par, par_dic, par_dic_pro, count, total_results ,matching_count ,count,pep
	if direction == 'ant' :
		outpep_dic = outpep_dic_ant
		dic = ant_dic
		dic_pro = ant_dic_pro
	else:
		outpep_dic = outpep_dic_par
		dic = par_dic
		dic_pro = par_dic_pro
	if pos ==0 :
		range_tail= range(len(tail))[::-1]
	else:
		range_tail = range(len(outpep[1:]))[::-1]
	#print(outpep,hbond)
	for tail_n in range_tail :
		tail_pos = pos + tail_n + 1
		if tail_pos in outpep_dic :
			for v, post_hbond in enumerate(set(list(outpep_dic[tail_pos]))) :
				post=post_hbond.split("_")[0]
				if pos == 0 :
					tail2 =tail[tail_n:]
				else:
					tail2 =outpep[1:][tail_n:]
				if post[:len(tail2)] == tail2 :
					#print(outpep, connected, post)
					connected_freg = connected + post[len(tail2):]
					set_post=[]
					set_post.append(outpep+'\t'+outpep_hbond)
					if len(connected_freg) == len(pep) :
						matching_count +=1
						print(matching_count)
						total_results.append(connected_freg)
						set_post.append(post+'\t'+post_hbond)
						r_li = pep+"\n"+ '\n'.join(set_post)
						total_3mer_list.append(r_li)
						only_3mer.append(connected_freg)
						print(r_li)
						r_li=''
						set_post=[]	
						
						del post, connected_freg
					else :
						#print(outpep,post,connected,connected_freg,outpep_hbond)
						merge (post, direction, tail2, tail_pos, connected_freg,outpep_hbond)
						

def cut(pep, direction) :
    global outpep_dic_ant, ant_dic, ant_dic_pro, outpep_dic_par, par_dic, par_dic_pro, start_dic_ant, start_dic_par
    if direction == 'ant' :
        outpep_dic = outpep_dic_ant
        dic = ant_dic
        dic_pro = ant_dic_pro
        start_dic = start_dic_ant
    elif direction == 'par' :
        outpep_dic = outpep_dic_par
        dic = par_dic
        dic_pro = par_dic_pro
        start_dic = start_dic_par
    for num in [0,1] :
        if num == 0 :
            normal_dic = dic
            pro_dic = dic_pro
        elif num == 1 :
            normal_dic = dic_pro
            pro_dic = dic
        for k in range(len(pep) - 2) :
           for b in range(len(pep) - 2 - k) :
              inpep = pep[k : k + 3 + b]
              if k == 0 :
                  if inpep in normal_dic :
                      outpep = normal_dic[inpep]
                      if k not in outpep_dic :
                          outpep_dic[k] = outpep
                      else :
                          outpep_dic[k] += outpep
                  start_dic[inpep] = ''
              else :
                  if len(inpep) >= 3 :
                      if inpep in normal_dic :
                          outpep = normal_dic[inpep]
                          if k not in outpep_dic :
                              outpep_dic[k] = outpep
                          else :
                              outpep_dic[k] += outpep

def run (direction) :
	global outpep_dic_ant, ant_dic, ant_dic_pro, outpep_dic_par, par_dic, par_dic_pro, start_dic_ant, start_dic_par,  matching_count, count
	if direction == 'ant' :
		outpep_dic = outpep_dic_ant
		dic = ant_dic
		dic_pro = ant_dic_pro
		start_dic = start_dic_ant
	elif direction == 'par' :
		outpep_dic = outpep_dic_par
		dic = par_dic
		dic_pro = par_dic_pro
		start_dic = start_dic_par
	for start_hbond in outpep_dic[0] :
		start=start_hbond.split("_")[0]
        #print(start_hbond,start)
        #print(hbond)
		tail = start[1:]
		if len(start) == len(pep) :
			matching_count +=1
			total_results.append(start)
			total_4mer_list.append(start_hbond)
			only_4mer.append(start)

		else : 
			merge (start, direction, tail, 0, start,start_hbond)
	del outpep_dic
	del dic
	del dic_pro
	del start_dic

if __name__ == "__main__":


	DB = sys.argv[1]
	pep = sys.argv[2]
	file_name = sys.argv[3]

	antigen_name = file_name.split("_")[1]

	db_dir= '/tools2/abscan_db'
	ant_db = '%s/ant_%s' %(db_dir,DB)
	ant_pro_db='%s/ant_%s' %(db_dir,DB.replace('.','_pro.'))
	par_db= '%s/par_%s' %(db_dir,DB)
	par_pro_db= '%s/par_%s' %(db_dir,DB.replace('.','_pro.'))

	with open(ant_db,'r') as indb:
		ant_dic = json.load(indb)
		indb.close()

	with open(ant_pro_db,'r') as indb:
		ant_dic_pro = json.load(indb)
		indb.close()
	
	with open(par_db,'r') as indb:
		par_dic = json.load(indb)
		indb.close()

	with open(par_pro_db,'r') as indb:
		par_dic_pro = json.load(indb)
		indb.close()

	start_dic_ant = {}
	start_dic_par = {}
	outpep_dic_ant = {}
	outpep_dic_par = {}
	total_results =[]

	count = 0
	matching_count=0

	total_4mer_list=[]
	total_3mer_list=[]
	only_4mer=[]
	only_3mer=[]

	sys.setrecursionlimit(10**6)

	cd = os.getcwd()
	result_dir= '%s/%s/%s/cdrgen_result' %(cd,antigen_name, file_name)
	tmp_result_dir = '%s/%s/%s/cdrgen_result/tmp' %(cd,antigen_name, file_name)

	if not os.path.exists(result_dir):
		os.makedirs(result_dir)
		os.makedirs(tmp_result_dir) # intermediate result 

	cut(pep, 'ant')
	cut(pep, 'par')
	print("cut_end")
	run('ant')
	print("ant_end")
	run('par')
	print("par_end")

	del ant_dic
	del ant_dic_pro
	del par_dic
	del par_dic_pro
	del start_dic_ant
	del start_dic_par
	del outpep_dic_ant
	del outpep_dic_par

	total_uniq=list(set(total_results))
	only_4mer_uniq=list(set(only_4mer))
	only_3mer_uniq=list(set(only_3mer))

	total_result= '%s/%s_%s_E_total_result.txt' %(result_dir,antigen_name,pep)
	frag_4mer_result= '%s/%s_%s_4mer_result.txt' %(tmp_result_dir,antigen_name,pep)
	frag_3mer_result= '%s/%s_%s_3mer_result.txt' %(tmp_result_dir,antigen_name,pep)
	frag_only_4mer_result = '%s/%s_%s_E_only_4mer_result.txt' %(tmp_result_dir,antigen_name,pep)
	frag_only_3mer_result = '%s/%s_%s_E_only_3mer_result.txt' %(tmp_result_dir,antigen_name,pep)
	total_count= '%s/%s_%s_E_total_result_count.txt' %(tmp_result_dir,antigen_name,pep)

	with open(total_result,'w') as t_u_w:
		t_u_w.write("\n".join(total_uniq))
		
	with open(frag_4mer_result,'w') as r_4:
		r_4.write("\n".join(total_4mer_list))

	with open(frag_3mer_result,'w') as r_3:
		r_3.write("\n".join(total_3mer_list))

	with open(frag_only_4mer_result,'w') as only_4:
		only_4.write("\n".join(only_4mer_uniq))
	
	with open(frag_only_3mer_result,'w') as only_3:
		only_3.write("\n".join(only_3mer_uniq))
	
	with open(total_count,'w' ) as co:
		co.write(str(matching_count))

	print(pep,"end")
