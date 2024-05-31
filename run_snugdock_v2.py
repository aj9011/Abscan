import os,sys
import shutil

file_name =sys.argv[1] 
antigen_name = file_name.split("_")[1]
sam = sys.argv[2] # pdb 
sam1 = sys.argv[3]  # nconf
pipeline = sys.argv[4] # pipeline

cd = os.getcwd()

snugdock_dir= '%s/%s/%s/%s/snugdock_input' %(cd,antigen_name,file_name,pipeline)
#sam_dir='%s/%s.LHA_sorted.pdb' %(snugdock_dir,sam)
tmp = '/lwork01/%s_abs_tmp_data'%(file_name)
os.makedirs(tmp,exist_ok =True)

ROSETTA_BIN = '/lwork01/rosetta_src_2019.40.60963_bundle/main/source/bin'
os.makedirs('%s/%s'%(tmp,sam),exist_ok =True)
os.makedirs('%s/%s/dock_res'%(tmp,sam),exist_ok=True)
os.chdir(tmp+'/'+sam)
shutil.copy('%s/%s.LHA_sorted.pdb'%(snugdock_dir,sam),'%s/%s'%(tmp,sam))

with open('%s_prepack_option' %(sam), 'w') as w:
	w.write('-in:file:s %s.LHA_sorted.pdb\n' %(sam))
	w.write('-ex1\n')
	w.write('-ex2aro\n')
	w.write('-partners LH_A\n')
	w.write('-docking:dock_rtmin\n')
	w.write('-docking:sc_min\n')
#	w.write('-out:path:pdb %s/pdb_data/%s' %(snugdock_dir,sam))

if os.path.exists(ROSETTA_BIN + '/docking_prepack_protocol.mpi.linuxgccrelease'):
	os.system(ROSETTA_BIN + '/docking_prepack_protocol.mpi.linuxgccrelease  @%s_prepack_option > %s_prepack_log' %(sam,sam))

with open('%s_snugdock_option' %(sam), 'w') as w2:
	w2.write('-s %s.LHA_sorted_0001.pdb\n' %(sam))
	w2.write('-antibody:auto_generate_kink_constraint\n')
	w2.write('-antibody:all_atom_mode_kink_constraint\n')
	w2.write('-partners LH_A\n')
	w2.write('-spin\n')
	w2.write('-dock_pert 3 8\n')
	w2.write('-ex1\n')
	w2.write('-ex2aro\n')
	w2.write('-nstruct %s\n'%(sam1))
	w2.write('-out:path:pdb dock_res\n')

if os.path.exists(ROSETTA_BIN + '/snugdock.mpi.linuxgccrelease'):
	os.system('mpirun -n 10 ' + ROSETTA_BIN + '/snugdock.mpi.linuxgccrelease @%s_snugdock_option > %s_run_log' %(sam,sam))
