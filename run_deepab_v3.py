import os,sys
import shutil
import subprocess as sp




fasta_file= sys.argv[1]
result_file =sys.argv[2]
sam1 = sys.argv[3] # cpu or gpu
sam2 = sys.argv[4] # if sam1 == 'gpu': 0,1,.... if sam1=='cpu': -

GEAR = '/lwork01/neoscan_gear'

cd =os.getcwd()
antigen_name = fasta_file.split("_")[0]
deepab_dir = '%s/%s/%s/DeepAb' %(cd,antigen_name,result_file)
ABS_model = '/KHIT1/ABS_model'

deepab_result_dir = '%s/result' %(deepab_dir)
tmp_deepab_dir = '/lwork01/DeepAb_tmp'

if not os.path.exists(deepab_result_dir):
	os.makedirs(deepab_result_dir,exist_ok=True)
if not os.path.exists(tmp_deepab_dir):
	os.makedirs(tmp_deepab_dir,exist_ok=True)

shutil.copy('%s/input/%s.fasta'%(deepab_dir,fasta_file),'%s/%s.fasta'%(tmp_deepab_dir,fasta_file))

if sam1 == 'gpu':
	os.environ['CUDA_VISIBLE_DEVICES'] = sam2	
	py_cmd="time python %s/predict.py %s/%s.fasta --model_dir %s/trained_models/ensemble_abresnet --target %s --pred_dir %s/%s --decoys 10 --renumber --use_gpu" %(GEAR,tmp_deepab_dir,fasta_file,ABS_model,fasta_file,tmp_deepab_dir,fasta_file)
#	print (py_cmd)
elif sam1 == 'cpu':
	py_cmd="time python %s/predict.py %s/%s.fasta --model_dir %s/trained_models/ensemble_abresnet --target %s --pred_dir %s/%s --decoys 10 --renumber" %(GEAR,tmp_deepab_dir,fasta_file,ABS_model,fasta_file,tmp_deepab_dir,fasta_file)
	
sp.call(py_cmd, shell= True)
shutil.copytree('%s/%s'%(tmp_deepab_dir,fasta_file),'%s/%s'%(deepab_result_dir,fasta_file))
shutil.rmtree('%s/%s'%(tmp_deepab_dir,fasta_file))
