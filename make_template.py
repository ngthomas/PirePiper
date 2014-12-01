#!/usr/bin/env python

import os, sys, re, argparse, string
import jtemplate

'''
This module reads the user's input config file and extracts information from it
'''
def read_config(fileName):
	configFile = open(fileName, "r")
	config=dict()
	for line in configFile.read().splitlines():
		if re.match("^#", line):
			next
		str_sep =re.split("=| ", line)
		#If user provides the sample label, the samples are separated out into a list
		if str_sep[0] == "sample":
			config[str_sep[0]]=str_sep[1].split(",")
		else:
			config[str_sep[0]]=str_sep[1]
		#config.setdefault(str_sep[0],[]).append(str_sep[1])
	configFile.close()
	return config
	
'''
This module makes sure all of the config's labels are present; if not, the system halts
'''
def check_config(config_param):
	for label in ['basePath', 'species', 'runN','sample']:
		if label in config_param:
			print "Provided ", label, " : ", config_param[label]
		else:
			print "Missing this label:", label
			sys.exit(1)

def check_or_make_dir(path):
	if not os.path.exists(path):
		os.makedirs(path)


'''
This module creates folders to hold fastq files, downstream results, and scripts. 
'''
def make_folders(config_param):

	path_pre = config_param['basePath'] +"/" + config_param['species']+"/" +config_param['runN']
			
	for sample in config_param['sample']:
		check_or_make_dir(path_pre+"/data/"+sample+"/raw") 
		check_or_make_dir(path_pre+"/data/"+sample+"/trim")
		check_or_make_dir(path_pre+"/data/"+sample+"/rna")
		check_or_make_dir(path_pre+"/data/"+sample+"/qc_report") 
		
		for i in ['pyRAD', 'Stacks', 'R', 'blat']:
			check_or_make_dir(path_pre+"/analysis/"+sample+"/"+i)
			check_or_make_dir(path_pre+"/scripts/"+sample+"/"+i)

def generate_scripts(config_param):

	path_pre = config_param['basePath'] +"/" + config_param['species']+"/" +config_param['runN']
			
	for sample in config_param['sample']:
		script_path = path_pre + "/scripts/"+sample
		rawfq_path = path_pre +"/data/"+sample+"/raw"
		trim_path = path_pre+"/data/"+sample+"/trim"
		barcode_path = path_pre +"/data/"+sample+"/barcode"
		
		print ("")
		print ("------------------------")
		print ("For "+sample+":")
		print ("3a. Transfer your fastq files to: " + rawfq_path)
		print ("3b. Edit barcode file (barcode w/ sample Name): vi " + barcode_path)
		print ("")
		print ("4a. To run quality summary testing: qsub "+script_path+"/run_qc.sh")
		print ("4b. To trim 5bp from the fastq read: qsub "+script_path+"/run_trim.sh")
		print ("")
		print ("5a. To demultiplex reads Using Stacks: qsub "+script_path+"/Stacks/run_processTags.sh")
		print ("5b. To call loci Using Stacks: qsub "+script_path+"/Stacks/run_stackCatalog.sh")
		print ("")
		print ("6a. To optimize stacks parameter: qsub ??")
		print ("6b. To generate loci quality caller: qsub ??")
		print ("6c. To align RNA to RAD: qsub "+script_path+"/blat/run_blat.sh")
		print ("")
		print ("7. To run pyRAD: qsub " +script_path+"/pyRAD/run_pyRAD.sh")
		print ("")
		print ("____________________________________________________________________")
		print ("")
		
		jtemplate.write_qc(sample, path_pre, config_param)
		jtemplate.write_trim(sample, path_pre, config_param)
		jtemplate.write_demulti(sample, path_pre, config_param)
		jtemplate.write_stack_core(sample, path_pre, config_param)
		jtemplate.write_blat(sample, path_pre, config_param)
		jtemplate.write_pyrad(sample, path_pre, config_param)
		jtemplate.write_align_summary(sample, path_pre, config_param)
	
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='create a scaffold to run PIRE RAD-seq pipeline')
	parser.add_argument('-c','--config',required=True, default=None,type=str,help='path of configuration file')
	opts = parser.parse_args()
	
	config_param=read_config(opts.config)
	check_config(config_param)
	make_folders(config_param)
	generate_scripts(config_param)

