#!/usr/bin/env python

import os, sys, re, argparse, string

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

'''
This module creates folders to hold fastq files, downstream results, and scripts. 
'''
def make_folders(config_param):

	path_pre = config_param['basePath'] +"/" + config_param['species']+"/" +config_param['runN']
			
	for sample in config_param['sample']:
		os.makedirs(path_pre+"/data/"+sample+"/raw") 
		os.makedirs(path_pre+"/data/"+sample+"/trim") 
		os.makedirs(path_pre+"/data/"+sample+"/qc_report") 
		
		for i in ['pyRAD', 'Stacks', 'R']:
			os.makedirs(path_pre+"/analysis/"+sample+"/"+i)
			os.makedirs(path_pre+"/scripts/"+sample+"/"+i)

def generate_scripts(config_param):

	script_path = config_param['basePath'] + "/" + config_param['species']+ "/" +config_param['runN'] + "/scripts" 
	path_pre = config_param['basePath'] +"/" + config_param['species']+"/" +config_param['runN']
			
	for sample in config_param['sample']:
		print ("")
		print ("------------------------")
		print ("For "+sample+":")
		print ("3a. Transfer your fastq files to: " + path_pre +"/data/"+sample+"/raw")
		print ("3b. Edit barcode file (barcode w/ sample Name): vi " + path_pre +"/data/"+sample+"/barcode")
		print ("")
		print ("4a. To run quality summary testing: qsub "+script_path+"/"+sample+"/run_qc.sh")
		print ("4b. To trim 5bp from the fastq read: qsub "+script_path+"/"+sample+"/run_trim.sh")
		print ("")
		print ("5a. To demultiplex reads Using Stacks: qsub "+script_path+"/"+sample+"/Stacks/run_processTags.sh")
		print ("5b. To call loci Using Stacks: qsub "+script_path+"/"+sample+"run_stackCatalog.sh")
		print ("")
		print ("6a. To optimize stacks parameter: qsub ??")
		print ("6b. To generate loci quality caller: qsub ??")
		print ("")
		print ("7. To run pyRAD: qsub " +script_path+"/"+sample+"/pyRAD/run_pyRAD.sh")
		print ("")
		print ("____________________________________________________________________")
		print ("")
	
	
	stack_script1_file = open (script_path+"/Stacks/run_processTags.sh", "w")
	stack_script1 = """#!/bin/bash
#$ -cwd
#$ -V 
#$ -N stack_process 
#$ -pe shared 1
#$ -l highp,h_data=4G,time=33:00:00 
#$ -m bea

## path: """+script_path+"""/Stacks
## Usage: Submit this script thru qsub: qsub run_processTags.sh

# base directory
WKDIR="""+config_param['basePath'] + "/" +config_param['species']+ "/" +config_param['runN']+"""
"""+string.join(["""mkdir -p ${WKDIR}/analysis/Stacks/processTags/"""+x for x in config_param['sample'] ],"\n")+"""

# making the stacks' barcode

"""+string.join(["""awk '{print $1}' ${WKDIR}/data/"""+x+"""/barcode > ${WKDIR}/data/"""+x+"""/RAD_barcode """ for x in config_param['sample'] ],"\n")+"""


# parameters for process_radtags
# -p: input directory - make sure you only have the fastq or fastq.gz files in the folder (no SampleSheet or anything else)
# -o: output directory 
# -b: barcode directory 
# -e: restriction enzyme 
# -r: rescue RAD-TAGS and barcode
# -c: clean data, remove uncalled base 
# -q: remove low quality reads 
# -i: input file type 

# -s (note: we need to figure out the best limit and add it below, or now we have set it to 20): set the score limit. If the average score within the sliding window drops below this value, the read is discarded (default 10).
"""+string.join(["""process_radtags -p ${WKDIR}/data/"""+x+""" -o ${WKDIR}/analysis/Stacks/processTags/"""+x+""" -b ${WKDIR}/data/"""+x+"""/RAD_barcode -e 'sbfI' -r -c -q -i gzfastq""" for x in config_param['sample'] ],"\n")+"""
# relabel fq barcode file to something more meaningful i.e. the sample name
"""+string.join(["""awk -v basePath=${WKDIR}/analysis/Stacks/processTags/"""+x+""" 'FS="," {if(NR>1) {print "mv "basePath"/sample_"$2".fq "basePath"/sample_"$1".fq"}}' ${WKDIR}/data/"""+x+"""/barcode |bash """ for x in config_param['sample'] ],"\n")

	stack_script1_file.write(stack_script1)
	
	stack_script2_file = open (script_path+"/Stacks/run_stackCatalog.sh", "w")
	stack_script2 = """#!/bin/bash
#$ -cwd
#$ -V 
#$ -N stackCatalog_process 
#$ -pe shared 6
#$ -l highp,h_data=4G,time=96:00:00 
#$ -m bea

## path: """+script_path+"""/Stacks
###Usage: Submit this script thru qsub: qsub run_stackCatalog.sh

# -m: min read depth 
# -M: max bp mismatch between loci for an indiv 
# -N: max bp mismatch allowed secondary reads to primary stack (default:M+2)
# -n: max bp mismatch between loci when building catalog
# -S: disable mysql recording 
# -b: batch id #
# -T: number of threads
# -t: remove highly rep seq 
# -s: path for fastq indiv samples
# -o: output path

WKDIR="""+config_param['basePath'] + "/" +config_param['species']+ "/" +config_param['runN']+"""

mkdir -p ${WKDIR}/analysis/Stacks/denovo

# make a list of input files
inputL=`ls ${WKDIR}/analysis/Stacks/processTags/*/*.fq | awk '{printf "-s "$1 " "}'`

echo "nohup denovo_map.pl -m 6 -M 2 -n 2 -S -b 1 -T 6 -t -o ${WKDIR}/analysis/Stacks/denovo" $inputL " >${WKDIR}/analysis/Stacks/denovo/nohup.out" | bash

"""	

	pyrad_script_file = open (script_path+"/pyRAD/run_pyRAD.sh", "w")
	pyrad_script = """#!/bin/bash
#$ -cwd
#$ -V 
#$ -N pyrad_process 
#$ -pe shared 4
#$ -l highp,h_data=4G,time=96:00:00 
#$ -m bea

## path: """+script_path+"""/pyRAD
###Usage: Submit this script thru qsub: qsub run_pyRAD.sh

#run all of the steps 1 to 7
WKDIR="""+config_param['basePath'] + "/" +config_param['species']+ "/" +config_param['runN']+"""

module load python/2.7
pyRAD -p $WKDIR/scripts/pyRAD/params.txt -s 234567

"""
	pyrad_script_file.write(pyrad_script)
	
	
	
	


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='create a scaffold to run PIRE RAD-seq pipeline')
	parser.add_argument('-c','--config',required=True, default=None,type=str,help='path of configuration file')
	opts = parser.parse_args()
	
	config_param=read_config(opts.config)
	check_config(config_param)
	make_folders(config_param)
	generate_scripts(config_param)

