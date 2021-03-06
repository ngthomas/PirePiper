def write_qc (sample, path_pre, config_param):
	script_path = path_pre + "/scripts/"+sample
	rawfq_path = path_pre +"/data/"+sample+"/raw"
		
	qc_FILE = open (script_path+"/run_qc.sh", "w")
	qc_script = """#!/bin/bash
#$ -cwd
#$ -V 
#$ -N fastq_QC 
#$ -pe shared 1
#$ -l highp,h_data=4G,time=33:00:00 
#$ -m bea

## path: """+script_path+"""
## Usage: Submit this script thru qsub: qsub run_qc.sh

cd """+rawfq_path+"""
cat *fastq* | fastqc /dev/stdin -o ../qc_report/
"""
	qc_FILE.write(qc_script)

def write_trim (sample, path_pre, config_param):
	script_path = path_pre + "/scripts/"+sample
	rawfq_path = path_pre +"/data/"+sample+"/raw"
	trim_path = path_pre+"/data/"+sample+"/trim"
	barcode_path = path_pre +"/data/"+sample+"/barcode"
	
	trim_FILE = open (script_path+"/run_trim.sh", "w")
	trim_script = """#!/bin/bash
#$ -cwd
#$ -V 
#$ -N trim_reads 
#$ -pe shared 1
#$ -l highp,h_data=4G,time=24:00:00 
#$ -m bea

## path: """+script_path+"""
## Usage: Submit this script thru qsub: qsub run_trim.sh

cd """+rawfq_path+"""

for fastq_f in *fastq*; do
	outfile=${fastq_f}
	if [[ ${outfile} == *.gz ]]; then
		outfile=`echo ${outfile} | rev| cut -c 4- |rev`; fi
		
	seqtk trimfq -e 5 ${fastq_f} > """+trim_path+"""/$outfile; 
done 

"""
	trim_FILE.write(trim_script)

def write_demulti (sample, path_pre, config_param):
	script_path = path_pre + "/scripts/"+sample
	trim_path = path_pre+"/data/"+sample+"/trim"
	barcode_path = path_pre +"/data/"+sample+"/barcode"
	
	stack_demulti_FILE = open (script_path+"/Stacks/run_processTags.sh", "w")
	stack_demulti = """#!/bin/bash
#$ -cwd
#$ -V 
#$ -N stack_process 
#$ -pe shared 1
#$ -l highp,h_data=4G,time=33:00:00 
#$ -m bea

## path: """+script_path+"""/Stacks
## Usage: Submit this script thru qsub: qsub run_processTags.sh

# base directory
WKDIR="""+path_pre+"""
mkdir -p ${WKDIR}/analysis/"""+sample+"""/Stacks/processTags

# making the stacks' barcode

awk '{print $1}' """+barcode_path+""" > ${WKDIR}/data/"""+sample+"""/RAD_barcode


# parameters for process_radtags
# -p: input directory - make sure you only have the fastq or fastq.gz files in the folder (no SampleSheet or anything else)
# -o: output directory 
# -b: barcode directory 
# -e: restriction enzyme 
# -r: rescue RAD-TAGS and barcode
# -c: clean data, remove uncalled base 
# -q: remove low quality reads 
# -i: input file type (gzfastq for gzip)

# -s (note: we need to figure out the best limit and add it below, or now we have set it to 20): set the score limit. If the average score within the sliding window drops below this value, the read is discarded (default 10).
process_radtags -p """+trim_path+""" -o ${WKDIR}/analysis/"""+sample+"""/Stacks/processTags -b ${WKDIR}/data/"""+sample+"""/RAD_barcode -e 'sbfI' -r -c -q 

# relabel fq barcode file to something more meaningful i.e. the sample name
awk -v basePath=${WKDIR}/analysis/"""+sample+"""/Stacks/processTags '{print "mv "basePath"/sample_"$1".fq "basePath"/sample_"$2".fq"}' ${WKDIR}/data/"""+sample+"""/barcode |bash """

	stack_demulti_FILE.write(stack_demulti)
	
def write_stack_core (sample, path_pre, config_param):
	script_path = path_pre + "/scripts/"+sample
	
	stack_catalog_FILE = open (script_path+"/Stacks/run_stackCatalog.sh", "w")
	stack_catalog = """#!/bin/bash
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

WKDIR="""+path_pre+"""

mkdir -p ${WKDIR}/analysis/"""+sample+"""/Stacks/denovo

# make a list of input files
inputL=`ls ${WKDIR}/analysis/"""+sample+"""/Stacks/processTags/*.fq | awk '{printf "-s "$1 " "}'`

echo "nohup denovo_map.pl -m 6 -M 2 -n 2 -S -b 1 -T 6 -t -o ${WKDIR}/analysis/"""+sample+"""/Stacks/denovo" $inputL " >${WKDIR}/analysis/"""+sample+"""/Stacks/denovo/nohup.out" | bash

"""	
	stack_catalog_FILE.write(stack_catalog)

def write_pyrad (sample, path_pre, config_param):
	script_path = path_pre + "/scripts/"+sample
	
	pyrad_FILE = open (script_path+"/pyRAD/run_pyRAD.sh", "w")
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
WKDIR=="""+path_pre+"""

module load python/2.7
pyRAD -p $WKDIR/scripts/"""+sample+"""/pyRAD/params.txt -s 234567
"""
	pyrad_FILE.write(pyrad_script)

def write_blat (sample, path_pre, config_param):
	script_path = path_pre + "/scripts/"+sample
	
	blat_FILE = open (script_path+"/blat/run_blat.sh", "w")
	blat_script = """#!/bin/bash
#$ -cwd
#$ -V 
#$ -N align_RNA_RAD
#$ -l highp,h_data=4G,time=46:00:00 
#$ -m bea

## path: """+script_path+"""/blat
###Usage: Submit this script thru qsub: qsub run_blat.sh

WKDIR="""+path_pre+"""

awk '{ split($8,con,","); for (i=1; i <= length(con); i++) {split(con[i], id, "_"); cons[id[1]]++} 
if (length(cons)>80) {a+=1; print ">" a "_" $3 "_" length(cons) "\\n" $9} delete cons}' ${WKDIR}/analysis/"""+sample+"""/Stacks/denovo/batch_1.catalog.tags.tsv > ${WKDIR}/analysis/"""+sample+"""/blat/RAD_consen.fa

cd ${WKDIR}/data/"""+sample+"""/rna
for i in *.fa; do
/u/local/apps/blat/34/bin/blat -q=rna -out=pslx ${WKDIR}/analysis/"""+sample+"""/blat/RAD_consen.fa $i ${WKDIR}/analysis/"""+sample+"""/blat/align_${i}.out; done
"""
	blat_FILE.write(blat_script)
	
def write_align_summary (sample, path_pre, config_param):
	script_path = path_pre + "/scripts/"+sample
	
	align_summary_FILE = open (script_path+"/blat/gen_tbl.sh", "w")
	script = """#!/bin/bash
	
## path: """+script_path+"""/blat
###Usage: Submit this script thru qsub: bash gen_tbl.sh

WKDIR="""+path_pre+"""
OUTDIR=${WKDIR}/analysis/"""+sample+"""/blat/tbl

mkdir -p ${OUTDIR}

cd ${WKDIR}/data/"""+sample+"""/rna
# grab read numbers
for i in *.fa; 
do numline=`wc -l $i|awk '{print $1/2}'`
echo $i $numline;
done > ${OUTDIR}/num_reads.data 

# grab number of reads aligned once and more than once
cd ${WKDIR}/analysis/"""+sample+"""/blat
for i in *.out; 
do numline=`awk ' {NR>5 && d[$10]++} 
    END{ 
    for(i in d) { 
        if(d[i]==1){ct_once++}
        else{ct_more++}}
    print ct_once "\\t" ct_more}' $i`;
    fileN=`echo $i | sed 's/align_//;s/.out//'`
echo $fileN $numline;
done > ${OUTDIR}/num_align.data 

cd ${OUTDIR}
join -1 1 -2 1 ${OUTDIR}/num_reads.data ${OUTDIR}/num_align.data |awk 'BEGIN{print "RNA File | # of reads | # of reads aligned once | # of reads aligned more than once"; print "---|:---:|:---:|:---:";} {OFS=" | "; CONVFMT="%0.3f"; print $1,$2, $3 " ( " $3*100/$2 " %)", $4 " (" $4*100/$2" % )" }'
> ${WKDIR}/analysis/"""+sample+"""/blat/alignment_report.md

cat ${WKDIR}/analysis/"""+sample+"""/blat/*.out > ${OUTDIR}/align.out.concat

for i in 30 50 70 90; do
awk '{if($5==0 && $2 < 3 && $1 > '$i'){if(e[$14"_"$16"_"$17]++==0) d[$14]++; }}
END{for (k in d) {print k"\\t"d[k] }}' ${OUTDIR}/align.out.concat > ${OUTDIR}/align_concat_${i}.tbl; done 

for i in 30 50 70 90; do
echo -n ">"$i" bp matches ";
for d in 1 5 10 50 100; do
awk '$2>'$d'' ${OUTDIR}/align_concat_${i}.tbl |wc -l |awk '{printf " | "$1 " (" "%3.1f" "%)",$1*100/38070}'; 
done
echo "";
done > ${WKDIR}/analysis/"""+sample+"""/blat/RAD_RNA_stat.md

"""

	align_summary_FILE.write(script)
