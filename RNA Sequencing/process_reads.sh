#!/bin/bash
echo -e "\n$(date +"%Y-%m-%d %T") DECLARING PARAMETERS"
echo -e "  Processing files:\n    $1\n    $2"
echo -e "  Using $3 threads and $4 MB RAM"
echo -e "  Working in $5"
echo -e "  Saving results to $6"
echo -e "  bowtie2 genome index prefix is $7"
echo -e "  Annotations (.gtf) are $8"
echo -e "  featureCounts will use $9"
echo -e "  Unique tag for this run is ${10}\n"
echo -e "$(date +"%Y-%m-%d %T") PREPARING ASSETS"

fr1=$(basename $1 .fastq.gz)
fr2=$(basename $2 .fastq.gz)

module load java
module load bowtie2
module load samtools

echo -e "$(date +"%Y-%m-%d %T") BEGINNING EXECUTION"

# copy compressed reads to node's local scratch space
command="cp $1 $5/${fr1}.fastq.gz"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command
command="cp $1 $5/${fr2}.fastq.gz"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command

# move to working directory
command="cd $5"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command

# unzip reads
command="gzip -d --force ${fr1}.fastq.gz"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command
command="gzip -d --force ${fr2}.fastq.gz"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command

# perform fastqc on raw reads
command="fastqc -t $3 ${fr1}.fastq ${fr2}.fastq"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command

# copy fastqc files to save directory

# remove read pairs w average qualtiy < 10
command="/bin/sh bbduk.sh in1=${fr1}.fastq in2=${fr2}.fastq out1=${fr1}c.fastq out2=${fr2}c.fastq maq=10 stats=${fr1}.bbdukstats"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command

# align reads to reference genome using bowtie2
command="bowtie2 --ff --threads $3 -x $7 -1 $5/${fr1}c.fastq -2 $5/${fr2}c.fastq -S ${10}.sam"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command >> ${10}.bowtie2stats 2>&1


# count aligned features of interest with featureCounts from the subread package
command="featureCounts -p -T $3 -t $9 -g gene_id -a $8 -o ${10}.counts.txt ${10}.sam"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command

# copy results to save folder (save to avoid file purge)
command="cp ${fr1}_fastqc.zip $6/${fr1}_fastqc.zip"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command
command="cp ${fr2}_fastqc.zip $6/${fr2}_fastqc.zip"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command
command="cp ${10}.bowtie2stats $6/${10}.bowtie2stats"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command
command="cp ${10}.counts.txt $6/${10}.counts.txt"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command
command="cp ${10}.counts.txt.summary $6/${10}.counts.txt.summary"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
$command
# diagnostic; save aligned reads
#command="cp ${10}.sam $6/${10}.sam"
#echo "$(date +"%Y-%m-%d %T") RUNNING $command"
#$command

# these copied files don't need to be deleted from the node; the node
# will purge these files once the job is complete


# NOT NECESSARY
# delete leftover read and alignment map files
command="rm -rf ${10}.sam"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
#$command
# wildcard .fastq* to remove the gzipped reads too
command="rm -rf ${fr1}.fastq*"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
#$command
command="rm -rf ${fr2}.fastq*"
echo "$(date +"%Y-%m-%d %T") RUNNING $command"
#$command


