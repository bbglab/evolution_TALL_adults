
#!/bin/bash


if [ ${#@} == 0 ] || [ $1 == "-h" ]; then
    echo "Usage: $0 param1 [param2, param3]"
    echo "* param1: <infile normal bam file>"
    echo "* param2: <infile tumor bam file>"
    echo "* param3: <file that snp-pileup creates>"
    echo "* param4: <output dir of FACETS results>"
    echo "* param5: <random seed number>"
else
    cd $4
    if [ -e "$3" ]; then
        echo "File exists"
    else 
        echo "File does not exist"
        snp-pileup -g ../../ext_files/run_FACETS/00-common_all.vcf -q15 -Q20 -P100 -r10,0 $3 $1 $2
    fi 
    /home/isentis/.conda/envs/copynumber/bin/Rscript --vanilla facets_script_wgs_hg19.R --input $3 --output $4 --seed $5
fi
 

