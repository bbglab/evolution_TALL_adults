#!/bin/bash


if [ ${#@} == 0 ] || [ $1 == "-h" ]; then
    echo "Usage: $0 param1 [param2, param3, param4]"
    echo "* param1: <VCF file to annotate>"
    echo "* param2: <Output tab file with the annotations>"
    echo "* param3: <path to the vep cache files to pass to the --dir vep command>"
    echo "* param4: <path to the Gnomad file to add population information>"
else
    source activate vep92
    vep -i $1 -o STDOUT -tab --assembly GRCh37 --no_stats --cache --symbol --protein --canonical --offline --af_1kg --dir $3 --custom $4,gnomADg,vcf,exact,0,AF,NFE > $2
fi
