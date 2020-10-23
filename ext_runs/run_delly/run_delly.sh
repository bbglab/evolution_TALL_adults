#!/bin/bash


if [ ${#@} == 0 ] || [ $1 == "-h" ]; then
    echo "Usage: $0 param1 [param2, param3, param4]"
    echo "* param1: <normal BAM path>"
    echo "* param2: <tumor BAM path>"
    echo "* param3: <comparison tumor_id_vs_normal_id>"
    echo "* param4: <reference genome fasta file>"
    echo "* param5: <excludable regions from the genome build of the fasta files>"
    echo "* param6: <enter type of SV. Options: BND, INV, DEL, DUP, INS>"
    echo "* param7: <output or working directory to write results>"
else
    source activate sv_delly_calling
    echo -e
    if [ "$6" == "BND" ]; then
        delly call -t BND -q 20 -x $5 -o $3'_bnd.bcf' -g $4 $1 $2
        delly filter -p -f somatic -m 0 -r 0.75 -a 0.1 -o $3'_bnd.bcf' -s $7