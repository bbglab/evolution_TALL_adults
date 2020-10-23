import sys, os
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
import click
import subprocess
import pandas as pd

@click.command()
@click.option('--output_path',
              '-out',
              type=click.Path(exists=True),
              required = True,
              help="Output path to write results")
@click.option('--normal_bam',
              '-nor_bam',
              type=click.Path(exists=True),
              required = True,
              help="BAM of the normal samples (or remission)",
              )
@click.option('--tumor_bam',
              '-tum_bam',
              type=click.Path(exists=True),
              required = True,
              help="BAM of the tumor samples (primary/relapse)",
              )
@click.option('--comparison',
              '-com',
              type=str,
              required = True,
              help="comparison as tumorid_vs_normalid",
              )
@click.option('--genome_reference_fasta',
              '-ref',
              type=click.Path(exists=True),
              required = True,
              help="FASTA file of the reference genome used in the aligments",
              )
@click.option('--excludable_regions',
              '-excl',
              type=click.Path(exists=True),
              required = True,
              help="External file with excludable regions such as telomere and centromere regions (provided by delly)",
              )


def cli(output_path, normal_bam,tumor_bam, comparison, genome_reference_fasta, excludable_regions):
    """
    It runs Delly and svtools to obtain bedpe files with the structural variant calls

    """
    tumor_id = comparison.split("_vs_")[0]
    normal_id = comparison.split("_vs_")[1]
    sample_tsv = pd.DataFrame.from_dict({tumor_id: 'tumor',
                               normal_id: 'control'}, orient='index')
    sample_tsv.to_csv(os.path.join(output_path, 'sample.tsv'), sep='\t',header=False)

    sv_list = ['bnd', 'del', 'dup', 'ins', 'inv']

    for sv in sv_list:
        call_file = os.path.join(output_path,comparison+"_"+sv+".bcf")
        subprocess.run("source activate sv_delly_calling && delly call -t "+sv.upper()+" -q 20 -x " + excludable_regions
                       +' -o ' + call_file + " -g " + genome_reference_fasta + ' ' + tumor_bam + ' ' + normal_bam,
                       shell=True, executable='/bin/bash')
        filter_file = os.path.join(output_path, comparison + "_" + sv + "_filtered.bcf")
        sample_tsv = os.path.join(output_path, 'sample.tsv')
        if sv == 'bnd':
            subprocess.run("source activate sv_delly_calling && delly filter -p -f somatic -m 0 -r 0.75 -a 0.1 -o "+
                           filter_file+" -s "+sample_tsv+" "+call_file,  shell=True, executable='/bin/bash')
        elif (sv == 'del') or (sv == 'ins'):
            subprocess.run("source activate sv_delly_calling && delly filter -p -f somatic -m 50 -r 0.75 -a 0.1 -o "+
                           filter_file+" -s "+sample_tsv+" "+call_file,  shell=True, executable='/bin/bash')
        elif sv == 'inv':
            subprocess.run("source activate sv_delly_calling && delly filter -p -f somatic -m 0 -r 0.75 -a 0.1 -o "+
                           filter_file+" -s "+sample_tsv+" "+call_file,  shell=True, executable='/bin/bash')
        elif sv == 'dup':
            subprocess.run("source activate sv_delly_calling && delly filter -p -f somatic -m 0 -r 0.75 -a 0.1 -o "+
                           filter_file+" -s "+sample_tsv+" "+call_file,  shell=True, executable='/bin/bash')
        else:
            "wrong sv, write: bnd, inv, del, ins, dup"
        call_vcf = os.path.join(output_path,comparison+"_"+sv+"_delly.vcf")
        subprocess.run("source activate sv_delly_calling && bcftools view " +filter_file+" > "+call_vcf,
                       shell=True, executable='/bin/bash')
        bedpe = os.path.join(output_path,comparison+"_"+sv+"_delly.bedpe")
        subprocess.run("source activate sv_delly_calling && svtools vcftobedpe -i "+call_vcf+" -o "+bedpe,
                       shell=True, executable='/bin/bash')

if __name__ == '__main__':
    cli()