import click
import os
from aux_functions import read_vcf, snvs_processing, indels_processing
from tqdm import tqdm
tqdm.pandas()


@click.command()
@click.option('--input_vcf',
              '-i',
              type=click.Path(exists=True),
              required = True,
              help="Input must be the path to the vcf")
@click.option('--vcf_type',
              '-vcf_t',
              type=str,
              required = True,
              help="VCF type: meaning whether it's indels or snvs vcf")
@click.option('--output_dire',
              '-o',
              type=click.Path(exists=True),
              required = True,
              help="Output directory to write the processed results")


def cli(input_vcf, vcf_type, output_dire):
    """
    It processes the results of Strelka and extracts information of the calls to return processed VCF in the form of
    TSV file
    """
    # read strelka results
    df = read_vcf(input_vcf)

    # first filters
    chroms = [str(x) for x in range(1, 23, 1)]
    chroms.extend(['X', 'Y'])
    df['#CHROM'].astype(str)
    df = df[df["#CHROM"].isin(chroms)]

    df = df[(df.FILTER == 'PASS') | (df.FILTER == 'DP')]

    # mut type processing
    if vcf_type == 'snvs':
        df_proc = df.progress_apply(lambda x: snvs_processing(x), axis=1)
    elif vcf_type == 'indels':
        df_proc = df.apply(lambda x: indels_processing(x), axis=1)
    else:
        print("Wrong vcf_type string. Write 'snvs' or 'indels' depending on the mut type")

    # write results
    df_proc[['n_alt_reads', 'n_ref_reads',
        't_alt_reads', 't_ref_reads', 'POS']] = df_proc[['n_alt_reads', 'n_ref_reads', 't_alt_reads',
                                                    't_ref_reads', 'POS']].astype(int)
    df_proc = df_proc[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'NORMAL', 'TUMOR', 'DP_tumor',
             't_alt_reads', 't_ref_reads', 'DP_normal', 'n_alt_reads', 'n_ref_reads', 'mut_type', 'GT_normal',
             'GT_tumor']]

    file = input_vcf.split("/")[-1]
    file = file.replace(".vcf.gz", "")

    print(len(df_proc))

    df_proc.to_csv(os.path.join(output_dire, file + ".maf"), sep='\t', index=False)


if __name__ == '__main__':
    cli()