import click
import os
import pandas as pd
from aux_functions import read_vcf, snvs_processing, indels_processing
from tqdm import tqdm
tqdm.pandas()

def sorting(df):
    """
    Sort DF by chromosome autosomal first then X and Y
    :param df: dataframe pandas
    :return: dataframe
    """

    chroms = [str(x) for x in range(1, 23, 1)]
    chroms.extend(['X', 'Y'])

    grps = df.groupby("#CHROM")

    df_to_return = pd.DataFrame()

    for c in chroms:
        try:
            dff = grps.get_group(c)
            df_to_return = df_to_return.append(dff, ignore_index=True)
        except KeyError as e:
            print(e)
    return df_to_return


@click.command()
@click.option('--output_path',
              '-out_path',
              type=click.Path(exists=True),
              required = True,
              help="Output path to write the filled mafs")
@click.option('--input_primary_maf',
              '-in_pry_maf',
              required=True,
              help="input primary maf to fill")
@click.option('--input_relapse_maf',
              '-in_rel_maf',
              required=True,
              help="input relapse maf to fill")
@click.option('--input_primary_vcf',
              '-pry_vcf',
              required=True,
              help="input primary vcf (raw calls)")
@click.option('--input_relapse_vcf',
              '-rel_vcf',
              required=True,
              help="input relapse vcf (raw calls)")
@click.option('--mutation_type',
              '-mut_type',
              required=True,
              help="options: snvs, indels")

def cli(output_path, input_primary_maf,input_relapse_maf,input_primary_vcf,input_relapse_vcf,mutation_type):
    """
    It gets processed maf files and recovers variants that are shared between primary and relapse samples but were
    filtered out in the previous processing step.
    """

    # read mafs
    maf_pry = pd.read_csv(input_primary_maf, sep='\t')
    maf_rel = pd.read_csv(input_relapse_maf, sep='\t')

    # read vcf
    vcf_pry = read_vcf(input_primary_vcf)
    vcf_rel = read_vcf(input_relapse_vcf)

    # create column 'Variant'
    list_df = [maf_pry, maf_rel, vcf_pry, vcf_rel]

    for i,df in enumerate(list_df):
        list_df[i]['Variant'] = list_df[i].apply(lambda x: str(x['#CHROM']) + '_' + str(x['POS']) + '_' + x['REF'] + '_' + x['ALT'],
                                       axis=1)

    # ensure no duplicates
    maf_pry.drop_duplicates(subset='Variant', inplace=True)
    maf_rel.drop_duplicates(subset='Variant', inplace=True)

    # get shared by intersection
    shared_variants = set(maf_pry['Variant'].tolist()).intersection(set(maf_rel['Variant'].tolist()))

    # define private variants
    private_pry_variants = set(maf_pry['Variant'].tolist()).difference(shared_variants)
    private_rel_variants = set(maf_rel['Variant'].tolist()).difference(shared_variants)

    # get the missed shared variants
    missing_pry_shared = vcf_pry[vcf_pry['Variant'].isin(private_rel_variants)]
    missing_rel_shared = vcf_rel[vcf_rel['Variant'].isin(private_pry_variants)]

    vcf_pry_missing = vcf_pry[vcf_pry['Variant'].isin(missing_pry_shared['Variant'])]
    vcf_rel_missing = vcf_rel[vcf_rel['Variant'].isin(missing_rel_shared['Variant'])]

    if mutation_type == 'snvs':
        df_pry = vcf_pry_missing.progress_apply(lambda x: snvs_processing(x), axis=1)
        df_rel = vcf_rel_missing.progress_apply(lambda x: snvs_processing(x), axis=1)
    elif mutation_type == 'indels':
        df_pry = vcf_pry_missing.apply(lambda x: indels_processing(x), axis=1)
        df_rel = vcf_rel_missing.apply(lambda x: indels_processing(x), axis=1)
    else:
        print("Write 'snvs' or 'indels' depending on the mut type")

    # append the recovered variants
    maf_pry = maf_pry.append(df_pry, ignore_index=True, sort=False)
    maf_rel = maf_rel.append(df_rel, ignore_index=True, sort=False)

    # sort them by chromosome
    maf_pry['#CHROM'] = maf_pry['#CHROM'].astype(str)
    maf_rel['#CHROM'] = maf_rel['#CHROM'].astype(str)

    maf_pry = sorting(maf_pry)
    maf_rel = sorting(maf_rel)

    # ensure data types and order of columns
    int_cols = ['n_alt_reads', 'n_ref_reads',
     't_alt_reads', 't_ref_reads', 'POS']

    order_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                     'NORMAL', 'TUMOR', 'DP_tumor', 't_alt_reads', 't_ref_reads', 'DP_normal',
                     'n_alt_reads', 'n_ref_reads', 'mut_type', 'GT_normal', 'GT_tumor', 'Variant']

    maf_pry[int_cols] = maf_pry[int_cols].astype(int)
    maf_rel[int_cols] = maf_rel[int_cols].astype(int)

    maf_pry = maf_pry[order_cols]
    maf_rel = maf_rel[order_cols]

    # create output file name
    file_pry = input_primary_maf.split("/")[-1]
    file_pry = file_pry.replace(".maf", "")
    comparison_pry = input_primary_maf.split("/")[-2]

    file_rel = input_relapse_maf.split("/")[-1]
    file_rel = file_rel.replace(".maf", "")
    comparison_rel = input_relapse_maf.split("/")[-2]

    print(len(maf_pry))
    print(len(maf_rel))

    # write results
    maf_pry.to_csv(os.path.join(output_path, comparison_pry, file_pry+"_sh.maf"), sep='\t', index=False)
    maf_rel.to_csv(os.path.join(output_path, comparison_rel, file_rel+"_sh.maf"), sep='\t', index=False)

if __name__ == '__main__':
    cli()
