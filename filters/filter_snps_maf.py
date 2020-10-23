import sys, os
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
import click
import pandas as pd
import numpy as np
import pybedtools

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from aux_data_in_pyvar import config_rcparams
from aux_functions import stage_mapping


def strelka_join_refilter(snvs_file, indels_file, mnvs_file):
    """
    Join all variants into one data frame and apply ultimate filters
    :param path: path to the input variant type files
    :param pat: patient
    :param com: alignment comparison (tumorid_vs_normalid) calling
    :return: Unique joined variant file
    """

    # read all variant type

    df_snvs = pd.read_csv(snvs_file, sep='\t')
    df_indels = pd.read_csv(indels_file, sep='\t')

    # join them
    df = df_snvs.append(df_indels, ignore_index=True)

    if mnvs_file != None: # no MNVs for some patients
        df_mnvs = pd.read_csv(mnvs_file, sep='\t')
        df = df.append(df_mnvs, ignore_index=True)

    # filter out rows that have weird variants. Those are variants that are either present in normal BAM or
    # that have 0 reads reported for the alternative allele (??)
    df = df[df['n_alt_reads'] < 2]
    df = df[df['t_alt_reads'] != 0]

    # finally after filtering multiple annotations for each variant randomly remove duplicates to ensure unique
    # row for each variant
    df['#CHROM'] = df['#CHROM'].astype(str)
    df['Variant'] = df.apply(lambda x: str(x['#CHROM'])+'_'+str(x['POS'])+'_'+x['REF']+'_'+x['ALT'], axis=1)
    df.drop_duplicates(subset='Variant', inplace=True, keep='first')
    return df


def check_chrom_m(rw):
    if (rw['#CHROM'] == 'X') and (rw['cnv_facets'] == 2):
        rw['cnv_facets'] = 1
    for chrom in range(1,23,1):
        if rw['#CHROM'] == str(chrom) and (rw['cnv_facets'] not in range(0,20,1)):
            rw['cnv_facets'] = 2
    return rw


def check_chrom_f(rw):
    for chrom in range(1,23,1):
        if rw['#CHROM'] == str(chrom) and (rw['cnv_facets'] not in range(0,10,1)):
            rw['cnv_facets'] = 2
        elif rw['#CHROM'] == "X" and (rw['cnv_facets'] not in range(0,10,1)):
            rw['cnv_facets'] = 2
        else:
            pass
    return rw


def add_cnv_facets(df_maf, com, clinical, cnv_file):

    """
    Add the copy number information to each variant
    :param df_maf: dataframe with variants
    :param com: comparison tumorid_vs_normalid that serves to organize data and folders
    :param clinical: clinical data with the sex of the patients
    :param cnv_file: copy number file with the results of FACETS
    :return: dataframe of variants with the copy number of each one added
    """

    # clinical subset
    sam = com.split("_vs_")[0]
    clinical_subset = clinical[clinical['SAMPLE'] == sam]
    clinical_subset.reset_index(inplace=True)
    # read facets results
    df_reg = pd.read_csv(cnv_file,  sep='\t')
    df_reg['chrom'] = df_reg['chrom'].apply(lambda x: "X" if x == 23 else str(x))

    col_list = list(df_maf.columns)
    df_maf['START'] = df_maf['POS'] - 1
    df_maf['END'] = df_maf['POS']

    df_maf[['#CHROM']] = df_maf[['#CHROM']].astype(str)
    df_reg[['chrom']] = df_reg[['chrom']].astype(str)

    df_maf[['START', 'END']] = df_maf[['START', 'END']].astype(int)
    df_reg[['start', 'end']] = df_reg[['start', 'end']].astype(int)

    # merge mutations with copy number segments
    maf = pybedtools.BedTool.from_dataframe(df_maf[['#CHROM', 'START', 'END']])
    reg = pybedtools.BedTool.from_dataframe(df_reg[['chrom', 'start', 'end']])
    result = maf.intersect(reg, loj=True)
    results = pd.read_table(result.fn, names=['#CHROM', 'START', 'END', 'chrom', 'start', 'end'])

    results[['#CHROM', 'chrom']] = results[['#CHROM', 'chrom']].astype(str)
    results[['START', 'END', 'start', 'end']] = results[['START', 'END', 'start', 'end']].astype(int)

    # restore the original columns and add them to the merged dataframe (results)
    df_maf = df_maf.merge(results, how='left', on=['#CHROM', 'START', 'END'])
    df_maf = df_maf.merge(df_reg, how='left', on=['chrom', 'start', 'end'])
    df_maf.drop_duplicates(inplace=True, keep='first')
    df_maf.rename(columns={'tcn.em': 'cnv_facets'}, inplace=True)
    col_list.append('cnv_facets')
    df_maf = df_maf[col_list]
    df_maf[['cnv_facets']] = df_maf[['cnv_facets']].astype(float)

    # correct minor CNV from sex chromosomes
    if clinical_subset.loc[0, 'SEX'] == 'M':
        df_maf.loc[df_maf['#CHROM'] == 'X', 'cnv_facets'] = df_maf.loc[df_maf['#CHROM'] == 'X', 'cnv_facets'].fillna(1)
        df_maf = df_maf.apply(lambda x: check_chrom_m(x), axis=1)
        df_maf.loc[df_maf['#CHROM'] == 'Y', 'cnv_facets'] = df_maf.loc[df_maf['#CHROM'] == 'Y', 'cnv_facets'].fillna(1)
    else:
        df_maf = df_maf.apply(lambda x: check_chrom_f(x), axis=1)
        df_maf = df_maf[df_maf['#CHROM'] != "Y"]
    df_maf['cnv_facets'] = df_maf['cnv_facets'].apply(lambda x: 1 if x == 0 else x)
    return df_maf


def categorize_muts(df, com, cnv_file, output_path_maf, clinical, pop_freq):

    """
    It classifies mutations by population frequencies and calls for function to add copy number to each variant
    :param df: dataframe with variants
    :param com: comparison tumorid_vs_normalid that serves to organize data and folders
    :param cnv_file: facets results with copy number segments
    :param file: file name with path of df to use it for writing other intermediate file versions of the input
    :param clinical: clinical data
    :return: writes intermediate data frames in the same path
    """

    # # less than 1% of the population
    # df['AF_less_001'] = df['gnomADg_AF'].apply(
    #     lambda x: 'yes' if x < 0.01 else ('no' if type(x) == float and x >= 0.01 else np.nan))

    df['AF_less_'+str(pop_freq)] = df['gnomADg_AF'].apply(
         lambda x: 'yes' if x < float(pop_freq) else ('no' if type(x) == float and x >= float(pop_freq) else np.nan))

    # get VAF
    df['VAF'] = df.apply(lambda x: x['t_alt_reads'] / (x['t_alt_reads'] + x['t_ref_reads']), axis=1)
    df.sort_values('VAF', ascending=True, inplace=True)
    df['bins'] = pd.cut(df['VAF'], bins=10).astype(str)

    # add copy number
    df = add_cnv_facets(df, com, clinical, cnv_file)

    # write

    out_f = os.path.join(output_path_maf , 'all_anno_vep92_categories.maf')
    out_f2 =  os.path.join(output_path_maf, 'all_anno_vep92_categories_filt_snps.maf')

    df.to_csv(out_f, sep='\t', index=False)
    df[df['AF_less_'+str(pop_freq)] == 'yes'].to_csv(out_f2, sep='\t', index=False)

    return df


def barplot_individual_filter(df, pat, com, out_path, pop_freq):
    """
    Makes stacked barplot with the filtering proportions
    :param df: dataframe with the categories
    :param pat: patient id
    :param com: comparison id (tumorid_vs_normalid) that serves to organize data and folders
    :param out_path: output path to create the plot
    :return: stacked barplot as pdf
    """

    # make list for the stacked bars
    j = 0

    bar_keep = list()
    bar_snps = list()
    r_pos = list()
    names = list()

    for b in df['bins'].unique():
        dff = df[df['bins'] == b]
        num_keep = len(set(dff[dff['AF_less_'+str(pop_freq)] == 'yes']['Variant'].unique()))
        bar_keep.append(num_keep)
        num_snps = len(set(dff[dff['AF_less_'+str(pop_freq)] == 'no']['Variant'].unique()))
        bar_snps.append(num_snps)

        names.append(b)
        r_pos.append(j)
        j = j + 1

    ## make plot

    config_rcparams()

    fig = plt.figure(figsize=(6, 6))

    outer = gridspec.GridSpec(1, 1, wspace=0, hspace=0)

    barWidth = 1

    ax = plt.subplot(outer[0, 0])

    ax.bar(r_pos, bar_keep, edgecolor='white', width=barWidth, label='keep', color='#bababa')
    ax.bar(r_pos, bar_snps, bottom=bar_keep, edgecolor='white', width=barWidth, label='population SNPs', color='#f4a582')

    # Custom axis
    plt.xticks(r_pos, names)
    plt.xlabel("VAF")
    plt.ylabel("Mutations")
    ax.tick_params(axis='both', which='major')

    plt.tight_layout()

    plt.xticks(rotation=90)
    ax.legend(prop={'size': 10})
    ax.set_title("Stacked Bar Plot mutations {} {}".format(pat, com), pad=30)

    fig.savefig(os.path.join(out_path, "barplot_{}_{}.pdf".format(pat, com)),
                dpi=300, bbox_inches='tight', orientation='portrait')
    plt.show()


@click.command()
@click.option('--input_snvs',
              '-in_snvs',
              required=True,
              help="path to the maf file with the SNVs")
@click.option('--input_indels',
              '-in_indels',
              required=True,
              help="path to the maf file with the InDels")
@click.option('--input_mnvs',
              '-in_mnvs',
              required=False,
              help="path to the maf file with the MNVs")
@click.option('--input_cnv',
              '-cnv',
              required=True,
              help="FACETS results file")
@click.option('--patient',
              '-pat',
              required=True,
              help="patient id")
@click.option('--output_path_plot',
              '-out_pp',
              type=click.Path(exists=True),
              required = True,
              help="Output path to make plot")
@click.option('--output_path_maf',
              '-out_pm',
              type=click.Path(exists=True),
              required = True,
              help="Output path to write the joined mafs. Give just the path because the script gives intermediate mafs")
@click.option('--comparison',
              '-com',
              required=True,
              help="Comparison structure showing somatic calls from tumor (normally TumorID_vs_normalID). This"
                   "is like a sample id in our project. In here is necessary for the clinical data")
@click.option('--clinical_data',
              '-clinic',
              required=True,
              help="TSV file with the clinical data")
@click.option('--population_frequency',
              '-pop_freq',
              required=True,
              help="Population frequency threshold to filter variants. Example: 0.01 represents that variants in at "
                   "least 1% in the population will be filter out"
              )
def cli(input_snvs, input_indels, input_mnvs, output_path_plot, output_path_maf, patient,comparison,input_cnv,
        clinical_data, population_frequency):

    """
    From VEP annotations it filters population SNPs and adds copy number information to the variants.

    """
    df = strelka_join_refilter(input_snvs, input_indels, input_mnvs)
    df.to_csv(os.path.join(output_path_maf,"all_anno_vep92.maf"),
              sep='\t', index=False) # write intermediate maf with final unique variant row

    df_clinical = pd.read_csv(clinical_data, sep='\t')
    df_clinical = stage_mapping(df_clinical)

    df = categorize_muts(df, comparison, input_cnv,output_path_maf,
                    df_clinical, population_frequency)
    barplot_individual_filter(df, patient, comparison, output_path_plot, float(population_frequency))


if __name__ == '__main__':
    cli()