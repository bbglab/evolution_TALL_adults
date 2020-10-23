import sys, os
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
import pandas as pd
import click
import subprocess
import collections


def identify_possible_mnvs(df):
    """
    Checks for consecutive variants in VCFs that might be dinucleotide variants. If this is the case they should be
    reported and annotated together.
    :param df: dataframe with the variants as rows (MAF file previosuly created from VCF)
    :return: dataframe of dinucleotide candidates
    """

    grps = df.groupby("#CHROM")
    dinucleotide_variants = pd.DataFrame()

    for g in grps.groups:
        df_chrom = grps.get_group(g)
        df_chrom['DIFF'] = df_chrom['POS'].diff()
        dnt_1 = df_chrom[df_chrom['DIFF'] < 2]
        indices_1 = dnt_1.index
        indices_2 = [x - 1 for x in indices_1]
        dnt_2 = df_chrom[df_chrom.index.isin(indices_2)]
        dnt_1.sort_values(by=['#CHROM', 'POS'], inplace=True)
        dnt_1.reset_index(drop=True, inplace=True)
        dnt_2.sort_values(by=['#CHROM', 'POS'], inplace=True)
        dnt_2.reset_index(drop=True, inplace=True)
        for i in range(0, len(dnt_1), 1):
            dinucleotide_variants = dinucleotide_variants.append(
                {'CHROM': g, 'START': dnt_2.loc[i, 'POS'],
                 'END': dnt_1.loc[i, 'POS'], 'CHR_START': g + ":" + str(int(dnt_2.loc[i, 'POS'])),
                 'DNT_ALT': dnt_2.loc[i, 'ALT'] +
                            dnt_1.loc[i, 'ALT'],
                 'DNT_REF': dnt_2.loc[i, 'REF'] +
                            dnt_1.loc[i, 'REF']}, ignore_index=True)

    print("num DNT candidates {}".format(len(dinucleotide_variants)))
    dinucleotide_variants[['START', 'END']] = dinucleotide_variants[['START', 'END']].astype(int)
    dinucleotide_variants = dinucleotide_variants[['CHROM', 'START', 'END', 'CHR_START',
                                                   'DNT_REF', 'DNT_ALT']]
    dinucleotide_variants.sort_values(by=['CHROM', 'START', 'END'], inplace=True)
    return dinucleotide_variants


def get_count_samtools(rw, path_bam):
    """
    For the positions determined in the rw (row) it gets the read counts supporting the alleles
    :param rw: row with the variant
    :param path_bam: path to the BAM file
    :return: defaultdict with counts for each pair of nucleotides
    """

    result = subprocess.getoutput('samtools tview -p ' + rw[
            'CHR_START'] + " -d T " + path_bam + " | cut -c1,2 | sort | uniq -c")
    result = result.split('\n')
    result = [item for item in result if not 'older' in item]

    result = [x.strip() for x in result]
    dnt_dict = collections.defaultdict(int)
    for info in result:
        if any(nt in info for nt in ['A', 'a', 'C', 'c', 'G', 'g', 'T', 't']):
            try:
                count, dnt = info.split(" ")
            except ValueError:
                try:
                    count, dnt = info.split("  ")
                except ValueError:
                    print(info)
            dnt_dict[dnt.upper()] += int(count)
    return dnt_dict


def count_pos_BAM(rw, path_bam_tumoral, path_bam_normal):

    """
    Adds read counts of dinucleotide positions calling another function that runs samtools
    :param rw: row with the variant info
    :param path_bam_tumoral: path to the tumoral BAM
    :param path_bam_normal: path to the normal BAM
    :return: row with columns with the read info
    """

    # TUMOR
    dnt_dict_t = get_count_samtools(rw, path_bam_tumoral)
    rw['SAMTOOLS_COUNT_DNT_ALT'] = dnt_dict_t[rw['DNT_ALT']]
    rw['SAMTOOLS_COUNT_DNT_REF'] = dnt_dict_t[rw['DNT_REF']]
    dnt_dict_writen = {rw['DNT_ALT']: dnt_dict_t[rw['DNT_ALT']], rw['DNT_REF']: dnt_dict_t[rw['DNT_REF']]}
    dnt_dict_left = {k: v for k, v in dnt_dict_t.items() if k not in dnt_dict_writen}
    rw['SAMTOOLS_COUNT_OTHERS'] = str(dnt_dict_left)
    rw['SAMTOOLS_COUNT_OTHERS'] = rw['SAMTOOLS_COUNT_OTHERS'].replace("{", "")
    rw['SAMTOOLS_COUNT_OTHERS'] = rw['SAMTOOLS_COUNT_OTHERS'].replace("}", "")

    # NORMAL
    dnt_dict_n = get_count_samtools(rw, path_bam_normal)
    rw['SAMTOOLS_COUNT_NOR_ALT'] = dnt_dict_n[rw['DNT_ALT']]
    rw['SAMTOOLS_COUNT_NOR_REF'] = dnt_dict_n[rw['DNT_REF']]
    dnt_dict_writen = {rw['DNT_ALT']: dnt_dict_n[rw['DNT_ALT']], rw['DNT_REF']: dnt_dict_n[rw['DNT_REF']]}
    dnt_dict_left = {k: v for k, v in dnt_dict_n.items() if k not in dnt_dict_writen}
    rw['SAMTOOLS_COUNT_NOR_OTHERS'] = str(dnt_dict_left)
    rw['SAMTOOLS_COUNT_NOR_OTHERS'] = rw['SAMTOOLS_COUNT_NOR_OTHERS'].replace("{", "")
    rw['SAMTOOLS_COUNT_NOR_OTHERS'] = rw['SAMTOOLS_COUNT_NOR_OTHERS'].replace("}", "")

    return rw


@click.command()
@click.option('--output_path',
              '-out_path',
              type=click.Path(exists=True),
              required = True,
              help="Output ")
@click.option('--input_maf',
              '-maf',
              required=True,
              help="maf to check for wrong MNVs")
@click.option('--input_tumor_bam',
              '-t_bam',
              required=True,
              help="tumor bam to check MNVs in reads")
@click.option('--input_normal_bam',
              '-n_bam',
              required=True,
              help="normal bam to check MNVs in reads")


def cli(output_path, input_maf,input_tumor_bam,input_normal_bam):

    # read input maf
    df = pd.read_csv(input_maf, sep='\t')
    df['chr_start'] = df.apply(lambda x: str(x['#CHROM'])+':'+str(x['POS']), axis=1)

    # ensure data types
    df[['#CHROM']] = df[['#CHROM']].astype(str)
    df[['POS']] = df[['POS']].astype(int)

    # ensure sorted positions
    df.sort_values(by=['#CHROM', 'POS'], inplace=True, ascending=True)
    df.reset_index(inplace=True)

    # get consecutive variants that are suspicious to be wrong annotated MNVs
    df_dnt = identify_possible_mnvs(df)

    # count reads with MNV in bam file with samtools
    dnt_var = df_dnt.apply(lambda x: count_pos_BAM(x, input_tumor_bam, input_normal_bam), axis=1)

    # filter out the individual mutations
    dnt_var = dnt_var[(dnt_var['SAMTOOLS_COUNT_DNT_REF'] != 0) & (dnt_var['SAMTOOLS_COUNT_DNT_ALT'] != 0)]

    dnt_var['CHR_END'] = dnt_var.apply(lambda x: x['CHROM'] + ':' + str(x['END']), axis=1)

    print("All variants in maf:{}".format(len(df)))
    df_filt = df[~df['chr_start'].isin(dnt_var['CHR_START'].tolist())]
    df_filt = df_filt[~df_filt['chr_start'].isin(dnt_var['CHR_END'].tolist())]
    df_filt.reset_index(drop=True, inplace=True)
    print("Filtered variants in maf:{}".format(len(df_filt)))

    df_mnvs = pd.DataFrame()

    for i, rw in dnt_var.iterrows():
        df_mnvs = df_mnvs.append({'#CHROM': rw['CHROM'], 'POS': rw['START'], 'ID': '.', 'REF': rw['DNT_REF'],
                                  'ALT': rw['DNT_ALT'], 'QUAL': '.', 'FILTER': 'PASS', 'INFO': '.', 'FORMAT': '.',
                                  'NORMAL': '.', 'TUMOR': '.',
                                  'DP_tumor': rw['SAMTOOLS_COUNT_DNT_REF'] + rw['SAMTOOLS_COUNT_DNT_ALT'],
                                  't_alt_reads': rw['SAMTOOLS_COUNT_DNT_ALT'],for_misss
                                  't_ref_reads': rw['SAMTOOLS_COUNT_DNT_REF'],
                                  'DP_normal': rw['SAMTOOLS_COUNT_DNT_REF'] + rw['SAMTOOLS_COUNT_DNT_ALT'],
                                  'n_alt_reads': rw['SAMTOOLS_COUNT_NOR_ALT'],
                                  'n_ref_reads': rw['SAMTOOLS_COUNT_NOR_REF'],
                                  'GT_normal': '.',
                                  'GT_tumor': '.',
                                  'mut_type': 'mnv'}, ignore_index=True)

    # ensure data types and order of columns

    int_cols = ['POS', 't_alt_reads', 't_ref_reads','n_alt_reads', 'n_ref_reads']
    order_cols = ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                       'NORMAL', 'TUMOR', 'DP_tumor', 't_alt_reads', 't_ref_reads', 'DP_normal',
                       'n_alt_reads', 'n_ref_reads', 'mut_type', 'GT_normal', 'GT_tumor']

    df_mnvs = df_mnvs[order_cols]
    df_mnvs[int_cols] = df_mnvs[int_cols].astype(int)

    df_filt = df_filt[order_cols]

    print(len(df_filt))
    print(len(df_mnvs))

    file_filt = input_maf.replace(".maf", "")
    df_filt.to_csv(os.path.join(output_path, file_filt + '_checked.maf'), sep='\t', index=False)

    file_mnvs = input_maf.replace("_snvs_sh.maf", "")
    df_mnvs.to_csv(os.path.join(output_path, file_mnvs + '_mnvs_sh_checked.maf'), sep='\t', index=False)


if __name__ == '__main__':
    cli()