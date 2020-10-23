import os
import pandas as pd
import numpy as np
from io import StringIO
import gzip
from bgreference import hg19
from aux_data_in_pyvar import CHANNELS

# COMMON PROCESSING FUNCTIONS
#
# def read_vcf(filename, comment='##', sep='\t'):
#     """
#     VCF has a long header starting with ##. The function reads the VCF and converts it into a Pandas dataframe
#     :param filename: input VCF file
#     :param comment: ## is the default but it can be any header preceding data in dataframe form
#     :param sep: tabulated space is the default but it can be anything
#     :return: Pandas DataFrame
#     """
#     if filename.endswith('vcf.gz'):
#         lines = ''.join([line for line in gzip.open(filename, 'rt') if not line.startswith(comment)])
#         return pd.read_csv(StringIO(lines), sep=sep, low_memory=False)
#     else:
#         lines = ''.join([line for line in open(filename) if not line.startswith(comment)])
#         return pd.read_csv(StringIO(lines), sep=sep, low_memory=False)
def read_vcf(filename, comment='##', sep='\t'):
    """
    #     VCF has a long header starting with ##. The function reads the VCF and converts it into a Pandas dataframe
    #     :param filename: input VCF file
    #     :param comment: ## is the default but it can be any header preceding data in dataframe form
    #     :param sep: tabulated space is the default but it can be anything
    #     :return: Pandas DataFrame
    #     """
    if filename.endswith('vcf.gz'):
        lines = ''.join([line for line in gzip.open(filename, 'rt') if not line.startswith(comment)])
        return pd.read_csv(StringIO(lines), sep=sep)
    else:
        lines = ''.join([line for line in open(filename) if not line.startswith(comment)])
        try:
            return pd.read_csv(StringIO(lines), sep=sep)
        except pd.io.common.EmptyDataError:
            return pd.DataFrame()


def get_reads_snvs(rw, status, var=None):
    """
    The name of the function defines the purpose of this function
    :param rw: row of the VCF with the variant (SNVs only)
    :param status: TUMOR or NORMAL columns with the read information
    :param var: ALT or REF depending on the info that one wants to retrieve
    :return: the same row pandas Series object but with added columns
    """
    '''
    rw: row and therefore variant of the dataframe
    '''
    if status == 'TUMOR':
        rw['DP_tumor'] = rw[status].split(':')[0]
    elif status == 'NORMAL':
        rw['DP_normal'] = rw[status].split(':')[0]
    else:
        print('Not a valid status, only TUMOR or NORMAL')

    if rw[var] == 'A':
        rw['reads'] = rw[status].split(':')[-4:-3][0].split(',')[0]
    elif rw[var] == 'C':
        rw['reads'] = rw[status].split(':')[-3:-2][0].split(',')[0]
    elif rw[var] == 'G':
        rw['reads'] = rw[status].split(':')[-2:-1][0].split(',')[0]
    elif rw[var] == 'T':
        rw['reads'] = rw[status].split(':')[-1:][0].split(',')[0]
    else:
        print('You should indicate whether is REF or ALT that you are referring at')
    return rw


def genotype_strelka(rw):
    """
    Gets the genotype estimates of each variant by Strelka. It is different in SNV than InDels
    :param rw: row of the VCF variant
    :return: row with the columns with the genotype information
    """
    if rw['mut_type'] == 'snv':
        info = rw['INFO'].replace('SOMATIC;', '')
        info = info.replace(';OVERLAP', '')
        info_things = dict(item.split("=") for item in info.split(";"))
        rw['GT'] = info_things['SGT']

        #rw['GT'] = rw['INFO'].split(';')[6].split('=')[1]
        if rw['GT'] == str(rw['REF']+rw['REF']+'->'+rw['REF']+rw['ALT']):
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '0/1'
        elif rw['GT'] == str(rw['REF']+rw['REF']+'->'+rw['ALT']+rw['REF']):
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '1/0'
        elif rw['GT'] == str(rw['REF']+rw['REF']+'->'+rw['ALT']+rw['ALT']):
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '1/1'
        else:
            print(rw['GT'])
    else:

        info = rw['INFO'].replace('SOMATIC;', '')
        info = info.replace(';OVERLAP', '')
        info_things = dict(item.split("=") for item in info.split(";"))
        rw['GT'] = info_things['SGT']

        #rw['GT'] = rw['INFO'].split(';')[6].split('=')[1]
        if rw['GT'] == 'ref->het':
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '0/1'
        elif rw['GT'] == 'ref->hom':
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '1/1'
        elif rw['GT'] == 'het->het':
            rw['GT_normal'] = '0/1'
            rw['GT_tumor'] = '1/0'
        elif rw['GT'] == 'ref->ref':
            rw['GT_normal'] = '0/1'
            rw['GT_tumor'] = '0/1'
        elif rw['GT'] == 'hom->hom':
            rw['GT_normal'] = '0/0'
            rw['GT_tumor'] = '1/1'
        else:
            print(rw['GT'])
    rw = rw.drop('GT', 0)
    return rw


def snvs_processing(rw):
    """
    process SNV lines from VCF
    :param rw: row of the tab delimited file (VCF)
    :return: same rw with more columns
    """
    rw['mut_type'] = 'snv'

    # Retrieve read info
    rw = get_reads_snvs(rw, 'TUMOR', 'ALT')
    rw.rename({'reads':'t_alt_reads'}, inplace=True)

    rw = get_reads_snvs(rw, 'TUMOR', 'REF')
    rw.rename({'reads':'t_ref_reads'}, inplace=True)

    rw = get_reads_snvs(rw, 'NORMAL', 'ALT')
    rw.rename({'reads':'n_alt_reads'}, inplace=True)

    rw = get_reads_snvs(rw, 'NORMAL', 'REF')
    rw.rename({'reads':'n_ref_reads'}, inplace=True)

    # Retrieve estimated genotype
    rw = genotype_strelka(rw)

    return rw


def reads_strelka_indels(rw, status):
    """
    The name of the function defines the purpose of this function.
    :param rw: row of the VCF with the variant (InDels only)
    :param status: TUMOR or NORMAL columns with the read information
    :return: the same row pandas Series object but with added columns
    """
    if status == 'TUMOR':
        rw['t_ref_reads'] = rw[status].split(':')[2:3][0].split(',')[0] # [2:3] corresponds to TAR
        rw['t_alt_reads'] = rw[status].split(':')[3:4][0].split(',')[0] # [3:4] corresponds to TIR
        rw['DP_tumor'] = rw[status].split(':')[0]
    elif status == 'NORMAL':
        rw['n_ref_reads'] = rw[status].split(':')[2:3][0].split(',')[0] # [2:3] corresponds to TAR
        rw['n_alt_reads'] = rw[status].split(':')[3:4][0].split(',')[0] # [3:4] corresponds to TIR
        rw['DP_normal'] = rw[status].split(':')[0]
    else:
        print('Not a valid status, only TUMOR or NORMAL')
    return rw


def indels_processing(rw):
    """
    process InDels lines from VCF
    :param rw: row of the tab delimited file (VCF)
    :return: same rw with more columns
    """
    rw['mut_type'] = "indels"

    # Retrieve read info
    rw = reads_strelka_indels(rw, 'TUMOR')
    rw = reads_strelka_indels(rw, 'NORMAL')

    # Retrieve estimated genotype
    rw = genotype_strelka(rw)

    return rw

# PIVOT CLINICAL DATA

def stage_mapping(df):

    """
    Changes the clinical input data to a more manageable data frame
    :param df: clinical data frame
    :return: transformed clinical data frame
    """
    df_pry = df[['Patient_id', 'Sex', 'Primary_seq_id', 'Remission_seq_id', 'Primary_diagnosis_age',
                'days_between_pry_rel','Primary_immunoclassification']]
    df_rel = df[['Patient_id', 'Sex', 'Relapse_seq_id', 'Remission_seq_id', 'Primary_diagnosis_age',
                'days_between_pry_rel', 'Primary_immunoclassification']]

    df_pry.rename(columns={'Patient_id':'PATIENT', 'Sex':'SEX', 'Primary_seq_id':'SAMPLE',
                   'Primary_immunoclassification':'SUBTYPE', 'Primary_diagnosis_age':'AGE'}, inplace=True)
    df_rel.rename(columns={'Patient_id': 'PATIENT', 'Sex': 'SEX', 'Relapse_seq_id': 'SAMPLE',
                   'Primary_immunoclassification': 'SUBTYPE','Primary_diagnosis_age':'AGE'}, inplace=True)


    df_pry['COMPARISON'] = df_pry['SAMPLE']+'_vs_'+df_pry['Remission_seq_id']
    df_rel['COMPARISON'] = df_rel['SAMPLE']+'_vs_'+ df_rel['Remission_seq_id']

    df_pry['STAGE'] = 'primary'
    df_rel['STAGE'] = 'relapse'

    dff = df_pry.copy()
    dff = dff.append(df_rel, ignore_index=True)

    dff['AGE_REL'] = dff.apply(lambda x: x['AGE'] + (x['days_between_pry_rel'] / 365) if pd.isnull(x['days_between_pry_rel']) == False else None, axis=1)

    return dff

# GET TRINUCLEOTIDE CONTEXT OF SNVS IN PYRIMIDE REF


def get_context_rev(rw):
    equival_nt = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    pos = rw['POS']
    left, ref, right = hg19(rw['#CHROM'], pos - 1, size=3)
    alt = rw['ALT']
    rw['TRIPLE'] = left + '[' + ref + '>' + alt + ']' + right
    rw['TRIPLE_COM'] = equival_nt[left] + '[' + equival_nt[ref] + '>' + equival_nt[alt] + ']' + equival_nt[right]
    rw['TRIPLE_COM_REV'] = equival_nt[right] + '[' + equival_nt[ref] + '>' + equival_nt[alt] + ']' + equival_nt[left]

    return rw


def add_pyrimidine_type(df):

    dff = pd.DataFrame()

    for i, rw in df.iterrows():
        if rw['TRIPLE'] in CHANNELS:
            rw['VARIANT_TYPE'] = rw['TRIPLE']
            df_rw_T = rw.to_frame().T
            dff = dff.append(df_rw_T, ignore_index=True)
        else:
            rw['VARIANT_TYPE'] = rw['TRIPLE_COM_REV']
            df_rw_T = rw.to_frame().T
            dff = dff.append(df_rw_T, ignore_index=True)
    dff.drop(columns=['TRIPLE','TRIPLE_COM', 'TRIPLE_COM_REV'], inplace=True)
    return dff


def count_variant_type(df):

    count = {ctxt: 0 for ctxt in CHANNELS}

    for i, rw in df.iterrows():
        count[rw['VARIANT_TYPE']] = count[rw['VARIANT_TYPE']] + 1

    return count


def df_to_dict(df):
    df.rename(columns={'CHROMOSOME': '#CHROM', 'POSITION': 'POS'}, inplace=True)
    df = df[['#CHROM', 'POS', 'REF', 'ALT']]
    df.drop_duplicates(keep='first', inplace=True)

    # make dictionary with counts centered in pyrimidines
    df = df.apply(lambda x: get_context_rev(x), axis=1)
    df = add_pyrimidine_type(df)
    count = count_variant_type(df)

    # from dict to pandas df
    df_96 = pd.DataFrame.from_dict(count, orient='index')
    df_96.reset_index(inplace=True)
    df_96.columns = ['change', 'count']

    # change to dictionary with special key annotation
    df_96[['count']] = df_96[['count']].astype(int)
    #     total = df_96['count'].sum()
    #     df_96['relative_count'] = df_96['count'].apply(lambda x: x/total)
    df_96 = df_96[['count', 'change']]
    df_96.set_index('change', inplace=True)
    df_96 = df_96.sort_index()
    dictionary = df_96.to_dict()['count']

    # fill missing context with 0 counts
    for k in CHANNELS:
        if k not in dictionary.keys():
            dictionary[k] = 0
    return dictionary


def get_muts_x_signature(sh, pp, pr, pat, sig, prob_file_path):
    # GET CONTEXTS
    dict_pry = df_to_dict(pp)
    dict_rel = df_to_dict(pr)
    dict_sh = df_to_dict(sh)

    # GET signature probabilities by context
    prob = pd.read_csv(os.path.join(prob_file_path, "mutation_sign_prob.tsv"), sep='\t')
    prob_pat = prob[prob['Sample'] == pat].set_index('Mutation_type')
    prob_pat.index.name = None

    count_pp = 0
    count_pr = 0
    count_sh = 0

    for cntxt, count in dict_pry.items():
        prob_sig = prob_pat.loc[cntxt, sig]

        count_pp = count_pp + count * prob_sig
        count_pr = count_pr + dict_rel[cntxt] * prob_sig
        count_sh = count_sh + dict_sh[cntxt] * prob_sig
    return count_pp, count_pr, count_sh


def check_rw(rw):
    if pd.isnull(rw['cf.em']) and pd.isnull(rw['lcn.em']):
        rw['tag'] = 'yes'
    else:
        rw['tag'] = 'no'
    return rw


def process_cnv(dire_cnv, clinical_data):
    facets_results = pd.DataFrame()

    for com in clinical_data['COMPARISON']:
        clinical_data_sample = clinical_data[clinical_data['COMPARISON'] == com].reset_index()
        facets_cnv = pd.read_csv(
            os.path.join(dire_cnv, clinical_data_sample.loc[0, 'PATIENT'], com, 'results_x_seg'), sep='\t')

        facets_cnv[['chrom']] = facets_cnv[['chrom']].astype(str)
        facets_cnv['chrom'] = facets_cnv['chrom'].str.replace("23", "X")
        facets_cnv = facets_cnv[~((facets_cnv['tcn.em'] == 2) & (facets_cnv['lcn.em'] == np.nan))]
        facets_cnv = facets_cnv.apply(lambda x: check_rw(x), axis=1)
        facets_cnv = facets_cnv[facets_cnv['tag'] != 'yes']
        facets_cnv['lcn.em'] = facets_cnv['lcn.em'].fillna(0)

        facets_cnv['stage'] = clinical_data_sample.loc[0, 'STAGE']
        facets_cnv['patient'] = clinical_data_sample.loc[0, 'PATIENT']
        facets_cnv['comparison'] = com
        facets_cnv['sample'] = com.split("_vs_")[0]

        facets_results = facets_results.append(facets_cnv, ignore_index=True)

    facets_results[['lcn.em', 'tcn.em']] = facets_results[['lcn.em', 'tcn.em']].astype(int)
    facets_results = facets_results[
        ['chrom', 'start', 'end', 'lcn.em', 'tcn.em', 'comparison', 'stage', 'sample', 'patient', 'cf.em']]
    return facets_results


def get_three_subsets(all_pry, all_rel):

    # CREATE SET OF VARIANTS CLONAL  SNVS

    all_pry_variants = set(all_pry['Variant'].unique())

    all_rel_variants = set(all_rel['Variant'].unique())

    shared_variants = all_pry_variants.intersection(all_rel_variants)

    private_pry_variants = all_pry_variants.difference(shared_variants)

    private_rel_variants = all_rel_variants.difference(shared_variants)

    return shared_variants, private_pry_variants, private_rel_variants