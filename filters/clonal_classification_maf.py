import sys, os
os.environ["PATH"] = os.path.dirname(sys.executable) + os.pathsep + os.environ["PATH"]
import click
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
from sklearn.neighbors.kde import KernelDensity
from sklearn.mixture import GaussianMixture
import scipy.stats as ss
from scipy.signal import find_peaks
from aux_data_in_pyvar import config_rcparams


def clustering_ccf(df):
    """
    Clusters CCF according to the peaks of their distribution
    :param df: dataframe of variants
    :return: returns cluster assignment to each variant as well as density values of the distribution
    """
    # Oriol Pich' piece of code to cluster ccf values

    # hardcoded!
    best_band = 0.09

    # remove extreme cases
    ccf_list = df['vaf*cnv']
    max_ccf = np.amax(df['vaf*cnv'])

    if max_ccf < 2.8:
        upbound = max_ccf
    else:
        print('there are ccf bigger than 2.8')
        upbound = 2.8

    # do the log2 of each of the ccf values
    ccf = [np.log2(x) for x in ccf_list]
    variant = df['Variant'].tolist()

    X = np.array(ccf).reshape(-1, 1)
    X_var = np.array(variant).reshape(-1, 1)

    kde = KernelDensity(kernel='gaussian', bandwidth=best_band).fit(X)

    grid2 = np.linspace(np.amin(ccf_list), upbound, num=150).reshape(-1, 1)
    grid2 = np.array([np.log2(x) for x in grid2])
    flat_array = grid2.flatten()

    log_density = kde.score_samples(grid2)
    density = np.exp(log_density)

    # find the maximum peaks
    number_components = len(find_peaks(density, height=0.1)[0])

    if number_components == 0:
        # at least 1 component which indicates one cluster
        print("peaks unfound")
        gmm = GaussianMixture(n_components=1, max_iter=2000).fit(X)
    else:
        gmm = GaussianMixture(n_components=number_components, max_iter=2000).fit(X)
    cluster_assign_val = defaultdict(list)
    cluster_assign_var = defaultdict(list)

    df_results = pd.DataFrame()
    for ix, prob in enumerate(np.argmax(gmm.predict_proba(X), axis=1)):
        cluster_assign_val[prob].append(X[ix])
        cluster_assign_var[prob].append(X_var[ix])
        df_results = df_results.append({'Variant': X_var[ix][0], 'ccf_log2': X[ix][0],
                                        'cluster': prob}, ignore_index=True)
    return df_results, cluster_assign_val, flat_array, density


def plot_ccf(cluster_assign, flat_array, density, outfile_plot, comparison, df_purities, cut_off):
    """
    Plot density plot of CCF variants and add peaks and cut -off of clonality categorization to have overview of the
    estimates
    """
    sam_purity = df_purities[df_purities['comparison'] == comparison].reset_index()

    ## make plot
    config_rcparams()

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    for clust, vals in cluster_assign.items():
        center = np.median(vals)
        ax.text(center + 0.1, 1, round(2 ** center, 3))
        ax.vlines(center, 0, 1, color='red')

    for cut, val in cut_off.items():
        ax.text(np.log2(val), 1, cut, fontsize=5, rotation=45)
        ax.vlines(np.log2(val), 0, 1, color='blue')

    ax.plot(flat_array, density, c='#31a354', lw=2)
    plt.xlabel('VAF*CNV')
    xtickslocs = ax.get_xticks()
    ax.set_xticklabels([2 ** i for i in xtickslocs], rotation=45)
    # Add purities estimated by other software to visually compare with the one estimated and used here
    ax.title.set_text(" Sample {} \n purity facets {} \n purity ascat {}".format(comparison,
                                                                                 sam_purity.loc[0, 'purity_facets'],
                                                                                 sam_purity.loc[0, 'purity_ascat']))

    plt.tight_layout()
    plt.savefig(os.path.join(outfile_plot, comparison + '_ccf_clustering.png'), dpi=200)
    plt.close()


def new_ccf(rw, purity):
    """
    Compute CCF with estimate of purity
    :param rw: row (variant)
    :param purity: purity
    :return: CCF of each variant
    """

    vaf = rw['t_alt_reads']/(rw['t_ref_reads']+rw['t_alt_reads'])

    # definition of CCF
    rw['ccf'] = vaf*(purity*rw['cnv_facets'] + (1-purity)*2)/purity
    return rw


def get_prob(rw, purity):
    """
    Assign probabilities to each variant from adjusted betabinomial distribution
    :param rw: row(variant)
    :param purity: purity
    :return:row with probability
    """

    max_CN = rw['cnv_facets']
    alt_count = rw['t_alt_reads']
    depth = rw['t_alt_reads']+rw['t_ref_reads']
    mu = purity/(purity*max_CN + (1-purity)*2)
    val = ss.betabinom.cdf(k=alt_count, n=depth, a=mu * (1 - 0.01) / 0.01, b=(1 - mu) * (1 - 0.01) / 0.01)
    rw['prob'] = val
    return rw


def take_closest(num,collection):
    return min(collection,key=lambda x:abs(x-num))



@click.command()
@click.option('--output_path_plot',
              '-out_pp',
              type=click.Path(exists=True),
              required = True,
              help="Output path to make plot")
@click.option('--output_path_maf',
              '-out_pp',
              type=click.Path(exists=True),
              required = True,
              help="Output path to write data frame")
@click.option('--input_path',
              '-in_path',
              required=True,
              help="input maf file")
@click.option('--input_purities',
              '-pur',
              required=True,
              help="Data frame with all purities estimated by sample by other software (ASCAT and FACETS)."
                   "you can find this in ../ext_files/purity_ploidy_TALL_adult.tsv")
@click.option('--comparison',
              '-com',
              required=True,
              help="Comparison structure showing somatic calls from tumor (normally TumorID_vs_normalID). This"
                   "is like a sample id in our project.")


def cli(output_path_plot, output_path_maf, input_path,comparison,input_purities):

    """
    Infer purity from distribution of CCF of variants as the maximum point and categorize mutations as clonal or
    subclonal
    """

    # read file
    df_in = pd.read_csv(input_path, sep='\t')
    # compute VAF
    df_in['vaf'] = df_in.apply(lambda x: x['t_alt_reads'] / (x['t_alt_reads'] + x['t_ref_reads']), axis=1)

    # compute depth and filter
    print(len(df_in))
    df_in['depth'] = df_in.apply(lambda x: x['t_alt_reads'] + x['t_ref_reads'], axis=1)
    df_in = df_in[df_in['depth'] > 5]
    print(len(df_in))

    # compute ccf assuming 100 of purity
    df_in['cnv_facets'] = df_in['cnv_facets'].astype(int)
    df_in['vaf*cnv'] = df_in['vaf'] * df_in['cnv_facets']

    # cluster ccf
    df_results, cluster_assign, flat_array, density = clustering_ccf(df_in)

    # parse results cluster assignment
    grps = df_results.groupby('cluster')
    df_results = pd.DataFrame()

    for g in grps.groups:
        cluster = grps.get_group(g)
        ccfs = cluster['ccf_log2']
        center = np.median(ccfs)
        name = round(2 ** center, 3)
        cluster['center_cluster'] = name
        df_results = df_results.append(cluster, ignore_index=True)

    df_out = df_in.merge(df_results, how='inner', on='Variant')

    df_pur_decision = df_out[['center_cluster', 'Variant']].groupby('center_cluster').count()
    suma = df_pur_decision['Variant'].sum()

    df_pur_decision['proportion'] = df_pur_decision['Variant'].apply(lambda x: x / suma)

    df_pur_decision.reset_index(inplace=True)

    df_pur_decision.rename(columns={'index': 'center_cluster'}, inplace=True)

    df_pur_decision.sort_values(by='center_cluster', ascending=False, inplace=True)

    df_pur_decision.reset_index(inplace=True)

    val = 0
    for i, rw in df_pur_decision.iterrows():
        if i == 0 and rw['proportion'] > 0.1:
            pur = rw['center_cluster']
            break
        elif i != 0 and val > 0.1:
            pur = rw['center_cluster']
            break
        else:
            val = val + rw['proportion']
            continue

    print(pur)
    # manually added after inspection, Hardcoded sorry
    if comparison == 'AE6526_vs_AE6525':
        pur = 0.253
    elif comparison == 'AE6533_vs_AE6534':
        pur = 0.65
    elif comparison == 'AE6518_vs_AE6519':
        pur = 0.345
    elif comparison == "AE6514_vs_AE6513":
        pur = 0.532
    elif comparison == 'AE6544_vs_AE6543':
        pur = 0.927
    elif comparison == 'SJTALL014_D_vs_SJTALL014_G':
        pur = 0.792
    else:
        pass
    df_out = df_out.apply(lambda x: get_prob(x, pur), axis=1)
    df_out.sort_values('prob', inplace=True, ascending=True)
    df_out['purity'] = pur
    cut_off = dict()

    cut_off['0.01'] = df_out[df_out['prob'] == take_closest(0.01, df_out['prob'].tolist())].reset_index().loc[0, 'vaf*cnv']
    df_out['clonal_classification'] = df_out.apply(lambda x: 'subclonal' if x['prob'] < 0.01 else 'clonal', axis=1)

    # density plot
    dff_purities = pd.read_csv(input_purities, sep='\t')
    plot_ccf(cluster_assign, flat_array, density, output_path_plot, comparison, dff_purities, cut_off)

    # compute ccf
    df_out = df_out.apply(lambda x: new_ccf(x, pur), axis=1)

    # drop unnecessary columns
    df_out.drop(labels=['VAF','ccf_log2'], axis='columns', inplace=True)

    # write results
    df_out.to_csv(output_path_maf, index=False, sep='\t')


if __name__ == '__main__':
    cli()