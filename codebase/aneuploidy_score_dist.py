
import matplotlib
matplotlib.use('Agg')

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.stats import hypergeom
from statsmodels.sandbox.stats.multicomp import fdrcorrection0
import constants
import seaborn as sns
from utils import calc_hg_enrichment_pval, calc_hg_depletion_pval
from scipy.stats import zscore, rankdata, wilcoxon, ks_2samp
import argparse
from utils import PVAL_METHODS

def calc_pval_dist(mat, aneu_types, chr_interaction, pval_method_name, label):
    hg_pvals=[]
    for a in np.arange(1,23):
        for arm_a in ['p','q']:
            for b in np.arange(a, 23):
                for arm_b in ['p', 'q']:
                    for aneu_type_a, aneu_type_b in aneu_types:

                        if (arm_a <= arm_b and b == a) or (chr_interaction and b == a) or (not chr_interaction and b != a): continue
                        if "{}{}".format(a, arm_a) not in mat.columns: continue
                        if "{}{}".format(b, arm_b) not in mat.columns: continue


                        pval=PVAL_METHODS[pval_method_name](mat, a, arm_a, aneu_type_a, b, arm_b, aneu_type_b)

                        hg_pvals.append({"{}{}{}-{}{}{}".format(a, arm_a, aneu_type_a, b, arm_b, aneu_type_b) : pval})

    qvals=fdrcorrection0([a.values()[0] for a in hg_pvals])[1]

    hg_qvals=[]
    for i, a in enumerate(hg_pvals):
        hg_qvals.append({a.keys()[0]:qvals[i]})

    df=pd.DataFrame(data=[[a.values()[0], hg_qvals[i].values()[0]] for i, a in enumerate(hg_pvals)], columns=['pval','qval'], index=[a.keys()[0] for a in hg_pvals])
    df["enrichment_score"]=-np.log10(df.loc[:,"pval"])
    df["zscore"]=zscore(df.loc[:,"enrichment_score"])
    df["rank"]=rankdata(df.loc[:,"pval"])
    df=df.sort_values(by=["pval"])
    df.to_csv(os.path.join(constants.OUTPUT_FOLDER,"hgs_emp_dist_{}_{}_{}.tsv".format(constants.CHR_INTERACTION_NAMES[chr_interaction], pval_method_name, label)), sep='\t')

    return df


def generate_dist_for_dataset(mats, file_names, pval_method_name):

    fig,axs=plt.subplots(2,2,figsize=(20,20))

    #  compare WGD within c
    chr_interaction=constants.INTRA
    dists=[]
    dfs_within=[]
    for mat, file_name in zip(mats, file_names):
        df=calc_pval_dist(mat, aneu_types=[(-1, 1), (1, -1)], chr_interaction=chr_interaction, pval_method_name=pval_method_name, label=file_name)
        dists.append(len(df.index))
        dfs_within.append(df)


    plot_dist(axs[0][0], dfs_within)

    s,p=ks_2samp(dfs_within[0].loc[:,"enrichment_score"], dfs_within[1].loc[:,"enrichment_score"])
    axs[0][0].set_title("{} chromosomal effect of WGD (p-value (KS-2samp)={})".format(constants.CHR_INTERACTION_NAMES[chr_interaction], "%.2E" % p))
    axs[0][0].legend(labels=["{} (n={})".format(a,b) for a,b in zip(file_names, dists)])

    df_summary_agg=pd.DataFrame()
    df_agg_pvals=pd.concat([dfs_within[0].loc[:,"pval"], dfs_within[1].loc[:,"pval"]] ,axis=1)
    df_agg_zscores=pd.concat([dfs_within[0].loc[:,"zscore"], dfs_within[1].loc[:,"zscore"]] ,axis=1)
    df_summary_agg["pval_ratio"]=(-np.log10(df_agg_pvals.iloc[:,0]))/(-np.log10(df_agg_pvals.iloc[:,1])+1)
    df_summary_agg["z_diff"]=df_agg_zscores.iloc[:,0]-df_agg_zscores.iloc[:,1]
    df_summary_agg["z_rank"]=rankdata(-df_summary_agg["z_diff"])
    df_summary_agg["pval_rank"]=rankdata(-df_summary_agg["pval_ratio"])
    df_summary_agg.to_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_summary_{}_{}.tsv".format(constants.CHR_INTERACTION_NAMES[chr_interaction], pval_method_name)), sep='\t')

    sns.distplot(df_summary_agg["z_diff"], norm_hist=False, kde=False, ax=axs[1][0])
    axs[1][0].set_yscale('log')
    axs[1][0].set_xlabel('z score difference')
    axs[1][0].set_ylabel('# of co-occurrences')
    axs[1][0].set_title('z-score difference in intra-chromosomal pairs between WGD plus and WGD minus datasets')
    axs[1][0].legend(labels=["z-score difference"])


    #  compare WGD between c
    chr_interaction=constants.INTER
    dists=[]
    dfs_between=[]
    for mat, label in zip(mats, file_names):
        df=calc_pval_dist(mat, aneu_types=[(-1, -1), (-1, 1), (1, -1), (1, 1)], chr_interaction=chr_interaction, pval_method_name=pval_method_name, label=label)
        dists.append(len(df.index))
        dfs_between.append(df)

    plot_dist(axs[0][1], dfs_between)

    s,p=ks_2samp(dfs_between[0].loc[:,"enrichment_score"].values, dfs_between[1].loc[:,"enrichment_score"].values)
    axs[0][1].set_title("{}-chromosomal effect of WGD (p-value KS-2samp)={})".format(constants.CHR_INTERACTION_NAMES[chr_interaction], "%.2E" % p))
    axs[0][1].legend(labels=["{} (n={})".format(a,b) for a,b in zip(file_names, dists)])

    df_summary_agg=pd.DataFrame()
    df_agg_pvals=pd.concat([dfs_between[0].loc[:,"pval"], dfs_between[1].loc[:,"pval"]] ,axis=1)
    df_agg_zscores=pd.concat([dfs_between[0].loc[:,"zscore"], dfs_between[1].loc[:,"zscore"]] ,axis=1)
    df_summary_agg["pval_ratio"]=(-np.log10(df_agg_pvals.iloc[:,0]))/(-np.log10(df_agg_pvals.iloc[:,1])+1)
    df_summary_agg["z_diff"]=df_agg_zscores.iloc[:,0]-df_agg_zscores.iloc[:,1]
    df_summary_agg["z_rank"]=rankdata(-df_summary_agg["z_diff"])
    df_summary_agg["pval_rank"]=rankdata(-df_summary_agg["pval_ratio"])
    df_summary_agg.to_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_summary_{}_{}.tsv".format(constants.CHR_INTERACTION_NAMES[chr_interaction], pval_method_name)), sep='\t')

    sns.distplot(df_summary_agg["z_diff"], norm_hist=False, kde=False, ax=axs[1][1])
    axs[1][1].set_yscale('log')
    axs[1][1].set_xlabel('z-score difference')
    axs[1][1].set_ylabel('# of co-occurrences')
    axs[1][1].set_title('z-score difference in inter-chromosomal pairs between WGD plus and WGD minus datasets')
    axs[1][1].legend(labels=["z-score difference"])


    plt.savefig(os.path.join(constants.OUTPUT_FOLDER, "figures", "hgs_agg_{}_{}_{}_{}.png".format(constants.CHR_INTERACTION_NAMES[chr_interaction], pval_method_name, *file_names)))

    # =============================== end fig 1 analysis ===============================================

    #  for each WGD type, compare directionality within c
    chr_interaction=constants.INTRA

    fig,axs=plt.subplots(3,2,figsize=(20,30))
    dists={}
    for mat, label, ax in zip(mats, file_names, axs.T):
        dists[label]=[]
        dists[label].append(calc_pval_dist(mat, aneu_types=[(-1, -1), (1, 1)], chr_interaction=chr_interaction, pval_method_name=pval_method_name, label="sm_{}".format(label)))
        dists[label].append(calc_pval_dist(mat, aneu_types=[(-1, 1), (1, -1)], chr_interaction=chr_interaction, pval_method_name=pval_method_name, label="op_{}".format(label)))
        plot_dist(ax[0], dists[label])
        ax[0].set_title("{}-chromosomal comparison of aberration orientations in {}".format(constants.CHR_INTERACTION_NAMES[chr_interaction], label))
        ax[0].legend(labels=["{} (n={})".format(a,b) for a,b in
                          zip(["same aberrations", "opposite aberrations"], [len(dists[label][0].index), len(dists[label][1].index)])])


        df_zscore_diff=pd.DataFrame()
        for a in np.arange(1,23):

            sm_indices = ["{}q1-{}p1".format(a,a), "{}q-1-{}p-1".format(a,a)]
            op_indices = ["{}q-1-{}p1".format(a,a), "{}q1-{}p-1".format(a,a)]

            try:
                df_zscore_diff.loc["chr{}".format(a),"zscore_diff"]=np.max(dists[label][0].loc[sm_indices,"zscore"]) - np.max(dists[label][1].loc[op_indices,"zscore"])
            except KeyError:
                continue
        z_diff=dists[label][0].loc[:, "enrichment_score"]-dists[label][1].loc[:, "enrichment_score"]

        sns.distplot(df_zscore_diff.loc[:,"zscore_diff"], norm_hist=False, kde=False, ax=ax[1], bins=20)
        ax[1].set_xlabel('z-score difference')
        ax[1].set_ylabel('# of co-occurrences')
        ax[1].set_title('z-score difference between intra and inter chromosomal pairs in {}'.format(label))
        ax[1].legend(labels=["z-score difference"])
        df_zscore_diff.to_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_intra_zscore_diff_sm_vs_op_{}.tsv".format(file_name)), sep='\t')


    df_zscore_diff_between_datasets={"sm" : pd.DataFrame(), "op": pd.DataFrame()}
    for a in np.arange(1,23):

        sm_indices = ["{}q1-{}p1".format(a,a), "{}q-1-{}p-1".format(a,a)]
        op_indices = ["{}q-1-{}p1".format(a,a), "{}q1-{}p-1".format(a,a)]

        try:
            df_zscore_diff_between_datasets["sm"].loc[sm_indices[0],"zscore_diff"]=dists[file_names[0]][0].loc[sm_indices[0],"zscore"] - dists[file_names[1]][0].loc[sm_indices[0],"zscore"]
            df_zscore_diff_between_datasets["sm"].loc[sm_indices[1],"zscore_diff"]=dists[file_names[0]][0].loc[sm_indices[1],"zscore"] - dists[file_names[1]][0].loc[sm_indices[1],"zscore"]
            df_zscore_diff_between_datasets["op"].loc[op_indices[0],"zscore_diff"]=dists[file_names[0]][1].loc[op_indices[0],"zscore"] - dists[file_names[1]][1].loc[op_indices[0],"zscore"]
            df_zscore_diff_between_datasets["op"].loc[op_indices[1],"zscore_diff"]=dists[file_names[0]][1].loc[op_indices[1],"zscore"] - dists[file_names[1]][1].loc[op_indices[1],"zscore"]
        except KeyError:
            continue

    sns.distplot(df_zscore_diff_between_datasets["sm"].loc[:,"zscore_diff"], norm_hist=False, kde=False, ax=axs[2,0], bins=20)
    df_zscore_diff_between_datasets["sm"].to_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_same_aneuploidy_orientation_{}_{}.tsv".format(file_names[0], file_names[1])))
    axs[2,0].set_xlabel('z-score difference')
    axs[2,0].set_ylabel('# of co-occurrences')
    axs[2,0].set_title('z-score difference between same aneuploidy orientations'.format(file_names[1]))
    axs[2,0].legend(labels=["z-score difference"])
    df_zscore_diff_between_datasets["sm"].to_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_intra_zscore_diff_sm_{}_{}_{}.tsv".format(pval_method_name, *file_names)), sep='\t')


    sns.distplot(df_zscore_diff_between_datasets["op"].loc[:,"zscore_diff"], norm_hist=False, kde=False, ax=axs[2,1], bins=20)
    axs[2,1].set_xlabel('z-score difference')
    axs[2,1].set_ylabel('# of co-occurrences')
    axs[2,1].set_title('z-score difference between opposite aneuploidy orientations'.format(file_names[1]))
    axs[2,1].legend(labels=["z-score difference"])
    df_zscore_diff_between_datasets["op"].to_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_intra_zscore_diff_op_{}_{}_{}.tsv".format(pval_method_name, *file_names)), sep='\t')

    df_zscore_diff.to_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_intra_zscore_diff_sm_vs_op_{}_{}_{}.tsv".format(pval_method_name, *file_name)), sep='\t')

    plt.savefig(os.path.join(constants.OUTPUT_FOLDER, "figures", "hgs_intra_sm_vs_op_{}_{}_{}.png".format(pval_method_name, *file_names)))

    chr_interaction=constants.INTER
    df_summary_agg=pd.DataFrame()
    df_agg_pvals=pd.concat([dfs_between[0].loc[:,"pval"], dfs_between[1].loc[:,"pval"]] ,axis=1)
    df_agg_zscores=pd.concat([dfs_between[0].loc[:,"zscore"], dfs_between[1].loc[:,"zscore"]] ,axis=1)
    df_summary_agg["pval_ratio"]=(-np.log10(df_agg_pvals.iloc[:,0]))/(-np.log10(df_agg_pvals.iloc[:,1])+1)
    df_summary_agg["z_diff"]=df_agg_zscores.iloc[:,0]-df_agg_zscores.iloc[:,1]
    df_summary_agg["z_rank"]=rankdata(-df_summary_agg["z_diff"])
    df_summary_agg["pval_rank"]=rankdata(-df_summary_agg["pval_ratio"])
    df_summary_agg.to_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_summary_{}_{}.tsv".format(constants.CHR_INTERACTION_NAMES[chr_interaction], pval_method_name)), sep='\t')


def plot_dist(ax, dfs):
    bins = np.histogram(np.hstack([df.loc[:, "enrichment_score"] for df in dfs]),
                        bins=30)[1]
    sns.distplot(dfs[0].loc[:, "enrichment_score"], norm_hist=False, kde=False, ax=ax, bins=bins)
    sns.distplot(dfs[1].loc[:, "enrichment_score"], norm_hist=False, kde=False, ax=ax, bins=bins)
    ax.set_yscale('log')
    ax.set_xlabel('enrichment score (-log10(pval))')
    ax.set_ylabel('# of events')


def main(file_names):

    mats=[]
    labels=[]
    for fname in file_names:
        mat=pd.read_csv(os.path.join(constants.DATASETS_FOLDER, "{}.tsv".format(fname)), sep='\t', index_col=0)
        mat[pd.isna(mat)] = 0
        mats.append(mat)

    for cur_pval_method_name in ["enr", "dpl"]:
        generate_dist_for_dataset(mats, file_names, pval_method_name=cur_pval_method_name)

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--file_names', dest='file_names', default="brca_wgd_plus,brca_wgd_minus")

    args = parser.parse_args()
    file_names=args.file_names.split(",")
    main(file_names)
