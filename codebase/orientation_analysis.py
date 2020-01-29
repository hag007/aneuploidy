import argparse
from bg_agg_common import *
from utils import PVAL_METHODS
from scipy.stats import chi2_contingency

def main(wgd_plus_file_name, wgd_minus_file_name):


    file_sm_format="emp_pvals_summary_intra_enr_sm_{}.tsv"
    file_op_format="emp_pvals_summary_intra_enr_op_{}.tsv"

    file_wgd_plus_sm=pd.read_csv(os.path.join(constants.OUTPUT_FOLDER, file_sm_format.format(wgd_plus_file_name)),sep='\t', index_col=0)
    file_wgd_minus_sm=pd.read_csv(os.path.join(constants.OUTPUT_FOLDER, file_sm_format.format(wgd_minus_file_name)),sep='\t', index_col=0)

    file_wgd_plus_op=pd.read_csv(os.path.join(constants.OUTPUT_FOLDER, file_op_format.format(wgd_plus_file_name)),sep='\t', index_col=0)
    file_wgd_minus_op=pd.read_csv(os.path.join(constants.OUTPUT_FOLDER, file_op_format.format(wgd_minus_file_name)),sep='\t', index_col=0)

    plot_zscore_diff(file_wgd_minus_op, file_wgd_minus_sm, file_wgd_plus_op, file_wgd_plus_sm)




def sm_op_chi_test(file_wgd_minus_op, file_wgd_minus_sm, file_wgd_plus_op, file_wgd_plus_sm):
    n_sig_sm_plus = len(file_wgd_plus_sm[file_wgd_plus_sm["qval"] < 0.05].index)
    n_sig_op_plus = len(file_wgd_plus_op[file_wgd_plus_op["qval"] < 0.05].index)
    n_sig_sm_minus = len(file_wgd_minus_sm[file_wgd_minus_sm["qval"] < 0.05].index)
    n_sig_op_minus = len(file_wgd_minus_op[file_wgd_minus_op["qval"] < 0.05].index)
    obs = [[n_sig_sm_plus, n_sig_op_plus], [n_sig_sm_minus, n_sig_op_minus]]
    print "obs: {}".format(obs)
    chi2, p, dof, ex = chi2_contingency(obs, correction=False)
    return p

def plot_zscore_diff(file_wgd_minus_op, file_wgd_minus_sm, file_wgd_plus_op, file_wgd_plus_sm):

    fig, axs = plt.subplots(2,2, figsize=(20, 20))
    bins=np.histogram(np.hstack((file_wgd_plus_sm.loc[:,"enrichment_score"],file_wgd_plus_op.loc[:,"enrichment_score"])), bins=20)[1]
    sns.distplot(file_wgd_plus_sm.loc[:,"enrichment_score"], norm_hist=False, kde=False, ax=axs[0][0], bins=bins)
    sns.distplot(file_wgd_plus_op.loc[:,"enrichment_score"], norm_hist=False, kde=False, ax=axs[0][0], bins=bins)
    axs[0][0].set_yticks(np.arange(6)*5)
    axs[0][0].set_ylim(0,30)
    axs[0][0].legend(["same", "opposite"])
    axs[0][0].set_title("enrichment score distribution of same and opposite orientation\nin WGD plus dataset")
    bins=np.histogram(np.hstack((file_wgd_minus_sm.loc[:,"enrichment_score"],file_wgd_minus_op.loc[:,"enrichment_score"])), bins=20)[1]
    sns.distplot(file_wgd_minus_sm.loc[:,"enrichment_score"], norm_hist=False, kde=False, ax=axs[0][1], bins=bins)
    sns.distplot(file_wgd_minus_op.loc[:,"enrichment_score"], norm_hist=False, kde=False, ax=axs[0][1], bins=bins)
    axs[0][1].set_yticks(np.arange(6)*5)
    axs[0][1].set_ylim(0,30)
    axs[0][1].legend(["same", "opposite"])
    axs[0][1].set_title("enrichment score distribution of same and opposite orientation\nin WGD minus dataset")
    df=pd.DataFrame()
    for a in np.arange(1,23):

        sm_indices = ["{}q1-{}p1".format(a,a), "{}q-1-{}p-1".format(a,a)]
        op_indices = ["{}q-1-{}p1".format(a,a), "{}q1-{}p-1".format(a,a)]

        try:
            plus_diff=np.max(file_wgd_plus_sm.loc[sm_indices,"zscore"]) - np.max(file_wgd_plus_op.loc[op_indices,"zscore"])
            minus_diff=np.max(file_wgd_minus_sm.loc[sm_indices,"zscore"]) - np.max(file_wgd_minus_op.loc[op_indices,"zscore"])
        except KeyError:
            continue

        agg_diff=plus_diff - minus_diff

        df.loc[a,"plus_diff"]=plus_diff
        df.loc[a,"minus_diff"]=minus_diff
        df.loc[a,"agg_diff"]=agg_diff

    bins=np.histogram(np.hstack((df.loc[:,"plus_diff"],df.loc[:,"minus_diff"])), bins=20)[1]
    sns.distplot(df.loc[:,"plus_diff"], norm_hist=False, kde=False, ax=axs[1][0], bins=bins)
    sns.distplot(df.loc[:,"minus_diff"], norm_hist=False, kde=False, ax=axs[1][0], bins=bins)
    axs[1][0].legend(["wgd plus", "wgd minus"])
    s,p_wilcoxon = wilcoxon(df.loc[:,"plus_diff"], df.loc[:,"minus_diff"])
    p_chi=sm_op_chi_test(file_wgd_minus_op, file_wgd_minus_sm, file_wgd_plus_op, file_wgd_plus_sm)
    axs[1][0].set_title("max-zscore differences between same and opposite orientation\n(chi p-value:{} wilcoxon p-value:{})".format("%.2E" % p_chi, "%.2E" % p_wilcoxon))

    sns.distplot(df.loc[:,"agg_diff"], norm_hist=False, kde=False, ax=axs[1][1], bins=20)
    axs[1][1].legend(["diff of zscore_diff"])
    axs[1][1].set_title("difference of max-zscore differences")
    df.to_csv(os.path.join(constants.OUTPUT_FOLDER, "orientation_analysis.tsv"), sep='\t')
    plt.savefig(os.path.join(constants.OUTPUT_FOLDER, "figures", "orientation_analysis.png"))



if __name__=="__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--wgd_plus_file_name', dest='wgd_plus_file_name', default="brca_wgd_plus")
    parser.add_argument('--wgd_minus_file_name', dest='wgd_minus_file_name', default="brca_wgd_minus")


    args = parser.parse_args()
    wgd_plus_file_name= args.wgd_plus_file_name
    wgd_minus_file_name= args.wgd_minus_file_name

    main(wgd_plus_file_name, wgd_minus_file_name)
