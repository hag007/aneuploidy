import argparse
from bg_agg_common import *
from utils import PVAL_METHODS

def main(file_names, n_start, n_end, chr_interaction, pval_method_name, th, p_factor, use_cache):

    original_summary=pd.read_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_summary_{}_{}.tsv".format(constants.CHR_INTERACTION_NAMES[chr_interaction], pval_method_name)), '\t', index_col=0)
    if chr_interaction==constants.INTER:
        if th <0:
            interactions=list((original_summary["z_diff"][original_summary["z_diff"]<th]).index)
        else:
            interactions=list((original_summary["z_diff"][original_summary["z_diff"]>th]).index)

    else:
        interactions=list(original_summary.index)

    print "total number of interactions: {}".format(len(original_summary.index))

    for file_name in file_names:
        file_suffix="{}_{}_{}".format(constants.CHR_INTERACTION_NAMES[chr_interaction], pval_method_name,th)

        original_mat=pd.read_csv(os.path.join(constants.DATASETS_FOLDER,file_name+".tsv"), sep='\t', index_col=0)
        bg=generate_bg_dists(aneu_types=[(-1, 1), (1, -1), (-1, -1), (1, 1)], chr_interaction=chr_interaction, file_name=file_name, n_start=n_start, n_end=n_end, original=original_mat, pval_method=PVAL_METHODS[pval_method_name], p_factor=p_factor, file_suffix=file_suffix, use_cache=use_cache)

        original_enrichment_scores=pd.read_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_emp_dist_{}_{}.tsv".format("_".join(file_suffix.split('_')[:-1]), file_name)), '\t', index_col=0)
        plot_dists(original_enrichment_scores, bg, interactions, file_name, p_factor, file_suffix, use_cache)


if __name__=="__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--p_factor', dest='p_factor', default="80")
    parser.add_argument('--n_start', dest='n_start', default="0")
    parser.add_argument('--n_end', dest='n_end', default="10000")
    parser.add_argument('--file_names', dest='file_names', default="brca_wgd_plus,brca_wgd_minus")
    parser.add_argument('--use_cache', dest='use_cache', default="1")
    parser.add_argument('--pval_method_name', dest='pval_method_name', default="enr")
    parser.add_argument('--chr_interaction', dest='chr_interaction', default=constants.INTER)
    parser.add_argument('--th', dest='th', default="-4")

    args = parser.parse_args()
    p_factor = int(args.p_factor)
    n_start= int(args.n_start)
    n_end= int(args.n_end)
    th=int(args.th)
    chr_interaction= int(args.chr_interaction)
    use_cache= int(args.use_cache)
    pval_method_name= args.pval_method_name
    file_names=args.file_names.split(",")

    main(file_names, n_start=n_start, n_end=n_end, chr_interaction=chr_interaction, pval_method_name=pval_method_name, th=th, p_factor=p_factor, use_cache=use_cache)
#