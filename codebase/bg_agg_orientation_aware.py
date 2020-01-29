import argparse
from bg_agg_common import *
from utils import PVAL_METHODS

def main(file_names, n_start, n_end, pval_method, p_factor, use_cache):


    for file_name in file_names:
        original_mat=pd.read_csv(os.path.join(constants.DATASETS_FOLDER,file_name+".tsv"), sep='\t', index_col=0)

        file_suffix="intra_enr_sm"
        original_enrichment_scores=pd.read_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_emp_dist_{}_{}.tsv".format(file_suffix, file_name)), '\t', index_col=0)
        interactions=sorted(list(original_enrichment_scores.index)) # ["19p1-20q1", "12q1-16p1", "2p-1-8q1", "12q1-20q1"]
        bg_similar=generate_bg_dists(aneu_types=[(-1, -1),(1, 1)], chr_interaction=INTRA, file_name=file_name, n_start=n_start, n_end=n_end, original=original_mat, pval_method=pval_method, p_factor=p_factor, file_suffix=file_suffix, use_cache=use_cache)
        intra_interactions=["{}{}{}1-{}{}{}1".format(a,b[0],c,a,b[1],c) for a in np.arange(1,24) for b in [("p","q"),("q","p")] for c in ["","-"]]
        intra_interactions=[a for a in intra_interactions if a in interactions]
        output_sim=plot_dists(original_enrichment_scores, bg_similar, intra_interactions, file_name, p_factor, file_suffix, use_cache)
        output_sim.rename(columns={"pval":"pval_sm", "qval":"qval_sm", "zscore":"zscore_sm"})

        file_suffix="intra_enr_op"
        original_enrichment_scores=pd.read_csv(os.path.join(constants.OUTPUT_FOLDER, "hgs_emp_dist_{}_{}.tsv".format(file_suffix, file_name)), '\t', index_col=0)
        interactions=sorted(list(original_enrichment_scores.index)) # ["19p1-20q1", "12q1-16p1", "2p-1-8q1", "12q1-20q1"]
        bg_opposite=generate_bg_dists(aneu_types=[(-1, 1),(1, -1)], chr_interaction=INTRA, file_name=file_name, n_start=n_start, n_end=n_end, original=original_mat, pval_method=pval_method, p_factor=p_factor, file_suffix=file_suffix, use_cache=use_cache)
        intra_interactions=["{}{}{}1-{}{}{}1".format(a,b[0],c[0],a,b[1],c[1]) for a in np.arange(1,24) for b in [("p","q"),("q","p")] for c in [("","-"), ("-","")]]
        intra_interactions=[a for a in intra_interactions if a in interactions]
        output_op=plot_dists(original_enrichment_scores, bg_opposite, intra_interactions, file_name, p_factor, file_suffix, use_cache)
        output_op.rename(columns={"pval":"pval_op", "qval":"qval_op", "zscore":"zscore_op"})

        output_agg=pd.concat([output_sim, output_op], axis=1)
        output_agg.to_csv(os.path.join(constants.OUTPUT_FOLDER, "emp_pvals_intra_agg_{}.tsv".format(file_name)), sep='\t')

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--p_factor', dest='p_factor', default="80")
    parser.add_argument('--n_start', dest='n_start', default="0")
    parser.add_argument('--n_end', dest='n_end', default="10000")
    parser.add_argument('--file_names', dest='file_names', default="brca_wgd_plus,brca_wgd_minus")
    parser.add_argument('--use_cache', dest='use_cache', default="1")
    parser.add_argument('--pval_method_name', dest='pval_method_name', default="enr")

    args = parser.parse_args()
    p_factor = int(args.p_factor)
    n_start= int(args.n_start)
    n_end= int(args.n_end)
    use_cache= int(args.use_cache)
    pval_method_name= args.pval_method_name
    file_names=args.file_names.split(",")

    main(file_names, n_start=n_start, n_end=n_end, pval_method=PVAL_METHODS[pval_method_name], p_factor=p_factor, use_cache=use_cache)
