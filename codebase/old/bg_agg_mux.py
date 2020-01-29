import argparse
from codebase.bg_agg_common import *

def main(file_names, n_start, n_end, p_factor, use_cache):

    interactions=None # ["19p1-20q1", "12q1-16p1", "2p-1-8q1", "12q1-20q1"]
    file_suffix="inter_mux"
    for file_name in file_names:
        original_enrichment_scores=pd.read_csv(os.path.join(constants.OUTPUT_FOLDER, ("hgs_{}_{}.tsv".format(file_suffix, file_name))), '\t', index_col=0)
        original_mat=pd.read_csv(os.path.join(constants.DATASETS_FOLDER,file_name+".tsv"), sep='\t', index_col=0)
        bg=generate_bg_dists(aneu_types=[(-1, 1),(1, -1), (-1, -1),(1, 1)], chr_interaction=INTER, file_name=file_name, n_start=n_start, n_end=n_end, original=original_mat, p_factor=p_factor, file_suffix=suffix_file_format, use_cache=use_cache)
        plot_dists(original_enrichment_scores, bg, interactions, file_name, p_factor, file_suffix, use_cache)

if __name__=="__main__":

    parser = argparse.ArgumentParser(description='args')
    parser.add_argument('--p_factor', dest='p_factor', default="80")
    parser.add_argument('--n_start', dest='n_start', default="0")
    parser.add_argument('--n_end', dest='n_end', default="10000")
    parser.add_argument('--file_names', dest='file_names', default="brca_wgd_plus")
    parser.add_argument('--use_cache', dest='use_cache', default="1")

    args = parser.parse_args()
    p_factor = int(args.p_factor)
    n_start= int(args.n_start)
    n_end= int(args.n_end)
    use_cache= int(args.use_cache)
    file_names=args.file_names.split(",")

    main(file_names, n_start=n_start, n_end=n_end, p_factor=p_factor, use_cache=use_cache)
