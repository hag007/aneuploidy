import matplotlib
matplotlib.use('Agg')

import seaborn as sns
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import codebase.constants
ANEU_DOUBLE_MINUS=0
ANEU_PLUS_MINUS=1
ANEU_DOUBLE_PLUS=2

WITHIN=0
BETWEEN=1

def get_count_dist(mat, aneu_type_a, aneu_type_b, chr_interaction, gd="wgd_plus"):
    co_occurrences=[]
    for a in np.arange(1,23):
        for arm_a in ['p','q']:
            for b in np.arange(a, 23):
                for arm_b in ['p', 'q']:
                    if (arm_a == arm_b and b == a) or (chr_interaction and b == a) or (not chr_interaction and b != a): continue
                    if "{}{}".format(a, arm_a) not in mat.columns: continue
                    if "{}{}".format(b, arm_b) not in mat.columns: continue

                    n_occurrence=np.sum(np.logical_or(np.logical_and(mat.loc[:,"{}{}".format(a,arm_a)].values==aneu_type_a,
                                            mat.loc[:,"{}{}".format(b,arm_b)].values==aneu_type_b),
                                            np.logical_and(mat.loc[:, "{}{}".format(a, arm_a)].values == aneu_type_b,
                                            mat.loc[:, "{}{}".format(b, arm_b)].values == aneu_type_a)))
                    co_occurrences.append(n_occurrence)

    sns.distplot(co_occurrences, norm_hist=False, kde=False)
    plt.title("n={}".format(len(co_occurrences)))
    chr_interaction_s = "between" if chr_interaction==0 else "within"
    aneu_type_a_s = "minus" if aneu_type_a == -1 else "plus"
    aneu_type_b_s = "minus" if aneu_type_a == -1 else "plus"
    plt.savefig(os.path.join(codebase.constants.OUTPUT_FOLDER, "figures", "counts_{}_{}_{}_{}.png".format(gd, chr_interaction_s, aneu_type_a_s, aneu_type_b_s)))
    plt.clf()







def main(mat, gd):
    #  within c + +
    get_count_dist(mat, aneu_type_a=1, aneu_type_b=1, chr_interaction=WITHIN, gd=gd)

    #  within c -
    get_count_dist(mat, aneu_type_a=-1, aneu_type_b=-1, chr_interaction=WITHIN, gd=gd)

    #  within c + -
    get_count_dist(mat, aneu_type_a=1, aneu_type_b=-1, chr_interaction=WITHIN, gd=gd)

    #  between c + +
    get_count_dist(mat, aneu_type_a=1, aneu_type_b=1, chr_interaction=BETWEEN, gd=gd)

    #  between c - -
    get_count_dist(mat, aneu_type_a=-1, aneu_type_b=-1, chr_interaction=BETWEEN, gd=gd)

    #  between c + -
    get_count_dist(mat, aneu_type_a=1, aneu_type_b=-1, chr_interaction=BETWEEN, gd=gd)


if __name__=="__main__":
    original_wgd_plus = pd.read_csv(os.path.join(codebase.constants.DATASETS_FOLDER, "brca_wgd_plus.tsv"), sep='\t', index_col=0)
    original_wgd_plus[pd.isna(original_wgd_plus)] = 0
    main(original_wgd_plus, gd="wgd_plus")

    original_wgd_minus = pd.read_csv(os.path.join(codebase.constants.DATASETS_FOLDER, "brca_wgd_minus.tsv"), sep='\t', index_col=0)
    original_wgd_minus[pd.isna(original_wgd_minus)] = 0
    main(original_wgd_minus, gd="wgd_minus")
