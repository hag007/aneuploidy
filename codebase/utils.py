import numpy as np
import pandas as pd
from scipy.stats import hypergeom, fisher_exact



def calc_hg_enrichment_pval(mat, a, arm_a, aneu_type_a, b, arm_b, aneu_type_b):
    n_overlap=np.sum(np.logical_and(mat.loc[:, "{}{}".format(a, arm_a)].values == aneu_type_a,
                                                     mat.loc[:, "{}{}".format(b, arm_b)].values == aneu_type_b))

    n_a=np.sum(np.logical_and(mat.loc[:, "{}{}".format(a, arm_a)].values == aneu_type_a,
                                 mat.loc[:, "{}{}".format(b, arm_b)].values != aneu_type_b))

    n_b=np.sum(np.logical_and(mat.loc[:, "{}{}".format(a, arm_a)].values != aneu_type_a,
                                       mat.loc[:, "{}{}".format(b, arm_b)].values == aneu_type_b))

    # pval=hypergeom.sf(n_overlap, mat.shape[0], n_overlap+n_a, n_overlap+n_b) \
    # + hypergeom.pmf(n_overlap, mat.shape[0], n_overlap+n_a, n_overlap+n_b)

    pval=hypergeom.sf(n_overlap, mat.shape[0], n_overlap+n_a, n_overlap+n_b) \
    + hypergeom.pmf(n_overlap, mat.shape[0], n_overlap+n_a, n_overlap+n_b) # n_overlap+n_a+n_overlap+n_b

    # tbl=[[n_overlap, n_b], [n_a, mat.shape[0]-(n_overlap+n_b+n_a)]]
    # pval_1=fisher_exact(tbl, 'greater')
    # if a==1 and arm_a=='p' and aneu_type_a==-1 and b==2 and  arm_b=='q' and aneu_type_b==-1:
    #     print (n_overlap, mat.shape[0], n_overlap+n_a, n_overlap+n_b)
    #     print pval, pval_1[1]


    return pval


def calc_hg_depletion_pval(mat, a, arm_a, aneu_type_a, b, arm_b, aneu_type_b):
    n_overlap=np.sum(np.logical_and(mat.loc[:, "{}{}".format(a, arm_a)].values == aneu_type_a,
                                                     mat.loc[:, "{}{}".format(b, arm_b)].values == aneu_type_b))

    n_a=np.sum(np.logical_and(mat.loc[:, "{}{}".format(a, arm_a)].values == aneu_type_a,
                                 mat.loc[:, "{}{}".format(b, arm_b)].values != aneu_type_b))

    n_b=np.sum(np.logical_and(mat.loc[:, "{}{}".format(a, arm_a)].values != aneu_type_a,
                                       mat.loc[:, "{}{}".format(b, arm_b)].values == aneu_type_b))

    pval=hypergeom.cdf(n_overlap, mat.shape[0], n_overlap+n_a, n_overlap+n_b)

    return pval

PVAL_METHODS={"enr" : calc_hg_enrichment_pval,
              "dpl" : calc_hg_depletion_pval}
