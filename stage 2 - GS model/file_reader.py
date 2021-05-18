import csv
import numpy as np
import consts as c
## Helps read Kcl 10 days



def get_gene_expression_matrix(file_path):
    '''
    loads gene expression matrix from esat file
    :param file_path: esat path
    :return: gene expression mat
    '''
    ## reader part
    reader = csv.reader(open(file_path,"r"),delimiter = ",")
    x = list(reader)
    iter_x = iter(x)
    repeats = np.array(next(iter_x)[1:]).astype("str")
    gene_exp_data = np.array([elem[1:len(repeats)+1] for elem in x])[1:,:].astype("float")

    return gene_exp_data



def get_gene_names(file_path):
    '''
    return all gene names in esat file  bbbb
    :param file_path: the input file
    :return: gene names list
    '''
    ## reader part
    reader = csv.reader(open(file_path,"r"),delimiter = ",")
    x = list(reader)
    gene_names = np.array([elem[0] for elem in x])[1:].astype("str")
    return gene_names


def set_rep_num(file_id):
    '''
    set number of repeats of specific experiment
    :param file_id: file id
    :return: number of repeats
    '''
    if file_id == c.FILE_ID_222:
        return 2
    if file_id in c.REP_3_LST:
        return 3
    if file_id == c.FILE_ID_222_21_10:
        return 5




