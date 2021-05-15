import csv
import numpy as np

## Helps read Kcl 10 days
REP_222  = 2
File_ID_222 = 222
IDX_222_1 = np.linspace(0, 12, 6, endpoint=False, dtype='int')
IDX_222_2 = np.linspace(1, 11, 6, endpoint=True, dtype='int')

## Helps read Pic 10 and 14 days
REP_21 = 3
File_ID_21_10 = 2110
IDX_21_10_1 = np.linspace(0, 18, 6, endpoint=False, dtype='int')
IDX_21_10_2 = np.linspace(1, 16, 6, endpoint=True, dtype='int')
IDX_21_10_3 = np.linspace(2, 17, 6, endpoint=True, dtype='int')

File_ID_21_14 = 2114
File_ID_21142 = 21142
IDX_21_14_1 = [0,3,6,9,10,13]
IDX_21_14_2 = [1,4,7,9,11,14]
IDX_21_14_3 = [2,5,8,9,12,15]

File_ID_90 = 90
IDX_21_14_1 = [0,3,6,9,10,13]
IDX_21_14_1 = [1,4,7,9,11,14]
IDX_21_14_1 = [2,5,8,9,12,15]

REP_222_21_10 = 5
File_ID_222_21_10 = 2222110
File_ID_89 = 89
File_ID_891 = 891

#"./proj/diff_eq/data_21_10_esat_t_normalized.csv","r"),delimiter = ","

def get_gene_expression_matrix(file_path):
    '''
    loads gene expression matrix from esat file
    :param file_path: esat path
    :return: gene expression mat
    '''
    ## reader part
    reader = csv.reader(open(file_path,"r"),delimiter = ",")
    #reader = pd.read_csv(file_path,encoding="ISO-8859-1")
    x = list(reader)
    gene_names = np.array([elem[0] for elem in x])[1:].astype("str")
    iter_x = iter(x)
    repeats = np.array(next(iter_x)[1:]).astype("str")
    gene_exp_data = np.array([elem[1:len(repeats)+1] for elem in x])[1:,:].astype("float")


    return gene_exp_data


def sort_by_repeats(file_id,expression_mat):
    '''
    sort the gene expression matrix by repeats(some experiment have, 2/3)
    :param file_id: the experiment id
    :param expression_mat: the input gene expression mat
    :return: reshaped gene expression mat
    '''

    if (file_id == File_ID_222):

        return np.concatenate((expression_mat[:,IDX_222_1],expression_mat[:,IDX_222_2]),0)

    if (file_id == File_ID_21_10):
        return np.concatenate((expression_mat[:, IDX_21_10_1], expression_mat[:, IDX_21_10_2],expression_mat[:, IDX_21_10_3]), 0)

    if(file_id == File_ID_21_14):
        return np.concatenate((expression_mat[:, IDX_21_14_1], expression_mat[:, IDX_21_14_2],expression_mat[:, IDX_21_14_3]), 0)
    else:
        return np.concatenate(
            (expression_mat[:, IDX_21_14_1], expression_mat[:, IDX_21_14_2], expression_mat[:, IDX_21_14_3]), 0)



def get_gene_names(file_path):
    '''
    return all gene names in esat file
    :param file_path: the input file
    :return: gene names list
    '''
    ## reader part
    reader = csv.reader(open(file_path,"r"),delimiter = ",")
    x = list(reader)
    gene_names = np.array([elem[0] for elem in x])[1:].astype("str")


    return gene_names

def get_repeats(file_path):
    '''

    :param file_path:
    :return:
    '''
    ## reader part
    reader = csv.reader(open(file_path,"r"),delimiter = ",")
    x = list(reader)
    iter_x = iter(x)

    repeats = np.array(next(iter_x)[1:]).astype("str")

    return repeats

def set_rep_num(file_id):
    '''
    set number of repeats of specific experiment
    :param file_id: file id
    :return: number of repeats
    '''
    if(file_id == File_ID_222):
        return REP_222
    if (file_id == File_ID_21_10  or file_id == File_ID_21_14 or file_id == File_ID_90
        or file_id == File_ID_89 or file_id == File_ID_891 or file_id == File_ID_21142):
        return REP_21
    if(file_id == File_ID_222_21_10):
        return REP_222_21_10

def get_genes_list(genes_file):
    reader = csv.reader(open(genes_file, "r"), delimiter=",")
    x = list(reader)
    return np.array([elem[1] for elem in x])[1:].astype("str")

def get_gauss_params(file_path):
    reader = csv.reader(open(file_path, "r"), delimiter=",")
    x = list(reader)
    return np.array([elem[1:7] for elem in x]).astype("float")

def experiment_name(exp_id):
    if (exp_id == 222):
        return 'KCl 10'
    elif (exp_id == 2110):
        return 'Pic 10'
    elif (exp_id == 2114):
        return 'Pic 14(21_14)'
    elif (exp_id == 21142):
        return 'Pic 14'
    elif (exp_id == 892):
        return 'Pic 14(89)'
    elif (exp_id == 89):
        return 'Pic 14(89)'
    elif (exp_id == 891):
        return 'Pic 14(89) one 480'
    elif (exp_id == 90):
        return 'KCl 14'
