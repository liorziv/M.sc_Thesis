import gs_model as gs
import sys
import file_reader
import consts as c

'''
Main part
'''
def write_res_to_csv(fitted_params_list, fitted_params_list_noised, gene_names, p_values, p_values_noised, start, end,
                     output_folder):
    '''
    saves all the results into csv
    :param fitted_params_list: a list of all model best params per gene
    :param fitted_params_list_noised: same but for noised y's
    :param gene_names: list of gene names
    :param p_values: p-values results for real y values
    :param p_values_noised:"" for noised y
    :param start: start gene
    :param end: end gene
    :param output_folder: output folder to save csv's
    '''
    df_1 = pd.DataFrame(gene_names[start:end], columns=['Genes'])

    df_2 = pd.DataFrame(fitted_params_list, columns=['c', 'h', 'b', 't0', 'sig'])
    df_1['p-val'] = p_values
    df_real = pd.concat((df_1, df_2), axis=1)

    df_noised = pd.DataFrame(fitted_params_list_noised, columns=['c', 'h', 'b', 't0', 'sig'])
    df_noised['p-val'] = p_values_noised
    df_noised = pd.concat((df_1, df_noised), axis=1)

    path_to_res = os.path.join(output_folder, "real_y_{s}_{e}.csv".format(s=start, e=end))
    path_to_model = os.path.join(output_folder, "noised_y_{s}_{e}.csv".format(s=start, e=end))
    df_noised.to_csv(path_to_model, index=None)
    df_real.to_csv(path_to_res, index=None)

def main():
    file_name, file_id, output_folder = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    start, end = int(sys.argv[4]), int(sys.argv[5])
    ##expression matrix from csv file
    exp_mat = file_reader.get_gene_expression_matrix(file_name)
    # number of time repeats
    rep_num = file_reader.set_rep_num(file_id)
    # vector of gene names
    gene_names = file_reader.get_gene_names(file_name)
    # in order to save p-values and fitted theta
    p_values = []
    p_values_model = []
    fitted_params_list = np.zeros((end - start, c.GS_MODEL_PARAM_NUM))
    fitted_params_list_noised = np.zeros((2 * (end - start), c.GS_MODEL_PARAM_NUM))
    str_lab = []

    plt_ctr = 0
    data_ctr = 0
    #go over each gene and find GS model fit result
    for i in range(start, end):
        curr_gene = gene_names[i]
        real_y = exp_mat[i, :]  # this is a gene from the expression matrix
        mean_of_real_y_reps = add_mean_point(file_id, real_y, rep_num, c.X_LABELS)
        # to avoid genes with low read number we add noise to the sample and see if the fit is still significant
        noised_y = real_y + np.random.normal(10, 2, len(real_y))
        mean_of_noised_y_reps = gs.add_mean_point(file_id, noised_y, rep_num, c.X_LABELS)
        # finds the fit params - try two runs of the model
        theta_hat_GS_real = gs.minimize_GS_model_fit(real_y, rep_num, file_id,mean_of_real_y_reps)[0]

        # continue on the best solution searching at its range
        theta_hat_GS_real = gs.minimiztion_optimiztion(real_y,rep_num,file_id,mean_of_real_y_reps,theta_hat_GS_real)
        theta_hat_GS_noised = gs.minimize_GS_model_fit(noised_y, rep_num, file_id,mean_of_noised_y_reps )[0]
        # in order the save into the csv
        fitted_params_list[plt_ctr, :] = theta_hat_GS_real
        fitted_params_list_noised[plt_ctr, :] = theta_hat_GS_noised

        # calculate the y's according to the two models(real + noised)
        reduced_y, y_fitted_GS, t_out, data_mean, _ = gs.calc_models_fitted_y(real_y, file_id, rep_num,
                                                                               theta_hat_GS_real, curr_gene, True)
        reduced_y_noised, y_fitted_GS_model, t_out_model, data_mean_model, _ = gs.calc_models_fitted_y(noised_y, file_id,
                                                                                                      rep_num, theta_hat_GS_noised,
                                                                                                       curr_gene, True)

        #calculate p-value using F-test
        p_val_real = gs.F_test(reduced_y, y_fitted_GS[t_out], data_mean)
        p_val_noised = gs.F_test(reduced_y_noised, y_fitted_GS_model[t_out_model], data_mean_model)

        gs.plot_fit(real_y, y_fitted_GS, y_fitted_GS_model, data_mean,  c.X_LABELS, rep_num,
                 gene_names[start + plt_ctr], p_val_real, p_val_noised,
                 theta_hat_GS_real, output_folder, file_id)

        p_values.append(p_val_real)
        p_values_model.append(p_val_real)
        p_values_model.append(p_val_noised)
        plt_ctr += 1
    write_res_to_csv(fitted_params_list, fitted_params_list_noised, gene_names, p_values, p_values_model, start, end,
                     output_folder)


main()


