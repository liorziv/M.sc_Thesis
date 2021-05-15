from gs_model import *
import sys
import file_reader


'''
Main part
'''



def main():
    # plotting_data_by_names("C:/Users/Lior/Desktop/GS_model_fit/RES/222/all_res_222.csv",
    #                        "C:/Users/Lior/Documents/high_ind.csv",
    #                        "C:/Users/Lior/Desktop/GS_model_fit/data_222_esat_t_normalized.csv",222,
    #                        'C:/Users/Lior/Desktop/GS_model_fit/High_ind')
    check_GS([0,0,0,480,20])
    FILE_NAME, FILE_ID, OUTPUT_FOLDER = sys.argv[1], sys.argv[2], sys.argv[3]
    start, end = int(sys.argv[4]), int(sys.argv[5])
    ## represents the file id
    file_id = int(FILE_ID)
    ## file full path
    file_name = FILE_NAME
    ##expression matrix from csv file
    exp_mat = file_reader.get_gene_expression_matrix(file_name)
    # number of time repeats
    rep_num = file_reader.set_rep_num(file_id)
    # vector of gene names
    gene_names = file_reader.get_gene_names(file_name)
    # in order to save p-values and fitted theta
    p_values = []
    p_values_model = []
    fitted_params_list = np.zeros((end - start, 5))
    fitted_params_list_model = np.zeros((2 * (end - start), 5))
    str_lab = []
    # labels for the plot
    x_labels = ['Cont', '30 Min', '60 Min', '120 Min', '240 Min', '480 Min']

    counter = 0
    conter_2 = 0
    for i in range(start, end):  # 10072
        curr_gene = gene_names[i]
        real_y = exp_mat[i, :]  # this is a line from the expression matrix
        new_real_y, new_rep_num = real_y, rep_num  #
        mean_of_real_y_reps = add_mean_point(file_id, real_y, rep_num, x_labels)
        model_y = real_y + np.random.normal(10, 2, len(real_y))  # to avoid genes with low read number
        mean_of_model_y_reps = add_mean_point(file_id, model_y, rep_num, x_labels)
        # finds the fit params - try two runs of the model
        theta_hat_GS_real_1,rmse_1 = minimize_GS_model_fit(new_real_y, new_rep_num, file_id,mean_of_real_y_reps)
        theta_hat_GS_real_2,rmse_2 = minimize_GS_model_fit(new_real_y, new_rep_num, file_id, mean_of_real_y_reps)

        # getting the real data fitted y
        idx = np.argmin([rmse_1,rmse_2])
        theta_hat_GS_real = minimiztion_optimiztion(real_y,new_rep_num,file_id,mean_of_real_y_reps,theta_hat_GS_real)[idx]

        # getting the model fitted y
        theta_hat_GS_model = minimize_GS_model_fit(model_y, new_rep_num, file_id,mean_of_model_y_reps )[0] #[1.1636571972164833, 0, 0, 139.8649874483564, 13.090487859211333]

        # in order the save into the csv
#        fitted_params_list[counter, :] = theta_hat_GS_real
        fitted_params_list_model[conter_2, :] = theta_hat_GS_real
        fitted_params_list_model[conter_2 + 1, :] = theta_hat_GS_model
        conter_2 += 2
        str_lab.append('real_y')
        str_lab.append('real_y + 10')

        # calculate the y's according to the two models
        print('real- fit')

        reduced_y, y_fitted_GS, t_out, data_mean, _ = calc_models_fitted_y(new_real_y, file_id, new_rep_num,
                                                                               theta_hat_GS_real, curr_gene, True)
        # print('model- fit')
        reduced_y_model, y_fitted_GS_model, t_out_model, data_mean_model, _ = calc_models_fitted_y(model_y, file_id,
                                                                                                      new_rep_num, theta_hat_GS_model,
                                                                                                       curr_gene, True)

        p_val_real = 0 #F_test(reduced_y, y_fitted_GS[t_out], data_mean)
        p_val_model = 0 #F_test(reduced_y_model, y_fitted_GS_model[t_out_model], data_mean_model)

        plot_fit(real_y, y_fitted_GS, y_fitted_GS_model, data_mean, x_labels, new_rep_num,
                 gene_names[start + counter], p_val_real, p_val_model,
                 theta_hat_GS_real, OUTPUT_FOLDER, file_id)

        p_values.append(p_val_real)
        p_values_model.append(p_val_real)
        p_values_model.append(p_val_model)
        counter += 1

    df_2 = pd.DataFrame(fitted_params_list, columns=['c', 'h', 'b', 't0', 'sig'])
    df_1 = pd.DataFrame(gene_names[start:end], columns=['Genes'])
    df_1['p-val'] = p_values
    #    df_1['p-val_model'] = p_values_model
    res_table = pd.concat((df_1, df_2), axis=1)
    path_to_res = os.path.join(OUTPUT_FOLDER, "res_table_{s}_{e}.csv".format(s=start, e=end))
    df_model = pd.DataFrame(fitted_params_list_model, columns=['c', 'h', 'b', 't0', 'sig'])
    df_model['p-val'] = p_values_model
    df_3 = pd.DataFrame([gene_names[i] for i in range(start, end) for j in range(2)], columns=['Genes'])
    df_3['Y-values'] = str_lab
    df_model = pd.concat((df_3, df_model), axis=1)
    path_to_model = os.path.join(OUTPUT_FOLDER, "model_{s}_{e}.csv".format(s=start, e=end))
    df_model.to_csv(path_to_model, index=None)
    res_table.to_csv(path_to_res, index=None)

main()


