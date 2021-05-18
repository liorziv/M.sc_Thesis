import numpy as np
from scipy.integrate import odeint
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import scipy.optimize as optimize
import scipy.stats
from scipy import signal
import pandas as pd
import csv
import os
import consts as c


def GS_model(theta, t_out_i, y0):
    '''
    Defines the model of gene expression levels in time
    using differential equations
    :param theta: the gs model params
    c       - basel transcription level
    h       - transcription level (coefficient for the GS)
    b       - degradation level
    t0      - peak place (mean)
    sig     - peak width (variance)
    :param t_out_i: time points of the result (to output)
    :param y0: the first value of y
    :return: fitted y by the model
    '''
    t_full = np.linspace(0, 483, 484, dtype='int')
    y_hat_1 = odeint(model, y0, t_full, args=(theta,))

    return y_hat_1[t_out_i]


def model(y, t, theta):
    '''
    Define the GS model differential equation
    :param y: amount of reads speard in a few repeats in time
    :param t: time
    :param theta:  (as in GS model)  y,t,c,h,b,mu,sig
    :return:
    '''
    c, h, b, t0, sig = theta[0], theta[1], theta[2], theta[3], theta[4]
    at = c + h * np.exp(-0.5 * (((t0 - t) / sig) ** 2))  # c
    dydt = at - b * y
    return dydt


def calc_rmse(real_y, fitted_y, peak_of_fit):
    '''
    calc the RMSE with a little help of the ratio between the peak itself
    usually when the fit is perfect for a gene we wouldn't believe it will
    happen with a small sigma and a really large peak at a place there is
    no data (meaning we wont use it in the computation of the RMSE)
    therefore is is added here small bug here is that according to this it will try to fit a really small peak
    to the model compared to the max of the real data - not good
    but luckily the entire score depends more on the fit of all the rest of the points
    therefore it it will only relay on the difference of the peak at real from the fitted y
    the RMSE will be bad. Therefore it probably tries to turn this ration close to 1.

    :param real_y: the gene expression data
    :param fitted_y: y we get from the model best fit
    :param peak_of_fit: maximal value of fit
    :return: rmse between those two with a small penalty to avoid mistakes of large gaps between real to fit peak
    '''
    return np.sqrt(np.mean((real_y[:, None] - fitted_y) ** 2)) + (peak_of_fit / np.max(real_y))



def calc_RSS(real_y, fitted_y):
    '''
    calculets the RSS of the fitted y values and the real gene expression levels
    :param real_y: the real gene expression data
    :param fitted_y: y we get from fitting the gs model
    :return: rss result
    '''
    return np.sum(np.square(fitted_y - real_y[:, None]))


def find_GS_model_fit(theta, real_y, rep_num, file_id):
    '''
    find the best parameters to fit the model
    to our gene expression data using fmin
    :param theta: the initial params guess
    :param real_y: the actual gene expression
    :param rep_num: number of repeats per time frame
    :param file_id: experiment id
    :return: rmse of the best guess found by f.min
    '''
    t_out = choose_time_line(file_id)[1]  # plot of real data timeline
    penalty = np.inf  # if c is negative

    # here the fitted y is calculated by the GS
    # differential equation model

    y_hat = GS_model(theta, np.linspace(0, 482, 483, dtype='int'), np.mean(real_y[:rep_num]))  # y0= c/b
    max_of_fitted_y = np.max(y_hat)
    y_hat = y_hat[t_out]  # get only the values at the real data index's

    # check that we dont have an outlier (fitted to the model but really show not affect
    # b- degradation rate is limited 0<= b <= 0.4
    # sig - the peak width sig >= 6 ( we get really thin peaks which can be fitted to even to
    # the smallest peaks [0,0,0,1,0,0] which of course we dont want.
    initial_cond = (
                theta[0] >= 0 and theta[1] >= 0.0 and theta[2] >= 0.0 and theta[2] <= 0.4 and theta[3] > 0 and theta[
            4] >= 6)
    sigma_sm_6_deg_gt_0003 = (theta[4] < 7 and theta[2] > 0.03)  # a condotion we dont want to have in our data
    ## by the second cond it seems that the sig is limited to >= 7 and 0<= b <=0.03
    if (initial_cond and not sigma_sm_6_deg_gt_0003):
        penalty = 0
    # if all the condition are met then the penalty is zeroed and the fit might work
    RMSE = calc_rmse(real_y, y_hat, max_of_fitted_y)

    # just for sanity checks
    # total_opt_output_list.append(RMSE + penalty)
    # opt_params.append(theta)
    # y_hat_list.append(y_hat)
    return RMSE + penalty  # this is what optimize try to minimize



def choose_time_line(exp_id):

    '''
    Define the time points for the data or the plot
    :param exp_id: a unique number which represents a specific experiment
    :return: time points array
    '''
    res = []
    if exp_id == c.FILE_ID_222:
        res = [c.TIME_LINE_2, c.TIME_LINE_2]
    elif exp_id == c.FILE_ID_21_10:
        res = [c.TIME_LINE_3, c.TIME_LINE_3]
    elif exp_id == c.FILE_ID_21_14:
        res = [c.TIME_LINE_3, np.concatenate((c.TIME_LINE_3[:10], c.TIME_LINE_3[12:]), axis=0)]
    elif exp_id == c.FILE_ID_89:
        res = [c.TIME_LINE_3, np.concatenate((c.TIME_LINE_3[:11], c.TIME_LINE_3[13:17]), axis=0)]
    elif exp_id == c.FILE_ID_90:
        res = [c.TIME_LINE_3,
               np.concatenate((c.TIME_LINE_3[:2],
                               c.TIME_LINE_3[3:10], c.TIME_LINE_3[12:17]),
                              axis=0)]
    return res




def minimize_GS_model_fit(real_y, rep_num, file_id, mean_values):
    '''
    Create 7 guesses for all 5 parameters
    :param real_y: real y data(source)
    :param rep_num: number of repeats per time frame
    :param file_id: unique id per experiment
    :param mean_values: the mean of repeats in each time point
    :return:
    '''
    # takes the first peak
    peaks = signal.find_peaks(mean_values)[0]
    peaks_val = mean_values[peaks]

    t0_idx = np.argmax(mean_values) * rep_num  # gets the peak index in the real_y_i array
    t0_1 = choose_time_line(file_id)[0][t0_idx]  # get the x(time) of the peak(for the plot)
    if len(peaks) == 0:
        t0_2 = choose_time_line(file_id)[0][len(real_y) - 1]  # the peak is at 480
    else:
        induction_start = np.argmax(peaks_val)
        t0_2 = choose_time_line(file_id)[0][max(peaks[induction_start] * 2 - 1, 0)]  # choose_time_line(rep_num)[t0_idx]

    b_1 = 0.001
    c_1 = np.mean(real_y[:rep_num])
    h_1 = 0.01
    sig = 10

    # a bit more induction
    b_2 = 0.2
    h_2 = 2
    h_3 = 0.2 * max(real_y)

    # high induction
    c_3 = np.mean(real_y[:rep_num]) / np.mean(real_y)
    h_3 = 10

    #high induction
    c_7 = 100
    h_7 = (np.max(real_y) - np.min(real_y)) / 2
    b_7 = 0.399

    initial_guess = np.zeros((7,5))
    initial_guess[0] = [c_1, h_1, b_1, t0_1, sig]  # basic - mean basal level
    initial_guess[1] = [c_1, h_2, b_2, t0_1, sig]  # higher degradation and induction
    initial_guess[2] = [c_3, h_3, b_2, t0_1, sig]  # mean of cont / mean of all as basel level, high induction
    initial_guess[3] = [c_1, h_3, b_1, choose_time_line(file_id)[0][t0_idx], 200]
    initial_guess[4] = [c_1, h_3, b_1, t0_2, 8]
    initial_guess[5] = [c_1, 0, 0, 480, sig]  # only going up
    initial_guess[6] = [c_7, h_7, b_7, t0_1, sig]  # genes which are usually not good(not affected by the treatment)

    # fmin try to minimize find GS model func which returns the RMSE of the model, real data.
    return find_best_fit(initial_guess, real_y,rep_num, file_id)


def find_best_fit(initial_guess, real_y, rep_num, file_id):
    '''
    get array of guesses as initialization for the fit, and find the best params using fmin
    :param initial_guess: number of initial guesses
    :param real_y: the actual gene expression
    :param rep_num: number of repaets per time frame
    :param file_id: id of experiment
    :return: the best params and their rmse
    '''
    num_of_guesses = initial_guess.shape[0]
    param_res = np.zeros((1, num_of_guesses))
    rmse_res = np.zeros((1, num_of_guesses))
    for i in range(num_of_guesses):
        param_res[i], rmse_res[i] = optimize.fmin(find_GS_model_fit, initial_guess[i], args=(real_y, rep_num, file_id),
                                                  full_output=True)[0:2]

    # find the guess with the smallest RMSE
    idx = np.argmin(rmse_res)
    return param_res[idx], rmse_res[idx]

def minimiztion_optimiztion(real_y, rep_num, file_id, theta):
    '''
    take the best guess(lowest RMSE) and search around it
    :param real_y: real y data(source)
    :param rep_num: number of repeats per time frame
    :param file_id: unique id per experiment
    :param theta: params of the best fit
    :return: new best params
    '''
    # trying to optimize the guess by moving a bit in each of the variables
    c, h, b, t0, sig = theta
    initial_guess[0] = np.zeros((7,5))
    initial_guess[1] = [c, h / 2, b, t0, sig + 10]
    initial_guess[2] = [c, h / 2, b, t0, sig + 20]
    initial_guess[3] = [c, h / 2, b, t0 - 10, sig + 10]
    initial_guess[4] = [c, h / 2, b, t0 - 20, sig + 20]
    initial_guess[5] = [c, h, b, t0, sig]

    return find_best_fit(initial_guess, real_y, rep_num, file_id)


def plot_fit(y_real, y_fitted, y_fitted_noised,
             data_mean, id, p_val,
             p_val_noised, theta, fd_name, file_id):
    '''
    plot  both of the models with the original data
    :param y_real: real gene expression values
    :param y_fitted: gs model fit gene expression values
    :param y_fitted_noised:  "" with noised y
    :param data_mean: mean of the real y
    :param id: experiment id
    :param p_val: real y p-values
    :param p_val_noised: noised
    :param theta: best fit params of gs model
    :param fd_name:
    :return:
    '''
    # plot the results
    plt.close("all")
    t_plot = choose_time_line(file_id)[1]
    ax1 = plt.subplot(111)
    ax1.plot(np.linspace(0, len(y_fitted), len(y_fitted)),
             10 + 5 * np.exp(-0.5 * (((60 - np.linspace(0, len(y_fitted), len(y_fitted))) / 10) ** 2)), label='fit')
    ax1.plot(t_plot, y_real, 'ro', label='data')
    ax1.plot(t_plot, data_mean, 'g+:', label='Mean')
    ax1.plot(np.linspace(0, len(y_fitted_noised), len(y_fitted_noised)), y_fitted_noised, 'b+:', label='fit+10')
    ax1.set_xlabel('Time')
    ax1.set_xticks(range(0, 481, 60))  # 30
    step = 10
    if (10 + max(np.max(y_fitted), np.max(y_real)) > 280 and max(np.max(y_fitted), np.max(y_real)) > 0):
        step = 20
    elif (10 + max(np.max(y_fitted), np.max(y_real)) > 480 and max(np.max(y_fitted), np.max(y_real)) > 0):
        step = 50

    ax1.set_yticks(range(0, int(10 + max(np.max(y_fitted), np.max(y_real))), step))
    ax1.set_xbound([-10, 490])
    ax1.set_ybound([0, 20])

    ax1.set_ylabel('Read Count')
    ax1.set_title('p-value: ' + str("%g" % p_val) + ' p-value +10: ' + str("%g" % p_val_noised) + '\n basal(c): ' + str(
        "%g" % round(theta[0], 3)) + ' ind(h): ' + str(
        "%g" % round(theta[1], 3)) + '\n ' + 'deg(b): ' + str("%g" % round(theta[2], 3)) + ' t0: ' + str(
        "%g" % round(theta[3], 3)) + ' sig: ' + str("%g" % round(int(theta[4]), 3)))
    plt.savefig(fd_name + str(id) + '.png', bbox_inches='tight')
    plt.close('all')


#################################3'
''' 
 Two more plot functions
'''
def plot_only_fit(y_real, y_fitted, data_mean, id, p_val, theta, fd_name, file_id, folder_name):
    '''
    plot only the fit data(no mean)
    '''
    # plot the results
    plt.close("all")
    t_plot = choose_time_line(file_id)[1]
    fig = plt.figure(figsize=(12, 3))

    ax1 = plt.subplot(111)
    ax1.plot(np.linspace(0, len(y_fitted), len(y_fitted)), y_fitted, label='fit')
    ax1.plot(t_plot, y_real, 'ro', label='data')
    ax1.plot(t_plot, data_mean, 'g+:', label='Mean')
    # ax1.plot(np.linspace(0, len(y_fitted_model), len(y_fitted_model)), y_fitted_model, 'b+:', label='fit+10')
    ax1.set_xlabel('Time')
    ax1.set_xticks(range(0, 481, 30))
    step = 10
    if (10 + max(np.max(y_fitted), np.max(y_real)) > 280 and max(np.max(y_fitted), np.max(y_real)) > 0):
        step = 20
    elif (10 + max(np.max(y_fitted), np.max(y_real)) > 480 and max(np.max(y_fitted), np.max(y_real)) > 0):
        step = 50

        ax1.set_yticks(range(0, int(10 + max(np.max(y_fitted), np.max(y_real))), step))
    ax1.set_xbound([-10, 490])
    ax1.set_ybound([0, 10 + max(np.max(y_fitted), np.max(y_real))])
    ax1.set_ylabel('Read Count')
    ax1.set_title('p-value: ' + str("%g" % p_val) + '\n basal(c): ' + str(
        "%g" % round(theta[0], 3)) + ' ind(h): ' + str(
        "%g" % round(theta[1], 3)) + '\n ' + 'deg(b): ' + str("%g" % round(theta[2], 3)) + ' t0: ' + str(
        "%g" % round(theta[3], 3)) + ' sig: ' + str("%g" % round(int(theta[4]), 3)))
    plt.savefig(folder_name + '/' + fd_name + str(id) + '.png', bbox_inches='tight')

    plt.close('all')


def plot_for_rep(y_fitted,  fd_name):
    # plot the fit results for figures
    plt.close("all")

    ax1 = plt.subplot(1, 1, 1)
    ax1.plot(np.linspace(0, len(y_fitted), len(y_fitted)), y_fitted, label='fit')
    ax1.axes.get_xaxis().set_visible(False)
    ax1.axes.get_yaxis().set_visible(False)
    plt.savefig(fd_name + str(id) + '.png', bbox_inches='tight')

    plt.close('all')



def F_test(real_y, y_fitted, to_print=False):
    '''
    F test for the two models :
    1. constant model (average)
    2. GS model
    :param real_y: real values from source file
    :param y_fitted: values from the GS model fit
    :param to_print: a boolean, if true print description
    of f test results
    :return: the p-value of the test
    '''
    # F test
    real_y = real_y + 1
    y_fitted = y_fitted + 1

    # H0 - const model
    df1 = 1
    data_mean = np.repeat(np.mean(real_y), len(real_y))
    model_1_RSS = calc_RSS(np.log(real_y), np.log(data_mean[:, None]))
    # H1 - differential equation
    model_2_RSS = calc_RSS(np.log(real_y), np.log(y_fitted))
    df2 = 5
    DFbet = df2 - df1

    DFwith = len(real_y) - DFbet - 1
    F = ((model_1_RSS - model_2_RSS) / DFbet) / (model_2_RSS / DFwith)
    '''
    the log on the data is since the F test assumes the data
    #  is normally distributed - but gene expression does not distribute normally.
    #  however log(gene_exp_data) does
    '''

    if (to_print):
        print("*******RSS(H0,H1)*******")
        print(model_1_RSS, model_2_RSS)
        print("************************")
        print("******F test value******")
        print(F)
        print("************************")
        print("********p-value*********")
        print(scipy.stats.f.ppf(0.01, DFwith, DFbet))
        print("*******numerator*******")
        print(((model_1_RSS - model_2_RSS) / DFbet))
        print("*******denumrator*******")
        print((model_2_RSS / DFwith))
        print(len(data_mean))
        print("*******RSS(H0,H1)*******")
        print(model_1_RSS, model_2_RSS)
        print("************************")


    print('F-critical: ' + str(scipy.stats.f.ppf(0.9, DFbet, DFwith)))
    print('p-value: ' + str(1 - scipy.stats.f.cdf(F, DFbet, DFwith)))
    print('F -statisti: ' + str(F))
    print('DFBet:' + str(DFbet))
    print('DFwith:' + str(DFwith))

    return 1 - scipy.stats.f.cdf(F, DFbet, DFwith)



def calc_models_fitted_y(real_y, file_id, rep_num,
                         theta_hat_GS, curr_gene, toPrint):
    '''
    recomputes the fitted values using the
    theta parameters and the GS model
    :param real_y: real gene expression values
    :param file_id: experiment id
    :param rep_num: number of repeats
    :param theta_hat_GS: best params using GS model fit
    :param curr_gene: the current gene name
    :param toPrint: if print results to stdout
    :return:
    '''
    # real data time line
    t_out = choose_time_line(file_id)[1]
    # expressive model getting the real fitted y
    y_fitted_GS = GS_model(theta_hat_GS, np.linspace(0, 483, 484, dtype='int'), np.mean(real_y[:rep_num]))
    # const model
    data_mean = np.repeat(np.mean(real_y), len(real_y))

    if (toPrint == True):
        print('#################')
        print(curr_gene)
        print('dif off peak : ' + str(np.max(y_fitted_GS) / np.max(real_y)))
        print('model-fit (RMSE): ' + str(calc_rmse(real_y, y_fitted_GS[t_out], np.max(y_fitted_GS))))
        print('const-fit (RMSE): ' + str(calc_rmse(real_y, data_mean, data_mean[0])))
        print('#################')

    return (new_real_y, y_fitted_GS, t_out, data_mean, calc_rmse(real_y, y_fitted_GS[t_out], np.max(y_fitted_GS)))


def add_mean_point(file_id, real_y, rep_num, x_labels):
    '''
     calculates the mean of each time frame
    :param file_id: unique id for each experiment
    :param real_y: data from source file
    :param rep_num: number of repeats per time frame
    :param x_labels: time frames
    :return: mean of repeats vec
    '''
    if file_id == c.FILE_ID_21_14:
        real_y = np.concatenate((real_y[:9], np.repeat(real_y[9], 3), real_y[10:]))
    if file_id == c.FILE_ID_90:
        real_y = np.concatenate(
            (real_y[:2], np.mean(real_y[:2])[None], real_y[3:9], np.repeat(real_y[9], 3), real_y[9:],
             np.mean(real_y[15:17])[None]))
    if file_id == c.FILE_ID_89:
        real_y = np.concatenate((
            real_y[:9], np.mean(real_y[9:11])[None], real_y[9:11], np.mean(real_y[11:13])[None], real_y[11:13],
            np.repeat(real_y[13:14], 3)))

    real_y_reshape = np.reshape(real_y, (int(len(real_y) / rep_num), rep_num))
    mean_of_rep = np.mean(real_y_reshape, axis=1)

    return mean_of_rep


def check_GS(thetha):
    '''
    for sanity check
    :param thetha: params fit
    '''
    x = np.linspace(0, 480, 481, dtype='int')
    y_fitted_GS = GS_model([0.2, 1, 0.2, 200, 5], np.linspace(0, 480, 481, dtype='int'), 0)
    plt.figure()
    plt.subplot(2, 1, 1)
    plt.plot(x, y_fitted_GS)
    plt.subplot(2, 1, 2)
    plt.plot(x, GS_model(thetha, np.linspace(0, 480, 481, dtype='int'), 0))
    plt.show()
    # plt.savefig('test.jpg')


def plotting_data_by_names(GS_path, names_path, data_path, file_id, folder_name):
    '''
    plots genes by name
    # plotting_data_by_names("C:/Users/Lior/Desktop/GS_model_fit/RES/222/all_res_222.csv",
    #                        "C:/Users/Lior/Documents/high_ind.csv",
    #                        "C:/Users/Lior/Desktop/GS_model_fit/data_222_esat_t_normalized.csv",222,
    #                        'C:/Users/Lior/Desktop/GS_model_fit/High_ind')
    :param GS_path: path to GS model results file ("C:/Users/Lior/Desktop/GS_model_fit/RES/222/all_res_222.csv")
    :param names_path: csv with all gene names ("gene_names.csv")
    :param data_path: normalized gene data
    :param file_id: experiment id
    :param folder_name: the output folder
    :return:
    '''
    ## reader part
    reader = csv.reader(open(names_path, "r"), delimiter=",")
    # reader = pd.read_csv(file_path,encoding="ISO-8859-1")
    x = list(reader)
    gene_names_ind = np.array([elem[1] for elem in x])[1:].astype("str")

    reader = csv.reader(open(GS_path, "r"), delimiter=",")
    x = list(reader)
    gene_names_res_file = list(np.array([elem[0] for elem in x])[1:].astype("str"))
    GS_fitted_params = np.array([elem[1:7] for elem in x]).astype("str")
    p_values = GS_fitted_params[:, 0]
    GS_fitted_params = GS_fitted_params[:, 1:6]
    # reader = pd.read_csv(file_path,encoding="ISO-8859-1")
    exp_mat = file_reader.get_gene_expression_matrix(data_path)
    gene_names_exp_mat = list(file_reader.get_gene_names(data_path))

    for i in range(0, len(gene_names_ind)):
        gene_GS_idx = gene_names_res_file.index(gene_names_ind[i])
        gene_exp_mat_idx = gene_names_exp_mat.index(gene_names_ind[i])
        real_y = exp_mat[gene_exp_mat_idx]
        theta_hat_GS = np.array(GS_fitted_params[gene_GS_idx + 1,], dtype='float')
        reduced_y, y_fitted_GS, t_out, data_mean, _ = calc_models_fitted_y(real_y, file_id, 3,
                                                                           theta_hat_GS, gene_names_ind[i], True)

        plot_only_fit(real_y, y_fitted_GS, data_mean, 3, np.float(p_values[gene_GS_idx + 1]), theta_hat_GS,
                      gene_names_ind[i], file_id, folder_name)

    return GS_fitted_params, gene_names_ind


def common_F_test_orig(real_y_1,real_y_2, y_1_fitted_common,y_2_fitted_common,y_1_fitted_old,y_2_fitted_old, data_mean_1,data_mean_2 ,to_print=False):

    real_y_1 = real_y_1 + 1
    real_y_2 = real_y_2 + 1

    y_1_fitted_common = y_1_fitted_common + 1
    y_2_fitted_common = y_2_fitted_common + 1

    y_1_fitted_old = y_1_fitted_old +1
    y_2_fitted_old = y_2_fitted_old +1

    # H0 - nine parameters model - one shared parameter
    df1 = 9
    model_1_RSS_1 = calc_RSS(np.log(real_y_1), np.log(y_1_fitted_common))
    model_1_RSS_2 = calc_RSS(np.log(real_y_2), np.log(y_2_fitted_common))
    model_1_RSS = model_1_RSS_1 + model_1_RSS_2
    # H1 - ten parameters model
    model_2_RSS_1 = calc_RSS(np.log(real_y_1), np.log(y_1_fitted_old))
    model_2_RSS_2 = calc_RSS(np.log(real_y_2), np.log(y_2_fitted_old))
    model_2_RSS = model_2_RSS_1 + model_2_RSS_2

    df2 = 10
    DFbet = df2 - df1  # number of free elemets  in first model 5 and const model 1
    DFwith = 2*len(real_y_1) - df2
    F = ((model_1_RSS - model_2_RSS) / DFbet) / (model_2_RSS / DFwith)


    if (to_print):
        print("*******RSS(9,10)*******")
        print(model_1_RSS, model_2_RSS)
        print("************************")
        print("******F test value******")
        print(F)
        print("************************")
        print("********p-value*********")
        print(scipy.stats.f.ppf(0.01, DFwith, DFbet))
        print("*******RMSE 9 params 1*******")
        print(calc_rmse(real_y_1,y_1_fitted_common))
        print("*******RMSE 9 params 2*******")
        print(calc_rmse(real_y_2, y_2_fitted_common))

        print("*******RMSE 10 params 1*******")
        print(calc_rmse(real_y_1, y_1_fitted_old, ))
        print("*******RMSE 9 params 2*******")
        print(calc_rmse(real_y_2, y_2_fitted_old, ))

        print("************************")
    if (np.isnan((1 - scipy.stats.f.cdf(F, DFbet, DFwith)))):
        print("yes")

    print('F-critical: ' + str(scipy.stats.f.ppf(0.9, DFbet, DFwith)))
    print('p-value: ' + str(1 - scipy.stats.f.cdf(F, DFbet, DFwith)))
    print('F -statisti: ' + str(F))
    print('DFBet:' + str(DFbet))
    print('DFwith:' + str(DFwith))

    return 1 - scipy.stats.f.cdf(F, DFbet, DFwith)

def common_F_test(real_y_1,real_y_2, y_1_fitted_common,y_2_fitted_common,y_1_fitted_old,y_2_fitted_old, data_mean_1,data_mean_2 ,to_print=False):

    real_y_1 = real_y_1 + 1
    real_y_2 = real_y_2 + 1

    y_1_fitted_common = y_1_fitted_common + 1
    y_2_fitted_common = y_2_fitted_common + 1

    y_1_fitted_old = y_1_fitted_old +1
    y_2_fitted_old = y_2_fitted_old +1

    # H0 - nine parameters model - one shared parameter
    df1 = 9
    model_1_RSS_1 = calc_RSS(np.log(real_y_1), np.log(y_1_fitted_common))
    model_1_RSS_2 = calc_RSS(np.log(real_y_2), np.log(y_2_fitted_common))
    model_1_RSS = model_1_RSS_1 + model_1_RSS_2
    # H1 - ten parameters model
    model_2_RSS_1 = calc_RSS(np.log(real_y_1), np.log(y_1_fitted_old))
    model_2_RSS_2 = calc_RSS(np.log(real_y_2), np.log(y_2_fitted_old))
    model_2_RSS = model_2_RSS_1 + model_2_RSS_2

    df2 = 10
    DFbet = df2 - df1  # number of free elemets  in first model 5 and const model 1
    DFwith = len(real_y_1)+ len(real_y_2) - df2
    F = ((model_1_RSS - model_2_RSS) / DFbet) / (model_2_RSS / DFwith)


    if (to_print):
        print("*******RSS(9,10)*******")
        print(model_1_RSS, model_2_RSS)
        print("************************")
        print("******F test value******")
        print(F)
        print("************************")
        print("********p-value*********")
        print(scipy.stats.f.ppf(0.01, DFwith, DFbet))
        print("*******RMSE 9 params 1*******")
        print(calc_rmse(real_y_1,y_1_fitted_common))
        print("*******RMSE 9 params 2*******")
        print(calc_rmse(real_y_2, y_2_fitted_common))

        print("*******RMSE 10 params 1*******")
        print(calc_rmse(real_y_1, y_1_fitted_old, ))
        print("*******RMSE 9 params 2*******")
        print(calc_rmse(real_y_2, y_2_fitted_old, ))

        print("************************")
    if (np.isnan((1 - scipy.stats.f.cdf(F, DFbet, DFwith)))):
        print("yes")

    print('F-critical: ' + str(scipy.stats.f.ppf(0.9, DFbet, DFwith)))
    print('p-value: ' + str(1 - scipy.stats.f.cdf(F, DFbet, DFwith)))
    print('F -statisti: ' + str(F))
    print('DFBet:' + str(DFbet))
    print('DFwith:' + str(DFwith))

    return 1 - scipy.stats.f.cdf(F, DFbet, DFwith)


