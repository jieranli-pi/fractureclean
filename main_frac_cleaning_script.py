"""
This module contains the main script for cleaning a raw fracture network

Author: Joey Herbold
"""
# TODO Cleanup and documentation!
# import modules ===================================================================================
import os
import math as m
import numpy as np
import matplotlib.pyplot as plt
import FractureCleaning
import matlab

FC = FractureCleaning.initialize()  # initialization of imported MATLAB functions
# It is advised to rewrite these functions in python but in order to speed up the project at this
# moment it is decided to use the functions directly from

from cleaning_python.calc_intersections_segm import calc_intersections_segm
from cleaning_python.closest_point_method import closest_point_method
# ==================================================================================================

# Optional Processing steps ========================================================================
apply_rotation = False
remove_duplicates = True
calc_init_intersections = False
straighten_fracs_pre = True
straighten_fracs_post = True
# ==================================================================================================

# Optional visualizations ==========================================================================
plot_frac_network_pre = True
plot_frac_network_post = True
plot_statistics = True
# ==================================================================================================

# Tolerances =======================================================================================
tolerance_frac = 0.5
tolerance_midpt = 0.45
tolerance_rad = 0.5
tolerance_zero = 1e-4
tolerance_intersect = 1e-10
tolerance_angle = 7.5
# ==================================================================================================

# Accuracy
num_decimals = 5
norm_order = 100

# Characteristic length scales
ratio_char_scales = [0.5, 1, 2, 4, 8]
char_len = 15
vec_char_len = [char_len * scale for scale in ratio_char_scales]

# Input/output file path ===========================================================================
BASE_DIR = 'tln_7_ef_high_filtered/'
IN_DIR = '../input/' + BASE_DIR
ith_real = 1  # File input

OUT_DIR = '../output/' + BASE_DIR
if not os.path.exists(OUT_DIR):
    os.mkdir(OUT_DIR)

for scale in vec_char_len:
    OUT_DIR_SCALE = OUT_DIR + 'char_scale_' + str(scale)
    if not os.path.exists(OUT_DIR_SCALE):
        os.mkdir(OUT_DIR_SCALE)
# ==================================================================================================

# Load data
frac_data_full = np.loadtxt(IN_DIR + 'real_' + str(ith_real) + '.txt')
num_main_segm = frac_data_full.shape[0]
if frac_data_full.shape[1] == 4:
    frac_data = np.zeros((num_main_segm, 6))
    frac_data[:, 2:] = np.round(frac_data_full, num_decimals)
else:
    frac_data = np.round(frac_data_full, num_decimals)
old_frac_data = frac_data
act_frac_sys = frac_data[:, 2:]

if apply_rotation:
    act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
    act_frac_sys_rot, rot_mat, bin_count, opt_angle = \
        FC.calc_angle_frac_sys(act_frac_sys, matlab.double([6], size=(1, 1)), nargout=4)
    act_frac_sys = np.array(act_frac_sys)
    act_frac_sys_rot = np.array(act_frac_sys_rot)
    if plot_frac_network_pre:
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.plot(act_frac_sys[:, [0, 2]].transpose(),
                 act_frac_sys[:, [1, 3]].transpose(), linewidth=1, color='black')
        plt.title('Original fracture network')

        plt.subplot(1, 2, 2)
        plt.plot(act_frac_sys_rot[:, [0, 2]].transpose(),
                 act_frac_sys_rot[:, [1, 3]].transpose(), linewidth=1, color='black')
        plt.title(f'Fractures after rotation of {opt_angle} degrees')
        plot_frac_network_pre = False

if plot_frac_network_pre:
    plt.figure()
    plt.plot(act_frac_sys[:, [0, 2]].transpose(),
             act_frac_sys[:, [1, 3]].transpose(), linewidth=1, color='black')
    plt.title('Original fracture network')

frac_set_vec = np.zeros((act_frac_sys.shape[0], 1))

if remove_duplicates:
    act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
    frac_set_vec = matlab.double(frac_set_vec.tolist(), size=frac_set_vec.shape)
    act_frac_sys, frac_set_vec = FC.extract_unique_segm(act_frac_sys, frac_set_vec,
                                                        matlab.double([tolerance_zero], size=(1, 1)),
                                                        nargout=2)
    act_frac_sys = np.round(np.array(act_frac_sys), num_decimals)
    frac_set_vec = np.array(frac_set_vec)

    num_main_segm = act_frac_sys.shape[0]
    frac_data = np.zeros((num_main_segm, 6))
    frac_data[:, 2:] = act_frac_sys

len_main_segm = np.sqrt((act_frac_sys[:, 0] - act_frac_sys[:, 2]) * (act_frac_sys[:, 0] - act_frac_sys[:, 2])
                        + (act_frac_sys[:, 1] - act_frac_sys[:, 3]) * (act_frac_sys[:, 1] - act_frac_sys[:, 3]))
sorted_len_main_segm = np.sort(len_main_segm, axis=0)[::-1]

if plot_statistics:
    plt.figure()
    plt.subplot(1, 2, 1)
    plt.hist(len_main_segm, bins=25, rwidth=0.8, edgecolor='black', orientation='horizontal')
    plt.ylim([0, 1.1 * np.max(len_main_segm)])
    plt.xlabel('Frequency')
    plt.ylabel('Length segments')
    plt.title('Histogram')

    plt.subplot(1, 2, 2)
    plt.plot(np.sort(len_main_segm, axis=0), linewidth=2, color='blue')
    plt.hlines(char_len, 1, len(len_main_segm), linewidth=2, color='red', label='l$_c$')
    plt.ylim([0, max(len_main_segm) * 1.1])
    plt.ylabel('Length segments')
    plt.title('Cumulative distr.')
    plt.legend(loc='upper left')

if straighten_fracs_pre:
    act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
    frac_set_vec = matlab.double(frac_set_vec.tolist(), size=frac_set_vec.shape)
    act_frac_sys, frac_set_vec = FC.straighten_frac_segm(act_frac_sys, frac_set_vec,
                                                         matlab.double([tolerance_zero], size=(1, 1)),
                                                         matlab.double([tolerance_angle], size=(1, 1)),
                                                         nargout=2)
    act_frac_sys = np.round(np.array(act_frac_sys), num_decimals)
    frac_set_vec = np.array(frac_set_vec)

if calc_init_intersections:
    act_frac_sys, frac_set_vec = calc_intersections_segm(act_frac_sys, frac_set_vec,
                                                         tolerance_intersect, tolerance_zero)
    act_frac_sys = np.round(act_frac_sys, num_decimals)

num_inner_iter = 2

for ith_iter, char_len_temp in enumerate(vec_char_len):
    for ith_inner_iter in range(num_inner_iter):
        act_frac_sys, frac_set_vec = closest_point_method(act_frac_sys, frac_set_vec,
                                                          char_len, norm_order,
                                                          tolerance_rad, tolerance_zero)
        act_frac_sys = np.round(act_frac_sys, num_decimals)
        if straighten_fracs_post:
            act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
            frac_set_vec = matlab.double(frac_set_vec.tolist(), size=frac_set_vec.shape)
            act_frac_sys, frac_set_vec = FC.straighten_frac_segm(act_frac_sys, frac_set_vec,
                                                                 matlab.double([tolerance_zero], size=(1, 1)),
                                                                 matlab.double([tolerance_angle], size=(1, 1)),
                                                                 nargout=2)
            act_frac_sys = np.round(np.array(act_frac_sys), num_decimals)
            frac_set_vec = np.array(frac_set_vec)

        act_frac_sys, frac_set_vec = calc_intersections_segm(act_frac_sys, frac_set_vec,
                                                             tolerance_intersect, tolerance_zero)
        act_frac_sys = np.round(act_frac_sys, num_decimals)

    act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
    act_frac_sys = FC.adjust_vertical_segm(act_frac_sys, matlab.double([tolerance_zero], size=(1, 1)))
    act_frac_sys = np.round(np.array(act_frac_sys), num_decimals)

    act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
    frac_set_vec = matlab.double(frac_set_vec.tolist(), size=frac_set_vec.shape)
    act_frac_sys, frac_set_vec = FC.find_partial_overlap_and_small_angles(act_frac_sys, frac_set_vec,
                                                         matlab.double([tolerance_zero], size=(1, 1)),
                                                         matlab.double([char_len_temp], size=(1, 1)),
                                                         nargout=2)
    act_frac_sys = np.round(np.array(act_frac_sys), num_decimals)
    frac_set_vec = np.array(frac_set_vec)

    act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
    frac_set_vec = matlab.double(frac_set_vec.tolist(), size=frac_set_vec.shape)
    act_frac_sys, frac_set_vec = FC.extract_unique_segm(act_frac_sys, frac_set_vec,
                                                        matlab.double([tolerance_zero], size=(1, 1)),
                                                        nargout=2)
    act_frac_sys = np.round(np.array(act_frac_sys), num_decimals)
    frac_set_vec = np.array(frac_set_vec)

    act_frac_sys, frac_set_vec = calc_intersections_segm(act_frac_sys, frac_set_vec,
                                                         tolerance_intersect, tolerance_zero)
    act_frac_sys = np.round(act_frac_sys, num_decimals)

    act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
    frac_set_vec = matlab.double(frac_set_vec.tolist(), size=frac_set_vec.shape)
    act_frac_sys, frac_set_vec = FC.extract_unique_segm(act_frac_sys, frac_set_vec,
                                                        matlab.double([tolerance_zero], size=(1, 1)),
                                                        nargout=2)
    act_frac_sys = np.round(np.array(act_frac_sys), num_decimals)
    frac_set_vec = np.array(frac_set_vec)

    act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
    frac_set_vec = matlab.double(frac_set_vec.tolist(), size=frac_set_vec.shape)
    act_frac_sys, frac_set_vec = FC.find_partial_overlap_and_small_angles(act_frac_sys, frac_set_vec,
                                                                          matlab.double([tolerance_zero], size=(1, 1)),
                                                                          matlab.double([char_len_temp], size=(1, 1)),
                                                                          nargout=2)
    act_frac_sys = np.round(np.array(act_frac_sys), num_decimals)
    frac_set_vec = np.array(frac_set_vec)

    act_frac_sys = matlab.double(act_frac_sys.tolist(), size=act_frac_sys.shape)
    frac_set_vec = matlab.double(frac_set_vec.tolist(), size=frac_set_vec.shape)
    act_frac_sys, frac_set_vec = FC.find_actual_overlap_segm(act_frac_sys, frac_set_vec,
                                                             matlab.double([tolerance_zero], size=(1, 1)),
                                                             nargout=2)
    act_frac_sys = np.round(np.array(act_frac_sys), num_decimals)
    frac_set_vec = np.array(frac_set_vec)

    if apply_rotation:
        act_frac_sys_rot = np.zeros_like(act_frac_sys)
        act_frac_sys_rot[:, [0, 1]] = act_frac_sys[:, [0, 1]] * rot_mat
        act_frac_sys_rot[:, [2, 3]] = act_frac_sys[:, [2, 3]] * rot_mat
        act_frac_sys = act_frac_sys_rot

    OUT_DIR_SCALE = OUT_DIR + 'char_scale_' + str(char_len_temp)
    np.savetxt(OUT_DIR_SCALE + '/real_' + str(ith_real) + '.txt', act_frac_sys, fmt='%8.5f', delimiter='\t')

if plot_frac_network_post:
    num_figures = len(vec_char_len) + 1
    num_columns = m.ceil(num_figures/2)

    plt.figure()
    plt.subplot(2, num_columns, 1)
    plt.plot(old_frac_data[:, [2, 4]].transpose(),
             old_frac_data[:, [3, 5]].transpose(), marker='.', markerfacecolor='red', color='black')
    plt.title('Raw network')

    for ith_iter, char_len_temp in enumerate(vec_char_len):
        char_len_str = str(char_len_temp)
        clean_frac_sys = np.loadtxt(OUT_DIR + 'char_scale_' + char_len_str + '/real_' + str(ith_real) + '.txt')
        plt.subplot(2, num_columns, ith_iter+1)
        plt.plot(clean_frac_sys[:, [0, 2]].transpose(),
                 clean_frac_sys[:, [1, 3]].transpose(), marker='.', markerfacecolor='red', color='black')
        plt.title('Clean network l$_c$= ' + char_len_str)

plt.show()
FC.terminate()
