"""
This module contains the closest point method

Author: Mieke Daemen
"""

import numpy as np
from cleaning_python.discr_segm import discr_segm


def closest_point_method(act_frac_sys, frac_set_vec, char_len, norm_ord, tolerance_rad, tolerance_zero):
    """
        This function takes a set of fractures an returns a set of merged fractures, using the ClosePoint method.
        :param act_frac_sys: A list where all rows describe a fracture in the system, and the four columns
                         describe the following: Column 0 = X value of the start of a fracture
                                                 Column 1 = Y value of the start of a fracture
                                                 Column 2 = X value of the end of a fracture
                                                 Column 3 = Y value of the end of a fracture
        :param frac_set_vec: This is a list where groups of fractures are defined that are close to each other.
        :param char_len: the length of the sub segments.
        :param norm_ord: the norm used to calculate distances
        :param tolerance_rad: This tolerance is used to check whether a segment is close to another segment
        :param tolerance_zero: This tolerance is used to check whether a segment is non-zero. In this case, the
                           length of the segment must be larger than this value.
        :return act_frac_sys: This returned matrix is similar to the input act_frac_sys, except that it now takes
                           The ClosePoints into account.
        :return frac_set_vec: This is a list where groups of fractures are defined that are close to each other.
        """

    # max number of segments
    num_tot_new_segm_reduced = len(act_frac_sys[:, 1])

    len_segm_new = np.sqrt((act_frac_sys[:, 2] - act_frac_sys[:, 0]) ** 2
                           + (act_frac_sys[:, 3] - act_frac_sys[:, 1]) ** 2)

    num_segm_per_frac = np.around(len_segm_new / char_len)

    tot_small_segm = int(np.sum(num_segm_per_frac) + np.sum(num_segm_per_frac == 0))

    # initialize matrices and variables
    new_frac_set_list = np.zeros(((tot_small_segm* 2), 1))
    new_frac_sys = np.ones(((tot_small_segm* 2), 4))
    unq_node_list = np.ones(((tot_small_segm * 2), 2)) * 9999999
    total_loop = 0
    total_nodes = 0

    # order in wich segments are discretized, and merged
    order_to_discr = np.argsort(len_segm_new)
    order_to_discr = order_to_discr[::-1]

    # the well coordinate method...still needs to be improved
    num_wells = 5
    well_coords = np.array([[100, 100], [900, 100], [500, 500], [100, 900], [900, 900]])
    method_ii = 1
    frac_length_threshold = 0.1

    well_frac_ids = np.zeros(num_wells)
    for well in range(num_wells):
        # find closest fracture to well coordinate
        if method_ii == 1:
            # method 1
            centroids = [np.mean(act_frac_sys[:, [0, 2]], axis=1), np.mean(act_frac_sys[:, [1, 3]], axis=1)]
            compare_set_2 = (
                    ((centroids[0] - well_coords[well, 0]) ** 2 + (centroids[1] - well_coords[well, 1]) ** 2) ** (
                    1 / 2))
            val_ii = np.sort(compare_set_2)
            id_ii = np.argsort(compare_set_2)
            segm_length_closest = len_segm_new[id_ii]
            id_ii_2 = np.where(segm_length_closest >= frac_length_threshold)
            id_ii_real = id_ii[id_ii_2]
            val_ii_real = val_ii[id_ii_2]
        else:
            # method 2
            compare_set_1 = ((act_frac_sys[:, 0] - well_coords[well:1]) ** 2 + (
                    act_frac_sys[:, 1] - well_coords[well, 1]) ** 2) ** (1 / 2)
            compare_set_2 = ((act_frac_sys[:, 2] - well_coords[well:1]) ** 2 + (
                    act_frac_sys[:, 3] - well_coords[well, 1]) ** 2) ** (1 / 2)
            val_1 = np.sort(compare_set_1)
            id_1 = np.argsort(compare_set_1)
            val_2 = np.sort(compare_set_2)
            id_2 = np.argsort(compare_set_2)
            segm_length_closest = len_segm_new[id_1]
            id_ii_1 = np.where(segm_length_closest >= frac_length_threshold)
            id_ii_real_1 = id_1[id_ii_1]
            val_ii_real_1 = val_1[id_ii_1]

            segm_length_closest = len_segm_new[id_2]
            id_ii_2 = np.where(segm_length_closest >= frac_length_threshold)
            id_ii_real_2 = id_2[id_ii_2]
            val_ii_real_2 = val_2[id_ii_2]

            if val_ii_real_1[0] < val_ii_real_2[0]:
                id_ii_real = id_ii_real_1
            else:
                id_ii_real = id_ii_real_2
        well_frac_ids[well] = id_ii_real[0]

    order_to_discr_new = np.zeros(num_tot_new_segm_reduced)
    # check if all the wells are in unique fractures
    well_frac_ids_unq = np.unique(well_frac_ids)
    num_unq_fracs = np.size(well_frac_ids_unq)
    indices = np.ones(num_tot_new_segm_reduced, dtype=bool)
    for ii in range(num_unq_fracs):
        indices[order_to_discr == well_frac_ids_unq[ii]] = False
    order_to_discr_new[0:num_unq_fracs] = well_frac_ids_unq
    order_to_discr_new[num_unq_fracs:] = order_to_discr[indices]
    order_to_discr = np.copy(order_to_discr_new)

    # Close Point method
    for i in range(num_tot_new_segm_reduced):
        # The ith fracture is discretized
        sorted_segm = int(order_to_discr[i])
        x_coor = [act_frac_sys[sorted_segm, 0], act_frac_sys[sorted_segm, 2]]
        y_coor = [act_frac_sys[sorted_segm, 1], act_frac_sys[sorted_segm, 3]]
        seg_mat = discr_segm(x_coor, y_coor, char_len)

        for j in range(np.size(seg_mat[:, 0])):
            # the j-th segment of the i-th fracture wether the 2 nodes are close to nodes in unique-list
            # if that's the case we ne
            xn_coor = [seg_mat[j, 0], seg_mat[j, 2]]
            yn_coor = [seg_mat[j, 1], seg_mat[j, 3]]

            dist_node_1 = (((xn_coor[0] - unq_node_list[0:total_nodes, 0]) ** norm_ord +
                            (yn_coor[0] - unq_node_list[0:total_nodes, 1]) ** norm_ord)) ** (1 / norm_ord)

            ids_1 = np.where(dist_node_1 < tolerance_rad * char_len)[
                0]

            if ids_1.size > 0:
                # Merge node with closest existing point in radius:
                ids_1_val = dist_node_1[ids_1]
                val, id_min = min(ids_1_val), np.argmin(ids_1_val)

                if val > tolerance_zero:
                    # check if node exists in left-column
                    dist_leftcol = (((act_frac_sys[:, 0] - seg_mat[j, 0]) ** norm_ord) + (
                            (act_frac_sys[:, 1] - seg_mat[j, 1]) ** norm_ord)) ** (1 / norm_ord)
                    ids_leftcol = np.where(dist_leftcol < tolerance_zero)[0]

                    if ids_leftcol.size > 0:
                        act_frac_sys[ids_leftcol, 0] = unq_node_list[ids_1[id_min], 0]
                        act_frac_sys[ids_leftcol, 1] = unq_node_list[ids_1[id_min], 1]
                    # check if node exists in right-column
                    dist_rightcol = ((act_frac_sys[:, 2] - seg_mat[j, 0]) ** norm_ord +
                                     (act_frac_sys[:, 3] - seg_mat[j, 1]) ** norm_ord) ** (1 / norm_ord)
                    ids_rightcol = np.where(dist_rightcol < tolerance_zero)[0]

                    # Do the same for discretized segments
                    if ids_rightcol.size > 0:
                        act_frac_sys[ids_rightcol, 2] = unq_node_list[ids_1[id_min], 0]
                        act_frac_sys[ids_rightcol, 3] = unq_node_list[ids_1[id_min], 1]

                    dist_lcol = ((seg_mat[:, 0] - seg_mat[j, 0]) ** norm_ord +
                                 (seg_mat[:, 1] - seg_mat[j, 1]) ** norm_ord) ** (1 / norm_ord)
                    ids_lcol = np.where(dist_lcol < tolerance_zero)[0]

                    if ids_lcol.size > 0:
                        seg_mat[ids_lcol, 0] = unq_node_list[ids_1[id_min], 0]
                        seg_mat[ids_lcol, 1] = unq_node_list[ids_1[id_min], 1]

                    dist_rcol = ((seg_mat[:, 2] - seg_mat[j, 0]) ** norm_ord +
                                 (seg_mat[:, 3] - seg_mat[j, 1]) ** norm_ord) ** (1 / norm_ord)

                    ids_rcol = np.where(dist_rcol < tolerance_zero)[0]

                    if ids_rcol.size > 0:
                        seg_mat[ids_rcol, 2] = unq_node_list[ids_1[id_min], 0]
                        seg_mat[ids_rcol, 3] = unq_node_list[ids_1[id_min], 1]

                seg_mat[j, 0:2] = unq_node_list[ids_1[id_min], :]
                new_frac_sys[total_loop, 0:2] = seg_mat[j, 0:2]
                new_frac_set_list[total_loop] = frac_set_vec[sorted_segm]

            else:
                # add new unique node
                unq_node_list[total_nodes, 0] = xn_coor[0]
                unq_node_list[total_nodes, 1] = yn_coor[0]
                new_frac_sys[total_loop, 0] = xn_coor[0]
                new_frac_sys[total_loop, 1] = yn_coor[0]
                new_frac_set_list[total_loop] = frac_set_vec[sorted_segm]
                total_nodes = total_nodes + 1

            dist_node_2 = (((xn_coor[1] - unq_node_list[0:total_nodes, 0]) ** norm_ord +
                            (yn_coor[1] - unq_node_list[0:total_nodes, 1]) ** norm_ord)) ** (1 / norm_ord)
            ids_2 = np.where(dist_node_2 < tolerance_rad * char_len)[0]

            if ids_2.size > 0:
                ids_2_val = dist_node_2[dist_node_2 < tolerance_rad * char_len]
                val, id_min = min(ids_2_val), np.argmin(ids_2_val)

                if val > tolerance_zero:
                    dist_leftcol = ((act_frac_sys[:, 0] - seg_mat[j, 2]) ** norm_ord +
                                    (act_frac_sys[:, 1] - seg_mat[j, 3]) ** norm_ord) ** (1 / norm_ord)
                    ids_leftcol = np.where(dist_leftcol < tolerance_zero)[0]

                    if ids_leftcol.size > 0:
                        act_frac_sys[ids_leftcol, 0] = unq_node_list[ids_2[id_min], 0]
                        act_frac_sys[ids_leftcol, 1] = unq_node_list[ids_2[id_min], 1]

                    dist_rightcol = ((act_frac_sys[:, 2] - seg_mat[j, 2]) ** norm_ord +
                                     (act_frac_sys[:, 3] - seg_mat[j, 3]) ** norm_ord) ** (1 / norm_ord)
                    ids_rightcol = np.where(dist_rightcol < tolerance_zero)[0]

                    if ids_rightcol.size > 0:
                        act_frac_sys[ids_rightcol, 2] = unq_node_list[ids_2[id_min], 0]
                        act_frac_sys[ids_rightcol, 3] = unq_node_list[ids_2[id_min], 1]

                    dist_lcol = ((seg_mat[:, 0] - seg_mat[j, 2]) ** norm_ord +
                                 (seg_mat[:, 1] - seg_mat[j, 3]) ** norm_ord) ** (1 / norm_ord)
                    ids_lcol = np.where(dist_lcol < tolerance_zero)[0]

                    if ids_lcol.size > 0:
                        seg_mat[ids_lcol, 0] = unq_node_list[ids_2[id_min], 0]
                        seg_mat[ids_lcol, 1] = unq_node_list[ids_2[id_min], 1]

                    dist_rcol = ((seg_mat[:, 2] - seg_mat[j, 2]) ** norm_ord +
                                 (seg_mat[:, 3] - seg_mat[j, 3]) ** norm_ord) ** (1 / norm_ord)
                    ids_rcol = np.where(dist_rcol < tolerance_zero)[0]

                    if ids_rcol.size > 0:
                        seg_mat[ids_rcol, 2] = unq_node_list[ids_2[id_min], 0]
                        seg_mat[ids_rcol, 3] = unq_node_list[ids_2[id_min], 1]

                seg_mat[j, 2:4] = unq_node_list[ids_2[id_min], :]
                new_frac_sys[total_loop, 2:4] = seg_mat[j, 2:4]
                new_frac_set_list[total_loop] = frac_set_vec[sorted_segm]

            else:

                unq_node_list[(total_nodes), 0] = xn_coor[1]
                unq_node_list[(total_nodes), 1] = yn_coor[1]
                new_frac_sys[total_loop, 2] = xn_coor[1]
                new_frac_sys[total_loop, 3] = yn_coor[1]
                new_frac_set_list[total_loop] = frac_set_vec[sorted_segm]
                total_nodes = total_nodes + 1
            total_loop = total_loop + 1

    new_frac_sys = new_frac_sys[~np.all(new_frac_sys == 1, axis=1)]

    return new_frac_sys, new_frac_set_list
