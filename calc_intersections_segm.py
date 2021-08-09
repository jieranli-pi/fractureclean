"""
This file describes the function calc_intersections_segm, which is called in the main cleaning function. This
function finds the intersections of the fractures in act_frac_sys, and splits these fractures at the intersection
points. The list act_frac_sys is then returned with the updated fractures.

Author: Ole de Koning

Last updated: 12/12/2020 by Ole de Koning
"""

import numpy as np
from cleaning_python.find_parametric_intersect import find_parametric_intersect


def calc_intersections_segm(act_frac_sys, frac_set_vec, tolerance_intersect, tolerance_zero):
    """
    :param act_frac_sys: A list where all rows describe a fracture in the system, and the four columns
                         describe the following: Column 0 = X value of the start of a fracture
                                                 Column 1 = Y value of the start of a fracture
                                                 Column 2 = X value of the end of a fracture
                                                 Column 3 = Y value of the end of a fracture

    :param frac_set_vec: This is a list where groups of fractures are defined that are close to each other.

    :param tolerance_intersect: A value that was set to determine whether two fractures intersect of not.
                                If t and s, the parametric distances along fracture 1 and 2 respectively are
                                smaller than this tolerance, it is counted as an intersection.

    :param tolerance_zero: This tolerance is used to check whether a segment is non-zero. In this case, the
                           length of the segment must be larger than this value.

    :return: act_frac_sys: This returned matrix is similar to the input act_frac_sys, except that it now takes
                           all intersection into account.

    :return: frac_set_vec: This returned matrix is similar to the input frac_set_vec, except that it now takes
                           all intersection into account. The groups formed in this list will thus differ from the
                           ones in the input.
    """

    n_fracs = np.size(act_frac_sys, 0)

    # TODO Segment list global (just a list of zeros, which can be changed to non-zero values later):
    #segment_id = frac_set_vec
    #new_segment_id = np.zeros(max_length_new_segm)

    # Counter for how many fractures ii have intersections.
    frac_int_counter = -1

    # Defining an array which will later hold the following information about the intersection:
    # [ x_int, y_int, ii_frac, jj_frac ]
    intersect_info = np.zeros((1,4))

    for ii in range(0, n_fracs):

        # Obtaining the x and y coords of the start and end of the frac
        ii_frac = act_frac_sys[ii - frac_int_counter - 1 , :]

        # Initialize the array where all new intersections will be added to
        current_ii_intersect = np.zeros((10,4))
        ii_int_counter = 0

        # Get intersection points for this frac jj previously calculated intersections
        # prev_inds = the row numbers of intersections calculated previously where jj was the current ii.
        prev_inds = np.where(intersect_info[:, 3] == ii)

        # prev_ii_intersect = x and y value of previously calculated intersections
        prev_ii_intersect = intersect_info[prev_inds, :2][0]

        for jj in range(ii+1, n_fracs):

            jj_frac = act_frac_sys[jj - frac_int_counter - 1, :]

            t, s, int_coord = find_parametric_intersect(ii_frac, jj_frac) # This function will be added later by Joey

            if (t >= (0 - tolerance_intersect) and t <= (1 + tolerance_intersect)) and \
               (s >= (0 - tolerance_intersect) and s <= (1 + tolerance_intersect)):

                # Only store intersections of segments that don't already share a node:
                if not (np.linalg.norm(ii_frac[:2] - jj_frac[:2]) < tolerance_intersect or
                        np.linalg.norm(ii_frac[:2] - jj_frac[2:]) < tolerance_intersect or
                        np.linalg.norm(ii_frac[2:] - jj_frac[:2]) < tolerance_intersect or
                        np.linalg.norm(ii_frac[2:] - jj_frac[2:]) < tolerance_intersect):

                    # Add the proper information to the correct row in the intersection information array
                    try:
                        current_ii_intersect[ii_int_counter, :] = np.array([int_coord[0], int_coord[1], ii, jj])
                        ii_int_counter += 1
                    except:
                        current_ii_intersect = np.vstack((current_ii_intersect, [int_coord[0], int_coord[1], ii, jj]))

        # Check if there even are fractures. If not, we go on to the next fracture ii.
        if np.size(prev_inds[0]) == 0 and ii_int_counter == 0:
            continue

        # Counter for how many fractures actually contain intersections
        frac_int_counter += 1

        # Get rid of the first row, because these were just zeros.
        current_ii_intersect = current_ii_intersect[:ii_int_counter]

        # Add this information to the general intersection information array for following fractures
        intersect_info = np.vstack((intersect_info, current_ii_intersect))

        # Clear memory of used intersection, as it will not be needed anymore
        intersect_info = np.delete(intersect_info, prev_inds, axis=0)

        # Here, the fracture indices (ii and jj) are deleted from current_ii_intersect and the x and y values
        # from previously found intersections prev_ii_intersect are added to this array.
        if ii != 0:
            current_ii_intersect = np.vstack((current_ii_intersect[:, :2], prev_ii_intersect))
        elif ii == 0:
            current_ii_intersect = current_ii_intersect[:, :2]

        # The x and y values of the intersections are sorted from low to high value of x, and the start and end
        # value of the original fracture are added at the top and bottom of this array, respectively.
        current_ii_intersect = np.vstack((np.vstack((ii_frac[:2], current_ii_intersect, ii_frac[2:]))))

        # Sorting
        current_ii_intersect = current_ii_intersect[np.lexsort((current_ii_intersect[:,1],current_ii_intersect[:,0]))]


        num_new_fracs = len(current_ii_intersect[:, 0]) - 1

        # Now, we have an array with all x and y value of the points where the fracture ii should be split at.
        new_fract_sys = np.zeros((num_new_fracs, 4))
        for mm in range(0, num_new_fracs):
            new_fract_sys[mm, :] = [current_ii_intersect[mm, 0], current_ii_intersect[mm, 1],
                                    current_ii_intersect[mm + 1, 0], current_ii_intersect[mm + 1, 1]]

        # The frac_int_counter counts how many rows have already been deleted. This matters for the index of the row
        # that we must delete now.
        act_frac_sys = np.vstack((np.delete(act_frac_sys, ii - frac_int_counter, axis=0), new_fract_sys))

        frac_set_vec = np.hstack((np.delete(frac_set_vec, ii - frac_int_counter),
                                  np.ones(num_new_fracs) * frac_set_vec[ii - frac_int_counter]))

    # Extract only actual segments
    #new_segment_id = new_segment_id[0:glob_segm_count-1] ?

    # Determine length of new "main" segments:
    len_segm_new = np.sqrt((act_frac_sys[:, 0] - act_frac_sys[:, 2])*(act_frac_sys[:, 0] - act_frac_sys[:, 2]) +
                           (act_frac_sys[:, 1] - act_frac_sys[:, 3])*(act_frac_sys[:, 1] - act_frac_sys[:, 3]))

    # Remove non-zero fracture segments: This comment does not make sense, probably meant the opposite
    nonzero_segm = np.where(len_segm_new > tolerance_zero)[0]
    act_frac_sys = act_frac_sys[nonzero_segm, :]

    return act_frac_sys, frac_set_vec