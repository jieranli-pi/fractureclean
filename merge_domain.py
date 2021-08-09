"""
Function that merges the domain partitions
"""

import numpy as np

##### IMPORTANT COMMENT #####
# The program is done and works, but must still discuss 4 things with Joey:
#       1. Estimate for the length of newly introduced arrays, which is estimated number of merged fractions.
#       2. How do we find out if fractions were cut at partition lines? Two options:
#               > We could use a threshold, but I am not sure how accurate this is.
#               > We could use ==, but you always have to be careful with that for floats. Right now tho, because
#                 these values are literally copied from eachother, it might be possible.
#       3. How are we going to build in a security in case that the threshold is not accurate enough and it
#          marks a fraction that actually was not partitioned.
#       4. I just used act_frac_sys as input now instead of act_frac_sys_list, so must change this.


def merge_domain(act_frac_sys_list, frac_set_vec_list, partition_lines, number_partitions_x):
    """

    :param act_frac_sys_list: Array of fractures that were once partitioned, and must now be merged back together.
    :param frac_set_vec_list: Array of grouping identifiers.
    :param partition_lines: Lines over which the domain was partitioned at the start of parallelizing.
    :param number_partitions_x: The numbers of partitions in x-direction.
    :return:
    """

    act_frac_sys = np.vstack(act_frac_sys_list)
    frac_set_vec = np.hstack(frac_set_vec_list)

    for i in range(len(partition_lines)):
        if i < number_partitions_x - 1:
            direction = 0
            other_dir = 1
        else:
            direction = 1
            other_dir = 0

        partition_line = partition_lines[i]

        # Each line will have one constant value, either on the x- or y-axis.
        const_value_par_line = partition_line[direction]

        merge_tolerance = 1e-10
        inds_left = np.where(np.abs(act_frac_sys[:, direction + 2] - const_value_par_line) <= merge_tolerance)
        inds_right = np.where(np.abs(act_frac_sys[:, direction] - const_value_par_line) <= merge_tolerance)

        #inds_left = np.where(act_frac_sys[:, direction + 2] == const_value_par_line)
        #inds_right = np.where(act_frac_sys[:, direction] == const_value_par_line)

        inds_afs = np.hstack((inds_left, inds_right))

        segments_left = act_frac_sys[inds_left]
        segments_right = act_frac_sys[inds_right]

        # Check if both arrays have the same length:
        if len(segments_right[:, 0]) != len(segments_left[:, 0]):
            raise Exception("Either a fraction was found that is not supposed to be merged \n "
                            "or a fraction was not found that was supposed to be merged.")

        # Sort the indices such that the segments can then be put in the right order.
        right_sort_inds = np.lexsort((segments_right[:, other_dir + 2], segments_right[:, other_dir]))
        left_sort_inds = np.lexsort((segments_left[:, other_dir], segments_left[:, other_dir + 2]))
        segments_right = segments_right[right_sort_inds]
        segments_left = segments_left[left_sort_inds]
        # Delete intersection information, it is not needed and stands in the way of merging segments.
        segments_right = segments_right[:, 2:]
        segments_left = segments_left[:, :2]

        # Grouping the left and right part or bottom and top part of the newly merged fracture together.
        restored_fracs = np.hstack((segments_left, segments_right))

        # Deleting the previous fractures that are now merged and adding the fractures to act_frac_sys.
        act_frac_sys = np.delete(act_frac_sys, inds_afs, axis=0)
        act_frac_sys = np.vstack((act_frac_sys, restored_fracs))

        # Update frac_set_vec
        frac_set_vec_new = frac_set_vec[inds_right][right_sort_inds]
        frac_set_vec = np.delete(frac_set_vec, inds_afs, axis=0)
        frac_set_vec = np.hstack((frac_set_vec, frac_set_vec_new))

    return act_frac_sys, frac_set_vec
