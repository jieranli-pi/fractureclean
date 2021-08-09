"""
Module that discretizes a fracture in to small sub fractures(nodes)
Rewritten from a matlab script by Stephan de Hoop

Author: Mieke Daemen/Jieran Li
Date: 22/11/2020
"""

import math as m
import numpy as np


def discr_segm(x_coor, y_coor, lencal):
    """
    This function takes a the x and y -coordinates of a fractures and returns a matrix with sub fractures.
    :param x_coor: an array with two elements that contains the x coordinates of the fracture.
                        the first element, x_coor[0], is the x-coordinate of one vertex
                        the second element, x_coor[1], is the x-coordinate of the other vertex
    :param y-coor: an array with two elements that contains the y coordinates of the fracture.
                        the first element, x_coor[0], is the x-coordinate of one vertex
                        the second element, x_coor[1], is the x-coordinate of the other vertex
    :param lencal: the length of the sub segments.
    :return: seg_mat: a matrix with the coordinates of the individual sub segments.

    """
    # length of a fracture
    len_frac = m.sqrt((x_coor[1]-x_coor[0])**2 + (y_coor[1]-y_coor[0])**2)
    # number of fractures in one fracture, round and make at least equal to 1
    num_segm = int(round(len_frac / lencal))
    if num_segm == 0:
        num_segm = 1

    seg_mat = np.zeros((num_segm, 4))
    for i in range(num_segm):
        seg_mat[i,0] = x_coor[0]+(i/num_segm)*(x_coor[1]-x_coor[0])
        seg_mat[i,1] = y_coor[0]+(i/num_segm)*(y_coor[1]-y_coor[0])
        seg_mat[i,2] = x_coor[0]+((i+1)/num_segm)*(x_coor[1]-x_coor[0])
        seg_mat[i,3] = y_coor[0]+((i+1)/num_segm)*(y_coor[1]-y_coor[0])
    return seg_mat
