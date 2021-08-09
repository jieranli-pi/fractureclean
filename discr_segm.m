function [segm_mat] = discr_segm(x_coor, y_coor, char_len)

% Calculate length of big fracture-segment:
len_frac = sqrt( (x_coor(2) - x_coor(1)).^2 + ...
                 (y_coor(2) - y_coor(1)).^2 );

% Find number of smaller segments (round and make at least 1 segm!):
num_segm = round(len_frac/char_len);
if num_segm == 0
    % Very small segment, will get merged anyway!
    num_segm = 1;
end

% Recalculate length of smaller segments (based on floored #segm):
len_segm = len_frac / num_segm;

% Find nodes-list of each smaller segment:
% [x_1, y_1, x_2, y_2]
segm_mat = zeros(num_segm, 4);

for ith_dis_segm = 1:num_segm
    % Add each small segment to list:
    % First segment, fix coordinate to original coordinate:
    segm_mat(ith_dis_segm, 1:2) = [x_coor(1)+((ith_dis_segm-1)/num_segm)*(x_coor(2)-x_coor(1)), ...
                                   y_coor(1)+((ith_dis_segm-1)/num_segm)*(y_coor(2)-y_coor(1))];
    segm_mat(ith_dis_segm, 3:4) = [x_coor(1)+(ith_dis_segm/num_segm)*(x_coor(2)-x_coor(1)), ...
                                   y_coor(1)+(ith_dis_segm/num_segm)*(y_coor(2)-y_coor(1))];
end

if num_segm == 0
    % Store normal segment:
    segm_mat = [x_coor(1), y_coor(1), x_coor(2), y_coor(2)];
end