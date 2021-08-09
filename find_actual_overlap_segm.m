function [act_frac_sys, frac_set_vec] = ...
    find_actual_overlap_segm(act_frac_sys, frac_set_vec, tolerance_zero)

% Script to remove segments that lie on top of each other but are different
% length:
norm_order = 2;
dist_crit_midpt = 0.35;

% Make sure segments are sorted from X_min to X_max:
act_frac_sys_dummy = act_frac_sys;
indices = act_frac_sys_dummy(:, 1) > act_frac_sys_dummy(:, 3);
act_frac_sys(indices, 1:2) = act_frac_sys_dummy(indices, 3:4);
act_frac_sys(indices, 3:4) = act_frac_sys_dummy(indices, 1:2);

% Extract unique segments:
act_frac_sys_dummy = act_frac_sys;
[act_frac_sys, id_frac_set, ~] = unique(act_frac_sys_dummy, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);

% Make sure collapsed segments are not present:
len_segm_new = sqrt( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                     (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 );
                 
nonzero_segm = len_segm_new>tolerance_zero;
act_frac_sys = act_frac_sys(nonzero_segm, :); 
frac_set_vec = frac_set_vec(nonzero_segm);

num_main_segm = size(act_frac_sys, 1);
removed_segm = zeros(num_main_segm, 1);
split_segment = zeros(num_main_segm, 1);
split_coord = zeros(num_main_segm, 2);
ith_split = 0;

% For later purposes:
len_segm_new = sqrt( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                     (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 );

% 1) Main loop over all segments:
for ith_segm = 1:num_main_segm
    if ~any(ith_segm == removed_segm)
    % Check every other segments that has either the same beginning 
    % or the same end node!
    same_left_node_x = find(abs(act_frac_sys(ith_segm, 1) - ...
                              act_frac_sys(:, 1)) < tolerance_zero);
    same_left_node_y = find(abs(act_frac_sys(ith_segm, 2) - ...
                              act_frac_sys(:, 2)) < tolerance_zero);
    same_left_node = intersect(same_left_node_x, same_left_node_y);
    
    same_right_node_x = find(abs(act_frac_sys(ith_segm, 3) - ...
                              act_frac_sys(:, 3)) < tolerance_zero);
    same_right_node_y = find(abs(act_frac_sys(ith_segm, 4) - ...
                              act_frac_sys(:, 4)) < tolerance_zero);
    same_right_node = intersect(same_right_node_x, same_right_node_y);
                           
    same_nodes_segm = union(same_left_node, same_right_node);
    same_nodes_segm(same_nodes_segm == ith_segm) = [];  % don't take segment itself
    
    % For the segments that have same node, calculate the angle:
    % If angle is exactly zero, then segment is either duplicate or
    if size(same_nodes_segm, 1) > 1
        same_nodes_segm = same_nodes_segm';
    end
    for ith_other_elem = same_nodes_segm
        % Find which nodes are shared:
        ithsegm_id_shared = [];
        ithsegm_id_other = [];
        other_id_shared = [];
        other_id_other = [];
        
        if abs(act_frac_sys(ith_segm, 1) - ...
               act_frac_sys(ith_other_elem, 1)) < tolerance_zero
           % Point (x2, y2) for both segments is non-shared is start
           ithsegm_id_shared = [1, 2];
           ithsegm_id_other = [3, 4];
           
           other_id_shared = [1, 2];
           other_id_other = [3, 4];

        elseif abs(act_frac_sys(ith_segm, 1) - ...
               act_frac_sys(ith_other_elem, 3)) < tolerance_zero 
           % Point (x2, y2) is non-shared for first segment, and 
           % (x1, y1) for other
           ithsegm_id_shared = [1, 2];
           ithsegm_id_other = [3, 4];
           
           other_id_shared = [3, 4];
           other_id_other = [1, 2];
           
        elseif abs(act_frac_sys(ith_segm, 3) - ...
               act_frac_sys(ith_other_elem, 1)) < tolerance_zero 
           % Point (x1, y1) is non-shared for first segment, and 
           % (x2, y2) for other
           ithsegm_id_shared = [3, 4];
           ithsegm_id_other = [1, 2];
           
           other_id_shared = [1, 2];
           other_id_other = [3, 4];
           
        elseif abs(act_frac_sys(ith_segm, 3) - ...
               act_frac_sys(ith_other_elem, 3)) < tolerance_zero 
           % Point (x1, y1) is non-shared for first segment, and 
           % (x2, y2) for other
           ithsegm_id_shared = [3, 4];
           ithsegm_id_other = [1, 2];
           
           other_id_shared = [3, 4];
           other_id_other = [1, 2];
        end
        
        vec_ith_segm = [act_frac_sys(ith_segm, ithsegm_id_other(1)) - act_frac_sys(ith_segm, ithsegm_id_shared(1));
                        act_frac_sys(ith_segm, ithsegm_id_other(2)) - act_frac_sys(ith_segm, ithsegm_id_shared(2))];
        vec_other_segm = [act_frac_sys(ith_other_elem, other_id_other(1)) - act_frac_sys(ith_other_elem, other_id_shared(1));
                          act_frac_sys(ith_other_elem, other_id_other(2)) - act_frac_sys(ith_other_elem, other_id_shared(2))];


        % Calculate actual angle:
        angle_deg = acos((vec_ith_segm' * vec_other_segm) ...
                    / (norm(vec_ith_segm)*norm(vec_other_segm)) ) * 180 / pi;

        % Check angle w.r.t. tolerance:
        if abs(angle_deg) < tolerance_zero
            % Segment has same start or end node AND has zero angle:
            % Either duplicate or partially overlapping segment, remove
            % smallest segment:
            len_ith_segm = norm(vec_ith_segm);
            len_other_segm = norm(vec_other_segm);
            
            if len_ith_segm < len_other_segm
                removed_segm(ith_segm) = ith_segm;
            else
                removed_segm(ith_other_elem) = ith_other_elem;
            end
        end
    end
    end
end

act_frac_sys = act_frac_sys(removed_segm == 0, :);
frac_set_vec = frac_set_vec(removed_segm == 0);

% Make sure segments are sorted from X_min to X_max:
act_frac_sys_dummy = act_frac_sys;
indices = act_frac_sys_dummy(:, 1) > act_frac_sys_dummy(:, 3);
act_frac_sys(indices, 1:2) = act_frac_sys_dummy(indices, 3:4);
act_frac_sys(indices, 3:4) = act_frac_sys_dummy(indices, 1:2);

% Extract unique segments:
act_frac_sys_dummy = act_frac_sys;
[act_frac_sys, id_frac_set, ~] = unique(act_frac_sys_dummy, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);