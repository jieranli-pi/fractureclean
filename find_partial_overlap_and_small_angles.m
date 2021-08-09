function [act_frac_sys, frac_set_vec] = ...
    find_partial_overlap_and_small_angles(act_frac_sys, frac_set_vec, tolerance_zero, char_len)

% Script to remove segments that lie on top of each other but are different
% length:
norm_order = 2;
dist_crit_midpt = 0.35;

% Make sure segments are sorted from X_min to X_max:
% % % act_frac_sys = act_frac_sys_dummy;
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
            else
                % Node has been split and no longer coincides with this node! 
                continue;
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

            elseif abs(angle_deg) < (atan(dist_crit_midpt) * 180 / pi)
                % Segments intersect at too low angle
                % Find which vector is the smallest (this point will be moved):
                col_ids_shared = [];
                col_ids_other = [];
                col_ids_move = [];

                len_ith_segm = norm(vec_ith_segm);
                len_other_segm = norm(vec_other_segm);

                if len_ith_segm < len_other_segm
                    id_to_keep = ith_other_elem;
                    id_to_move = ith_segm;
                    col_ids_shared = other_id_shared;
                    col_ids_other = other_id_other;
                    col_ids_move = ithsegm_id_other;
                else
                    id_to_keep = ith_segm;
                    id_to_move = ith_other_elem;
                    col_ids_shared = ithsegm_id_shared;
                    col_ids_other = ithsegm_id_other;
                    col_ids_move = other_id_other; 
                end

                % Set vector to keep and move:
                if id_to_keep == ith_other_elem
                    vec_to_keep = vec_other_segm;
                    vec_point_move = vec_ith_segm;
                else
                    vec_to_keep = vec_ith_segm;
                    vec_point_move = vec_other_segm;
                end

                t = vec_point_move' * vec_to_keep / (norm(vec_to_keep)^2);

                % Intersection on the left of right:
                new_coord_p = [act_frac_sys(id_to_keep, col_ids_shared(1)) + vec_to_keep(1)*t, ...
                               act_frac_sys(id_to_keep, col_ids_shared(2)) + vec_to_keep(2)*t];


                % Check if new node is not too close to original node:
                old_node_coord = act_frac_sys(id_to_move, col_ids_move);
                dist_new_node = ( (act_frac_sys(id_to_keep, col_ids_other(1)) - new_coord_p(1)).^2 + ...
                                  (act_frac_sys(id_to_keep, col_ids_other(2)) - new_coord_p(2)).^2).^(1/2);

                % Move or merge:
                merge = dist_new_node < (char_len / 4);
                if merge
                    new_coord_p = act_frac_sys(id_to_keep, col_ids_other);
                    removed_segm(id_to_move) = id_to_move;
                end

               % Find all segments that have the node that is required to
               % move:
               dist_leftcol = ( (act_frac_sys(:, 1) - old_node_coord(1)).^norm_order + ...
                                (act_frac_sys(:, 2) - old_node_coord(2)).^norm_order).^(1/norm_order);
               ids_leftcol = find(dist_leftcol < tolerance_zero);

               dist_rightcol = ( (act_frac_sys(:, 3) - old_node_coord(1)).^norm_order + ...
                                 (act_frac_sys(:, 4) - old_node_coord(2)).^norm_order).^(1/norm_order);
               ids_rightcol = find(dist_rightcol < tolerance_zero);

               ids_incl = [id_to_keep, id_to_move, ids_leftcol', ids_rightcol'];

               % Move/merge nodes:
               nr_move_leftnodes = length(ids_leftcol);
               if ~isempty(ids_leftcol)
                   act_frac_sys(ids_leftcol, [1, 2]) = ones(nr_move_leftnodes, 1) * new_coord_p;
               end

               nr_move_rightnodes = length(ids_rightcol);
               if ~isempty(ids_rightcol)
                   act_frac_sys(ids_rightcol, [3, 4]) = ones(nr_move_rightnodes, 1) * new_coord_p;
               end

               % Re-sort move/merged segments on x1 < x2:
               move_merge_segm = union(ids_leftcol', ids_rightcol');
               indices_to_swap = act_frac_sys(move_merge_segm, 1) > act_frac_sys(move_merge_segm, 3);

               dummy_temp = act_frac_sys(move_merge_segm, :);
               if any(indices_to_swap)
                   act_frac_sys(move_merge_segm(indices_to_swap), [1, 2]) = dummy_temp(indices_to_swap, [3, 4]);
                   act_frac_sys(move_merge_segm(indices_to_swap), [3, 4]) = dummy_temp(indices_to_swap, [1, 2]);
               end

               if ~merge
                   % Try to split segments already during loop and add at the end
                   % to act_frac_sys:
                   if new_coord_p(1) > act_frac_sys(id_to_keep, 1)
                       old_end_coord = act_frac_sys(id_to_keep, [3, 4]);

                       act_frac_sys(id_to_keep, :) = [act_frac_sys(id_to_keep, 1), ...
                                                      act_frac_sys(id_to_keep, 2), ...
                                                      new_coord_p(1), ...
                                                      new_coord_p(2)];

                       % Create entry to extra segment:
                       extra_frac_sys = [new_coord_p(1), new_coord_p(2), old_end_coord];                          
                   else
                       % Nodes are swapped
                       old_end_coord = act_frac_sys(id_to_keep, [1, 2]);

                       act_frac_sys(id_to_keep, :) = [act_frac_sys(id_to_keep, 3), ...
                                                      act_frac_sys(id_to_keep, 4), ...
                                                      new_coord_p(1), ...
                                                      new_coord_p(2)];

                       % Create entry to extra segment:
                       extra_frac_sys = [new_coord_p(1), new_coord_p(2), old_end_coord];
                   end

                   % Add extra segment (after splitting to 
                   act_frac_sys = [act_frac_sys; extra_frac_sys];
                   frac_set_vec = [frac_set_vec; frac_set_vec(id_to_keep, 1);];
                   removed_segm = [removed_segm; 0];
               end

               % Check for unique/collapsed segments:
               for kk = move_merge_segm
                   % Compare with first segment:
                   % (happens always, even when segments are merged instead of move)
                   % NOTE: ONLY WORKS WHEN SORTED!!!
                   if norm(act_frac_sys(id_to_keep, :) - act_frac_sys(kk, :)) < tolerance_zero
                       % Segments are the saim:
                       % Remove moved segment:
                       removed_segm(kk) = kk;
                   end

                   if ~merge
                       % Check also split segment:
                       if norm(act_frac_sys(end, :) - act_frac_sys(kk, :)) < tolerance_zero
                           % Segments are the saim:
                           % Remove moved segment:
                           removed_segm(kk) = kk;
                       end
                   end    
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

% Extract only 
len_segm_new = sqrt( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                     (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 );
nonzero_segm = len_segm_new>tolerance_zero;
act_frac_sys = act_frac_sys(nonzero_segm, :); 
frac_set_vec = frac_set_vec(nonzero_segm);