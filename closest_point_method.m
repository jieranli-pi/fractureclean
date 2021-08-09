function [act_frac_sys, frac_set_vec] = ...
    closest_point_method (act_frac_sys, frac_set_vec, ...
    norm_order, char_len, tolerance_rad, tolerance_zero)

% Maximum number of segments:
num_tot_new_segm_reduced = size(act_frac_sys, 1);
len_segm_new = sqrt( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                     (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 );
                      
num_segm_per_frac = round(len_segm_new/char_len);
tot_small_segm = sum(num_segm_per_frac) + sum(num_segm_per_frac==0);

% Allocate some matrices:
new_frac_set_list = zeros(tot_small_segm, 1); % Allocate memory for fracture set list
new_frac_sys = zeros(tot_small_segm, 4); % Allocate max size new frac sys
unq_node_list = ones(tot_small_segm*2, 2) * 9999999; % Max number of nodes if 2*nr_segm
total_loop = 0; % Counter for total number of segments
total_nodes = 0; % Counter for total number of (unique) nodes

% Create order in which to discretize and add segments:
[~, order_to_discr] = sort(len_segm_new, 'descend');

% Change order slightly based on segments that penetrate the well
num_wells = 5;
well_coords = zeros(num_wells, 2);
well_coords(1, :) = [100, 100];
well_coords(2, :) = [900, 100]; 
well_coords(3, :) = [500, 500];
well_coords(4, :) = [100, 900]; 
well_coords(5, :) = [900, 900];
method_ii = 1;
frac_length_threshold = 15;

[order_to_discr] = change_order_discr(act_frac_sys, len_segm_new, ...
    order_to_discr, num_wells, well_coords, frac_length_threshold, method_ii);


for ith_main_segm = 1:num_tot_new_segm_reduced
    % For each main segment, populate the discretized (smaller) segments in
    % the model and find conflicts with currently existing (smaller)
    % fracture segments in the model:
    % Get x- and y-coordinate of main segment (to be discretized):
    sorted_ith_main_segm = order_to_discr(ith_main_segm);
    x_coor = [act_frac_sys(sorted_ith_main_segm, 1); act_frac_sys(sorted_ith_main_segm, 3)];
    y_coor = [act_frac_sys(sorted_ith_main_segm, 2); act_frac_sys(sorted_ith_main_segm, 4)];
    
    % Try fix 1: always sort from x+ and y+:
    % Other fix is to modify and node from the frac_data matrix which is
    % merged such that any other segment that contains these nodes will 
    % merge accordingly.
    
    % Discretize main segment:
    [segm_mat] = discr_segm(x_coor, y_coor, char_len);
    
    % Check if conflicts with existing discrete segments:
    for ith_new_segm = 1:size(segm_mat, 1)
        total_loop = total_loop + 1;
        
        % Check for each node to be added if within radius of other node:
        % Store coordinate of new (to be added) discretized segment:
        x_new_coor = segm_mat(ith_new_segm, [1, 3]);
        y_new_coor = segm_mat(ith_new_segm, [2, 4]);
        
        % Check node 1:
        dist_node_1 = ( (x_new_coor(1) - unq_node_list(1:total_nodes, 1)).^norm_order + ...
                        (y_new_coor(1) - unq_node_list(1:total_nodes, 2)).^norm_order ).^(1/norm_order);
        ids_1 = find(dist_node_1 < tolerance_rad*char_len);
        
        if ~isempty(ids_1)
            % Merge node with closest existing point in radius:
            [val, id_min] = min(dist_node_1(ids_1));
            
            if val > tolerance_zero
                % Find if node exists in left column:
                dist_leftcol = ( (act_frac_sys(:, 1) - segm_mat(ith_new_segm, 1)).^norm_order + ...
                                 (act_frac_sys(:, 2) - segm_mat(ith_new_segm, 2)).^norm_order).^(1/norm_order);
                ids_leftcol = find(dist_leftcol < tolerance_zero);
                if ~isempty(ids_leftcol)
                    act_frac_sys(ids_leftcol, 1) = unq_node_list(ids_1(id_min), 1);
                    act_frac_sys(ids_leftcol, 2) = unq_node_list(ids_1(id_min), 2);
                end

                % Find if node exists in left column:
                dist_rightcol = ( (act_frac_sys(:, 3) - segm_mat(ith_new_segm, 1)).^norm_order + ...
                                  (act_frac_sys(:, 4) - segm_mat(ith_new_segm, 2)).^norm_order).^(1/norm_order);
                ids_rightcol = find(dist_rightcol < tolerance_zero);
                if ~isempty(ids_rightcol)
                    act_frac_sys(ids_rightcol, 3) = unq_node_list(ids_1(id_min), 1);
                    act_frac_sys(ids_rightcol, 4) = unq_node_list(ids_1(id_min), 2);
                end
                % Do the same of discretized segments:
                dist_leftcol = ( (segm_mat(:, 1) - segm_mat(ith_new_segm, 1)).^norm_order + ...
                                 (segm_mat(:, 2) - segm_mat(ith_new_segm, 2)).^norm_order).^(1/norm_order);
                ids_leftcol = find(dist_leftcol < tolerance_zero);
                if ~isempty(ids_leftcol)
                    segm_mat(ids_leftcol, 1) = unq_node_list(ids_1(id_min), 1);
                    segm_mat(ids_leftcol, 2) = unq_node_list(ids_1(id_min), 2);
                end

                dist_rightcol = ( (segm_mat(:, 3) - segm_mat(ith_new_segm, 1)).^norm_order + ...
                                  (segm_mat(:, 4) - segm_mat(ith_new_segm, 2)).^norm_order).^(1/norm_order);
                ids_rightcol = find(dist_rightcol < tolerance_zero);
                if ~isempty(ids_rightcol)
                    segm_mat(ids_rightcol, 3) = unq_node_list(ids_1(id_min), 1);
                    segm_mat(ids_rightcol, 4) = unq_node_list(ids_1(id_min), 2);
                end
            end
            
            segm_mat(ith_new_segm, 1:2) = unq_node_list(ids_1(id_min), :);
            new_frac_sys(total_loop, 1:2) = segm_mat(ith_new_segm, 1:2);
            new_frac_set_list(total_loop) = frac_set_vec(sorted_ith_main_segm);
        else
            % Add new non-existing node:
            total_nodes = total_nodes + 1;
            unq_node_list(total_nodes, :) = [x_new_coor(1), y_new_coor(1)];
            new_frac_sys(total_loop, 1:2) = [x_new_coor(1), y_new_coor(1)];
            new_frac_set_list(total_loop) = frac_set_vec(sorted_ith_main_segm);
        end
        
        % Check node 2:
        dist_node_2 = ( (x_new_coor(2) - unq_node_list(1:total_nodes, 1)).^norm_order + ...
                        (y_new_coor(2) - unq_node_list(1:total_nodes, 2)).^norm_order ).^(1/norm_order);
        ids_2 = find(dist_node_2 < tolerance_rad*char_len);
        
        if ~isempty(ids_2)
            % Merge node with closest existing point in radius:
            [val, id_min] = min(dist_node_2(ids_2));
            
            if val > tolerance_zero
                % Find if node exists in left column:
                dist_leftcol = ( (act_frac_sys(:, 1) - segm_mat(ith_new_segm, 3)).^norm_order + ...
                                 (act_frac_sys(:, 2) - segm_mat(ith_new_segm, 4)).^norm_order).^(1/norm_order);
                ids_leftcol = find(dist_leftcol < tolerance_zero);
                if ~isempty(ids_leftcol)
                    act_frac_sys(ids_leftcol, 1) = unq_node_list(ids_2(id_min), 1);
                    act_frac_sys(ids_leftcol, 2) = unq_node_list(ids_2(id_min), 2);
                end

                % Find if node exists in left column:
                dist_rightcol = ( (act_frac_sys(:, 3) - segm_mat(ith_new_segm, 3)).^norm_order + ...
                                  (act_frac_sys(:, 4) - segm_mat(ith_new_segm, 4)).^norm_order).^(1/norm_order);
                ids_rightcol = find(dist_rightcol < tolerance_zero);
                if ~isempty(ids_rightcol)
                    act_frac_sys(ids_rightcol, 3) = unq_node_list(ids_2(id_min), 1);
                    act_frac_sys(ids_rightcol, 4) = unq_node_list(ids_2(id_min), 2);
                end

                % Do the same of discretized segments:
                dist_leftcol = ( (segm_mat(:, 1) - segm_mat(ith_new_segm, 3)).^norm_order + ...
                                 (segm_mat(:, 2) - segm_mat(ith_new_segm, 4)).^norm_order).^(1/norm_order);
                ids_leftcol = find(dist_leftcol < tolerance_zero);
                if ~isempty(ids_leftcol)
                    segm_mat(ids_leftcol, 1) = unq_node_list(ids_2(id_min), 1);
                    segm_mat(ids_leftcol, 2) = unq_node_list(ids_2(id_min), 2);
                end

                dist_rightcol = ( (segm_mat(:, 3) - segm_mat(ith_new_segm, 3)).^norm_order + ...
                                  (segm_mat(:, 4) - segm_mat(ith_new_segm, 4)).^norm_order).^(1/norm_order);
                ids_rightcol = find(dist_rightcol < tolerance_zero);
                if ~isempty(ids_rightcol)
                    segm_mat(ids_rightcol, 3) = unq_node_list(ids_2(id_min), 1); 
                    segm_mat(ids_rightcol, 4) = unq_node_list(ids_2(id_min), 2);
                end
            end
            
            segm_mat(ith_new_segm, 3:4) = unq_node_list(ids_2(id_min), :);
            new_frac_sys(total_loop, 3:4) = segm_mat(ith_new_segm, 3:4);
            new_frac_set_list(total_loop) = frac_set_vec(sorted_ith_main_segm);
        else
            % Add new non-existing node:
            total_nodes = total_nodes + 1;
            unq_node_list(total_nodes, :) = [x_new_coor(2), y_new_coor(2)];
            new_frac_sys(total_loop, 3:4) = [x_new_coor(2), y_new_coor(2)];
            new_frac_set_list(total_loop) = frac_set_vec(sorted_ith_main_segm);
        end
    end
end

act_frac_sys = new_frac_sys(1:total_loop, :);
frac_set_vec = new_frac_set_list(1:total_loop);