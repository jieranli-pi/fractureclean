%% Main code for fracture network processing 
clear all; close all; clc;

% Define several input parameters:
% ratio_char_scales = [0.25, 0.5, 1, 2, 4, 8];  % relative char. scales
ratio_char_scales = [0.5, 1, 2, 4, 8];  % relative char. scales
apply_rotation = false;

% Directory with raw fracture network files
ith_real = 1;
BASE_DIR = 'tln_7_ef_high_filtered\\';
IN_DIR = ['..\\input\\' BASE_DIR];

% Output directory with cleaned (processed) fracture networks
OUT_DIR = ['..\\output\\' BASE_DIR];
if ~exist(OUT_DIR, 'dir')
   mkdir(OUT_DIR)
end

char_len = 15; % Absolute char. length
vec_char_len = ratio_char_scales * char_len;
num_iters = length(vec_char_len);
for ith_scale = 1:num_iters
    OUT_DIR_SCALE = [OUT_DIR 'char_scale_' num2str(vec_char_len(ith_scale))];
    if ~exist(OUT_DIR_SCALE, 'dir')
        mkdir(OUT_DIR_SCALE)
    end
end

% Additional parameters for script:
tolerance_frac = 0.5; % Relative cut-off value (fraction from intersection) 
tolerance_midpt = 0.45; % Relative cut-off value (radius around midpoint) 
tolerance_rad = 0.5; % Relative cut-off value (radius around end nodes) 
tolerance_zero = 1e-4; % Cut-off value for duplicated nodes
tolerance_intersect = 1e-10; % Cut-off value for duplicated nodes
num_decimals = 5; % Number of decimals accuracy (applicable for small datasets)
remove_duplicates = true; % Set true if want to remove duplicate fracture-segment
norm_order = 100; % L_n where n is the order of the norm: (sum(x - x0)^L_n)^(1/L_n)

% If no small segments occur for straight fractures, so usually this only
% happens when you have image interpreted fractures (automatically), then
% turn both on. Otherwise, only straight the fractures post cleaning! 
calc_init_intersections = 1;
straighten_fracs_pre = 1; % Straight fracture pre-cleaning
straighten_fracs_post = 1; % Straight fracture post-cleaning
tol_angle = 7.5; % Angle tolerance for straightness fractures

% Visualization related parameters:
plot_frac_network_pre = false;
plot_statistics = false;
plot_frac_network_post = true;

% Start main algorith:   
% NOTE: here only segment list is present, sometimes data contains two 
%   additional columns with information of fracture sets in the form
%   | large_set_id | sub_set_id | x_1 | y_1 | x_2 | y_2 | 
%   here large_set_id and sub_set_id are assumed zero (no info on frac sets
frac_data_full = load([IN_DIR 'real_' num2str(ith_real) '.txt']);
num_main_segm = size(frac_data_full, 1);
if size(frac_data_full, 2) == 4
    frac_data = zeros(num_main_segm, 6);
    frac_data(:, [3:6]) = round(frac_data_full(:, [1:4]) * 10^num_decimals)/(10^num_decimals);
else
    frac_data = round(frac_data_full * 10^num_decimals)/(10^num_decimals);
end
old_frac_data = frac_data;
act_frac_sys = frac_data(:, [3:6]);

if apply_rotation
    % Call script to optimally rotate fractures:
    [act_frac_sys_rot, rot_mat, bin_count, opt_angle] = ...
        calc_angle_frac_sys(act_frac_sys, 6);
    
    if plot_frac_network_pre
        figure();
        subplot(1, 2, 1);
        plot(act_frac_sys(:, [1, 3])', act_frac_sys(:, [2, 4])', ...
            'LineWidth', 1, 'color', [0, 0, 0])
        title('Original fracture network')
        
        subplot(1, 2, 2);
        plot(act_frac_sys_rot(:, [1, 3])', act_frac_sys_rot(:, [2, 4])', ...
            'LineWidth', 1, 'color', [0, 0, 0])
        title(['Fractures after rotation of ' num2str(opt_angle) ' degrees'])
    end
    %act_frac_sys = act_frac_sys_rot;
    %frac_data(:, [3:6]) = act_frac_sys;  % why? 
end

% Vector containing all to which fracture set each segment belong:
% (assumed zero here, see explanation above)
frac_set_vec = zeros(size(frac_data, 1), 1);

if remove_duplicates
    % Remove any collapsed (when both nodes of the segment are mapped to the
    % same point):
    [act_frac_sys, frac_set_vec] = extract_unique_segm(act_frac_sys, frac_set_vec, tolerance_zero);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
    
    % Calculate new main segments length (if duplicated nodes exists):
    len_main_segm = sqrt( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                          (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 ); 
    
    % Calculate new number of nodes:
    num_main_segm = length(len_main_segm);
    
    % Update fracture dataset:
    frac_data = zeros(num_main_segm, 6);
    frac_data(:, [3:6]) = act_frac_sys;
end

% Create order in which to discretize and add segments:
[sorted_len_main_segm, order_to_discr] = sort(len_main_segm, 'descend');

if plot_statistics
    figure();
    subplot(1,2,1);
    [counts,bins] = hist(len_main_segm, 25); %# get counts and bin locations
    barh(bins,counts)
    % xticks([])
    ylim([0, max(len_main_segm)*1.1])
    xlabel('Frequency', 'Interpreter', 'latex')
    ylabel('Length segments', 'Interpreter', 'latex')
    title('Histogram', 'Interpreter', 'latex')
    set(gca, 'FontSize', 14)

    subplot(1,2,2);
    plot(sort(len_main_segm), 'LineWidth', 2, 'color', [0, 0, 1])
    hold on
    cut_off = plot([1, length(len_main_segm)], [char_len, char_len], 'LineWidth', 2, 'color', [1, 0, 0]);
    hold off
    legend(cut_off, {'l_c'}, 'location', 'northwest')
    ylim([0, max(len_main_segm)*1.1])
    ylabel('Length segments', 'Interpreter', 'latex')
    title('Cumulative distr.', 'Interpreter', 'latex')
    set(gca, 'FontSize', 14)
end

%% Pre-processing algorithm:
if straighten_fracs_pre == 1
    [act_frac_sys, frac_set_vec] = ...
        straighten_frac_segm(act_frac_sys, frac_set_vec, tolerance_zero, tol_angle);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
end
if calc_init_intersections
    [act_frac_sys, frac_set_vec] = ...
        calc_intersections_segm(act_frac_sys, frac_set_vec, tolerance_intersect, tolerance_zero);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
end

[act_frac_sys] = adjust_vertical_segm(act_frac_sys, tolerance_zero);
act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);

num_inner_iter = 2; % determines how many times segment points are checked 

figure(); 
plot(act_frac_sys(:, [1, 3])', act_frac_sys(:, [2, 4])', 'color', [0, 0, 0])
hold on
plot(act_frac_sys(:, [1, 3])', act_frac_sys(:, [2, 4])', '.', 'color', [1, 0, 0])
hold off

% Try iterative scheme:
for ith_iter = 1:num_iters
    % Update new char_len:
    char_len_temp = vec_char_len(ith_iter);
    
    % Perform 1 inner iteration: 
    for ith_inner_iter = 1:num_inner_iter
        % Run main sequential fracture discretization and cleaning script:
        [act_frac_sys, frac_set_vec] = ...
            closest_point_method(act_frac_sys, frac_set_vec, ...
            norm_order, char_len_temp, tolerance_rad, tolerance_zero); 
        act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);

        % Perform an extra step of straightening too much refined segments:
        if straighten_fracs_post == 1
            % Execute script for straigthening fractures (pre-processing step):     
            [act_frac_sys, frac_set_vec] = ...
                straighten_frac_segm(act_frac_sys, frac_set_vec, tolerance_zero, tol_angle);
            act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
        end

        % Make sure intersections are caught:
        [act_frac_sys, frac_set_vec] = ...
            calc_intersections_segm(act_frac_sys, frac_set_vec, tolerance_intersect, tolerance_zero);
        act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
    end
    
    % Change vertical fractures:
    [act_frac_sys] = adjust_vertical_segm(act_frac_sys, tolerance_zero);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
    
    % Call new script:
    [act_frac_sys, frac_set_vec] = ...
        find_partial_overlap_and_small_angles(act_frac_sys, frac_set_vec, tolerance_zero, char_len_temp);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);

    % Remove any collapsed (when both nodes of the segment are mapped to the
    % same point):
    [act_frac_sys, frac_set_vec] = extract_unique_segm(act_frac_sys, frac_set_vec, tolerance_zero);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
    
    % Calculate intersections of all fractures:
    [act_frac_sys, frac_set_vec] = ...
        calc_intersections_segm(act_frac_sys, frac_set_vec, tolerance_intersect, tolerance_zero);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
    
    % Remove any collapsed (when both nodes of the segment are mapped to the
    % same point):
    [act_frac_sys, frac_set_vec] = extract_unique_segm(act_frac_sys, frac_set_vec, tolerance_zero);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
    
    % Call new script:
    [act_frac_sys, frac_set_vec] = ...
        find_partial_overlap_and_small_angles(act_frac_sys, frac_set_vec, tolerance_zero, char_len_temp);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
    
    % Determine total amount of fracture segments and lengths of each segment:
    num_tot_new_segm = size(act_frac_sys, 1);
    len_segm_new = sqrt( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                         (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 );
                     
    [act_frac_sys, frac_set_vec] = ...
        find_actual_overlap_segm(act_frac_sys, frac_set_vec, tolerance_zero);
    act_frac_sys = round(act_frac_sys * 10^num_decimals) * 10^(-num_decimals);
    
    % Finally re-rotate the fracture network back to original position:
    if apply_rotation
        % Rotate frac_sys back to original orientation:
        act_frac_sys_rot = zeros(size(act_frac_sys));
        act_frac_sys_rot(:, [1, 2]) = act_frac_sys(:, [1, 2]) * rot_mat;
        act_frac_sys_rot(:, [3, 4]) = act_frac_sys(:, [3, 4]) * rot_mat;
        act_frac_sys = act_frac_sys_rot;
    end

    % Save cleaned act_frac_sys to file
    OUT_DIR_SCALE = [OUT_DIR 'char_scale_' num2str(vec_char_len(ith_iter))];
    fid3 = fopen([OUT_DIR_SCALE '\\real_' num2str(ith_real) '.txt'],'w+');
    for ii = 1:size(act_frac_sys, 1)
        fprintf(fid3,'%8.5f \t %8.5f \t %8.5f \t %8.5f\n', act_frac_sys(ii, :));
    end
    fclose(fid3);
end

if plot_frac_network_post
    % Plot results:
    figure();
    sizefig = length(vec_char_len) + 1;
    subplot(2, ceil(sizefig/2), 1);
    plot(old_frac_data(:, [3, 5])', old_frac_data(:, [4, 6])', 'color', [0, 0, 0])
    hold on
    plot(old_frac_data(:, [3, 5])', old_frac_data(:, [4, 6])', '.', 'color', [1, 0, 0])
    hold off
    title('Raw network')
    
    for ii = 1:length(vec_char_len)
        char_len_str = num2str(vec_char_len(ii));
        clean_frac_sys = load([OUT_DIR 'char_scale_' char_len_str '\\real_' num2str(ith_real) '.txt']);
        subplot(2, ceil(sizefig/2), ii + 1);
        plot(clean_frac_sys(:, [1, 3])', clean_frac_sys(:, [2, 4])', 'color', [0, 0, 0])
        hold on
        plot(clean_frac_sys(:, [1, 3])', clean_frac_sys(:, [2, 4])', '.', 'color', [1, 0, 0])
        hold off
        title(['Clean Network lc = ' char_len_str])
    end
end