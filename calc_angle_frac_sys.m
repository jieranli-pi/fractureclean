function [act_frac_sys_rot, rot_mat, bin_count, opt_angle] = ...
    calc_angle_frac_sys(act_frac_sys, num_bins)

% Calculate angle of all fractures in system weighted by their length:
num_frac = size(act_frac_sys, 1);
length_fracs = sqrt( (act_frac_sys(:, 1) - act_frac_sys(:, 3)).^2 + ...
                     (act_frac_sys(:, 2) - act_frac_sys(:, 4)).^2 );
tot_length = sum(length_fracs);
angle_fracs = zeros(num_frac, 1);
weighted_angle = length_fracs ./ tot_length;
bin_size = 90 / num_bins;
bin_count = zeros(num_bins, 1);

% Set reference vector w.r.t. x-axis
ref_vec = [1; 0];  

for ith_frac = 1:num_frac
    % Calculate vector of fracture segment:
    frac_vec = [act_frac_sys(ith_frac, 3) - act_frac_sys(ith_frac, 1);
                act_frac_sys(ith_frac, 4) - act_frac_sys(ith_frac, 2)];
    frac_vec = frac_vec / norm(frac_vec);
    
    % Angle:
    angle = acos( ref_vec' * frac_vec );
    angle = angle * 180 / pi;
    
    while angle > 90 || angle < 0
       if angle > 90
           angle = angle - 90;
       else
          angle = angle + 90; 
       end
    end
    
    angle_fracs(ith_frac) = angle;
    
    % Find bin (bin 1 == 0 degrees and loops around to 90 - bin_size/2, 
    % since rotating by 0 or 90 degrees in this case is equivalent):
    if (angle < bin_size/2) || (angle > (90 - bin_size/2))
        bin_ii = 1;
    else
        bin_ii = ceil((angle + bin_size/2)/bin_size);
    end
    bin_count(bin_ii) = bin_count(bin_ii) + weighted_angle(ith_frac);
end

% Plot figure:
% % % figure();
bin_center = (bin_size/2):bin_size:(90 - bin_size/2);
bin_center = 0:bin_size:(90 - bin_size);
% % % plot(bin_center, bin_count, 'LineWidth', 2, 'color', [0, 0, 0])

[val, id] = max(bin_count);
opt_angle = bin_center(id) * pi / 180;

% Apply optimal rotation:
rot_mat = [cos(opt_angle), -sin(opt_angle); 
           sin(opt_angle), cos(opt_angle)];

% Rotate frac_sys:
act_frac_sys_rot = zeros(size(act_frac_sys));
act_frac_sys_rot(:, [1, 2]) = act_frac_sys(:, [1, 2]) * rot_mat';
act_frac_sys_rot(:, [3, 4]) = act_frac_sys(:, [3, 4]) * rot_mat';