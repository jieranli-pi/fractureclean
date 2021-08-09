function [order_to_discr] = change_order_discr(act_frac_sys, len_segm_new, ...
    order_to_discr, num_wells, well_coords, frac_length_threshold, method_ii)

% Top left and bottom right corners:
% Step 1) Find fracture that is closest to the well
%   Based on centroid? or based on one of the coordinates?
%   How to make sure that in sim. always same fracture is penetrated?
% Step 2) Set these two as the starting fractures
%   Make new empty list --> fill first two as top and bot frac --> fill
%   rest as ~= top and bot frac :)
% Continue normally
if num_wells ~= size(well_coords, 1)
    error('Please specify correctly the well-coordinates!')
end
num_segm = length(len_segm_new);
well_frac_ids = zeros(num_wells, 1);

% Step 1):
for ith_well = 1:num_wells
   % Find closest fracture to well coord:
   if method_ii == 1
       % Method 1: centroid
       centroids = [mean(act_frac_sys(:, [1, 3]), 2), mean(act_frac_sys(:, [2, 4]), 2)];
       [val_ii, id_ii] = sort( ( (centroids(:, 1) - well_coords(ith_well, 1)).^2 + ...
                                 (centroids(:, 2) - well_coords(ith_well, 2)).^2).^(1/2) );
                             
       % Find closest fracture of certain length treshold:
       segm_length_closest = len_segm_new(id_ii);
       id_ii_2 = find(segm_length_closest >= frac_length_threshold, 1);
       id_ii_real = id_ii(id_ii_2);
       val_ii_real = val_ii(id_ii_2);
       
   else
       % Method 2: either of segment coordinates
       [val_1, id_1] = sort( ( (act_frac_sys(:, 1) - well_coords(ith_well, 1)).^2 + ...
                               (act_frac_sys(:, 2) - well_coords(ith_well, 2)).^2).^(1/2) );

       [val_2, id_2] = sort( ( (act_frac_sys(:, 3) - well_coords(ith_well, 1)).^2 + ...
                               (act_frac_sys(:, 4) - well_coords(ith_well, 2)).^2).^(1/2) );
        
       % Find closest fracture of certain length treshold:
       segm_length_closest = len_segm_new(id_1);
       id_ii_1 = find(segm_length_closest >= frac_length_threshold, 1);
       id_ii_real_1 = id_1(id_ii_1);
       val_ii_real_1 = val_1(id_ii_1);       
       
       segm_length_closest = len_segm_new(id_2);
       id_ii_2 = find(segm_length_closest >= frac_length_threshold, 1);
       id_ii_real_2 = id_2(id_ii_2);
       val_ii_real_2 = val_2(id_ii_2);
                           
       if val_ii_real_1 < val_ii_real_2 
           id_ii_real = id_ii_real_1;
       else
           id_ii_real = id_ii_real_2;
       end
   end
   
   well_frac_ids(ith_well) = id_ii_real;
end

% Step 2):
order_to_discr_new = zeros(num_segm, 1);

% Check if all wells are in unique fractures:
[well_frac_ids_unq, ~] = unique(well_frac_ids);
num_unq_fracs = length(well_frac_ids_unq);
indices = true(num_segm, 1);
for ii = 1:num_unq_fracs
    indices(order_to_discr==well_frac_ids_unq(ii)) = false;
end
order_to_discr_new(1:num_unq_fracs) = well_frac_ids_unq;
order_to_discr_new((num_unq_fracs + 1):end) = order_to_discr(indices);

if sum(order_to_discr_new) == sum(1:num_segm)
    order_to_discr = order_to_discr_new;
else
    error('BUG ALERT!')
end