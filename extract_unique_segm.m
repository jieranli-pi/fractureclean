function [act_frac_sys, frac_set_vec] = extract_unique_segm(act_frac_sys, frac_set_vec, tolerance_zero)

% Extract unique rows:
[act_frac_sys_dummy, id_frac_set, ~] = unique(act_frac_sys, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);

% Determine length of new "main" segments:
len_segm_new = sqrt( (act_frac_sys_dummy(:, 1) - act_frac_sys_dummy(:, 3)).^2 + ...
                     (act_frac_sys_dummy(:, 2) - act_frac_sys_dummy(:, 4)).^2 );

% Remove non-zero fracture segments:                      
nonzero_segm = len_segm_new>tolerance_zero;
act_frac_sys_dummy_2 = act_frac_sys_dummy(nonzero_segm, :); 
frac_set_vec = frac_set_vec(nonzero_segm);

% Filter on y-coordinate and check for unique:
act_frac_sys_dummy_3 = act_frac_sys_dummy_2;
indices = act_frac_sys_dummy_2(:, 2) > act_frac_sys_dummy_2(:, 4);
act_frac_sys_dummy_3(indices, 1:2) = act_frac_sys_dummy_2(indices, 3:4);
act_frac_sys_dummy_3(indices, 3:4) = act_frac_sys_dummy_2(indices, 1:2);

% Extract once more the unique rows:
[act_frac_sys_dummy_4, id_frac_set, ~] = unique(act_frac_sys_dummy_3, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);

% Filter on x-coordinate and check again for unique:
act_frac_sys_dummy_5 = act_frac_sys_dummy_4;
indices = act_frac_sys_dummy_4(:, 1) > act_frac_sys_dummy_4(:, 3);
act_frac_sys_dummy_5(indices, 1:2) = act_frac_sys_dummy_4(indices, 3:4);
act_frac_sys_dummy_5(indices, 3:4) = act_frac_sys_dummy_4(indices, 1:2);

% Extract once more the unique rows:
[act_frac_sys, id_frac_set, ~] = unique(act_frac_sys_dummy_5, 'rows');
frac_set_vec = frac_set_vec(id_frac_set);