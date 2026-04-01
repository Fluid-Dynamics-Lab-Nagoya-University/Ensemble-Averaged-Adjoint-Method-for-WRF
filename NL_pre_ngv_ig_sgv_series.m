% Generates perturbed WRF initial conditions for nonlinear runs,
% combining adjoint sensitivity fields with background noise
% across a range of input energy rates and noise rates.
%
% by Shan Jiang, FDL, Nagoya University

clear
close all

fontsize = 20; fontname = "helvetica";

if isempty(gcp('nocreate'))
    parpool('local', 6);
end

input_mode      = 2; %1:all 2:positive only 3: negative only
input_pos       = 1; %1:all 2:actuators 3:actuators with constant input
s_mode          = 1; %1: ground 2: volume
n_mode          = 1; %1: ground 2: volume
QV_mode         = 2; %1: consider original QVAPOR 2: not consider QVAPOR

ngEQsg = 1;  % 1: ng equals sg; 0: ng fixed at 0

allTogether_flag = 0;  % 0: single run with above parameters; 1: run three combinations: inputMode1pos1,inputMode3pos1,inputMode3pos2


%==============mode tags setting============
if s_mode == 1
    s_tag = "sg"; % ground sensitivity noise
elseif s_mode == 2
    s_tag = "sv"; % volume sensitivity noise
else
    error('s_mode must be 1 (ground) or 2 (volume).');
end

if n_mode == 1
    n_tag = "ng"; % ground valid noise
elseif n_mode == 2
    n_tag = "nv"; % volume valid noise
else
    error('n_mode must be 1 (ng) or 2 (nv).');
end

if QV_mode == 1
    QV_tag = ""; % original QVAPOR considered
elseif QV_mode == 2
    QV_tag = "NOV0.02_"; % QVAPOR not considered, reference value 0.02
else
    error('QV_mode must be 1 (original QVAPOR considered) or 2 (not considered).');
end

pert_dir = "./pert";

%===========================================

input_eng_rate_list = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];
noise_rate_list = [0, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5];

% Pre-load all perturbation files (up to max_ens_member)
nx=50;
ny=50;
max_ens_member = 100;
noise_all = zeros(nx, ny, max_ens_member);

for j = 1:max_ens_member
    pert_file = fullfile(pert_dir, sprintf("pert%d.dat", j));
    if ~isfile(pert_file)
        error("Perturbation file not found: %s", pert_file);
    end
    fid_noise = fopen(pert_file, 'r');
    if fid_noise < 0
        error("Cannot open file: %s", pert_file);
    end
    noise_1d = fscanf(fid_noise, '%f', [nx*ny, 1]);
    fclose(fid_noise);
    if numel(noise_1d) ~= nx*ny
        error("Number of data in %s is %d, expected %d", ...
              pert_file, numel(noise_1d), nx*ny);
    end
    noise_all(:,:,j) = reshape(noise_1d, [nx, ny]);
end

% Read QVAPOR from wrfinput_d01_wide_6h_org once (reused inside parfor)
val        = ncread("./clean/wrfinput_d01_clean","QVAPOR");
init_QVAPOR = val(:,:,1);

% Determine (input_mode, input_pos) combinations to run
if allTogether_flag == 1
    % Three fixed combinations:
    % (1) input_mode=1, input_pos=1
    % (2) input_mode=3, input_pos=1
    % (3) input_mode=3, input_pos=2
    combo_list = [
        1 1;
        3 1;
        3 2;
    ];
else
    % Single run with current settings
    combo_list = [input_mode, input_pos];
end

for combo_i = 1:size(combo_list,1)
    input_mode = combo_list(combo_i,1);
    input_pos  = combo_list(combo_i,2);

    % Recompute i_tag and ipos_tag for each combo (ensures filename and logic consistency)
    if input_mode == 1
        i_tag = "";    % both
    elseif input_mode == 2
        i_tag = "P";   % positive only
    elseif input_mode == 3
        i_tag = "N";   % negative only
    else
        error('input_mode must be 1 (both) or 2(positive only) or 3(negative only)');
    end

    if input_pos == 1
        ipos_tag = "";       % all
    elseif input_pos == 2
        ipos_tag = "actr";   % actuators
    elseif input_pos == 3
        ipos_tag = "actrC";  % actuators with constant input
    else
        error('input_pos must be 1 (all) or 2 (actuators).');
    end

    for rate_idx = 1:length(input_eng_rate_list)
        input_eng_rate = input_eng_rate_list(rate_idx);


        for noise_idx = 1:length(noise_rate_list)
            sens_noise_rate = noise_rate_list(noise_idx);
            if ngEQsg == 1
                bgnoise_amp_valid_rate = noise_rate_list(noise_idx);
            else
                bgnoise_amp_valid_rate = 0;
            end

            % Determine ens_member based on bgnoise_amp_valid_rate
            if bgnoise_amp_valid_rate == 0
                ens_member = 1;
            else
                ens_member = 100;
            end

            sens_file = "./dats/A_QVAPOR_mean_absG_i0M_" + QV_tag + s_tag + sens_noise_rate + ".dat";
            data_flat = load(sens_file);

            s_absG_i0M = reshape(data_flat, [31, 50, 50]);
            s_absG_i0M = permute(s_absG_i0M, [3 2 1]);
            s_absG_i0M = s_absG_i0M(:,:,1);

            varname="QVAPOR";   % amount of steam

            dir_dataout_woinput = "./prepare/absG_i0M"+ "_" + QV_tag + n_tag +bgnoise_amp_valid_rate...
                + "_ig"+input_eng_rate+i_tag+ipos_tag+"_" + s_tag+sens_noise_rate;
            mkdir(dir_dataout_woinput)

    %===============actuator selection==============

            figure(16)

            imagesc(flipud(s_absG_i0M(:,:)'));
            colormap(jet)
            colorbar
            clim([-6000, 6000]);
            set(gca, 'FontName', fontname, 'FontSize', fontsize);
            title(sprintf('Actuators placement'));

            axis equal
            axis image

            hold on;

            if input_mode == 1
                [A_RAINNC_mean_sort, linear_idx] = sort(abs(s_absG_i0M(:)),'descend');
            elseif input_mode == 2
                [A_RAINNC_mean_sort, linear_idx] = sort(s_absG_i0M(:),'descend');
            elseif input_mode == 3 %default sorting in matlab is ascend
                [A_RAINNC_mean_sort, linear_idx] = sort(s_absG_i0M(:));
            end

            orders = [1, 2, 3, 4, 5, 6, 7, 8];
            selected_values = A_RAINNC_mean_sort(orders);
            selected_indices = linear_idx(orders);

            [row_idx, col_idx] = ind2sub(size(s_absG_i0M), selected_indices);

            disp('Situations of the 1-8 smallest value');
            for i = 1:length(orders)
                fprintf('%dth: %.4f, pos: (%d, %d)\n', orders(i), selected_values(i), row_idx(i), col_idx(i));
            end

            x_plot = row_idx;
            y_plot = size(s_absG_i0M,2) - col_idx + 1;

            scatter(x_plot, y_plot, 30, 'magenta', 'filled', 'MarkerEdgeColor','k');

            hold off;
    %=======================================================

            S_max = max(abs(s_absG_i0M(:)));                 % for normalization
            s_norm_abs = abs(s_absG_i0M) / S_max;            % normalize
            s_norm = s_norm_abs .* sign(s_absG_i0M);         % sign recover

            switch input_mode
                case 1 % both
                    % no-op
                case 2 % positive-only
                    s_norm(s_norm < 0) = 0;
                case 3 % negative-only
                    s_norm(s_norm > 0) = 0;
            end

            if QV_mode ==1  % consider original QVAPOR
                input_eng_g = s_norm.* input_eng_rate .* init_QVAPOR;
            else
                input_eng_g = s_norm.* input_eng_rate .* 0.02;
            end

            if input_pos == 3
                input_eng_g = -0.5 .* input_eng_rate .* init_QVAPOR; %constant input magnitude by actuator
            end

            % Adjust Frobenius norm to reach same energy
            if QV_mode == 1
                tmpFroNorm = s_norm .* init_QVAPOR;
            else
                tmpFroNorm = s_norm .* 0.02;
            end
            switch input_mode
                case 1 %both
                    target_froNorm = 0.095222;   % sg0.1 igNOV0.02 absG_i0M:0.095222
                case 2 %positive
                    target_froNorm = 0.076545;   % sg0.1 igNOV0.02 absG_i0M:0.076545
                case 3 %negative
                    target_froNorm = 0.056641;   % sg0.1 igNOV0.02 absG_i0M:0.052743
                    %note: the L2-norm set here makes maximum QVAPOR change is roughly 10% when sg=0.1
            end

            % Compute Frobenius norm
            current_froNorm = norm(tmpFroNorm,'fro');

            % Print case info and norm
            fprintf("Case: s_tag=%s, i_tag=%s, n_tag=%s, ig=%g, sens_noise=%g\n", ...
                    char(s_tag), char(i_tag), char(n_tag), input_eng_rate, sens_noise_rate);

            fprintf("    target_norm  = %.6f\n", target_froNorm);
            fprintf("    current_norm = %.6f\n\n", current_froNorm);

            % Scale by Frobenius norm
            input_eng_g = input_eng_g .* (target_froNorm./current_froNorm);

            actatr_mark = zeros(size(s_absG_i0M));
            if input_pos == 1
                actatr_mark = input_eng_g;
            else
                for i = 1:length(row_idx)
                    actatr_mark(row_idx(i), col_idx(i)) = input_eng_g(row_idx(i), col_idx(i));
                end
            end

            parfor j=1:ens_member
                valmod = val;

                % Inside parfor: access noise by index, no file I/O
                noise2d = noise_all(:,:,j);

                switch n_mode
                    case 1
                        if QV_mode == 1
                            valmod(:,:,1)=valmod(:,:,1)+noise2d.*bgnoise_amp_valid_rate.*valmod(:,:,1);
                        else
                            valmod(:,:,1)=valmod(:,:,1)+noise2d.*bgnoise_amp_valid_rate.*0.02;
                        end
                        valmod(:,:,1) = valmod(:,:,1) + actatr_mark;
                        valmod(valmod < 0) = 0;
                    case 2
                        valmod=valmod+randn(size(val)).*bgnoise_amp_valid_rate.*valmod;
                        valmod(:,:,1) = valmod(:,:,1) + actatr_mark;
                        valmod(valmod < 0) = 0;
                end

                filename = dir_dataout_woinput + "/" + "wrfinput_d01_woinput_" + n_tag + "_" + j;
                copyfile('./clean/wrfinput_d01_clean', filename);
                ncwrite(filename,varname,valmod);
                disp(['finished preparing noise pattern ', num2str(j)])
            end

        end
    end
end  % end of combo_i loop
