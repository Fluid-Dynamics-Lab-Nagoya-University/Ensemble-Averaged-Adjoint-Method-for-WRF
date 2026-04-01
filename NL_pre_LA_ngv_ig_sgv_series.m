% Generates perturbed WRF initial conditions for nonlinear runs using the
% linear-approximation sensitivity field (Smat from linear_approx_pre.m).
% Must be run in the same MATLAB session as linear_approx_pre.m (Smat is
% taken from the workspace; clear/close all are intentionally omitted).
%
% by Shan Jiang, FDL, Nagoya University

%clear
%close all
if isempty(gcp('nocreate'))
    parpool('local', 6);
end
input_mode = 1; %1:all 2:positive only 3: negative only
s_mode = 1; %1: ground 2: volume
n_mode = 1; %1: ground 2: volume
valid_mode = 2; %1: ensemble validation (with noise) 2: deterministic validation (no noise)

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
input_eng_rate_list = [0.3];
noise_rate_list = [0.3];

% ==== preload external noise files ====
nx = 50;
ny = 50;
max_ens_member = 100;
pert_dir = "./pert";
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

% ==== read QVAPOR from clean wrfinput ====
basefile = "./clean/wrfinput_d01_clean";
val = ncread(basefile, "QVAPOR");
varname = "QVAPOR";
for rate_idx = 1:length(input_eng_rate_list)
    input_eng_rate = input_eng_rate_list(rate_idx);
    for noise_idx = 1:length(noise_rate_list)
        sens_noise_rate = noise_rate_list(noise_idx);

        % determine ensemble size and noise amplitude based on valid_mode
        switch valid_mode
            case 1
                bgnoise_amp_valid_rate = sens_noise_rate;
                ens_member = 100;
                valid_ng = sens_noise_rate;
            case 2
                bgnoise_amp_valid_rate = 0;
                ens_member = 1;
                valid_ng = 0;
            otherwise
                error('valid_mode must be 1 or 2.');
        end

        s_LA = Smat;  % Smat already negated: gradient descent to reduce precipitation
        withSgm = "";
        dir_dataout_woinput = "LA" + "_" + n_tag + valid_ng ...
            + "_ig" + input_eng_rate + "_" + s_tag + sens_noise_rate + withSgm;
        mkdir(dir_dataout_woinput)
        S_max = max(abs(s_LA(:)));
        s_norm_abs = abs(s_LA) / S_max;
        s_norm = s_norm_abs .* sign(s_LA);
        % input_mode=1 (both): no sign filtering; QV_mode=2: reference QVAPOR = 0.02
        input_eng_g = s_norm .* input_eng_rate .* 0.02;

        % scale to target Frobenius norm
        tmpFroNorm = s_norm .* 0.02;
        target_froNorm = 0.095222;
        current_froNorm = norm(tmpFroNorm, 'fro');
        fprintf("Case: s_tag=%s, n_tag=%s, ig=%g, sens_noise=%g, valid_mode=%d\n", ...
                char(s_tag), char(n_tag), input_eng_rate, sens_noise_rate, valid_mode);
        fprintf("    valid_ng     = %g\n", valid_ng);
        fprintf("    ens_member   = %d\n", ens_member);
        fprintf("    target_norm  = %.6f\n", target_froNorm);
        fprintf("    current_norm = %.6f\n\n", current_froNorm);
        input_eng_g = input_eng_g .* (target_froNorm ./ current_froNorm);
        parfor j = 1:ens_member
            valmod = val;
            noise2d = noise_all(:,:,j);
            switch n_mode
                case 1
                    % QV_mode=2: reference QVAPOR = 0.02
                    valmod(:,:,1) = valmod(:,:,1) + noise2d .* bgnoise_amp_valid_rate .* 0.02;
                    valmod(:,:,1) = valmod(:,:,1) + input_eng_g;
                    valmod(valmod < 0) = 0;
                case 2
                    valmod = valmod + randn(size(val)) .* bgnoise_amp_valid_rate .* valmod;
                    valmod(:,:,1) = valmod(:,:,1) + input_eng_g;
                    valmod(valmod < 0) = 0;
            end
            filename = dir_dataout_woinput + "/" + "wrfinput_d01_woinput_" + n_tag + "_" + j;
            copyfile('./clean/wrfinput_d01_clean', filename);
            ncwrite(filename, varname, valmod);
            disp(['finished preparing noise pattern ', num2str(j)])
        end

        % ==== verify control input Frobenius norm ====
        verify_file = dir_dataout_woinput + "/" + "wrfinput_d01_woinput_" + n_tag + "_1";
        val_verify = ncread(verify_file, varname);
        diff_field = val_verify(:,:,1) - val(:,:,1);
        % subtract background noise to isolate the control input
        diff_field = diff_field - noise_all(:,:,1) .* bgnoise_amp_valid_rate .* 0.02;
        fprintf("    Verification: Frobenius norm of control input = %.6f\n", norm(diff_field, 'fro'));
        fprintf("    Expected: input_eng_rate * target_froNorm = %.6f\n\n", input_eng_rate * target_froNorm);
    end
end
