% Generates perturbed WRF initial conditions for adjoint runs.
% Adds Gaussian noise to QVAPOR across a range of noise rates.
%
% by Shan Jiang, FDL, Nagoya University

clear
close all

ens_member = 100;
noise_rate_list = [0, 0.001, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5];

varname = "QVAPOR";   % amount of water vapor

val = ncread("./clean/wrfinput_d01_clean", "QVAPOR");

for rate = noise_rate_list

    dir_dataout_woinput = ['./prepare/woinput_absG_i0M_NOV0.02_sg_', num2str(rate)];
    mkdir(dir_dataout_woinput)
    copyfile('./prepare/A_RAINNC_absG_i0M.dat', [dir_dataout_woinput, '/A_RAINNC_absG_i0M.dat']);

    for j = 1:ens_member
        valmod = val;
        rng(j)
        valmod(:,:,1,1) = valmod(:,:,1,1) + randn(size(val(:,:,1,1))) .* rate .* 0.02;
        valmod(valmod < 0) = 0;

        filename = [dir_dataout_woinput, '/wrfinput_d01_woinput_', num2str(j)];
        copyfile('./clean/wrfinput_d01_clean', filename);
        ncwrite(filename, varname, valmod);
        disp(['finished preparing noise pattern ', num2str(j)])
    end

end

figure(12)
imagesc(flipud(valmod(:,:,1)'));
colormap(jet);
colorbar;
