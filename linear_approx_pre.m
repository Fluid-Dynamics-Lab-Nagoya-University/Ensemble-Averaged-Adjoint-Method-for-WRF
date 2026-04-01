% Estimates the sensitivity field via linear approximation (pseudo-inverse).
% Constructs a perturbation matrix D from ensemble noise files, computes
% the precipitation response vector r using a 3x3 Gaussian-weighted kernel,
% and recovers the sensitivity as s = r' * pinv(D).
%
% by Shan Jiang, FDL, Nagoya University

clear
close all
fontsize = 20;
fontname = "times";
ens_member = 100;
bgnoise_amp=0.3;
row = 17;
col = 35;   %position of reducing precipitation
dir = './wrfout/woinput_absG_i0M_NOV0.02_sg_0.3/';
dir_data = './data/';

% ==== preload external noise files ====
nx = 50;
ny = 50;
pert_dir = "./pert";
noise_all = zeros(nx, ny, ens_member);
for j = 1:ens_member
    pert_file = fullfile(pert_dir, sprintf("pert%d.dat", j));
    fid_noise = fopen(pert_file, 'r');
    noise_1d = fscanf(fid_noise, '%f', [nx*ny, 1]);
    fclose(fid_noise);
    noise_all(:,:,j) = reshape(noise_1d, [nx, ny]);
end

% read clean data
filename_clean = "./clean/wrfinput_d01_clean";
RAINNC_clean_terminal       = ncread(filename_clean,"RAINNC");
QVAPORclean                 = ncread(filename_clean,"QVAPOR");

% ==== construct Dmat: pert file * bgnoise_amp * 0.02 ====
Dmat = zeros(2500, 100);
for j=1:ens_member
    noise2d = noise_all(:,:,j);
    Dtmp = noise2d(:) .* bgnoise_amp .* 0.02;
    Dmat(:, j) = Dtmp;
end

figure(15)
noise_last = noise_all(:,:,ens_member) .* bgnoise_amp .* 0.02;
imagesc(flipud(noise_last(:,:)'));
colormap(jet)
colorbar
set(gca, 'FontName', fontname, 'FontSize', fontsize);
title(sprintf('Random noise (pert file)'));
axis equal
axis image

%read noisy data
val_diff = zeros(50,50,ens_member);
for i=1:ens_member
    filename = sprintf('%swrfout_d01_2018-07-05_120000_woinput_AD%d', dir, i);
    XLONG   = ncread(filename,"XLONG"); % longitude
    XLAT    = ncread(filename,"XLAT"); % latitude
    Times   = ncread(filename,"Times"); % time
    val       = ncread(filename,"RAINNC");
    val_diff(:,:,i) = val(:,:,7)-RAINNC_clean_terminal(:,:,7);
end

figure(16)
landmask = ncread(filename, 'LANDMASK');
LM = landmask(:,:,1);
lon_min = min(XLONG(:));
lon_max = max(XLONG(:));
lat_min = min(XLAT(:));
lat_max = max(XLAT(:));
imagesc(flipud(RAINNC_clean_terminal(:,:,7)'));
colormap(jet)
colorbar
set(gca, 'FontName', fontname, 'FontSize', fontsize);
title(sprintf('RAINNC_{clean}'));
axis equal
axis image
hold on;
contour(flipud(LM'), [0.5 0.5], 'k', 'LineWidth', 2);
hold off;

figure(17)
landmask = ncread(filename, 'LANDMASK');
LM = landmask(:,:,1);
imagesc(flipud(QVAPORclean(:,:,1,1)'));
colormap(jet)
colorbar
clim([0.0138, 0.0211]);
set(gca, 'FontName', fontname, 'FontSize', fontsize);
title(sprintf('QVAPOR_{clean}'));
axis equal
axis image
hold on;
contour(flipud(LM'), [0.5 0.5], 'k', 'LineWidth', 2);
hold off;

% Gaussian filter and construct r
G = [0.25 0.5 0.25;
     0.5  1.0 0.5;
     0.25 0.5 0.25];
A = zeros(3,3);
rvec = zeros(100,1);
for i=1:ens_member
    A = val_diff(row-1:row+1, col-1:col+1, i);
    rtmp = sum(sum(A .* G)) / 4;
    rvec(i,1) = rtmp;
end

Dp = pinv(Dmat);
svec = rvec' * Dp;
Smat = -reshape(svec, 50, 50);  % negative sign: gradient descent to reduce precipitation

figure(18)
imagesc(flipud(Smat(:,:)'));
colormap(jet)
colorbar
set(gca, 'FontName', fontname, 'FontSize', fontsize);
title(sprintf('Sensitivity (negated for control)'));
axis equal
axis image

