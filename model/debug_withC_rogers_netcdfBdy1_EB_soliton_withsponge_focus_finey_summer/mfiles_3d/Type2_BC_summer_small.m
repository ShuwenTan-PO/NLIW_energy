% This code generate a netcdf file that contains type-2 boundary condition 
% the boundary condition contains a solitary wave solved using DJLES solver
% source from: https://github.com/mdunphy/DJLES
% =========================================================================
% vars:
% Nt Nk, Ntype2
% time: [Nt]
% z: [Nk]
% xe, ye, edgep: [Ntype2]
% boundary_u, boundary_v, boundary_w, boundary_T, boundary_S: [Nt Nk, Ntype2]
% =========================================================================
% update by S.Tan, 2023-03-06
% method_double = 1 double wave amplitude
% 
close all 
addpath(genpath('/Users/stan/Desktop/Github/DJLES'))
addpath('../../solitary_wave_boundary')
datadir='../rundata';
verbose=1;
method_double = 0;

%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define APE for wave 
% A  = 1.88e+7; % APE for wave (kg m/s^2)
A  = 1.88e+7; % APE for wave (kg m/s^2)
% Alist  = A.*[.1 .5 1 2 4 8]; 
% define domin width
L = 8e+3;
% Llist = [1e+3 5e+3 10e+3 15e+3 20e+3]; % domain width (m)
filepath = ' ./';
eval(strcat('load ', filepath, 'SUNTANS_grid.mat INPUT Nx Ny Nk'))
H  = -INPUT.z0; % domain depth (m)

%% Solve the DJL equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DJLES_solution = wave_generate_Rogers2022_summer(L, H, A, datadir);
% DJLES_solution = wave_generate_Vitousek2011(L, H, A, datadir);

%% double wave amplitude %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method_double
    den_double = nan(size(DJLES_solution.density));
    for i = 1:length(DJLES_solution.xc)
        zc = DJLES_solution.zc;
        den_ori = DJLES_solution.density(:,i);
        zc_new = zc*2;
        den_double(:,i) = interp1(zc_new,den_ori,zc);
    end
    figure
    plot(den_ori, zc)
    hold on
    plot(den_double(:,i), zc)
    hold on
    plot(den_ori, zc_new, '-.')
    ylim([min(zc), max(zc)])
    DJLES_solution.density = den_double;
end

save soliton_file_summer_smallAPE.mat DJLES_solution

