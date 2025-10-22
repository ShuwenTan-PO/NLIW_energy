%% compute c from DJLEs using summer stratification for different APE
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
clear all
close all 
addpath(genpath('/Users/stan/Desktop/Github/DJLES'))
datadir='../debug_withC_rogers_netcdfBdy1_EB_soliton_withsponge_focus_finey_summer/rundata';

%% Specify the parameters of the problem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define APE for wave 
% A  = 1.88e+7; % APE for wave (kg m/s^2)
% A  = 13e+7; % APE for wave (kg m/s^2)
Alist  = [1.88e+7, 13e+7]; 
% define domin width
L = 8e+3;
% Llist = [1e+3 5e+3 10e+3 15e+3 20e+3]; % domain width (m)
Hlist  = [100:1:600]; % domain depth (m)
c = nan(length(Alist), length(Hlist));
A_w = nan(length(Alist), length(Hlist));
L_w = nan(length(Alist), length(Hlist));
%% Solve the DJL equation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Alist)
    for j=1:length(Hlist)
        DJLES_solution = wave_generate_Rogers2022_summer(L, Hlist(j), Alist(i), datadir);
        c(i,j) = DJLES_solution.c;
        A_w(i,j) = DJLES_solution.wave_ampl;
        L_w(i,j) = DJLES_solution.wavelength;
    end
end

save c_vs_APE_summer.mat Alist Hlist c A_w L_w

% i=1;j=310;