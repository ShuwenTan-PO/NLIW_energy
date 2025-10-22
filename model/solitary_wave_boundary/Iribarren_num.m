close all
clear all
clc
S = .05;
load soliton_file_summer.mat
Sw_summer = -DJLES_solution.wave_ampl/(DJLES_solution.wavelength/2);
Ir_summer = S/sqrt(-DJLES_solution.wave_ampl/(DJLES_solution.wavelength/2));
a_summer = S/Sw_summer;
load soliton_file.mat
Sw_winter = -DJLES_solution.wave_ampl/(DJLES_solution.wavelength/2);
Ir_winter = S/sqrt(-DJLES_solution.wave_ampl/(DJLES_solution.wavelength/2));
a_winter = S/Sw_winter;

disp(strcat('summer Ir = ',num2str(Ir_summer), ' Sw = ',num2str(Sw_summer), 'a = ', num2str(a_summer)));
disp(strcat('winter Ir = ',num2str(Ir_winter), ' Sw = ',num2str(Sw_winter), 'a = ', num2str(a_winter)));

% summer Ir =0.54987 Sw =0.0082684a =6.0471
% winter Ir =0.37611 Sw =0.017673a =2.8292