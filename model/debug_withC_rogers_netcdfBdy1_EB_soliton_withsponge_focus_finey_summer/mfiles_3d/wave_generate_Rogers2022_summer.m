function DJLES_solution = wave_generate_Rogers2022_summer(L, H, A, datadir)
% ======= Solve DJL equation using the DJLES by Dunphy et al. (2011) ======
% source from: https://github.com/mdunphy/DJLES
% the density profile is specified according to Rogers et al. (2022)
% background flow is set to zero for now
% A is the available potential energy [kg m s-2]
% L is the domain width
% H is the domain depth
%
% created by S.Tan 2023/01/10
% 1) evenly spaced x and z for now, might worth to adjust to the stretched 
% SUNTANS grid
% 2) zero background flow for now, when adding background tides might be
% good to include such
% S.Tan 2023/03/20
% !!!!!!fixed a bug: rho_b = -gamma*rho_0*T instead of -gamma*T+rho0!!!!!!!
% TO UPDATE
% 1) no meaningful djle solutions from the current density profile in
% Rogers 2022, need careful revisit
%
%%% rho and N2 profiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho0 = 1000;
% gamma taken from suntans.dat
gamma = getvalue([datadir,'/suntans.dat'],'gamma');

% % Vitousek 2011 density profile
% a1 = .5*.001;%./2.1e-4; 
% a2 = 2/tanh(.99)/200; a3 = 250; %a4 = .001;
% frho=@(z,rho0,a1,a2,a3) -a1*tanh(a2*(z+a3))+rho0;%-a4*z;
% frhoz=@(z,a1,a2,a3) -a1*a2*(1-tanh(a2*(z+a3)).^2);%+a4;
% rho  = @(z) frho(z,rho0,a1,a2,a3);
% rhoz = @(z) frhoz(z,a1,a2,a3);

% Rogers 2022 density profile
% parameters taken from initialization.c
a1 = 23.36; a2 = 3.13; a3 = -44.12; a4 = 293.12;
% idealized rho profile taken from state.c
frho=@(z,gamma,rho0,a1,a2,a3,a4) -gamma*rho0*(a1*exp(-(-z+a3)/a4)+a2)+rho0;
frhoz=@(z,gamma,rho0,a1,a3,a4) -gamma*rho0*a1*exp(-(-z+a3)/a4)/a4;
% below is wrong, but somehow it gives DJLE solutions
% frho=@(z,gamma,rho0,a1,a2,a3,a4) -gamma*rho0*(a1*exp(-(-z+a3)/a4)+a2)+rho0;
% frhoz=@(z,gamma,rho0,a1,a3,a4) -gamma*rho0*a1*exp(-(-z+a3)/a4)/a4;
rho  = @(z) frho(z, gamma,rho0,a1,a2,a3,a4);
rhoz = @(z) frhoz(z, gamma,rho0,a1,a3,a4);

%%% background velocity profiles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The velocity profile (zero for this case) (m/s)
Ubg=@(z) 0*z; Ubgz=@(z) 0*z; Ubgzz=@(z) 0*z;

%%% Find the solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
start_time = clock;
% Specify resolution according to this schedule and run iteratively
% NXlist=[  1 2 4 8 16].*Ny;
% NZlist=[  1 2 4 8 16].*Nk;
NXlist=[  64   128    256   512   1024   1024   1024];
NZlist=[  32    64    128   256    512   1024   1024*2];
for Nindex=1:length(NXlist)
    % Resolution for this wave
    NX = NXlist(Nindex);
    NZ = NZlist(Nindex);
    if Nindex==length(NXlist)
        epsilon=1e-5;
    end
    % Iterate the DJL solution
    djles_refine_solution
    djles_diagnostics; djles_plot; % uncomment to view progress at each step
end

end_time=clock;
fprintf('Total wall clock time: %f seconds\n',etime(end_time, start_time));

%%% make plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
subplot(221)
contourf(xc, zc, density)
title('density')
colorbar
subplot(222)
contourf(xc, zc, eta)
title('eta')
colorbar
subplot(223)
contourf(xc, zc, u)
title('u')
colorbar
subplot(224)
contourf(xc, zc, w)
title('w')
colorbar
print -djpeg -r300 DJLE_solution
% 
figure
subplot(121)
plot(rho(z),z)
subplot(122)
plot(rhoz(z),z)
print -djpeg -r300 DJLE_solution_rho

ftemp=@(z,a1,a2,a3,a4) a1*exp(-(-z+a3)/a4)+a2;
temperature  = @(z) ftemp(z, a1,a2,a3,a4);
% ftemp=@(z,gamma,a1,a2,a3) a1*tanh(a2*(z+a3))/gamma + 20;
% temperature  = @(z) ftemp(z,gamma, a1,a2,a3);
figure
plot(temperature(z),z)
print -djpeg -r300 DJLE_solution_temp

figure
subplot(121)
plot(rho(zc), zc)
title('rho')
subplot(122)
plot(Ubg(zc), zc)
title('Ubg')
%%% save results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DJLES_solution.A=A; DJLES_solution.c=c; 
DJLES_solution.wave_ampl=wave_ampl; DJLES_solution.wavelength=wavelength;

DJLES_solution.L=L; DJLES_solution.H=H;
DJLES_solution.ze=ze; DJLES_solution.zc=zc;
DJLES_solution.xe=xe; DJLES_solution.xc=xc;
DJLES_solution.density=density; DJLES_solution.eta=eta;
DJLES_solution.u=u; DJLES_solution.w=w;

end