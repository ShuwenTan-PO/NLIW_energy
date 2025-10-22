function [] = SUNTANS_reduce(tsave)
disp('loading data to reduce')
% tsave = 7.5*12.42*3600;

load('./SUNTANS_results.mat',...
't','D','Nkmax','z','dz','rho','N2','rho0','xv','yv','xe','ye',...
'eta','z_pyn','Ubar','Utilde','rho_prime','u_tildeprime','u_barprime',...
'Vbar','Vtilde','v_tildeprime','v_barprime',...
'x_sponge','x_wavemaker','depth','Nout_keep',...
'Fx_0bar','Fx_primebar','Fx_0','Fx_prime','Fy_0bar','Fy_primebar',...
'T','S','bottom_cells',...
'u','v','etabar','indx_pycnocline','Tavg','w','wbar');

AVG = load('./SUNTANS_results_average.mat');

%% save only nearest time step

[~,indx] = min(abs(tsave-t));
tsave = t(indx);

N2_init = N2(:,:,1);
rho_init = rho(:,:,1);

eta = eta(:,:,indx);
N2 = N2(:,:,indx);
rho = rho(:,:,indx);
rho_prime = rho_prime(:,:,indx);
u = u(:,:,indx);
v = v(:,:,indx);
w = w(:,:,indx);
% wbar = wbar(:,:,indx);
u_tildeprime = u_tildeprime(:,:,indx);
Utilde = Utilde(:,indx);
v_tildeprime = v_tildeprime(:,:,indx);
Vtilde = Vtilde(:,indx);
z_pyn = z_pyn(:,indx);
S = S(:,:,indx);
for i=1:length(xv)
   Tbot(i,:) = squeeze(T(i,bottom_cells(i),:));
end
T = T(:,:,indx);
indx = AVG.t>(AVG.t(end)-Tavg);

rho_b = squeeze(mean(AVG.rho_b(:,:,indx),3,'omitnan'));
rho_prime_bar = squeeze(mean(AVG.rho_prime(:,:,indx),3,'omitnan'));
etabar = squeeze(mean(AVG.eta(:,:,indx),3,'omitnan'));
Tbar = squeeze(mean(AVG.T(:,:,indx),3,'omitnan'));
Sbar = squeeze(mean(AVG.S(:,:,indx),3,'omitnan'));
% nuTbar = squeeze(nanmean(AVG.nuT(:,:,indx),3));
% kappaTbar = squeeze(nanmean(AVG.kappaT(:,:,indx),3));

clear AVG
disp('saving reduced data...')
save ./SUNTANS_results_reduce
clear;
disp('done!')
