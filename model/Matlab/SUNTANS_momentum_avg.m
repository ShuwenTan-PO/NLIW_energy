% clear
close all
set(0,'defaulttextinterpreter','latex')
datadir='./rundata'; % data directory
% do_momentum_xy = 1;
% do_momentum_cyl = 1;
% % SUNTANS_momentum_avg
%%
disp('load data and process grid...')

plt = {'-','--','-.','-','--'};

load('SUNTANS_results_average.mat')
GRID = load('SUNTANS_grid.mat');

Nx = GRID.Nx;
Ny = GRID.Ny;
Nt = length(t);
Nk = length(GRID.z);
ylims = (GRID.INPUT.y0 +GRID.INPUT.r*[-1 1]);
xlims = (GRID.INPUT.x0 +GRID.INPUT.r*[-1 2]);
zlims = [GRID.INPUT.z0/3 0];
xv = reshape(GRID.xv,Nx,Ny);
yv = reshape(GRID.yv,Nx,Ny);
z = GRID.z;
dz = GRID.dz;
depth = reshape(GRID.Depth,Nx,Ny);

% get some relevant parameters from the input file
f =  getvalue([datadir,'/suntans.dat'],'Coriolis_f');
period =  getvalue([datadir,'/suntans.dat'],'TM2');
nu =  getvalue([datadir,'/suntans.dat'],'nu');
nu_H =  getvalue([datadir,'/suntans.dat'],'nu_H');
Cd =  getvalue([datadir,'/suntans.dat'],'CdB');
% low frequency parameters
lowfreq_nudging =  getvalue([datadir,'/suntans.dat'],'lowfreq_nudging');
TauL =  3600*getvalue([datadir,'/suntans.dat'],'TauL');
ULm =  getvalue([datadir,'/suntans.dat'],'ULm');
VLm =  getvalue([datadir,'/suntans.dat'],'VLm');
grav =  getvalue([datadir,'/suntans.dat'],'grav');
% grav = 9.81;

% back out lat and b factor
% f = 2*omega*sin(lat)
% b = 2*omega*cos(lat)
omega = 7.27e-5; %rad/s
lat = asind(f/2/omega);
b = 2*omega*cosd(lat);

% get important variables, all variables here are time avg so drop bar notation
disp('reshape variables...')

eta = reshape(eta,Nx,Ny);

U = reshape(U,Nx,Ny);
V = reshape(V,Nx,Ny);
u_prime = reshape(u_prime,Nx,Ny,Nk);
v_prime = reshape(v_prime,Nx,Ny,Nk);
u = reshape(u,Nx,Ny,Nk);
v = reshape(v,Nx,Ny,Nk);
w = reshape(w,Nx,Ny,Nk+1);
w_center = w;
w_center(isnan(w_center))=0; %do this so w=0 at boundary
w_center = 0.5*(w_center(:,:,2:end) + w_center(:,:,1:end-1));
w_center = w_center + 0*u; %now nans below depth

S = reshape(S,Nx,Ny,Nk);
T = reshape(T,Nx,Ny,Nk);
rho = reshape(rho,Nx,Ny,Nk);
rho_b = reshape(rho_b,Nx,Ny,Nk);
rho_prime = reshape(rho_prime,Nx,Ny,Nk);
p0 = reshape(p0,Nx,Ny,Nk);
p_b = reshape(p_b,Nx,Ny,Nk);
p_prime = reshape(p_prime,Nx,Ny,Nk);
nuT = reshape(nuT,Nx,Ny,Nk);
kappaT = reshape(kappaT,Nx,Ny,Nk);

empties = reshape(empties,Nx,Ny,Nk);

if exist('upup','var') %reynolds stresses
    upup = reshape(upup,Nx,Ny,Nk);
    vpvp = reshape(vpvp,Nx,Ny,Nk);
    wpwp = reshape(wpwp,Nx,Ny,Nk);
    upvp = reshape(upvp,Nx,Ny,Nk);
    upwp = reshape(upwp,Nx,Ny,Nk);
    vpwp = reshape(vpwp,Nx,Ny,Nk);
else
    upup = u + nan;
    vpvp = u + nan;
    wpwp = u + nan;
    upvp = u + nan;
    upwp = u + nan;
    vpwp = u + nan;
end
if exist('p_nhyd','var')
    p_nhyd = reshape(p_nhyd,Nx,Ny,Nk);
else
    p_nhyd = u + nan;
end
% compute pressure from reference density
p_hyd_rho0 = 0*p0;
for k=1:Nk
   p_hyd_rho0(:,:,k) = grav * rho0 * eta; 
end
% add nans to below depth
p_hyd_rho0 = p_hyd_rho0 + 0*u;
% compute pressure from background/deviation density
p_hyd_rhob = p_b + p_prime;

% create 3D coordinates for vectorized compute
disp('create vectorized coordinates...')

X = 0*u;
Y = 0*u;
Z = 0*u;
for i=1:length(z)
    X(:,:,i) = xv;
    Y(:,:,i) = yv;
    Z(:,:,i) = z(i);
end

% cut cells
disp('create cut cell dz')
zb = z-dz/2;
zt = z+dz/2;
dZ = nan+Z;

DEPTH = 0*u;
ZB = 0*u;
ZT = 0*u;
for i=1:length(z)
    ZB(:,:,i) = zb(i);
    ZT(:,:,i) = zt(i);
    DEPTH(:,:,i) = depth;
end
indx = ZB>=-DEPTH;
dZ(indx)= ZT(indx) -ZB(indx);
indx = ZB < -DEPTH & ZT > -DEPTH;
dZ(indx)= ZT(indx) + DEPTH(indx);

% clear some variables we don't need to keep memory low
clear p0 pfs p_b p_prime rho_b rho_prime kappaT nuT w
clear DEPTH ZB ZT zb zt indx
%% cartesian momentum
if do_momentum_xy
    
    % momentum in X direction, put everything on LHS
    % du/dt + duu/dx + dvu/dy + dwu/dz  (Du/Dt = US + NL1 + NL2 + NL3)
    % + dupup/dx + dupvp/dy + dupwp/dz  (Reynolds = Re1 + Re2 + Re3)
    % -fv + bw (coriolis = Cor1 + Cor2)
    % +1/rho0*dp_hyd/dx + 1/rho0*dp_nhyd/dx (pressure = PGh + PGnh)
    % -nu_H*d2u/dx2 -nu_H*d2u/dy2 - nu*d2u/dz2 (visc = Diss1 + Diss2 + Diss3)
    % -Fl (low freq nudging)
    % = 0
    % derivatives, use central diff where possible (in reality suntans uses
    % upwind)
    % dir 1 = x, dir 2 = y, dir 3 = z
    disp('compute average x mometum terms...')

    % +duu/dx
    Mx_Nl1 = u +nan;
    Mx_Nl1(2:end-1,:,:) = (u(3:end,:,:).^2 - u(1:end-2,:,:).^2)...
                          ./ (X(3:end,:,:) - X(1:end-2,:,:));
    % +dvu/dy
    Mx_Nl2 = u + nan;
    Mx_Nl2(:,2:end-1,:) = (v(:,3:end,:).*u(:,3:end,:) - v(:,1:end-2,:).*u(:,1:end-2,:))...
                          ./ (Y(:,3:end,:) - Y(:,1:end-2,:));
    % +dwu/dz, only problem with this approach is top/bottom are nan,
    % but will be consistent with Reynolds stress calc
    Mx_Nl3 = u + nan;
    Mx_Nl3(:,:,2:end-1) = (w_center(:,:,3:end).*u(:,:,3:end) - w_center(:,:,1:end-2).*u(:,:,1:end-2))...
                          ./ (Z(:,:,3:end) - Z(:,:,1:end-2));
    % +dupup/dx
    Mx_Re1 = u + nan;
    Mx_Re1(2:end-1,:,:) = (upup(3:end,:,:) - upup(1:end-2,:,:))...
                          ./ (X(3:end,:,:) - X(1:end-2,:,:));
    % +dupvp/dy 
    Mx_Re2 = u + nan;
    Mx_Re2(:,2:end-1,:) = (upvp(:,3:end,:) - upvp(:,1:end-2,:))...
                          ./ (Y(:,3:end,:) - Y(:,1:end-2,:));
    % +dupwp/dz                  
    Mx_Re3 = u + nan;
    Mx_Re3(:,:,2:end-1) = (upwp(:,:,3:end) - upwp(:,:,1:end-2))...
                          ./ (Z(:,:,3:end) - Z(:,:,1:end-2));
    % -fv + bw                 
    Mx_Cor1 = -f * v;                  
    % + bw 
    Mx_Cor2 = b * w_center;

    % +1/rho0*dp_hyd_rho0/dx 
    Mx_PGh_rho0 = u + nan;
    Mx_PGh_rho0(2:end-1,:,:) = (1/rho0) * (p_hyd_rho0(3:end,:,:) - p_hyd_rho0(1:end-2,:,:))...
                          ./ (X(3:end,:,:) - X(1:end-2,:,:));
    % +1/rho0*dp_hyd_rhob/dx 
    Mx_PGh_rhob = u + nan;
    Mx_PGh_rhob(2:end-1,:,:) = (1/rho0) * (p_hyd_rhob(3:end,:,:) - p_hyd_rhob(1:end-2,:,:))...
                          ./ (X(3:end,:,:) - X(1:end-2,:,:));
    % +1/rho0*dp_nhyd/dx
    Mx_PGnh = u + nan;
    Mx_PGnh(2:end-1,:,:) = (1/rho0) * (p_nhyd(3:end,:,:) - p_nhyd(1:end-2,:,:))...
                          ./ (X(3:end,:,:) - X(1:end-2,:,:));
    % -nu_H*d2u/dx2
    Mx_Diss1 = u + nan;
    dxtemp = diff(0.5*(X(1:end-1,:,:)+X(2:end,:,:)), 1, 1);
    Mx_Diss1(2:end-1,:,:) = -nu_H * diff(u,2,1) ./ dxtemp.^2;

    % -nu_H*d2u/dy2
    Mx_Diss2 = u + nan;
    dxtemp = diff(0.5*(Y(:,1:end-1,:)+Y(:,2:end,:)), 1, 2);
    Mx_Diss2(:,2:end-1,:) = -nu_H * diff(u,2,2) ./ dxtemp.^2;
    Mx_DissH = Mx_Diss1 + Mx_Diss2;
    clear Mx_Diss1 Mx_Diss2

    % - nu*d2u/dz2
    Mx_DissV = u + nan;
    dxtemp = diff(0.5*(Z(:,:,1:end-1)+Z(:,:,2:end)), 1, 3);
    Mx_DissV(:,:,2:end-1) = -nu * diff(u,2,3) ./ dxtemp.^2;
    clear dxtemp
    
    % -Fl = (u - ULM)/TauL
    Mx_Fl = lowfreq_nudging*(u - ULm)/TauL;
    
    % take depth avg
    disp('take depth average of x momentum...')
    Mx_Nl1_2d = get_depth_avg(Mx_Nl1,dZ,3);
    Mx_Nl2_2d = get_depth_avg(Mx_Nl2,dZ,3);
    Mx_Nl3_2d = get_depth_avg(Mx_Nl3,dZ,3);
    Mx_Re1_2d = get_depth_avg(Mx_Re1,dZ,3);
    Mx_Re2_2d = get_depth_avg(Mx_Re2,dZ,3);
    Mx_Re3_2d = get_depth_avg(Mx_Re3,dZ,3);
    Mx_Cor1_2d = get_depth_avg(Mx_Cor1,dZ,3);
    Mx_Cor2_2d = get_depth_avg(Mx_Cor2,dZ,3);
    Mx_PGh_rho0_2d = get_depth_avg(Mx_PGh_rho0,dZ,3);
    Mx_PGh_rhob_2d = get_depth_avg(Mx_PGh_rhob,dZ,3);
    Mx_PGnh_2d = get_depth_avg(Mx_PGnh,dZ,3);
    Mx_DissH_2d = get_depth_avg(Mx_DissH,dZ,3);
    Mx_DissV_2d = get_depth_avg(Mx_DissV,dZ,3);
    Mx_Fl_2d = get_depth_avg(Mx_Fl,dZ,3);
    
    disp('saving avg x momentum...')
    save('SUNTANS_results_momentum_avg_x.mat','-v7.3','Mx*','b','f','Cd','nu','nu_H','lat','xv','yv','z','mtime')
    clear Mx* %keep memory use low
    
    %% momentum in Y direction, put everything on LHS
    % dv/dt + dvv/dy + dvu/dx + dwv/dz  (Dv/Dt = US + NL1 + NL2 + NL3)
    % + dvpvp/dy + dupvp/dx + dvpwp/dz  (Reynolds = Re1 + Re2 + Re3)
    % +fu (coriolis = Cor1 )
    % +1/rho0*dp_hyd/dy + 1/rho0*dp_nhyd/dy (pressure = PGh + PGnh)
    % -nu_H*d2v/dx2 -nu_H*d2v/dy2 - nu*d2v/dz2 (visc = Diss1 + Diss2 + Diss3)
    % - Fl (low frequency nudging)
    % = 0
    % derivatives, use central diff where possible (in reality suntans uses
    % upwind)
    % dir 1 = x, dir 2 = y, dir 3 = z
    disp('compute average y mometum terms...')

    % +dvv/dy
    My_Nl1 = u +nan;
    My_Nl1(:,2:end-1,:) = (v(:,3:end,:).^2 - v(:,1:end-2,:).^2)...
                          ./ (Y(:,3:end,:) - Y(:,1:end-2,:));
    % +dvu/dx
    My_Nl2 = u +nan;
    My_Nl2(2:end-1,:,:) = (v(3:end,:,:).*u(3:end,:,:) - v(1:end-2,:,:).*u(1:end-2,:,:))...
                          ./ (X(3:end,:,:) - X(1:end-2,:,:));           
    % +dwv/dz, only problem with this approach is top/bottom are nan,
    % but will be consistent with Reynolds stress calc
    My_Nl3 = u + nan;
    My_Nl3(:,:,2:end-1) = (w_center(:,:,3:end).*v(:,:,3:end) - w_center(:,:,1:end-2).*v(:,:,1:end-2))...
                          ./ (Z(:,:,3:end) - Z(:,:,1:end-2));
    % +dupvp/dy 
    My_Re1 = u + nan;
    My_Re1(:,2:end-1,:) = (vpvp(:,3:end,:) - vpvp(:,1:end-2,:))...
                          ./ (Y(:,3:end,:) - Y(:,1:end-2,:));   
    % +dupvp/dx
    My_Re2 = u + nan;
    My_Re2(2:end-1,:,:) = (upvp(3:end,:,:) - upvp(1:end-2,:,:))...
                          ./ (X(3:end,:,:) - X(1:end-2,:,:));
    % +dvpwp/dz                  
    My_Re3 = u + nan;
    My_Re3(:,:,2:end-1) = (vpwp(:,:,3:end) - vpwp(:,:,1:end-2))...
                          ./ (Z(:,:,3:end) - Z(:,:,1:end-2));
    % +fu                  
    My_Cor = f * u;
    % there is no second Cor term, but for consistency with x, add dummy
    % variable
%     My_Cor2 = 0*u;
    % +1/rho0*dp_hyd_rho0/dy 
    My_PGh_rho0 = u + nan;
    My_PGh_rho0(:,2:end-1,:) = (1/rho0) * (p_hyd_rho0(:,3:end,:) - p_hyd_rho0(:,1:end-2,:))...
                          ./ (Y(:,3:end,:) - Y(:,1:end-2,:));
    % +1/rho0*dp_hyd_rhob/dy 
    My_PGh_rhob = u + nan;
    My_PGh_rhob(:,2:end-1,:) = (1/rho0) * (p_hyd_rhob(:,3:end,:) - p_hyd_rhob(:,1:end-2,:))...
                          ./ (Y(:,3:end,:) - Y(:,1:end-2,:));
    % +1/rho0*dp_nhyd/dy
    My_PGnh = u + nan;
    My_PGnh(:,2:end-1,:) = (1/rho0) * (p_nhyd(:,3:end,:) - p_nhyd(:,1:end-2,:))...
                          ./ (Y(:,3:end,:) - Y(:,1:end-2,:));
    % -nu_H*d2v/dx2
    My_Diss1 = u + nan;
    dxtemp = diff(0.5*(X(1:end-1,:,:)+X(2:end,:,:)), 1, 1);
    My_Diss1(2:end-1,:,:) = -nu_H * diff(v,2,1) ./ dxtemp.^2;

    % -nu_H*d2v/dy2
    My_Diss2 = u + nan;
    dxtemp = diff(0.5*(Y(:,1:end-1,:)+Y(:,2:end,:)), 1, 2);
    My_Diss2(:,2:end-1,:) = -nu_H * diff(v,2,2) ./ dxtemp.^2;
    My_DissH = My_Diss1 + My_Diss2;
    clear My_Diss1 My_Diss2

    % - nu*d2v/dz2
    My_DissV = u + nan;
    dxtemp = diff(0.5*(Z(:,:,1:end-1)+Z(:,:,2:end)), 1, 3);
    My_DissV(:,:,2:end-1) = -nu * diff(v,2,3) ./ dxtemp.^2;
    clear dxtemp
    
    % -Fl = (v - VLM)/TauL
    My_Fl = lowfreq_nudging*(v - VLm)/TauL;
    
    % take depth avg
    disp('take depth average of y momentum...')
    My_Nl1_2d = get_depth_avg(My_Nl1,dZ,3);
    My_Nl2_2d = get_depth_avg(My_Nl2,dZ,3);
    My_Nl3_2d = get_depth_avg(My_Nl3,dZ,3);
    My_Re1_2d = get_depth_avg(My_Re1,dZ,3);
    My_Re2_2d = get_depth_avg(My_Re2,dZ,3);
    My_Re3_2d = get_depth_avg(My_Re3,dZ,3);
    My_Cor_2d = get_depth_avg(My_Cor,dZ,3);
%     My_Cor2_2d = get_depth_avg(My_Cor2,dZ,3);
    My_PGh_rho0_2d = get_depth_avg(My_PGh_rho0,dZ,3);
    My_PGh_rhob_2d = get_depth_avg(My_PGh_rhob,dZ,3);
    My_PGnh_2d = get_depth_avg(My_PGnh,dZ,3);
    My_DissH_2d = get_depth_avg(My_DissH,dZ,3);
    My_DissV_2d = get_depth_avg(My_DissV,dZ,3);
    My_Fl_2d = get_depth_avg(My_Fl,dZ,3);
    
    disp('saving avg y momentum...')
    save('SUNTANS_results_momentum_avg_y.mat','-v7.3','My*','b','f','Cd','nu','nu_H','lat','xv','yv','z','mtime')
    clear My* % keep memory use low
    %% momentum in Z direction, put everything on LHS
    % dw/dt + duw/dx + dvw/dy + dww/dz  (Dw/Dt = US + NL1 + NL2 + NL3)
    % + dupwp/dx + dvpwp/dy + dwpwp/dz  (Reynolds = Re1 + Re2 + Re3)
    % -bu (coriolis = Cor1)
    % +1/rho0*dp_nhyd/dz (pressure = PGnh)
    % -nu_H*d2w/dx2 -nu_H*d2w/dy2 - nu*d2w/dz2 (visc = Diss1 + Diss2 + Diss3)
    % = 0
    % derivatives, use central diff where possible (in reality suntans uses
    % upwind)
    % dir 1 = x, dir 2 = y, dir 3 = z
    disp('compute average z mometum terms...')

    % +dwu/dx
    Mz_Nl1 = u +nan;
    Mz_Nl1(2:end-1,:,:) = (w_center(3:end,:,:).*u(3:end,:,:) - w_center(1:end-2,:,:).*u(1:end-2,:,:))...
                          ./ (X(3:end,:,:) - X(1:end-2,:,:));   
    % +dvw/dy
    Mz_Nl2 = u + nan;
    Mz_Nl2(:,2:end-1,:) = (v(:,3:end,:).*w_center(:,3:end,:) - v(:,1:end-2,:).*w_center(:,1:end-2,:))...
                          ./ (Y(:,3:end,:) - Y(:,1:end-2,:));
    % +dww/dz, only problem with this approach is top/bottom are nan,
    % but will be consistent with Reynolds stress calc
    Mz_Nl3 = u + nan;
    Mz_Nl3(:,:,2:end-1) = (w_center(:,:,3:end).^2 - w_center(:,:,1:end-2).^2)...
                          ./ (Z(:,:,3:end) - Z(:,:,1:end-2)); 
                      
    % ** need to add RE terms here **
    % +dupwp/dx
    Mz_Re1 = u + nan;
    Mz_Re1(2:end-1,:,:) = (upwp(3:end,:,:) - upwp(1:end-2,:,:))...
                          ./ (X(3:end,:,:) - X(1:end-2,:,:));
    % +dvpwp/dy 
    Mz_Re2 = u + nan;
    Mz_Re2(:,2:end-1,:) = (vpwp(:,3:end,:) - vpwp(:,1:end-2,:))...
                          ./ (Y(:,3:end,:) - Y(:,1:end-2,:));
    % +dwpwp/dz                  
    Mz_Re3 = u + nan;
    Mz_Re3(:,:,2:end-1) = (wpwp(:,:,3:end) - wpwp(:,:,1:end-2))...
                          ./ (Z(:,:,3:end) - Z(:,:,1:end-2));
    % -bu                  
    Mz_Cor = -b * u;
    % +1/rho0*dp_nhyd/dz 
    Mz_PGnh = u + nan;
    Mz_PGnh(:,:,2:end-1) = (1/rho0) * (p_nhyd(:,:,3:end) - p_nhyd(:,:,1:end-2))...
                          ./ (Z(:,:,3:end) - Z(:,:,1:end-2));
    % -nu_H*d2w/dx2
    Mz_Diss1 = u + nan;
    dxtemp = diff(0.5*(X(1:end-1,:,:)+X(2:end,:,:)), 1, 1);
    Mz_Diss1(2:end-1,:,:) = -nu_H * diff(w_center,2,1) ./ dxtemp.^2;

    % -nu_H*d2w/dy2
    Mz_Diss2 = u + nan;
    dxtemp = diff(0.5*(Y(:,1:end-1,:)+Y(:,2:end,:)), 1, 2);
    Mz_Diss2(:,2:end-1,:) = -nu_H * diff(w_center,2,2) ./ dxtemp.^2;
    Mz_DissH = Mz_Diss1 + Mz_Diss2;
    clear Mz_Diss1 Mz_Diss2

    % -nu*d2w/dz2
    Mz_DissV = u + nan;
    dxtemp = diff(0.5*(Z(:,:,1:end-1)+Z(:,:,2:end)), 1, 3);
    Mz_DissV(:,:,2:end-1) = -nu * diff(w_center,2,3) ./ dxtemp.^2;
    clear dxtemp
    
    % take depth avg
    disp('take depth average of z momentum...')
    Mz_Nl1_2d = get_depth_avg(Mz_Nl1,dZ,3);
    Mz_Nl2_2d = get_depth_avg(Mz_Nl2,dZ,3);
    Mz_Nl3_2d = get_depth_avg(Mz_Nl3,dZ,3);
    Mz_Re1_2d = get_depth_avg(Mz_Re1,dZ,3);
    Mz_Re2_2d = get_depth_avg(Mz_Re2,dZ,3);
    Mz_Re3_2d = get_depth_avg(Mz_Re3,dZ,3);
    Mz_Cor_2d = get_depth_avg(Mz_Cor,dZ,3);
    Mz_PGnh_2d = get_depth_avg(Mz_PGnh,dZ,3);
    Mz_DissH_2d = get_depth_avg(Mz_DissH,dZ,3);
%     Mz_Diss2_2d = get_depth_avg(Mz_Diss2,dZ,3);
    Mz_DissV_2d = get_depth_avg(Mz_DissV,dZ,3);

    disp('saving avg z momentum...')
    save('SUNTANS_results_momentum_avg_z.mat','-v7.3','Mz*','b','f','Cd','nu','nu_H','lat','xv','yv','z','mtime')
    clear Mz* % keep memory use low
end
%% remove some varibles
% clear up* vp* wp* rho* p* nuT kappaT T S

%%
if do_momentum_cyl
% transform coordinates
    disp('create cylindrical coordinates...')

    r_x = sqrt((X - GRID.INPUT.x0).^2+(Y- GRID.INPUT.y0).^2);
    theta_x = atan2(Y- GRID.INPUT.y0,X- GRID.INPUT.x0);
    
    % create new coordinate system in r,theta space
    % grid is ~90m resolution at edge
    dr = 4.5*GRID.INPUT.dxmin;
    % dr/radius = tan(dtheta); keep grid square at edge of slope
    dtheta = atan2(dr,GRID.INPUT.width_top/2);
    % don't start at r=0 to avoid 1/r blowup at origin
    r = dr:dr:GRID.INPUT.width_base*0.65;
    theta = -pi:dtheta:pi;
    [theta,r] = meshgrid(theta,r);
    [Nr,Ntheta] = size(r);

    xr = r.*cos(theta) + GRID.INPUT.x0;
    yr = r.*sin(theta) + GRID.INPUT.y0;
%     X_r = zeros(Nr,Ntheta,Nk);
%     Y_r = zeros(Nr,Ntheta,Nk);
    Z_r = zeros(Nr,Ntheta,Nk);
    R = zeros(Nr,Ntheta,Nk);
    THETA = zeros(Nr,Ntheta,Nk);
    for i=1:length(z)
%         X_r(:,:,i) = xv_r;
%         Y_r(:,:,i) = yv_r;
        Z_r(:,:,i) = z(i);
        R(:,:,i) = r;
        THETA(:,:,i) = theta;
    end
 
    % compute reynolds stress in cyl coord
    disp('compute variables in cylindrical coordinates')
    dZ_r = interp_grid(xv,yv,dZ,xr,yr);
    depth = interp_grid(xv,yv,depth,xr,yr);
    % velocities, rotate, interp to new grid
    [u_r,u_theta] = cart_to_cyl(u,v,theta_x);
    % mean forcing velocities
    [ULm_r,ULm_theta] = cart_to_cyl(ULm + 0*u,VLm + 0*v,theta_x);
    u_r = interp_grid(xv,yv,u_r,xr,yr);
    u_theta = interp_grid(xv,yv,u_theta,xr,yr);
    ULm_r = interp_grid(xv,yv,ULm_r,xr,yr);
    ULm_theta = interp_grid(xv,yv,ULm_theta,xr,yr);
    w_center = interp_grid(xv,yv,w_center,xr,yr);
    U_r = get_depth_avg(u_r,dZ_r,3);
    U_theta = get_depth_avg(u_theta,dZ_r,3);
    
    
    
    
    % reynolds stresses, rotate, interp to new grid
    urpurp = upup.*cos(theta_x).^2 + 2*upvp.*cos(theta_x).*sin(theta_x) + vpvp.*sin(theta_x).^2;
    urpurp = interp_grid(xv,yv,urpurp,xr,yr);
    
    utputp = upup.*sin(theta_x).^2 - 2*upvp.*cos(theta_x).*sin(theta_x) + vpvp.*cos(theta_x).^2;
    utputp = interp_grid(xv,yv,utputp,xr,yr);
    
    urputp = -upup.*cos(theta_x).*sin(theta_x) ...
        + upvp.*(cos(theta_x).^2 - sin(theta_x).^2)+...
        vpvp.*cos(theta_x).*sin(theta_x);
    urputp = interp_grid(xv,yv,urputp,xr,yr);
    
    urpwp = upwp.*cos(theta_x) + vpwp.*sin(theta_x);
    urpwp = interp_grid(xv,yv,urpwp,xr,yr);
    
    utpwp = -upwp.*sin(theta_x) + vpwp.*cos(theta_x);
    utpwp = interp_grid(xv,yv,utpwp,xr,yr);
    
    wpwp = interp_grid(xv,yv,wpwp,xr,yr);
    
    p_hyd_rho0 = interp_grid(xv,yv,p_hyd_rho0,xr,yr);
    p_hyd_rhob = interp_grid(xv,yv,p_hyd_rhob,xr,yr);
    p_nhyd = interp_grid(xv,yv,p_nhyd,xr,yr);
    
    eta = interp_grid(xv,yv,eta,xr,yr);
    
    % get tracers on new grid
    S_r = interp_grid(xv,yv,S,xr,yr);
    T_r = interp_grid(xv,yv,T,xr,yr);
    rho_r = interp_grid(xv,yv,rho,xr,yr);
    
    clear up* vp* u v S T rho
    %% momentum in theta direction, put everything on LHS
    % dutheta/dt + durutheta/dr + 1/r*dutheta^2/dtheta + dwutheta/dz - ur*utheta/r  (Du/Dt = US + NL1 + NL2 + NL3 + NL4)
    % + durputhetap/dr + 1/r*duptheta^2/dtheta + dupthetawp/dz - upruptheta/r  (Reynolds = Re1 + Re2 + Re3 + Re4)
    % -fur (coriolis = Cor1 + Cor2)
    % +1/rho0*dp_hyd/dtheta + 1/rho0*dp_nhyd/dtheta (pressure = PGh + PGnh)
    % -nu_H*d2u/dx2 -nu_H*d2u/dy2 - nu*d2u/dz2 (visc = Diss1 + Diss2 + Diss3)
    % -Fl (low freq nudging)
    % = 0
    % derivatives, use central diff where possible (in reality suntans uses
    % upwind)
    % dir 1 = r, dir 2 = theta, dir 3 = z
    disp('compute average theta mometum terms...')
    
%     disp('NL1 theta')
    % +durut/dr 
    temp = u_r.*u_theta;
    Mtheta_Nl1 = u_theta +nan;
    Mtheta_Nl1(2:end-1,:,:) = (temp(3:end,:,:) - temp(1:end-2,:,:))...
                          ./ (R(3:end,:,:) - R(1:end-2,:,:));                   
%     disp('NL2 theta')
%     1/r * d u_theta^2/dtheta
    temp = u_theta.^2;
    Mtheta_Nl2 = u_theta +nan;
    Mtheta_Nl2(:,2:end-1,:) = (temp(:,3:end,:) - temp(:,1:end-2,:))...
                          ./ (THETA(:,3:end,:) - THETA(:,1:end-2,:)) ...
                          ./R(:,2:end-1,:);
%     disp('NL3 theta')
%     d w*u_theta/dz
    temp = u_theta.*w_center;
    Mtheta_Nl3 = u_theta + nan;
    Mtheta_Nl3(:,:,2:end-1) = (temp(:,:,3:end) - temp(:,:,1:end-2))...
                          ./ (Z_r(:,:,3:end) - Z_r(:,:,1:end-2));
%     disp('NL4 theta')
%     ur * u_theta / r
    Mtheta_Nl4 = -u_r.*u_theta./R;
    
%     disp('Re1 theta')
%    +durputp/dr 
    temp = urputp;
    Mtheta_Re1 = u_theta +nan;
    Mtheta_Re1(2:end-1,:,:) = (temp(3:end,:,:) - temp(1:end-2,:,:))...
                          ./ (R(3:end,:,:) - R(1:end-2,:,:));    
%     disp('Re2 theta')
%     1/r * d ututp^2/dtheta 
    temp = utputp;
    Mtheta_Re2 = u_theta +nan;
    Mtheta_Re2(:,2:end-1,:) = (temp(:,3:end,:) - temp(:,1:end-2,:))...
                          ./ (THETA(:,3:end,:) - THETA(:,1:end-2,:)) ...
                          ./R(:,2:end-1,:);
%     disp('Re3 theta')
%     d utpwp/dz
    temp = utpwp;
    Mtheta_Re3 = u_theta + nan;
    Mtheta_Re3(:,:,2:end-1) = (temp(:,:,3:end) - temp(:,:,1:end-2))...
                          ./ (Z_r(:,:,3:end) - Z_r(:,:,1:end-2));
%     disp('Re4 theta')
%     urputp / r
    Mtheta_Re4 = -urputp./R;
    
%     disp('Cor1 theta')
%     f * u_r
    Mtheta_Cor1 = f*u_r;
    
    % -bw sin(theta)
    Mtheta_Cor2 = -b*w_center.*sin(THETA);
    
    % +1/rho0/r * dp_hyd_rho0/dtheta
    temp = p_hyd_rho0;
    Mtheta_PGh_rho0 = u_r + nan;
    Mtheta_PGh_rho0(:,2:end-1,:) = (temp(:,3:end,:) - temp(:,1:end-2,:))...
                          ./ (THETA(:,3:end,:) - THETA(:,1:end-2,:))...
                          ./ (rho0 .* R(:,2:end-1,:));
    % +1/rho0/r * dp_hyd_rhob/dtheta 
    temp = p_hyd_rhob;
    Mtheta_PGh_rhob = u_r + nan;
    Mtheta_PGh_rhob(:,2:end-1,:) = (temp(:,3:end,:) - temp(:,1:end-2,:))...
                          ./ (THETA(:,3:end,:) - THETA(:,1:end-2,:))...
                          ./ (rho0 .* R(:,2:end-1,:));
    % +1/rho0/r * dp_nhyd/dtheta
    temp = p_nhyd;
    Mtheta_PGnh = u_r + nan;
    Mtheta_PGnh(:,2:end-1,:) = (temp(:,3:end,:) - temp(:,1:end-2,:))...
                          ./ (THETA(:,3:end,:) - THETA(:,1:end-2,:))...
                          ./ (rho0 .* R(:,2:end-1,:));
    % Horiz Dissipation terms from cartesian Mx,My, rotate, then interp
    if do_momentum_xy
        load('SUNTANS_results_momentum_avg_x.mat', 'Mx_DissH');
        load('SUNTANS_results_momentum_avg_y.mat', 'My_DissH');
        [Mr_DissH,Mtheta_DissH] = cart_to_cyl(Mx_DissH, My_DissH,theta_x);
        clear Mx_DissH My_DissH                     
        Mr_DissH = interp_grid(xv,yv,Mr_DissH,xr,yr);
        Mtheta_DissH = interp_grid(xv,yv,Mtheta_DissH,xr,yr);
    else
        Mr_DissH = 0.*Mtheta_Cor1;
        Mtheta_DissH = 0.*Mtheta_Cor1;
    end
    
    % - nuV * d2ut/dz2
    Mtheta_DissV = u_r + nan;
    dxtemp = diff(0.5*(Z_r(:,:,1:end-1)+Z_r(:,:,2:end)), 1, 3);
    Mtheta_DissV(:,:,2:end-1) = -nu * diff(u_theta,2,3) ./ dxtemp.^2;
    clear dxtemp
    % (u - ULm)/TauL
    Mtheta_Fl = lowfreq_nudging*(u_theta - ULm_theta)/TauL;
          
     % take depth avg
    disp('take depth average of theta momentum...')
    Mtheta_Nl1_2d = get_depth_avg(Mtheta_Nl1,dZ_r,3);
    Mtheta_Nl2_2d = get_depth_avg(Mtheta_Nl2,dZ_r,3);
    Mtheta_Nl3_2d = get_depth_avg(Mtheta_Nl3,dZ_r,3);
    Mtheta_Nl4_2d = get_depth_avg(Mtheta_Nl4,dZ_r,3);
    Mtheta_Re1_2d = get_depth_avg(Mtheta_Re1,dZ_r,3);
    Mtheta_Re2_2d = get_depth_avg(Mtheta_Re2,dZ_r,3);
    Mtheta_Re3_2d = get_depth_avg(Mtheta_Re3,dZ_r,3);
    Mtheta_Re4_2d = get_depth_avg(Mtheta_Re4,dZ_r,3);
    Mtheta_Cor1_2d = get_depth_avg(Mtheta_Cor1,dZ_r,3);
    Mtheta_Cor2_2d = get_depth_avg(Mtheta_Cor2,dZ_r,3);
    Mtheta_PGh_rho0_2d = get_depth_avg(Mtheta_PGh_rho0,dZ_r,3);
    Mtheta_PGh_rhob_2d = get_depth_avg(Mtheta_PGh_rhob,dZ_r,3);
    Mtheta_PGnh_2d = get_depth_avg(Mtheta_PGnh,dZ_r,3);
    Mtheta_DissH_2d = get_depth_avg(Mtheta_DissH,dZ_r,3);
    Mtheta_DissV_2d = get_depth_avg(Mtheta_DissV,dZ_r,3);
    Mtheta_Fl_2d = get_depth_avg(Mtheta_Fl,dZ_r,3);
    
    
    disp('saving avg theta momentum...')   
    save('SUNTANS_results_momentum_avg_theta.mat','-v7.3','u_theta','U_theta','eta','depth',...
        'b','f','Cd','nu','nu_H','lat','xr','yr','z','r','theta','mtime','Mtheta*')
    
    
    %% momentum in r direction, put everything on LHS
    % dur/dt + durur/dr + 1/r*duthetaur/dtheta + dwur/dz - utheta^2/r  (Du/Dt = US + NL1 + NL2 + NL3 + NL4)
    % + d upr upr/dr + 1/r*d uptheta upr/dtheta + d upr wp/dz - uptheta^2/r  (Reynolds = Re1 + Re2 + Re3 + Re4)
    % -futheta + bw (coriolis = Cor1 + Cor2)
    % +1/rho0*dp_hyd/dr + 1/rho0*dp_nhyd/dr (pressure = PGh + PGnh)
    % -nu_H*d2ur/dx2 -nu_H*d2ur/dy2 - nu*d2ur/dz2 (visc = Diss1 + Diss2 + Diss3)
    % -Fl (low freq nudging)
    % = 0
    % derivatives, use central diff where possible (in reality suntans uses
    % upwind)
    % dir 1 = r, dir 2 = theta, dir 3 = z
    disp('compute average R mometum terms...')
    
%     disp('NL1 r')
    % +durur/dr 
    temp = u_r.*u_r;
    Mr_Nl1 = u_theta +nan;
    Mr_Nl1(2:end-1,:,:) = (temp(3:end,:,:) - temp(1:end-2,:,:))...
                          ./ (R(3:end,:,:) - R(1:end-2,:,:)); 
                      
%     disp('NL2 r')
%     1/r * d u_theta ur /dtheta
    temp = u_theta.* u_r;
    Mr_Nl2 = u_theta +nan;
    Mr_Nl2(:,2:end-1,:) = (temp(:,3:end,:) - temp(:,1:end-2,:))...
                          ./ (THETA(:,3:end,:) - THETA(:,1:end-2,:)) ...
                          ./R(:,2:end-1,:);
%     disp('NL3 r')
%     d w*u_r/dz
    temp = u_r.*w_center;
    Mr_Nl3 = u_theta + nan;
    Mr_Nl3(:,:,2:end-1) = (temp(:,:,3:end) - temp(:,:,1:end-2))...
                          ./ (Z_r(:,:,3:end) - Z_r(:,:,1:end-2));
%     disp('NL4 r')
%     u_theta^2 / r
    Mr_Nl4 = -u_theta.*u_theta./R;
        
%     disp('Re1 r')
%    +durpurp/dr 
    temp = urpurp;
    Mr_Re1 = u_theta +nan;
    Mr_Re1(2:end-1,:,:) = (temp(3:end,:,:) - temp(1:end-2,:,:))...
                          ./ (R(3:end,:,:) - R(1:end-2,:,:));    
%     disp('Re2 r')
%     1/r * d urputp^2/dtheta 
    temp = urputp;
    Mr_Re2 = u_theta +nan;
    Mr_Re2(:,2:end-1,:) = (temp(:,3:end,:) - temp(:,1:end-2,:))...
                          ./ (THETA(:,3:end,:) - THETA(:,1:end-2,:)) ...
                          ./R(:,2:end-1,:);
%     disp('Re3 r')
%     d urpwp/dz
    temp = urpwp;
    Mr_Re3 = u_theta + nan;
    Mr_Re3(:,:,2:end-1) = (temp(:,:,3:end) - temp(:,:,1:end-2))...
                          ./ (Z_r(:,:,3:end) - Z_r(:,:,1:end-2));
%     disp('Re4 r')
%     utputp / r
    Mr_Re4 = -utputp./R;
    
%     disp('Cor1 r')
%     -f * u_theta
    Mr_Cor1 = -f*u_theta ;
    
    % +bw cos(theta)
    Mr_Cor2 = b*w_center.*cos(THETA);
    
    % +1/rho0 * dp_hyd_rho0/dr
    temp = p_hyd_rho0;
    Mr_PGh_rho0 = u_r + nan;
    Mr_PGh_rho0(2:end-1,:,:) = (temp(3:end,:,:) - temp(1:end-2,:,:))...
                          ./ (R(3:end,:,:) - R(1:end-2,:,:))...
                          ./ (rho0);
    % +1/rho0 * dp_hyd_rhob/dr
    temp = p_hyd_rhob;
    Mr_PGh_rhob = u_r + nan;
    Mr_PGh_rhob(2:end-1,:,:) = (temp(3:end,:,:) - temp(1:end-2,:,:))...
                          ./ (R(3:end,:,:) - R(1:end-2,:,:))...
                          ./ (rho0);
    % +1/rho0 * dp_nhyd/dr
    temp = p_nhyd;
    Mr_PGnh = u_r + nan;
    Mr_PGnh(2:end-1,:,:) = (temp(3:end,:,:) - temp(1:end-2,:,:))...
                          ./ (R(3:end,:,:) - R(1:end-2,:,:))...
                          ./ (rho0);
    
    % Mr_dissH -> we already got this above, rotating x,y components
    
    % - nuV * d2ur/dz2
    Mr_DissV = u_r + nan;
    dxtemp = diff(0.5*(Z_r(:,:,1:end-1)+Z_r(:,:,2:end)), 1, 3);
    Mr_DissV(:,:,2:end-1) = -nu * diff(u_r,2,3) ./ dxtemp.^2;
    clear dxtemp
    
    % (u - ULm)/TauL
    Mr_Fl = lowfreq_nudging*(u_r - ULm_r)/TauL;
    
     % take depth avg
    disp('take depth average of R momentum...')
    Mr_Nl1_2d = get_depth_avg(Mr_Nl1,dZ_r,3);
    Mr_Nl2_2d = get_depth_avg(Mr_Nl2,dZ_r,3);
    Mr_Nl3_2d = get_depth_avg(Mr_Nl3,dZ_r,3);
    Mr_Nl4_2d = get_depth_avg(Mr_Nl4,dZ_r,3);
    Mr_Re1_2d = get_depth_avg(Mr_Re1,dZ_r,3);
    Mr_Re2_2d = get_depth_avg(Mr_Re2,dZ_r,3);
    Mr_Re3_2d = get_depth_avg(Mr_Re3,dZ_r,3);
    Mr_Re4_2d = get_depth_avg(Mr_Re4,dZ_r,3);
    Mr_Cor1_2d = get_depth_avg(Mr_Cor1,dZ_r,3);
    Mr_Cor2_2d = get_depth_avg(Mr_Cor2,dZ_r,3);
    Mr_PGh_rho0_2d = get_depth_avg(Mr_PGh_rho0,dZ_r,3);
    Mr_PGh_rhob_2d = get_depth_avg(Mr_PGh_rhob,dZ_r,3);
    Mr_PGnh_2d = get_depth_avg(Mr_PGnh,dZ_r,3);
    Mr_DissH_2d = get_depth_avg(Mr_DissH,dZ_r,3);
    Mr_DissV_2d = get_depth_avg(Mr_DissV,dZ_r,3);
    Mr_Fl_2d = get_depth_avg(Mr_Fl,dZ_r,3);
    
    disp('saving avg R momentum...')     
    save('SUNTANS_results_momentum_avg_r.mat','-v7.3','u_r','U_r','eta','depth',...
        'b','f','Cd','nu','nu_H','lat','xr','yr','z','r','theta','mtime','Mr*',...
        'S_r','T_r','rho_r')
                      
    %% go term by term to save memory
%     terms = {'Nl1','Nl2','Nl3','Re1','Re2','Re3','Cor1','Cor2','PGh_rho0',...
%         'PGh_rhob','PGnh','Diss1','Diss2','Diss3'};
%     % this works for vectors but not for tensors
%     for i=1:length(terms)
%         disp(['processing ' terms{i}])
%         % read 3D terms    
%         load('SUNTANS_results_momentum_avg_x.mat', ['Mx_' terms{i}]);
%         load('SUNTANS_results_momentum_avg_y.mat', ['My_' terms{i}]);
%         %rename variable
%         eval(['Fx=Mx_'  terms{i} ';'])
%         eval(['Fy=My_'  terms{i} ';'])
%         % transform coordinates
%         [Fr,Ftheta] = cart_to_cyl(Fx, Fy, theta);
%         %rename variable
%         eval(['Mr_'  terms{i} '=Fr;'])
%         eval(['Mtheta_'  terms{i} '=Ftheta;'])
% 
%         % read 2D terms    
%         load('SUNTANS_results_momentum_avg_x.mat', ['Mx_' terms{i} '_2d']);
%         load('SUNTANS_results_momentum_avg_y.mat', ['My_' terms{i} '_2d']);
%         %rename variable
%         eval(['Fx=Mx_'  terms{i} '_2d;'])
%         eval(['Fy=My_'  terms{i} '_2d;'])
%         % transform coordinates
%         [Fr,Ftheta] = cart_to_cyl(Fx, Fy, squeeze(theta(:,:,1)));
%         %rename variable
%         eval(['Mr_'  terms{i} '_2d=Fr;'])
%         eval(['Mtheta_'  terms{i} '_2d=Ftheta;'])
% 
%         % save results
%         save('SUNTANS_results_momentum_avg_r.mat','-append', ['Mr_'  terms{i}],['Mr_'  terms{i} '_2d'])
%         save('SUNTANS_results_momentum_avg_theta.mat','-append', ['Mtheta_'  terms{i}], ['Mtheta_'  terms{i} '_2d'])
%         clear Mr* Mtheta* Mx* My* Fx Fy
%     end

end
disp('done!')
% clear
