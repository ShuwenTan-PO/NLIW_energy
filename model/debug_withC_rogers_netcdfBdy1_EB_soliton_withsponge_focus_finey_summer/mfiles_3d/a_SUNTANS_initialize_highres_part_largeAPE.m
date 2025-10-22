close all
clear all
clc

diary a_SUNTANS_initialize_log.txt
addpath(genpath('../../../Matlab/'))

igrid = 1;
idepth = 1;
iprofile = 1;
iinit_ic = 0;
iinit_bc = 1;

lambdaM2=70.1e3;
% C = 1.57;
D0=600;
He = 313.1;
datadir='../rundata';

l_sponge = getvalue([datadir,'/suntans.dat'],'sponge_distance');
L=1.5*lambdaM2+2*l_sponge;
W=1.5*lambdaM2;
INPUT.slope = 0.05;%300/4.5e3;%
INPUT.width_he2 = 26000;%width at He/2
INPUT.width_top = INPUT.width_he2 - 2*(He/2)/INPUT.slope; %width of top
INPUT.r = INPUT.width_he2/2+(D0-He/2)/INPUT.slope;% radius of base
INPUT.width_base = 2*INPUT.r; % width of base
INPUT.x0 = 2*l_sponge+1.5*lambdaM2/2;%L/2;
INPUT.y0 = W/2;
INPUT.z0 = -D0; % level of base
INPUT.zmax = -4;

if igrid
    dx=80; %550;%550/1.5; % average dx in side regions if stretching
    dxr=20; % 30;%30/1.5; % dx in central refined region if stretching
    BC = [2 5 2 5]; %E,N,W,S Boundary conditions, 
    % 1 solid free-slip, 2 velocity, 3 free-surface, 5 periodic, 6 mixed 2/3
    STRETCHING=false;
    FOCUS = true; %false;
    CHEBYCHEV=false;
    rmax=1.05;   
    x0=0;
    y0=0;
    theta=0;
    if ~STRETCHING
        Nx = round(L/dx);% total cells no stretching
        Ny = round(W/dx);
        Lr =0;
        Nr=0;
    else
        Ls = 0.2*W*[1 1]; % size of stretching region
        Lr = [L W] - 2*Ls;
        Nr = ceil(Lr./dxr); % number of cells in refined region
        Ns = ceil(2*Ls./dx); % number of cells in side regions    
        Nx = Nr(1)+Ns(1); % total cells with stretching
        Ny = Nr(2)+Ns(2);    
    end 
    INPUT.refine_radius =2*INPUT.r; % radius of refinement for FOCUS
    K=(Nr./Lr).*([L W]./[Nx Ny]); % K is refinement relative to avg dx
    INPUT.rad = INPUT.refine_radius*Nx/L;%10;% radius/dx
    INPUT.xc = INPUT.x0*Nx/L;%100; % xc/dx
    INPUT.yc = INPUT.y0*Ny/W;%40; % yc/dy
    INPUT.dxmax = dx;%100;
    INPUT.dxmin = dxr;%10;
    save SUNTANS_grid.mat -v7.3
    quadgrid_periodic(datadir,L,W,Nx,Ny,BC,STRETCHING,CHEBYCHEV,Lr,K,rmax,x0,y0,theta,FOCUS,INPUT)
    

end

if idepth
    load SUNTANS_grid.mat
    shape='cone';

    x0 = 0;
    y0 = 0;
    theta=0;
    dtrend = 'yes';
    bathy_dir = [];
    dmin=[];
    extreme=[];
    
    save SUNTANS_grid.mat -v7.3
    SUNTANS_depth_ideal(datadir,D0,shape,bathy_dir,x0,y0,theta,dtrend,dmin,extreme,INPUT)
end


if iprofile
    load SUNTANS_grid.mat
    depth_prof=75;
    xc = INPUT.x0;
    yc = INPUT.y0;
    r_bot = INPUT.r;
    r_top = INPUT.r-(D0-depth_prof)/INPUT.slope;
    r_mid = INPUT.r-(D0-He/2)/INPUT.slope;
    r_mid_crnr = sind(45)*r_mid;
    name = {'Center','Front_slope_top','Front_slope_mid','Front',...
        'Back','Top','Bottom',...
        'Back_slope_top','Top_slope_top','Bottom_slope_top',...
        'Back_slope_mid','Top_slope_mid','Bottom_slope_mid',...
        'NE_slope_mid','NW_slope_mid',...
        'SW_slope_mid','SE_slope_mid',...
        'Mid'};
    N = [yc, yc,yc,yc,...
        yc,yc+r_bot,yc-r_bot,...
        yc,yc+r_top,yc-r_top,...
        yc,yc+r_mid,yc-r_mid,...
        yc+r_mid_crnr,yc+r_mid_crnr,...
        yc-r_mid_crnr,yc-r_mid_crnr,...
        W/2];
    E = [xc, xc+r_top,xc+r_mid, xc+r_bot,...
        xc-r_bot,xc,xc,...
        xc-r_top,xc,xc,...
        xc-r_mid,xc,xc,...
        xc+r_mid_crnr,xc-r_mid_crnr,...
        xc-r_mid_crnr,xc+r_mid_crnr,...
        L/2];

    % reef top
    r_reeftop = INPUT.r-D0/INPUT.slope;
    % around island
    n_r = 24*2;
    n_theta = 48*2;
    r_n = INPUT.r-(1-linspace(0,1,n_r))*D0/INPUT.slope;
    dr = diff(r_n); dr = dr(1);
    r_n = [r_reeftop-dr, r_n];
%     ratio_n = 10;
%     r_n = [r_reeftop-100, r_reeftop, INPUT.r-(1-ratio_n.^(linspace(-1,0,n_r)))*D0/INPUT.slope];
    n_r = length(r_n);
%     r_n = [r_top, r_n];
%     theta_n = [linspace(-180,-120,3), linspace(-90,90,n_theta), linspace(120, 150, 2)];   
    theta_n = linspace(-180,180-360/n_theta,n_theta);
    name_ = cell(1,n_r*n_theta);
    N_ = nan(1,n_r*n_theta);
    E_ = nan(1,n_r*n_theta);
    for i=1:n_r
        for j=1:n_theta
            name_{1,(i-1)*n_theta+j} = strcat('r_',num2str(i),'_theta_',num2str(j));
            N_(1,(i-1)*n_theta+j) = yc+sind(theta_n(j))*r_n(i);
            E_(1,(i-1)*n_theta+j) = xc+cosd(theta_n(j))*r_n(i);
        end
    end

    % upstream and downstream sections
%     n_x = 4;
    n_x = 3;
    n_y = 100;%20;
    dx_ = 25000; %13819;
    y_n = linspace(0,2*INPUT.y0,n_y);
%     x_n = [xv(1)+dx_*[1:n_x/2], xv(end)+dx_*[-n_x/2:-1]];
%     x_n = [xv(1)+dx_*[0:n_x/2-1], xv(end)+dx_*[-n_x/2+1:0]];
    x_n = [xv(1)+dx_*[1:2], xv(end)-dx_];
    y_n(1) = yv(1); y_n(end) = yv(end);
    name__ = cell(1,n_x*n_y);
    N__ = nan(1,n_x*n_y);
    E__ = nan(1,n_x*n_y);
    for i = 1:n_x
        for j = 1:n_y
            name__{1,(i-1)*n_y+j} = strcat('x_',num2str(i),'_y_',num2str(j));
            N__(1,(i-1)*n_y+j) = y_n(j);
            E__(1,(i-1)*n_y+j) = x_n(i);
        end
    end

    figure
    scatter(xv,yv,5,-Depth)
    hold on
    plot(E,N,'k*');
    hold on
    plot(E_,N_,'r.');
    hold on
    plot(E__,N__,'r*');
    xlim([0, max(xv)]); ylim([0, max(xv)]);
    pbaspect([1 1 1])
    print -djpeg -r300 suntans_prof_loc

    name=[name, name_, name__, {'EBC','NBC','WBC','SBC'}];
    N = [N, N_, N__];
    E = [E, E_, E__];

    profile_points_ideal(datadir,name,E,N)
   
end

if iinit_ic
   initial_ideal = 1;
   init_filename = '/sun_IC.nc';
   SUNTANS_initial_conditions(initial_ideal,datadir,...
                 init_filename,'');
end

% %% create DJL solution
% % compute alpha and write it in suntans.dat
% % data from Davis
% depth_OR3 = [0 10 20 30 50 80 100 150 250 500];
% density_OR3 = [1020.94 1021.10 1021.59 1022.37 1023.38 1024.55 1025.10 1025.96 1027.04 1029.09];
% temp_OR3 = [30 29.63 28.33 26.37 23.85 20.83 19.24 16.62 13.32 8.00];
% salt_OR3 = [33.95 33.95 33.96 34.10 34.31 34.56 34.61 34.60 34.48 34.41];
% 
% lat = 23;
% long = 117;
% p = gsw_p_from_z(-depth_OR3,lat);
% 
% [SA, ] = gsw_SA_from_SP(salt_OR3,p,long,lat);
% CT = gsw_CT_from_t(SA,temp_OR3,p);
% 
% alpha = gsw_alpha(SA,CT,p);
% disp(strcat('mean gamma = ',string(mean(alpha)),' based on mooring data'))
% disp(strcat('while the gamma in suntans.dat = ',string(getvalue([datadir,'/suntans.dat'],'gamma'))))
% 
% run Type2_BC
% x_offset = 4*getvalue([datadir,'/suntans.dat'],'thetaramptime')*DJLES_solution.c;
% DJLES_solution.xc = DJLES_solution.xc+x_offset;
% DJLES_solution.xe = DJLES_solution.xe+x_offset;
% % feed wave from eastern boundary
% DJLES_solution.u = -DJLES_solution.u;
% figure
% contourf(DJLES_solution.xc, DJLES_solution.zc, DJLES_solution.u)
% print -djpeg -r300 DJLES_solution_u
% save soliton_file.mat DJLES_solution
% 
% dist_run = getvalue([datadir,'/suntans.dat'],'nsteps')*getvalue([datadir,'/suntans.dat'],'dt')/DJLES_solution.c/lambdaM2;
% disp(strcat('Distance that solitary wave propagates during the model run = ',string(dist_run), ' lambda(M2)'));
% disp(strcat('Wavelength of solitary wave = ',string(DJLES_solution.wavelength), ' m'));
% disp(strcat('dx at boundary = ',string(max(diff(xv))), ' m is ', string(max(diff(xv))/DJLES_solution.wavelength), 'wave length'));
% disp(strcat('Total run time = ',string(getvalue([datadir,'/suntans.dat'],'nsteps')*getvalue([datadir,'/suntans.dat'],'dt')/3600), ' hr'));
% disp(strcat('Time for solitary wave propagates across the boundary = ',string(DJLES_solution.wavelength/DJLES_solution.c/3600), ' hr'));
% disp(strcat('Time for solitary wave propagates across the sponge layer = ',string(getvalue([datadir,'/suntans.dat'],'sponge_distance')/DJLES_solution.c/3600), ' hr'));
% disp(strcat('While the ramping time scale = ', string(getvalue([datadir,'/suntans.dat'],'thetaramptime')/3600), ' hr'));

%% 
if iinit_bc   
   n_soliton = 1;
   bc_filename = '/sun_BC.nc';
   soliton_file = 'soliton_file_summer_largeAPE.mat';
%    bc_dt = getvalue([datadir,'/suntans.dat'],'dt')*10;  % ideal bc dt step [-1/24:1/24:4]*24*3600; % seconds since startime
   bc_dt = getvalue([datadir,'/suntans.dat'],'bc_dt');  % ideal bc dt step [-1/24:1/24:4]*24*3600; % seconds since startime
   soliton_interval = 0;
   method_sponge = 1;
   SUNTANS_boundary_conditions_soliton_short(n_soliton,datadir,...
       bc_filename,bc_dt,soliton_file,soliton_interval,method_sponge);   
%    SUNTANS_boundary_conditions_soliton_nopar_short_new(n_soliton,datadir,...
%        bc_filename,bc_dt,soliton_file,soliton_interval,method_sponge);   
%    SUNTANS_boundary_conditions_rogers(bc_ideal,datadir,...
%        bc_filename,bc_dt);   
end

% ======================  make plots on test grid =========================
load('SUNTANS_grid.mat')
load SUNTANS_grid_quadgrid.mat xp yp
Xv = reshape(xv,Nx,Ny); Yv = reshape(yv,Nx,Ny);
Xp = reshape(xp,Nx+1,Ny+1); Yp = reshape(yp,Nx+1,Ny+1);

figure(1)
plot(Xp,Yp,'k.-',markersize=2)
hold on
plot(Xp',Yp','k.-',markersize=2)
hold on
plot(Xv,Yv,'b.-',markersize=2)
hold on
plot(Xv',Yv','b.-',markersize=2)
% mesh(Xv, Yv, D)
pbaspect([2 1 1])

print -djpeg -r300 figure_edge_cell

figure
plot(Xp(:,1),'k.')
hold on
plot(Xv(:,1),'r.')

figure
plot(Yp(1,:),'k.')
hold on
plot(Yv(1,:),'r.')

figure
plot(Xp(1:end-1,1)-Xv(:,1),'c.')
hold on
plot(Yp(1,1:end-1)-Yv(1,:),'b.')

figure
plot((Xp(1:end-1,1)-Xv(:,1))./diff(Xp(:,1)),'c.')
hold on
plot((Yp(1,1:end-1)-Yv(1,:))./diff(Yp(1,:)),'b.')

diary off