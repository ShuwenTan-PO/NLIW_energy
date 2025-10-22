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
INPUT.slope = 0.05;%300/4.5e3;%
INPUT.width_he2 = 26000;%width at He/2
INPUT.width_top = INPUT.width_he2 - 2*(He/2)/INPUT.slope; %width of top
INPUT.r = INPUT.width_he2/2+(D0-He/2)/INPUT.slope;% radius of base
INPUT.width_base = 2*INPUT.r; % width of base
L=1.5*lambdaM2-INPUT.width_top-(INPUT.width_base-INPUT.width_top)/2;
W=l_sponge/2;
INPUT.x0 = L/2;
INPUT.y0 = W/2;
INPUT.z0 = -D0; % level of base
INPUT.zmax = -4;

if igrid
    dx=550; %80;   % average dx in side regions if stretching
    dxr=30; %20;  % dx in central refined region if stretching
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
%         Ls = 0.2*W*[1 1]; % size of stretching region
        Ls = 0.2*L*[1 1]; % size of stretching region
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
    shape='leftslope_2d_island';

    INPUT.way_point1 = L/2-(INPUT.width_base-INPUT.width_top)/4;
    INPUT.way_point2 = L/2+(INPUT.width_base-INPUT.width_top)/4;
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

    % upstream and downstream sections
    n_x = 4;
    n_y = 20;
    dx_ = 14000; %13819;
    y_n = linspace(0,2*INPUT.y0,n_y);
    x_n = [xv(1)+dx_*[1:4]];
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
    plot(E__,N__,'r*');
    xlim([0, max(xv)]); ylim([0, max(xv)]);
    pbaspect([1 1 1])
%     print -djpeg -r300 suntans_prof_loc

    name=name__;
    N = N__;
    E = E__;

    profile_points_ideal(datadir,name,E,N)
   
end

if iinit_ic
   initial_ideal = 1;
   init_filename = '/sun_IC.nc';
   SUNTANS_initial_conditions(initial_ideal,datadir,...
                 init_filename,'');
end

%% 
if iinit_bc   
   n_soliton = 1;
   bc_filename = '/sun_BC_2d.nc';
   soliton_file = 'soliton_file_summer_smallAPE.mat';
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