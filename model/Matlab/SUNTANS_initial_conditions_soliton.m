function [] = SUNTANS_initial_conditions_soliton(datadir,...
                 init_filename,soliton_file,x_soliton1,n_soliton,method_init)
%%
% SUNTANS Initial Conditions
% Justin Rogers
% Stanford Univerity
%% %%%

% updates
% by S.Tan 2023/02/24
% 1) initialize the model with only temperature
% 2) n_soliton: default = 1 meaning one solitary wave
% by S.Tan 2023/02/27
% 1) method_input = 0: initial condition as temp field
%    method_input = 1: initial condition as uf field
% TO UPDATE: 
% 1) parameters needed for DJL_solution.density getting back to temperature
% , need to update function input, like param_T
% NOTE:
% 1) suntans.dat currently turned on readinitialnc=1, gamma \= 0, so read 
% temp from the nc file; to read in velocity, must turn on initialUNC=1 
% in suntans.dat
delete(gcp('nocreate'))
poolobj = parpool('local'); 

if nargin<5
    n_soliton = 1; 
end
if nargin<6
    method_init = 0;
end

if method_init == 1
    if length(getvalue([datadir,'/suntans.dat'],'initialUNC'))==0
        error('Error: to turn on method_init = 1, initialUNC in suntans.dat must set to 1')
    end
end

%% load grid data
c = load([datadir,'/cells.dat']);
Nc = size(c,1);
if size(c,2)<9
    xv = c(:,1);
    yv = c(:,2);
else
    xv = c(:,2);
    yv = c(:,3);
end
cells = c(:,4:7);
neigh = c(:,8:11);
nfaces = c(:,1);
numsides =nfaces(1,1);

points = load([datadir,'/points.dat']);
Np = size(points,1);
xp = points(:,1);
yp = points(:,2);

ed = load([datadir,'/edges.dat']);
Ne = size(ed,1);
edges = ed(:,1:2);
mark = ed(:,3);
grad =ed(:,4:5);
edgep = 0:(Ne-1); % this is the edge index in C format
for j=1:Ne
    xe(j) = 0.5*(xp(edges(j,1)+1,1)+xp(edges(j,2)+1,1));
    ye(j) = 0.5*(yp(edges(j,1)+1,1)+yp(edges(j,2)+1,1));
    
    % get outward unit normal vectors for computational points only
    n1(j) = yp(edges(j,1)+1,1)-yp(edges(j,2)+1,1);
    n2(j) = xp(edges(j,1)+1,1)-xp(edges(j,2)+1,1);
    n = sqrt(n1(j).^2+n2(j).^2);
    n1(j) = n1(j)/n;
    n2(j) = n2(j)/n;
end
n2=-n2; % set correct direction
% we don't care about boundary points for sponge fluxes, set to 0
n1(mark>0)=0;
n2(mark>0)=0;


depth=load([datadir,'/depth.dat']);
dv = depth(:,3);

Nkmax = getvalue([datadir,'/suntans.dat'],'Nkmax');
Nk = Nkmax+zeros(Nc,1);
Nkw=Nkmax+1;
basetime = sprintf('%14.6f',getvalue([datadir,'/suntans.dat'],'basetime'));
mtime_base = datenum(basetime,'yyyymmdd.HHMMSS');
starttime = sprintf('%14.6f',getvalue([datadir,'/suntans.dat'],'starttime'));
mtime_start = datenum(starttime,'yyyymmdd.HHMMSS');
Toffset = mtime_start-mtime_base; % days between basetime and starttime

files = dir(['../data/vertspace.dat']);
if ~isempty(files) % load in suntans file
    dz = load(['../data/vertspace.dat']);
else % figure out myself
    rstretch = getvalue([datadir,'/suntans.dat'],'rstretch');
    dz=1;
    for i=1:Nkmax-1
        dz(i+1) = rstretch*dz(i);        
    end
    dz = dz'*max(dv)/sum(dz);
end
z_r = getz(dz); %depth
% get 2d normal vectors for easy computation
[~,N1]=meshgrid(z_r,n1);
[~,N2]=meshgrid(z_r,n2);


save('SUNTANS_grid.mat','z_r','dz','-append')


%% create netcdf file
disp('creating netcdf IC file');

clear mex
delete([datadir '/*sun_IC.nc']);

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
ncid = netcdf.create([datadir init_filename],cmode);

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Created', ['Created on ' datestr(now)]);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author', '');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description', 'SUNTANS Initial Conditions File');

% define dimensions
Nc_dimid = netcdf.defDim(ncid,'Nc',Nc);
Np_dimid = netcdf.defDim(ncid,'Np',Np);
Ne_dimid = netcdf.defDim(ncid,'Ne',Ne);
Nk_dimid = netcdf.defDim(ncid,'Nk',Nkmax);
Nkw_dimid = netcdf.defDim(ncid,'Nkw',Nkw);
numsides_dimid = netcdf.defDim(ncid,'numsides',numsides);
two_dimid = netcdf.defDim(ncid,'two',2);
time_dimid = netcdf.defDim(ncid,'time',netcdf.getConstant('NC_UNLIMITED'));

% define variables
varid = netcdf.defVar(ncid,'Nk','NC_INT',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','Number of layers at face');

varid = netcdf.defVar(ncid,'neigh','NC_INT',[numsides_dimid Nc_dimid ]);
netcdf.putAtt(ncid,varid,'long_name','Maps every face to its neighbouring faces');
netcdf.putAtt(ncid,varid,'cf_role','face_face_connectivity');

varid = netcdf.defVar(ncid,'xp','NC_DOUBLE',Np_dimid);
netcdf.putAtt(ncid,varid,'long_name','Easting of 2D mesh node');
netcdf.putAtt(ncid,varid,'standard_name','Easting');

varid = netcdf.defVar(ncid,'xv','NC_DOUBLE',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','Easting of 2D mesh face');
netcdf.putAtt(ncid,varid,'standard_name','Easting');

varid = netcdf.defVar(ncid,'z_r','NC_DOUBLE',Nk_dimid);
netcdf.putAtt(ncid,varid,'long_name','depth at layer mid points');
netcdf.putAtt(ncid,varid,'standard_name','ocean_z_coordinate');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'positive','up');

varid = netcdf.defVar(ncid,'mark','NC_INT',Ne_dimid);
netcdf.putAtt(ncid,varid,'long_name','Edge marker type');
netcdf.putAtt(ncid,varid,'units','0 - comp, 1 - boundary, 2, 3');
netcdf.putAtt(ncid,varid,'coordinates','xe, ye');

varid = netcdf.defVar(ncid,'edges','NC_INT',[two_dimid Ne_dimid]);
netcdf.putAtt(ncid,varid,'long_name','Maps every edge to the two nodes it connects');
netcdf.putAtt(ncid,varid,'cf_role','edge_node_connectivity');

varid = netcdf.defVar(ncid,'dz','NC_DOUBLE',Nk_dimid);
netcdf.putAtt(ncid,varid,'long_name','z layer spacing');
netcdf.putAtt(ncid,varid,'units','m');

varid = netcdf.defVar(ncid,'dv','NC_DOUBLE',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','sea floor depth');
netcdf.putAtt(ncid,varid,'standard_name','sea_floor_depth_below_geoid');
netcdf.putAtt(ncid,varid,'units','m');
netcdf.putAtt(ncid,varid,'positive','down');
netcdf.putAtt(ncid,varid,'coordinates','xv yv');
netcdf.putAtt(ncid,varid,'mesh','suntans_mesh');
netcdf.putAtt(ncid,varid,'location','face');

varid = netcdf.defVar(ncid,'yp','NC_DOUBLE',Np_dimid);
netcdf.putAtt(ncid,varid,'long_name','Northing of 2D mesh node');
netcdf.putAtt(ncid,varid,'standard_name','Northing');

varid = netcdf.defVar(ncid,'yv','NC_DOUBLE',Nc_dimid);
netcdf.putAtt(ncid,varid,'long_name','Northing of 2D mesh face');
netcdf.putAtt(ncid,varid,'standard_name','Northing');

varid = netcdf.defVar(ncid,'cells','NC_INT',[numsides_dimid Nc_dimid ]);
netcdf.putAtt(ncid,varid,'long_name','Maps every face to its corner nodes');
netcdf.putAtt(ncid,varid,'standard_name','face_node_connectivity');

varid = netcdf.defVar(ncid,'nfaces','NC_INT',[Nc_dimid]);
netcdf.putAtt(ncid,varid,'long_name','Number of cell faces');

varid = netcdf.defVar(ncid,'grad','NC_INT',[two_dimid Ne_dimid ]);
netcdf.putAtt(ncid,varid,'long_name','Maps every edge to the two faces it connects');
netcdf.putAtt(ncid,varid,'standard_name','edge_face_connectivity');

varid = netcdf.defVar(ncid,'time','NC_DOUBLE',time_dimid);
netcdf.putAtt(ncid,varid,'units','seconds since 1990-01-01 00:00:00');
netcdf.putAtt(ncid,varid,'long_name','time');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'eta','NC_DOUBLE',[Nc_dimid time_dimid ]);
netcdf.putAtt(ncid,varid,'units','meters');
netcdf.putAtt(ncid,varid,'long_name','sea surface elevation');
netcdf.putAtt(ncid,varid,'coordinates','time xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'uf','NC_DOUBLE',[Ne_dimid Nk_dimid time_dimid]);
netcdf.putAtt(ncid,varid,'units','meters second-1');
netcdf.putAtt(ncid,varid,'long_name','Easward water velocity component');
netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

% varid = netcdf.defVar(ncid,'uc','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
% netcdf.putAtt(ncid,varid,'units','meters second-1');
% netcdf.putAtt(ncid,varid,'long_name','Easward water velocity component');
% netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
% netcdf.defVarFill(ncid,varid,false,999999);
% 
% varid = netcdf.defVar(ncid,'vc','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
% netcdf.putAtt(ncid,varid,'units','meters second-1');
% netcdf.putAtt(ncid,varid,'long_name','Northward water velocity component');
% netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
% netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'salt','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
netcdf.putAtt(ncid,varid,'units','ppt');
netcdf.putAtt(ncid,varid,'long_name','Salinity');
netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'temp','NC_DOUBLE',[Nc_dimid Nk_dimid time_dimid]);
netcdf.putAtt(ncid,varid,'units','degrees C');
netcdf.putAtt(ncid,varid,'long_name','Water temperature');
netcdf.putAtt(ncid,varid,'coordinates','time z_r xv yv');
netcdf.defVarFill(ncid,varid,false,999999);

% end define variables
netcdf.endDef(ncid);

% write grid variables
varid = netcdf.inqVarID(ncid,'Nk');
netcdf.putVar(ncid,varid,Nk);

varid = netcdf.inqVarID(ncid,'neigh');
netcdf.putVar(ncid,varid,neigh);

varid = netcdf.inqVarID(ncid,'xp');
netcdf.putVar(ncid,varid,xp);

varid = netcdf.inqVarID(ncid,'xv');
netcdf.putVar(ncid,varid,xv);

varid = netcdf.inqVarID(ncid,'z_r');
netcdf.putVar(ncid,varid,z_r);

varid = netcdf.inqVarID(ncid,'mark');
netcdf.putVar(ncid,varid,mark);

varid = netcdf.inqVarID(ncid,'edges');
netcdf.putVar(ncid,varid,edges);

varid = netcdf.inqVarID(ncid,'dz');
netcdf.putVar(ncid,varid,dz);

varid = netcdf.inqVarID(ncid,'dv');
netcdf.putVar(ncid,varid,dv);

varid = netcdf.inqVarID(ncid,'yp');
netcdf.putVar(ncid,varid,yp);

varid = netcdf.inqVarID(ncid,'yv');
netcdf.putVar(ncid,varid,yv);

varid = netcdf.inqVarID(ncid,'cells');
netcdf.putVar(ncid,varid,cells);

varid = netcdf.inqVarID(ncid,'nfaces');
netcdf.putVar(ncid,varid,nfaces);

varid = netcdf.inqVarID(ncid,'grad');
netcdf.putVar(ncid,varid,grad);



% close file
netcdf.close(ncid);

%% compute fields
disp('computing IC variables')

time = Toffset*86400; % seconds since basetime

% free surface
eta = zeros(Nc,1);

% velocity
uc = zeros(Nc,Nkmax,1);
vc = uc;
ue = zeros(Ne,Nkmax,1);
ve = ue;

% salt
salt = uc;

% temperature
temp = uc;

% age parameters
agec = uc;
agealpha = uc;
agesource = uc;   

eval("load " + soliton_file)
gamma = getvalue([datadir,'/suntans.dat'],'gamma');

rho0 = 1000;

x_range1 = [x_soliton1-max(DJLES_solution.xc)/2, x_soliton1+max(DJLES_solution.xc)/2];
% if x_range1(1)<0
%     x_range1 = x_range1 - x_range1(1);
%     disp('given center location for the first soliton too close to the boundary, location adjusted')        
% end

if method_init == 0
    fprintf('T...')
    if n_soliton == 1
        idx_soliton = find((xv<=x_range1(2))&(xv>=x_range1(1)));
        for i = 1:length(idx_soliton)
            for j = 1:Nkmax
                temp(idx_soliton(i),j,1) = -(interp2(DJLES_solution.zc, ...
                    DJLES_solution.xc+x_range1(1), ...
                    DJLES_solution.density', ...
                    z_r(j), ... 
                    xv(idx_soliton(i)))-rho0)/gamma;
            end
        end
    else
        disp('need to fill this up') 
    end
    temp(find(xv<x_range1(1)),:,1) = ones(length(find(xv<x_range1(1))),1)*temp(idx_soliton(1),:,1);
    temp(find(xv>x_range1(2)),:,1) = ones(length(find(xv>x_range1(2))),1)*temp(idx_soliton(end),:,1);
elseif method_init == 1
    fprintf('T...')
    if n_soliton == 1
        idx_soliton = find((xv<=x_range1(2))&(xv>=x_range1(1)));
        for i = 1:length(idx_soliton)
            for j = 1:Nkmax
                temp(idx_soliton(i),j,1) = -(interp2(DJLES_solution.zc, ...
                    DJLES_solution.xc+x_range1(1), ...
                    DJLES_solution.density', ...
                    z_r(j), ... 
                    xv(idx_soliton(i)))-rho0)/gamma;
            end
        end
    else
        disp('need to fill this up') 
    end
    temp(find(xv<x_range1(1)),:,1) = ones(length(find(xv<x_range1(1))),1)*temp(idx_soliton(1),:,1);
    temp(find(xv>x_range1(2)),:,1) = ones(length(find(xv>x_range1(2))),1)*temp(idx_soliton(end),:,1);
    
%     temp_ = -(interp1(DJLES_solution.zc, ...
%                 DJLES_solution.density(:,1), ...
%                 z_r)-rho0)/gamma;
%     temp = ones(length(xv),1)*temp_';
    if n_soliton == 1
        idx_soliton = find((xe<=x_range1(2))&(xe>=x_range1(1)));
        for i = 1:length(idx_soliton)
            for j = 1:Nkmax
                ue(idx_soliton(i),j,1) = interp2(DJLES_solution.zc, ...
                    (DJLES_solution.xc+x_range1(1)), ...
                    DJLES_solution.u', ...
                    z_r(j), ... 
                    xe(idx_soliton(i)));
            end
        end        
    else
        disp('need to fill this up') 
    end
else
    error('Error, method_init must be 0 (initial condition temp field) or 1 (initial condition uf field)')
end

% compute u at cell faces
fprintf('uface...')
uf = N1.*ue + N2.*ve;

%% write variables

disp('writing variables to netcdf IC file')

ncid = netcdf.open([datadir init_filename],'WRITE');

varid = netcdf.inqVarID(ncid,'time');
netcdf.putVar(ncid, varid,0,1, time);

varid = netcdf.inqVarID(ncid,'eta');
netcdf.putVar(ncid, varid,eta);

varid = netcdf.inqVarID(ncid,'uf');
netcdf.putVar(ncid, varid,uf);

% varid = netcdf.inqVarID(ncid,'uc');
% netcdf.putVar(ncid, varid,uc);
% 
% varid = netcdf.inqVarID(ncid,'vc');
% netcdf.putVar(ncid, varid,vc);

varid = netcdf.inqVarID(ncid,'salt');
netcdf.putVar(ncid, varid,salt);

varid = netcdf.inqVarID(ncid,'temp');
netcdf.putVar(ncid, varid,temp);

netcdf.close(ncid);


%% plot IC conditions
load SUNTANS_grid.mat Nx Ny
GRID.Nx = Nx; GRID.Ny=Ny;
xplot = reshape(xv,GRID.Nx,GRID.Ny)/1000;
yplot = reshape(yv,GRID.Nx,GRID.Ny)/1000;
xeplot = xe/1000;
yeplot = ye/1000;
etaplot = reshape(eta,GRID.Nx,GRID.Ny);
Tplot = reshape(temp(:,1),GRID.Nx,GRID.Ny);
Splot = reshape(salt(:,1),GRID.Nx,GRID.Ny);
if method_init == 0
    Uplot = reshape(uc(:,1),GRID.Nx,GRID.Ny);
    Vplot = reshape(vc(:,1),GRID.Nx,GRID.Ny);
    UplotM = sqrt(Uplot.^2+Vplot.^2);
else
    UplotM = squeeze(sqrt(uf(:,1).^2));
end

figure(3452423)
subplot(2,2,1)
if GRID.Ny>2
pcolor(xplot,yplot,etaplot)
cb=colorbar;
ylabel(cb,'$\eta$ (m)','interpreter','latex')
ylabel('y (km)')
else
    plot(xplot,etaplot,'.')
    ylabel('$\eta$ (m)','interpreter','latex')
end

if exist('NC')
    title(['IC netcdf ' datestr(temp_mtime)])
else
    title(['IC ideal ' datestr(mtime_start)])
end
    
if method_init == 0
    subplot(2,2,2)
    if GRID.Ny>2
    pcolor(xplot,yplot,UplotM)
    cb=colorbar;
    ylabel(cb,'$|u|_{surf}$ (m/s)','interpreter','latex')
    else
        plot(xplot,UplotM,'.')
        ylabel('$|u|_{surf}$ (m/s)','interpreter','latex')
    end
else
    subplot(2,2,2)
    if GRID.Ny>2
    scatter(xeplot,yeplot,5,UplotM)
    cb=colorbar;
    ylabel(cb,'$|u|_{surf}$ (m/s)','interpreter','latex')
    else
        plot(xeplot,UplotM,'.')
        ylabel('$|u|_{surf}$ (m/s)','interpreter','latex')
    end
end

subplot(2,2,3)
if GRID.Ny>2
pcolor(xplot,yplot,Tplot)
cb=colorbar;
ylabel(cb,'$T_{surf}$ (C)','interpreter','latex')
ylabel('y (km)')
else
    plot(xplot,Tplot)
    ylabel('$T_{surf}$ (C)','interpreter','latex');
end

subplot(2,2,4)
if GRID.Ny>2
pcolor(xplot,yplot,Splot)
cb=colorbar;
ylabel(cb,'$S_{surf}$ (psu)','interpreter','latex')
else
    plot(xplot,Splot)
    ylabel('$S_{surf}$ (psu)','interpreter','latex')
end
xlabel('x (km)')

print -djpeg -r300 figure_IC_soliton
close


delete(poolobj)






