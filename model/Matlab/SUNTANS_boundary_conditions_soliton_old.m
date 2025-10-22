function [] = SUNTANS_boundary_conditions_soliton(n_soliton,datadir,...
                   bc_filename,bc_dt,soliton_file,soliton_interval,method_sponge)


%% %%%%%%%%
% SUNTANS Boundary Conditions
% Justin Rogers
% Stanford University
%% %%%%%%%

% clear
% n_soliton = 1;
% datadir='../rundata';
% bcfilename = '/sun_BC.nc';
% bc_time = [0:100:1000]; % seconds since 1990-1-1
% UTM_zone = 50;
% filter option: 'all', 'low','high'

% updates
% by S.Tan 2023/01/10
% 1) mask out the step for converting to lat lon
% 2) change the input NCOM_file to soliton_file, to read in solitary wave
% boundary conditions
% by S.Tan 2023/02/24
% 1) if n_soliton==1: one solitary wave
% 2) if n_soliton>1, then multiple solitary wave, time_interval is an input
% S.Tan 2023/03/20
% !!!!!!fixed a bug: rho_b = -gamma*rho_0*T instead of -gamma*T+rho0!!!!!!!

%% load grid data
tic
if nargin<5
    soliton_interval = 0;
end
if nargin<6
    method_sponge = 1;
end
disp('creating boundary condition files...')

% clear mex
% delete(gcp('nocreate'))
% poolobj = parpool('local'); 
% warning('on')
%%

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
cellp=0:(Nc-1);

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

basetime = sprintf('%14.6f',getvalue([datadir,'/suntans.dat'],'basetime'));
mtime_base = datenum(basetime,'yyyymmdd.HHMMSS');
starttime = sprintf('%14.6f',getvalue([datadir,'/suntans.dat'],'starttime'));
mtime_start = datenum(starttime,'yyyymmdd.HHMMSS');
Toffset = mtime_start-mtime_base; % days between basetime and starttime
dt = getvalue([datadir,'/suntans.dat'],'dt');
nsteps = getvalue([datadir,'/suntans.dat'],'nsteps');
% nsteps = getvalue([datadir,'/suntans.dat'],'nbc');
mtime_end = mtime_start + dt*nsteps/86400; % ending time of simulation
ntaverage =  getvalue([datadir,'/suntans.dat'],'ntaverage');
ntaveragestore =  getvalue([datadir,'/suntans.dat'],'ntaveragestore');
sponge_distance =  getvalue([datadir,'/suntans.dat'],'sponge_distance');
wave_nesting =  getvalue([datadir,'/suntans.dat'],'wave_nesting');
lowfreq_nudging =  getvalue([datadir,'/suntans.dat'],'lowfreq_nudging');
Tfilt = dt*ntaverage*ntaveragestore; % filter time window, seconds

beta =  getvalue([datadir,'/suntans.dat'],'beta');
gamma =  getvalue([datadir,'/suntans.dat'],'gamma');
rho0=1000;

%% find indices of type2 and type3 bc

type2 = mark==2; % edges with type 2
type3_edge = grad(mark==3,1); % index of cells neighboring type 3 edge
type3 = logical(0*xv);
type3(type3_edge+1)=1;

Ntype2 = sum(type2);
Ntype3 = sum(type3);

%% get outward normal vector on edges
disp('computing outward normal vectors')

r_v = 0*xv';
rn1_v = 0*xv';
rn2_v = 0*xv';

for i=1:Nc    
  r2 = (xe-xv(i)).^2+(ye-yv(i)).^2;
  [r2,indx]=min(r2(type2));
  
  xb = xe(type2);
  xb = xb(indx);
  
  yb = ye(type2);
  yb = yb(indx);
  
  rn1_v(i) = (xb-xv(i))/sqrt(r2);
  rn2_v(i) = (yb-yv(i))/sqrt(r2);
  
  r_v(i) = sqrt((xb-xv(i)).^2+(yb-yv(i)).^2);

end

% interp to edges
F = scatteredInterpolant(xv,yv,rn1_v','linear','linear');
rn1_e = F(xe,ye);

F = scatteredInterpolant(xv,yv,rn2_v','linear','linear');
rn2_e = F(xe,ye);

F = scatteredInterpolant(xv,yv,r_v','linear','linear');
r_e = F(xe,ye);

% normalize so unit vector
r2 = sqrt(rn1_e.^2+rn2_e.^2);
if isempty(r2)
    disp('r2 is empty, not interpolating')
    return
end
rn1_e = rn1_e./r2;
rn2_e = rn2_e./r2;

D_hat = exp(-4.0*r_e/sponge_distance);

clear F r2

% get index of sponge variables
% include a few extra points
spongei = D_hat >= exp(-4.1);
Ns = sum(spongei);

%% get time information

dtlow = Tfilt;    
Nbc_end = ceil(dt*nsteps/bc_dt)+2; % add extra pts on end
time = bc_dt*[-1:Nbc_end]+Toffset*86400;% time in s since basetime
Nt = length(time);

% get low time variable   
nskip =round(dtlow/bc_dt);
indx_tl = [1:nskip:(length(time)-1) length(time)];
if length(indx_tl)<4 % minimum length
    indx_tl=[1 round(length(time)/2) (length(time)-1) length(time)];
    disp('time window too small for low freq u, NTl<4, adjusted points')        
end
time_low = time(indx_tl);
Ntl =length(indx_tl);


%%
disp('writing netcdf BC file')

% ensure no open/locked files
clear mex
ncid = netcdf.create([datadir bc_filename],'CLOBBER');
netcdf.close(ncid);
% delete(['..\rundata\sun_BC.nc'])

cmode = netcdf.getConstant('NETCDF4');
cmode = bitor(cmode,netcdf.getConstant('CLASSIC_MODEL'));
ncid = netcdf.create([datadir bc_filename],cmode);

netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Created', ['Created on ' datestr(now)]);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author', '');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description', 'SUNTANS Boundary Conditions File');

% define dimensions
Nk_dimid = netcdf.defDim(ncid,'Nk',Nkmax);
Nt_dimid = netcdf.defDim(ncid,'Nt',netcdf.getConstant('NC_UNLIMITED'));
Ntl_dimid = netcdf.defDim(ncid,'Ntl',Ntl);
Ne_dimid = netcdf.defDim(ncid,'Ne',Ne);
Ns_dimid = netcdf.defDim(ncid,'Ns',Ns);
if Ntype2
    Ntype2_dimid = netcdf.defDim(ncid,'Ntype2',Ntype2);
end
if Ntype3
    Ntype3_dimid = netcdf.defDim(ncid,'Ntype3',Ntype3);
end

% define variables
varid = netcdf.defVar(ncid,'z','NC_DOUBLE',Nk_dimid);
netcdf.putAtt(ncid,varid,'long_name','Vertical grid mid-layer depth');
netcdf.putAtt(ncid,varid,'units','meters');

varid = netcdf.defVar(ncid,'time','NC_DOUBLE',Nt_dimid);
netcdf.putAtt(ncid,varid,'units','seconds since 1990-01-01 00:00:00');
netcdf.putAtt(ncid,varid,'long_name','Boundary time');
netcdf.defVarFill(ncid,varid,false,999999);

varid = netcdf.defVar(ncid,'time_low','NC_DOUBLE',Ntl_dimid);
netcdf.putAtt(ncid,varid,'units','seconds since 1990-01-01 00:00:00');
netcdf.putAtt(ncid,varid,'long_name','Boundary time');
netcdf.defVarFill(ncid,varid,false,999999);

if Ntype2
    varid = netcdf.defVar(ncid,'xe','NC_DOUBLE',Ntype2_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Easting of type-2 boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'ye','NC_DOUBLE',Ntype2_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Northing of type-2 boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'xe_all','NC_DOUBLE',Ne_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Easting of edge boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'ye_all','NC_DOUBLE',Ne_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Northing of edge boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'edgep','NC_INT',Ntype2_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Index of suntans grid edge corresponding to type-2 boundary');
    netcdf.putAtt(ncid,varid,'units','');
    
    varid = netcdf.defVar(ncid,'edgep_spg','NC_INT',Ns_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Index of suntans grid edge for all points');
    netcdf.putAtt(ncid,varid,'units','');
    
    varid = netcdf.defVar(ncid,'edgep_all','NC_INT',Ne_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Index of suntans grid edge for all points');
    netcdf.putAtt(ncid,varid,'units','');
    
    varid = netcdf.defVar(ncid,'boundary_h','NC_DOUBLE',[Ntype2_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Free-surface elevation at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'boundary_u','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Eastward velocity at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'boundary_v','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Northward velocity at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'sponge_uf','NC_DOUBLE',[Ns_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Eastward sponge velocity at edges');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
%     
%     varid = netcdf.defVar(ncid,'lowfreq_uf','NC_DOUBLE',[Ne_dimid Nk_dimid Ntl_dimid]);
%     netcdf.putAtt(ncid,varid,'long_name','Eastward low frequency velocity at edges');
%     netcdf.putAtt(ncid,varid,'units','meters second-1');
    
%     varid = netcdf.defVar(ncid,'sponge_u','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
%     netcdf.putAtt(ncid,varid,'long_name','Eastward sponge velocity at edges');
%     netcdf.putAtt(ncid,varid,'units','meters second-1');
    
%     varid = netcdf.defVar(ncid,'sponge_v','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
%     netcdf.putAtt(ncid,varid,'long_name','Northward sponge velocity at edges');
%     netcdf.putAtt(ncid,varid,'units','meters second-1');
%     
%     varid = netcdf.defVar(ncid,'lowfreq_u','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
%     netcdf.putAtt(ncid,varid,'long_name','Eastward low frequency velocity at edges');
%     netcdf.putAtt(ncid,varid,'units','meters second-1');
%     
%     varid = netcdf.defVar(ncid,'lowfreq_v','NC_DOUBLE',[Ne_dimid Nk_dimid Nt_dimid]);
%     netcdf.putAtt(ncid,varid,'long_name','Northward low frequency velocity at edges');
%     netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'boundary_w','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Vertical velocity at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'boundary_T','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Water temperature at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','degrees C');
    
    varid = netcdf.defVar(ncid,'boundary_S','NC_DOUBLE',[Ntype2_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Salinity at type-2 boundary point');
    netcdf.putAtt(ncid,varid,'units','psu');
end

if Ntype3
    varid = netcdf.defVar(ncid,'xv','NC_DOUBLE',Ntype3_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Easting of type-3 boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'yv','NC_DOUBLE',Ntype3_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Northing of type-3 boundary points');
    netcdf.putAtt(ncid,varid,'units','meters');
    
    varid = netcdf.defVar(ncid,'cellp','NC_INT',Ntype3_dimid);
    netcdf.putAtt(ncid,varid,'long_name','Index of suntans grid cell corresponding to type-3 boundary');
    netcdf.putAtt(ncid,varid,'units','');
    
    varid = netcdf.defVar(ncid,'uc','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Eastward velocity at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'vc','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Northward velocity at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'wc','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Vertical velocity at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters second-1');
    
    varid = netcdf.defVar(ncid,'T','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Temperature at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','degrees C');
    
    varid = netcdf.defVar(ncid,'S','NC_DOUBLE',[Ntype3_dimid Nk_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Salinity at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','psu');
    
    varid = netcdf.defVar(ncid,'h','NC_DOUBLE',[Ntype3_dimid Nt_dimid]);
    netcdf.putAtt(ncid,varid,'long_name','Water surface elevation at type-3 boundary point');
    netcdf.putAtt(ncid,varid,'units','meters');
    
   
end
% end define variables
netcdf.endDef(ncid);

% write grid variables
varid = netcdf.inqVarID(ncid,'z');
netcdf.putVar(ncid,varid,z_r);

if Ntype2
varid = netcdf.inqVarID(ncid,'xe');
netcdf.putVar(ncid,varid,xe(type2));

varid = netcdf.inqVarID(ncid,'ye');
netcdf.putVar(ncid,varid,ye(type2));

varid = netcdf.inqVarID(ncid,'xe_all');
netcdf.putVar(ncid,varid,xe);

varid = netcdf.inqVarID(ncid,'ye_all');
netcdf.putVar(ncid,varid,ye);

varid = netcdf.inqVarID(ncid,'edgep');
netcdf.putVar(ncid,varid,edgep(type2));

varid = netcdf.inqVarID(ncid,'edgep_spg');
netcdf.putVar(ncid,varid,edgep(spongei));

varid = netcdf.inqVarID(ncid,'edgep_all');
netcdf.putVar(ncid,varid,edgep);

varid = netcdf.inqVarID(ncid,'time_low');
netcdf.putVar(ncid, varid,0,length(time_low), time_low);
    
end

if Ntype3
varid = netcdf.inqVarID(ncid,'xv');
netcdf.putVar(ncid,varid,xv(type3));

varid = netcdf.inqVarID(ncid,'yv');
netcdf.putVar(ncid,varid,yv(type3));

varid = netcdf.inqVarID(ncid,'cellp');
netcdf.putVar(ncid,varid,cellp(type3));
end

% close file
netcdf.close(ncid);


%% compute fields
disp('computing BC variables')

if n_soliton==1 %one solitary wave bc

    if Ntype2
        % type2 boundary and sponge layer locations
        xe_type2 = xe(type2);
        xe_spongei = xe(spongei);
        ye_type2 = ye(type2);
        ye_spongei = ye(spongei);

        % free surface
        boundary_h = zeros(Ntype2,Nt);
        % velocity
        boundary_u = zeros(Ntype2,Nkmax,Nt);
        boundary_v = 0*boundary_u;
        boundary_w = 0*boundary_u;
        % salt
        boundary_S = 0*boundary_u;
        % temp
        boundary_T = 0*boundary_u;

        % sponge layer: now is no sponge layer, need to code up
        Fhat_x = getvalue([datadir '/suntans.dat'],'Fhat_x');
        Fhat_y = getvalue([datadir '/suntans.dat'],'Fhat_y');      
        F_hat = 0*xe; 
        F_hat(Fhat_x.*rn1_e + Fhat_y.*rn2_e<0)=1;
        F_hat_type2 = F_hat(type2);
        F_hat_spongei = F_hat(spongei);

        sponge_u = zeros(length(xe(spongei)),Nkmax,Nt);
        sponge_v = zeros(length(xe(spongei)),Nkmax,Nt);
        sponge_uf = zeros(length(xe(spongei)),Nkmax,Nt);

        % interpret DJL solution onto c*cumsum(dt)
        disp('loading soliton data')
        load(soliton_file);
        if DJLES_solution.wavelength/DJLES_solution.c<bc_dt*10
            warning('WARNING: time step for read boundary condition is too coarse to resolve the solitary wave!')
        end
        
        % background T at boundary
        T_bound = -(interp1(DJLES_solution.zc, DJLES_solution.density(:,1), z_r)-rho0)/gamma/rho0;
        for i=1:Ntype2
            boundary_T(i,:,:) = T_bound*ones(1,Nt);
        end   
        
        time_ = time-time(1);
        time_onesoliton = DJLES_solution.xc(end)/DJLES_solution.c;
        idx_onesoliton = find(time_<=time_onesoliton);

        % only force incoming wave at F_hat == 1
        idx_type2 = find(F_hat_type2 == 1);
        idx_spongei = find(F_hat_spongei == 1);
        xe_spongei_elements = unique(xe_spongei(idx_spongei));
        
        for i=1:Nt %length(idx_onesoliton)
            if mod(i,round(Nt/10)) == 0
                fprintf('At time iteration %d...\n',i);
            end
            x_=time_(i)*DJLES_solution.c;
            if method_sponge == 1
                % sponge layer: incoming wave
                for j=1:length(xe_spongei_elements)
                    idx_spongei_elements = find(xe_spongei(idx_spongei)==xe_spongei_elements(j));
                    x_s=xe_type2(idx_type2) - xe_spongei_elements(j);
                    x_s=x_s(1);
                    x_s_=-x_s+x_;
                    if x_s_<DJLES_solution.xc(1) || x_s_>DJLES_solution.xc(end)    
                        sponge_u(idx_spongei(idx_spongei_elements),:,i) = 0; %interp1(DJLES_solution.zc, DJLES_solution.u(:,1), z_r)';
                    elseif x_s_>=DJLES_solution.xc(1) && x_s_<=DJLES_solution.xc(end) 
                        for k=1:Nkmax
                            sponge_u(idx_spongei(idx_spongei_elements),k,i) = interp2(DJLES_solution.zc, ...
                                DJLES_solution.xc, ...
                                DJLES_solution.u', ...
                                z_r(k), ...
                                x_s_)';
                        end
                    end
                end
            end
            % apply F_hat
            sponge_u(:,:,i) = F_hat_spongei'.*sponge_u(:,:,i);
            sponge_v(:,:,i) = F_hat_spongei'.*sponge_v(:,:,i);

            % type2 boundary: incoming wave
            if x_<DJLES_solution.xc(1) || x_>DJLES_solution.xc(end)
                boundary_u(idx_type2,:,i) = 0; %ones(length(idx_type2),1)*interp1(DJLES_solution.zc, DJLES_solution.u(:,1), z_r)';
%                 boundary_w(idx_type2,:,i) = 0; %ones(length(idx_type2),1)*interp1(DJLES_solution.zc, DJLES_solution.w(:,1), z_r)';
                boundary_T(idx_type2,:,i) = -(ones(length(idx_type2),1)*interp1(DJLES_solution.zc, DJLES_solution.density(:,1), z_r)'-rho0)./gamma/rho0;           
            elseif x_>=DJLES_solution.xc(1) && x_<=DJLES_solution.xc(end) 
                for k=1:Nkmax
                    boundary_u(idx_type2,k,i) = interp2(DJLES_solution.zc, ...
                        DJLES_solution.xc, ...
                        DJLES_solution.u', ...
                        z_r(k), ...
                        x_);
%                     boundary_w(idx_type2,k,i) = interp2(DJLES_solution.zc, ...
%                         DJLES_solution.xc, ...
%                         DJLES_solution.w', ...
%                         z_r(k), ...
%                         x_);          
                    boundary_T(idx_type2,k,i) = -(interp2(DJLES_solution.zc, ...
                        DJLES_solution.xc, ...
                        DJLES_solution.density', ...
                        z_r(k), ...
                        x_)-rho0)/gamma/rho0;
                end
            end

        end

        % face fluxes
        for i=1:Nt
            sponge_uf(:,:,i) = N1(spongei,:).*sponge_u(:,:,i) + N2(spongei,:).*sponge_v(:,:,i);
        end

    end

    
    if Ntype3
        % eta
        h = zeros(Ntype3,Nt);
        
        % velocity
        uc = zeros(Ntype3,Nkmax,Nt);
        vc = uc;
        wc = uc;
        
        % salt
        S = uc;
        
        % temp
        disp('loading soliton data')
        load(soliton_file);

        T_bound = -(interp1(DJLES_solution.zc, DJLES_solution.density(:,1), z_r)-rho0)/gamma/rho0;
        for i=1:Ntype3
            T(i,:,:) = T_bound*ones(1,Nt);
        end   
        
    end

else % a train of solitary waves bc

    if Ntype2
        % type2 boundary and sponge layer locations
        xe_type2 = xe(type2);
        xe_spongei = xe(spongei);
        ye_type2 = ye(type2);
        ye_spongei = ye(spongei);

        % free surface
        boundary_h = zeros(Ntype2,Nt);
        % velocity
        boundary_u = zeros(Ntype2,Nkmax,Nt);
        boundary_v = 0*boundary_u;
        boundary_w = 0*boundary_u;
        % salt
        boundary_S = 0*boundary_u;
        % temp
        boundary_T = 0*boundary_u;

        % sponge layer: now is no sponge layer, need to code up
        Fhat_x = getvalue([datadir '/suntans.dat'],'Fhat_x');
        Fhat_y = getvalue([datadir '/suntans.dat'],'Fhat_y');      
        F_hat = 0*xe; 
        F_hat(Fhat_x.*rn1_e + Fhat_y.*rn2_e<0)=1;
        F_hat_type2 = F_hat(type2);
        F_hat_spongei = F_hat(spongei);

        sponge_u = zeros(length(xe(spongei)),Nkmax,Nt);
        sponge_v = zeros(length(xe(spongei)),Nkmax,Nt);
        sponge_uf = zeros(length(xe(spongei)),Nkmax,Nt);

        % interpret DJL solution onto c*cumsum(dt)
        disp('loading soliton data')
        load(soliton_file);
        if DJLES_solution.wavelength/DJLES_solution.c<bc_dt*10
            warning('WARNING: time step for read boundary condition is too coarse to resolve the solitary wave!')
        end
        
        % background T at boundary
        T_bound = -(interp1(DJLES_solution.zc, DJLES_solution.density(:,1), z_r)-rho0)/gamma/rho0;
        for i=1:Ntype2
            boundary_T(i,:,:) = T_bound*ones(1,Nt);
        end   
        
        time_ = time-time(1);
        time_nsoliton = n_soliton*DJLES_solution.xc(end)/DJLES_solution.c+...
                        (n_soliton-1)*soliton_interval;
        idx_nsoliton = find(time_<=time_nsoliton);

        % only force incoming wave at F_hat == 1
        idx_type2 = find(F_hat_type2 == 1);
        idx_spongei = find(F_hat_spongei == 1);
        xe_spongei_elements = unique(xe_spongei(idx_spongei));

        % propagation range of the train of solitons
        xc_starts = zeros(n_soliton,1);
        xc_ends = zeros(n_soliton,1);
        for ii=1:n_soliton
            xc_starts(ii) = DJLES_solution.xc(1)+(ii-1)*soliton_interval*DJLES_solution.c;
            xc_ends(ii) = DJLES_solution.xc(end)+(ii-1)*soliton_interval*DJLES_solution.c;
        end
        if xc_ends(end)>max(xv)
            warning('WARNING: the domain is too short to put in all the sequencial solitary waves')
        end

        for i=1:Nt %length(idx_onesoliton)
            if mod(i,round(Nt/10)) == 0
                fprintf('At time iteration %d...\n',i);
                toc
            end
            x_=time_(i)*DJLES_solution.c;

            if method_sponge == 1
                % sponge layer: incoming wave
                for j=1:length(xe_spongei_elements)
                    idx_spongei_elements = find(xe_spongei(idx_spongei)==xe_spongei_elements(j));
                    x_s=xe_type2(idx_type2) - xe_spongei_elements(j);
                    x_s=x_s(1);
                    x_s_=-x_s+x_;
                    if x_s_<xc_starts(1) || x_s_>xc_ends(end)    
                        sponge_u(idx_spongei(idx_spongei_elements),:,i) = 0; %interp1(DJLES_solution.zc, DJLES_solution.u(:,1), z_r)';
                    end
                    for ii=1:n_soliton-1
                        if x_s_>xc_ends(ii) && x_s_<xc_starts(ii+1)
                            sponge_u(idx_spongei(idx_spongei_elements),:,i) = 0;
                        end
                    end
                    for ii=1:n_soliton
                        if x_s_>=xc_starts(ii) && x_s_<=xc_ends(ii) 
                            for k=1:Nkmax
                                sponge_u(idx_spongei(idx_spongei_elements),k,i) = interp2(DJLES_solution.zc, ...
                                    DJLES_solution.xc, ...
                                    DJLES_solution.u', ...
                                    z_r(k), ...
                                    x_s_-xc_starts(ii)+xc_starts(1))';
                            end   
                        end
                    end
                end
            end
            % apply F_hat
            sponge_u(:,:,i) = F_hat_spongei'.*sponge_u(:,:,i);
            sponge_v(:,:,i) = F_hat_spongei'.*sponge_v(:,:,i);

            % type2 boundary: incoming wave
            if x_<xc_starts(1) || x_>xc_ends(end)
                boundary_u(idx_type2,:,i) = 0; %ones(length(idx_type2),1)*interp1(DJLES_solution.zc, DJLES_solution.u(:,1), z_r)';
%                 boundary_w(idx_type2,:,i) = 0; %ones(length(idx_type2),1)*interp1(DJLES_solution.zc, DJLES_solution.w(:,1), z_r)';
                boundary_T(idx_type2,:,i) = -(ones(length(idx_type2),1)*interp1(DJLES_solution.zc, DJLES_solution.density(:,1), z_r)'-rho0)./gamma/rho0;           
            end
            for ii=1:n_soliton-1
                if x_>xc_ends(ii) && x_<xc_starts(ii+1)
                    boundary_u(idx_type2,:,i) = 0;
                    boundary_T(idx_type2,:,i) = -(ones(length(idx_type2),1)*interp1(DJLES_solution.zc, DJLES_solution.density(:,1), z_r)'-rho0)./gamma/rho0;    
                end
            end   
            for ii=1:n_soliton
                if x_>=xc_starts(ii) && x_<=xc_ends(ii) 
                   for k=1:Nkmax
                        boundary_u(idx_type2,k,i) = interp2(DJLES_solution.zc, ...
                            DJLES_solution.xc, ...
                            DJLES_solution.u', ...
                            z_r(k), ...
                            x_-xc_starts(ii)+xc_starts(1));
                        boundary_T(idx_type2,k,i) = -(interp2(DJLES_solution.zc, ...
                            DJLES_solution.xc, ...
                            DJLES_solution.density', ...
                            z_r(k), ...
                            x_-xc_starts(ii)+xc_starts(1))-rho0)/gamma/rho0;
                    end   
                end
            end
        end

        % face fluxes
        for i=1:Nt
            sponge_uf(:,:,i) = N1(spongei,:).*sponge_u(:,:,i) + N2(spongei,:).*sponge_v(:,:,i);
        end

    end

    
    if Ntype3
        % eta
        h = zeros(Ntype3,Nt);
        
        % velocity
        uc = zeros(Ntype3,Nkmax,Nt);
        vc = uc;
        wc = uc;
        
        % salt
        S = uc;
        
        % temp
        disp('need to code this')
        
    end

end
    
%% write variables
disp('writing variables to netcdf BC file')

ncid = netcdf.open([datadir bc_filename],'WRITE');

varid = netcdf.inqVarID(ncid,'time');
netcdf.putVar(ncid, varid,0,length(time), time);

if Ntype2

    varid = netcdf.inqVarID(ncid,'boundary_h');
    netcdf.putVar(ncid,varid,boundary_h);

    varid = netcdf.inqVarID(ncid,'boundary_u');
    netcdf.putVar(ncid,varid,boundary_u);

    varid = netcdf.inqVarID(ncid,'boundary_v');
    netcdf.putVar(ncid,varid,boundary_v);

    varid = netcdf.inqVarID(ncid,'sponge_uf');
    netcdf.putVar(ncid,varid,sponge_uf);
% 
%         varid = netcdf.inqVarID(ncid,'lowfreq_uf');
%         netcdf.putVar(ncid,varid,lowfreq_uf);

    varid = netcdf.inqVarID(ncid,'boundary_w');
    netcdf.putVar(ncid,varid,boundary_w);

    varid = netcdf.inqVarID(ncid,'boundary_T');
    netcdf.putVar(ncid,varid,boundary_T);

    varid = netcdf.inqVarID(ncid,'boundary_S');
    netcdf.putVar(ncid,varid,boundary_S);    
end

if Ntype3 
    varid = netcdf.inqVarID(ncid,'h');
    netcdf.putVar(ncid,varid,h);

    varid = netcdf.inqVarID(ncid,'uc');
    netcdf.putVar(ncid,varid,uc);

    varid = netcdf.inqVarID(ncid,'vc');
    netcdf.putVar(ncid,varid,vc);

    varid = netcdf.inqVarID(ncid,'wc');
    netcdf.putVar(ncid,varid,wc);

    varid = netcdf.inqVarID(ncid,'T');
    netcdf.putVar(ncid,varid,T);

    varid = netcdf.inqVarID(ncid,'S');
    netcdf.putVar(ncid,varid,S);    
end

netcdf.close(ncid);




%% plot bc

figure
plot(xe, ye, 'k.');
hold on
plot(xe_spongei, ye_spongei, 'r.');
hold on
plot(xe_type2, ye_type2, 'b.')        
hold on
plot(xe_spongei(idx_spongei), ye_spongei(idx_spongei), 'r+');
hold on
plot(xe_type2(idx_type2), ye_type2(idx_type2), 'b+')
hold on
plot(xe,1e+5*D_hat,'.')
print -djpeg -r300 figure_BC

figure
subplot(211)
contourf(time_,z_r,squeeze(boundary_T(2,:,:)));
title('temp')
xlabel('t (s)')
ylabel('z (m)')

subplot(212)
contourf(time_,z_r,squeeze(boundary_u(2,:,:)));
title('u')
xlabel('t (s)')
ylabel('z (m)')

print -djpeg -r300 figure_BC_input

figure
subplot(211)
j=1;
contourf(time_,z_r,squeeze(sponge_uf(j,:,:)));
title(strcat('sponge_uf at x=', string(xe(idx_spongei(j))), ...
   'y=', string(ye(idx_spongei(j)))))
xlabel('t (s)')
ylabel('z (m)')

subplot(212)
j=12;
contourf(time_,z_r,squeeze(sponge_uf(j,:,:)));
title(strcat('sponge_uf at x=', string(xe(idx_spongei(j))), ...
   'y=', string(ye(idx_spongei(j)))))
xlabel('t (s)')
ylabel('z (m)')

print -djpeg -r300 figure_sponge_input
   

    
end

