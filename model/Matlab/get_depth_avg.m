function [U] = get_depth_avg(u,dz,dir,mask)
% dz is vectorized cut cell (nan below bottom)
% dir is direction of depth
% mask is 2d matrix, 0= open, nan=mask
% U = nansum(u.*dz,dir)./nansum(dz,dir);
N = size(u);
if nargin < 4
    mask = zeros(N(1),N(2));
end
U = nansum(u.*dz,dir)./nansum(dz,dir) + mask;
% mask = nansum(dz,dir) == 0;

end