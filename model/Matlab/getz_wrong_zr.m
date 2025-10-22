function z = getz(dz, rogers)
  
  N = length(dz);
  if (nargin==2)
      dzh = 0.5*(dz(1:end-1)+dz(2:end));
      % Rogers' original way to reconstruct z from dz
  else
      dzh = dz(2:end);
      % S. Tan's adaption to that (this is now the default), z is the same
      % as z_r in the netcdf file, the depth of layer centers
  end

  z = zeros(N,1);
  z(1) = -dz(1)/2;
  z(2:N) = z(1)-cumsum(dzh);