function [unew] = interp_grid(x,y,u,xnew,ynew)
% x(2d),y(2d),u(2d,3d,4d) are the old coordinates
% xnew(2d), ynew(2d) are the new coordinates

Nx = size(ynew,1);
Ny = size(ynew,2);
if length(size(u))==2 %2D
    F = griddedInterpolant(x,y,u);
    unew = F(xnew,ynew);
elseif length(size(u))==3 %3D
    Nk = size(u,3);
    unew = zeros(Nx,Ny,Nk);
    for i = 1:Nk
        F = griddedInterpolant(x,y,u(:,:,i));
        unew(:,:,i) = F(xnew,ynew);
    end
elseif length(size(u))==4 %2D
    Nk = size(u,3);
    Nt = size(u,4);
    unew = zeros(Nx,Ny,Nk,Nt);
    for i = 1:Nk
        for j=1:Nt
            F = griddedInterpolant(x,y,u(:,:,i,j));
            unew(:,:,i,j) = F(xnew,ynew);
        end
    end
else
    disp('error: variable is not 2D, 3D, or 4D')
    return
end

end