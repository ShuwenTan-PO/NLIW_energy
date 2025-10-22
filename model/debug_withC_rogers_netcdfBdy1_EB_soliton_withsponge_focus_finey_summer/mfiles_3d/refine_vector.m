%%
function [ xnew ] =  refine_vector(xv,dxmax,dxmin,xc,rad,tol)
% Nx = round(L/dxmax);
% x0 = 1;
% xv = linspace(x0,x0+Nx*dxmax,Nx);
% rad = 0.4*L;
% xc = 0.3*L;
    xnew=xv(1);
    dxmaxi = (max(xv)-min(xv))/length(xv);
    dxratio = dxmin/dxmax;
    k=1;
    while abs(xnew(end)-xv(end))>tol && k<1000 % iterate until length agrees

        disp(['k=' num2str(k) ', dx error =' num2str(xnew(end)-xv(end))])
        xnew=xv(1);
        for kk=1:length(xv)-1
            dx = dxmaxi -(dxmaxi-dxratio*dxmaxi)* exp(-(xnew-xc).^2/(2*rad^2));
            xnew(kk+1) = xnew(kk)+dx(kk);
        end
        dxmaxi = dxmaxi * (1 + 0.10*(xv(end)-xv(1) - xnew(end)+xnew(1))/(xv(end)-xv(1)));
        k=k+1;
    end

end
