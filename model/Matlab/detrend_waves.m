function [eta_bar,eta_tide,eta_prime] = detrend_waves(eta,mtime,Tfilt,xv,yv,dx,dxFilt,vardx,method)

if vardx
   xl =min( min(xv)):dx:max(max(xv));
   yl = min( min(yv)):dx:max(max(yv));
   [Yl,Xl]=meshgrid(yl,xl);
end

eta_bar = nan+eta;
eta_prime = nan+eta;
eta_tide = nan+eta;
for i=1:length(mtime)
    
    indx = mtime >= mtime(i)-Tfilt/2 &...
           mtime(i) + Tfilt/2 >= mtime;
       
%     indx = mtime >= mtime(i)-Tfilt &...
%            mtime(i) >= mtime;
    eta_bar(:,:,i) = nanmean(eta(:,:,indx),3); 
    
    if strcmp(method,'poly')
        %     polynomial fit
        f=eta(:,:,i)-eta_bar(:,:,i);
        indx_good = ~isnan(eta(:,:,i));
        f = reshape(f(indx_good),[],1);
        x=reshape(xv(indx_good),[],1);
        y=reshape(yv(indx_good),[],1);
        p = polyFit2D(f,x,y,1,1);  
        eta_tide(:,:,i) = polyVal2D(p,xv,yv,1,1);
    
    elseif strcmp(method,'linear')
    
    % linear fit
        ydata=eta(:,:,i)-eta_bar(:,:,i);
        indx_good = ~isnan(eta(:,:,i));
        ydata = reshape(ydata(indx_good),[],1);
        xdata=reshape(xv(indx_good),[],1);
        xdata2=reshape(yv(indx_good),[],1);
        [px]=polyfit(xdata,ydata,1);
        [py]=polyfit(xdata2,ydata,1);
        eta_tide(:,:,i) = px(1)*xv + py(1)*yv + px(2);

    elseif strcmp(method,'imfilter')
    
    %  imfilter   
        filtersize = round(dxFilt/dx);
        h = ones(filtersize,filtersize)/filtersize^2;
        ydata=eta(:,:,i)-eta_bar(:,:,i);
        ydata =naninterp(ydata);

        if vardx % stretched grid
            % interp to regular grid
           F = scatteredInterpolant(reshape(xv,[],1),...
               reshape(yv,[],1),reshape(ydata,[],1));
           ydata=F(Xl,Yl);
           % filter
           ydata = imfilter(ydata,h);
           % transform back
           F = scatteredInterpolant(reshape(Xl,[],1),...
                reshape(Yl,[],1),reshape(ydata,[],1));
           eta_tide(:,:,i) = F(xv,yv);
           
           trims = round(filtersize/4);
           eta_tide(1:trims,:,i)=nan;
           eta_tide((end-trims):end,:,i)=nan;
           eta_tide(:,1:trims,i)=nan;
           eta_tide(:,(end-trims):end,i)=nan;

        else %regular grid
           eta_tide(:,:,i) = imfilter(ydata,h);
           trims = round(filtersize/4);
           eta_tide(1:trims,:,i)=nan;
           eta_tide((end-trims):end,:,i)=nan;
           eta_tide(:,1:trims,i)=nan;
           eta_tide(:,(end-trims):end,i)=nan;
        end

    elseif strcmp(method,'gaussian')
    % %   gaussian filter  
%         eta_tide(:,:,i) =imgaussfilt(ydata, 5);
        filtersize = round(dxFilt/dx);
        h = ones(filtersize,filtersize)/filtersize^2;
        ydata=eta(:,:,i);
%         ydata=eta(:,:,i)-eta_bar(:,:,i);
        ydata =naninterp(ydata);

        if vardx % stretched grid
            % interp to regular grid
           F = scatteredInterpolant(reshape(xv,[],1),...
               reshape(yv,[],1),reshape(ydata,[],1));
           ydata=F(Xl,Yl);
           % filter
           ydata = imgaussfilt(ydata,filtersize);
           % transform back
           F = scatteredInterpolant(reshape(Xl,[],1),...
                reshape(Yl,[],1),reshape(ydata,[],1));
           eta_tide(:,:,i) = F(xv,yv);
           
%            trims = round(filtersize/4);
%            eta_tide(1:trims,:,i)=nan;
%            eta_tide((end-trims):end,:,i)=nan;
%            eta_tide(:,1:trims,i)=nan;
%            eta_tide(:,(end-trims):end,i)=nan;

        else %regular grid
           eta_tide(:,:,i) = imgaussfilt(ydata,filtersize);
%            trims = round(filtersize/4);
%            eta_tide(1:trims,:,i)=nan;
%            eta_tide((end-trims):end,:,i)=nan;
%            eta_tide(:,1:trims,i)=nan;
%            eta_tide(:,(end-trims):end,i)=nan;
        end
    end
    eta_prime(:,:,i) = eta(:,:,i) -eta_tide(:,:,i);
%     eta_prime(:,:,i) = eta(:,:,i) - eta_bar(:,:,i)-eta_tide(:,:,i);
    % demean eta_prime and add to eta_tide
%     resid = nanmean(nanmean(eta_prime(:,:,i)));
%     eta_tide(:,:,i) = eta_tide(:,:,i) + resid;
%     eta_prime(:,:,i) = eta_prime(:,:,i)-resid;
%     eta_prime(:,:,i) = eta_prime(:,:,i) - nanmean(nanmean(eta_prime(:,:,i)));
    
end
   


end