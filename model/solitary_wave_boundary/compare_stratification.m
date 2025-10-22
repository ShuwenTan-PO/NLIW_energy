rho0 = 1000;
% gamma taken from suntans.dat
gamma = 2.1e-4;

%% winter
% Rogers 2022 density profile
% parameters taken from initialization.c
a1 = 14.7098; a2 = 7.4166; a3 = -59.1951; a4 = 173.7317;
% idealized rho profile taken from state.c
frho=@(z,gamma,rho0,a1,a2,a3,a4) -gamma*rho0*(a1*exp(-(-z+a3)/a4)+a2)+rho0;
frhoz=@(z,gamma,rho0,a1,a3,a4) -gamma*rho0*a1*exp(-(-z+a3)/a4)/a4;
% below is wrong, but somehow it gives DJLE solutions
% frho=@(z,gamma,rho0,a1,a2,a3,a4) -gamma*rho0*(a1*exp(-(-z+a3)/a4)+a2)+rho0;
% frhoz=@(z,gamma,rho0,a1,a3,a4) -gamma*rho0*a1*exp(-(-z+a3)/a4)/a4;
rho  = @(z) frho(z, gamma,rho0,a1,a2,a3,a4);
rhoz = @(z) frhoz(z, gamma,rho0,a1,a3,a4);

figure(1)
subplot(121)
plot(rho([0:-1:-600]),[0:-1:-600],'k','linewidth',3);
ax = gca;
ax.FontSize = 14; 
xlabel('\rho [kgm^{-3}]','FontSize',14)
ylabel('z [m]','FontSize',14)

subplot(122)
plot(10*rhoz([0:-1:-600])/1000,[0:-1:-600],'k','linewidth',3)
ax = gca;
ax.FontSize = 14; 
xlabel('N^2 [s^{-2}]','FontSize',14)
ylabel('z [m]','FontSize',14)

%% winter
% Rogers 2022 density profile
% parameters taken from initialization.c
a1 = 23.36; a2 = 3.13; a3 = -44.12; a4 = 293.12;
% idealized rho profile taken from state.c
frho=@(z,gamma,rho0,a1,a2,a3,a4) -gamma*rho0*(a1*exp(-(-z+a3)/a4)+a2)+rho0;
frhoz=@(z,gamma,rho0,a1,a3,a4) -gamma*rho0*a1*exp(-(-z+a3)/a4)/a4;
% below is wrong, but somehow it gives DJLE solutions
% frho=@(z,gamma,rho0,a1,a2,a3,a4) -gamma*rho0*(a1*exp(-(-z+a3)/a4)+a2)+rho0;
% frhoz=@(z,gamma,rho0,a1,a3,a4) -gamma*rho0*a1*exp(-(-z+a3)/a4)/a4;
rho  = @(z) frho(z, gamma,rho0,a1,a2,a3,a4);
rhoz = @(z) frhoz(z, gamma,rho0,a1,a3,a4);

figure(1)
subplot(121)
hold on
plot(rho([0:-1:-600]),[0:-1:-600],'r','linewidth',3);
legend('Winter','Summer')

subplot(122)
hold on
plot(10*rhoz([0:-1:-600])/1000,[0:-1:-600],'r','linewidth',3)

print -djpeg -r300 den_winter_summer
