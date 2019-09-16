%Test of the multilayer function. The script calculates and plots the
%reflectance of a Mo/Si multilayer with a Ru capping layer and 2
%intermixing MoSi2 layers and a Si substrate.

parameters = struct;
parameters.n0 = [1;0.886358-1i*0.01707];

parameters.nsubs = 0.99932-1i*0.00183;
parameters.n = [0.99932-1i*0.00183,...
    0.9693-1i*0.00433,...
    0.92108-1i*0.00644,...
    0.9693-1i*0.00433];
parameters.h = [2.506,0.802,1.904,1.844];
parameters.h0 = [0,2];
parameters.hsubs = 1000;
parameters.N = 40;
parameters.theta0 = (0:0.1:25)*pi/180;
parameters.wavelength = 13.5;

close
hh = parameters.h;
reflectivity = Multilayer(parameters);

plot(parameters.theta0*180/pi,reflectivity.rs,...
    parameters.theta0*180/pi,reflectivity.rp,'linewidth',2)
grid;
legend({'Polarization s','Polarization p'})
xlabel('Angle of incidence [deg]')
ylabel('Reflectance')
title('Reflectivity of standard EUV multilayer stack at \lambda = 13.5 nm')
set(gca,'fontsize',14)