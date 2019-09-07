function reflectivity = Multilayer(parameters)
%Multilayer02 Reflectivity of a multilayer mirror.
%   Reflectivity = Multilayer02(parameters) calculates the reflectivity of
%   a multilayer mirror defined by the values specified in the
%   structure 'parameters'.
%   The structure parameters needs the following fields:
%   wavelength: wavelength of the incoming light [nm]
%   n0: array of the complex index of refraction of the propagating medium 
%   and the capping layer(s).
%   n: array of the complex refraction index of each of the layers
%   constituting the multilayer period unit.
%   nsubs: complex index of refraction of the mirror substrate.
%   h: array of the heigth values of each of the layers constituting the
%   multilayer period unit.
%   h0: array of the heigth values of the propagation medium (ingnored) 
%   and of each of the capping layers.
%   N: number of repetitions of the periodic layers.
%   theta0: array of the incidence angles [rad].
%   OUTPUT
%   reflectivity.rp: reflectance for p-polarized light
%   reflectivity.tp: transmittance for p-polarized light
%   reflectivity.rs: reflectance for s-polarized light
%   reflectivity.ts: transmittance for s-polarized light

wavelength = parameters.wavelength;
n0 = parameters.n0;
n = parameters.n;
nsubs = parameters.nsubs;
n = n(:);
n0 = n0(:);
nsubs = nsubs(:);
h = parameters.h;
h0 = parameters.h0;
hsubs = parameters.hsubs;
N = parameters.N;
nv = zeros(length(n0)+length(n)*N+length(nsubs),1);
hv = nv;

nv(1:length(n0),1) = n0;
hv(1:length(h0),1) = h0;
for m = 1:N
    nv((length(n0)+1+(m-1)*length(n)):(length(n0)+m*length(n))) = n;
    hv((length(h0)+1+(m-1)*length(h)):(length(h0)+m*length(h))) = h;
end
nv(end) = nsubs;
hv(end) = hsubs;

theta0 = parameters.theta0;
r = zeros(size(theta0));
t = zeros(size(theta0));
nr = (nv);
for nn = 1:length(theta0)
    
    %%Calculate refraction angles and phase shifts (polarization S)
    k0 = 2*pi/wavelength;
    Mn = eye(2);
    
    theta = zeros(size(hv));
    theta(1) = theta0(nn);
    for cnt = 2:length(nv)
        theta(cnt)=asin(sin(theta(cnt-1))*nr(cnt-1)/nr(cnt));
        hh=nv(cnt)*hv(cnt)*cos(theta(cnt));
        Y = nv(cnt-1)*nv(cnt)*cos(theta(cnt));
        Mn = Mn*[cos(k0*hh) 1i*sin(k0*hh)/Y;...
            1i*sin(k0*hh)*Y cos(k0*hh)];
    end
    p1 = nv(1)*cos(theta(1));
    pl = nv(end)*cos(theta(end));
    t(nn) = 2*p1/((Mn(1,1) + Mn(1,2)*pl)*p1 + (Mn(2,1) + Mn(2,2)*pl));
    r(nn) = ((Mn(1,1) + Mn(1,2)*pl)*p1 - (Mn(2,1) + Mn(2,2)*pl))/...
        ((Mn(1,1) + Mn(1,2)*pl)*p1 + (Mn(2,1) + Mn(2,2)*pl));
end
reflectivity.rs = abs(r).^2;
reflectivity.ts = abs(t).^2;

for nn = 1:length(theta0)
    
    %%Calculate refraction angles and phase shifts (polarization P)
    k0 = 2*pi/wavelength;
    Mn = eye(2);
    
    theta = zeros(size(hv));
    theta(1) = theta0(nn);
    for cnt = 2:length(nv)
        theta(cnt)=asin(sin(theta(cnt-1))*nr(cnt-1)/nr(cnt));
        hh=nv(cnt)*hv(cnt)*cos(theta(cnt));
        Y = nv(cnt-1)*nv(cnt)/cos(theta(cnt));
        Mn = Mn*[cos(k0*hh) 1i*sin(k0*hh)/Y;...
            1i*sin(k0*hh)*Y cos(k0*hh)];
    end
    p1 = nv(1)*cos(theta(1));
    pl = nv(end)*cos(theta(end));
    t(nn) = 2*p1/((Mn(1,1) + Mn(1,2)*pl)*p1 + (Mn(2,1) + Mn(2,2)*pl));
    r(nn) = ((Mn(1,1) + Mn(1,2)*pl)*p1 - (Mn(2,1) + Mn(2,2)*pl))/...
        ((Mn(1,1) + Mn(1,2)*pl)*p1 + (Mn(2,1) + Mn(2,2)*pl));
end
reflectivity.rp = abs(r).^2;
reflectivity.tp = abs(t).^2;


