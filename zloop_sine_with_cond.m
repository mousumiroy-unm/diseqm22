% Script to solve nondimensional Eqns 5 and 6 in the Solid Earth paper:
%
% Roy, M., Assessing the role of thermal disequilibrium in the evolution of 
% the lithosphere-asthenosphere boundary: An idealized model of heat exchange during 
% channelized melt-transport, Solid Earth, se-2022-13, 2022
% 
% This code loops through z,d, etc., using the with_cond_all_zloop_sinefn.m code, 
% to solve the coupled system, storing the results in different folders with 
% parameters indicated in the folder names
%
% Author: Mousumi Roy
% University of New Mexico
% May 2022
%
% Here we are using a sinusoidal pulse in the inlet temperature, but one
% can change that in the with_cond_all_zloop_sinefn.m -- and rename, if
% needed
%
% Notes: the stability criterion is modified after adding in (axial) diffusion.
% When the axial diffusion terms are large, the fluid and solid are nearly 
% always in eqm… we will get oscillations because of small differences 
% in Tf-Ts from the linear driving term.  Things behave nicely when 
% Df is <= 1.7 and Ds <=80 with the canned matlab gradient and del2 
% functions --> for my paper, this implies channel spacing d>=100 m and 
% z>=0.049 or so...So, for stability, an if-statement finds when we have 
% too large Df and Ds (threshold depends on d, z) … simpler than an 
% analytic derivation of the CFL-type of condition using eigenmodes.
% 
%
clear all 
close all

n     = 5000;    % grid resolution
cs    = 4.125e6; % in J/(K m^3) % from values in Table 1, with pure phases
cf    = 3.920e6; % in J/(K m^3)
beta = 6;
Nu = 12; 
kf = 1; % in W/mK
ks = 2.5; % in W/mK
condratio = ks/kf;
v = 1./(60*60*24*365); % 1 m/yr = 1/(60*60*24*365 im m/s
fac1 = kf/(v*v*cf*cf);

darr = [1000]; %[500]; % in m
 
zarr = [0.129];
ztophi = @(z)z*cs/(cf+z*cs)
phi = arrayfun(ztophi, zarr);

for i = 1:length(zarr)
 for j = 1:length(darr)    

     Lx   = 100;
     Lx   = 5000;
     tmax = 15000;
     dx   = Lx/n;
     dt   = 1.0e-4*dx;          % small, for stability 

     denom = 1/(kf*Nu) + 1/(beta*ks);
     term1 = 2*(1-phi(i))/denom;
     term2 = 1./(darr(j)*darr(j));
     bigK  = term1.*term2; % heat transfer coefficient
     Df    = fac1*bigK./phi(i)
     Ds    = condratio*Df*(1-phi(i))/phi(i)
     startind = 40;

     period   = [20];% period for sine, or width of pulse for tanh, in dimensionless time
      
      close all
      for jj = 1:length(period)
         Lx   = 1000;
         tmax = 150;
         %tmax = 2.1*period(jj)+Lx*zarr(i);
         w0   = period(jj)*0.1; %sharpness of pulse edges if using tanh
         %Df=0; % no diffusion (debugging)
         %Ds=0;
         with_cond_all_zloop_sinefn(zarr(i), Lx, tmax, dx, dt, startind, period(jj),w0, Df, Ds, darr(j));
         close all
      end
end
end