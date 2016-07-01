%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Calculate 2d Barotropic Rossby ray paths over specified background 
%%  fields.
%
%   Based on Karoly, 1983
%   Varables:
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all; clear all

addpath matlab

rad     = 6.371e6; % radius of sphere having same volume as Earth (m)
e_omega = 7.292e-5; % rotation rate of Earth (rad/s)
dtr     = pi/180;
rtd     = 180/pi

day = 24*60*60; % in s
dt  = 60*60;  %
integration_time=10*day;
Nsteps = round(integration_time/dt);

% Periods ('Inf' for quasi-stationary waves) 
Periods=[ Inf 50 20 ]*day;
freq=2*pi./Periods;
Nfrequencies=length(freq);

% Initial k wave number:
k_wavenumbers=[1 2 3 4 5 6];

% Starting point of ray (location)

lon0 = [ 0   0]  ; %deg.E
lat0 = [30 -30] ; %deg.N

disp('-----------------------------------')

disp(['k wavenumbers: ' num2str(k_wavenumbers(:)') '']) ;
disp(['Periods: ' num2str(Periods(:)'/day) ' days']) ;
disp(['Integration time: ', num2str(integration_time'/day),' days']) ;
disp(['dt: ' num2str(dt'/60) ' min']) ;
disp(['Stating point: ' num2str(lon0) ' E   ' num2str(lat0) ' N']) ;

disp('-----------------------------------')

% Read data



%------------------------------------------------
%%%  !!! add time from the start of inegration (s) to the output !!!
%fout = sprintf('../output/ray_loc%dE_%dN_period%d_k%d_root%d'...
%                         ,lon0(),lat0(),period,kk,RR);
