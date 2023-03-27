clear all; close all; clc

conditions = "earth";

if conditions == "earth"
% Constants For Earth
mu = 398600; % mu in km^3/s^2
radius = 6370; % in km
j2 = 0.00108;
%Given Initial Conditions
peri_altitude = 600; % km 
period_planet = 23*60*60 + 56*60 + 4; %(seconds) for 1 day
m = 3; % number of sattelite orbits per Earth day

elseif conditions == "mars"
% Constants For mars
mu = 42820; % mu in km^3/s^2
radius = 3390; % in km
j2 = 0.00196;
peri_altitude = 400; % km
period_planet = 24*60*60 + 39*60 + 35;
m = 1; % number of sattelite orbits per Mars day
end

deg = 180/pi;

% Solution 
i = asin(0.8^0.5);
a = (mu*(period_planet/(2*pi))^2)^(1/3); % semi-major axis (km) ....from Period = 2π*sqrt(a^3 / mu)
rp = peri_altitude + radius; %(km) perigee

e = 0;
itr = 0;
g = 1;
while (e < 1) && (e >= 0) && abs(g) > 1e-9
e = 1 - rp/a; 
n = sqrt(mu / (a)^3);
period_current = 2*pi*sqrt(a^3/ mu); %(seconds) sattelite orbits earth m times per day

omega_star = 2*pi / period_planet; %rad/s Earth's rotation rate
omega_dot = -3/2 * (n * j2) * (radius/a)^2 * cos(i)/(1 - e^2)^2; %rad/s
n_s = m*(omega_star - omega_dot); % rad/s 
g = n - n_s;    
ra = a*(1+e); %(km) apogee

% fprintf('itr = %d \n',itr)
% fprintf('e = %d \n',e)
% fprintf('a = %d \n',a)
% fprintf('Nodal precession d\x03A9/dt  = %E deg/s or %.4f deg/day\n',omega_dot*deg, omega_dot*deg*60*60*24);
% Update
a = a - 0.1;
itr = itr + 1;
end
a = a + 0.1;

fprintf('Semi-major axis a\t = %.4f km\n',a);
fprintf('Eccentricity e\t\t = %.4f \n',e);
fprintf('Inclination i\t\t = %.4f deg\n',i*deg);
fprintf('Nodal precession d\x03A9/dt  = %E deg/s or %.4f deg/day\n',omega_dot*deg, omega_dot*deg*60*60*24);
