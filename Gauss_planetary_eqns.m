
clear all; clear;clc
% This function solves Example 10.6 by using MATLAB’s ode45 to numerically
% integrate Equations 10.84 (the Gauss planetary equations) to determine
% the J2 perturbation of the orbital elements.
%
% User M-functions required: None
% User subfunctions required: rates
% ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
%...Preliminaries:

%...Conversion factors:
hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds
deg = pi/180; %Degrees to radians
%...Constants:
mu = 398600; %Gravitational parameter (km^3/s^2)
RE = 6370; %Earth’s radius (km)
J2 = 0.00108; %Earth’s J2
%...Initial orbital parameters (given):
RA0 = 90*deg; %Right ascencion of the node (radians)
i0 = 1.10654; %Inclination (radians)
w0 = 5*deg; %Argument of perigee (radians)
e0 = 0.74; %eccentricity
a0 = 26600; %Semimajor axis (km)
M = 10*deg; %Mean anomaly (radians)
%...Initial orbital parameters (inferred):
T0 = 2*pi/sqrt(mu)*a0^1.5; %Period (s)
h0 = sqrt(mu*a0*(1 - e0^2)); %angular momentrum (km^2/s)


% Findnig Eccentric anomaly (initial)
E = M; %Initial guess for eccentric anomaly
g = 1;
itr = 0;
while abs(g) > 1e-15  
g = E-e0*sin(E) - M; % rad/s (modified Kepler's equation)
dgdE = 1-e0*cos(E);
E_new = E - g/dgdE;
% Update
E = E_new;
itr = itr + 1;
end

% Find True Anomaly from Eccentric Anomaly
 f1 = (2*atan(sqrt((1 + e0) / (1 - e0)) * tan(E/2))); 
 kk = E/(2*pi);              %    Case 1       or   Case 2
 k_round = round(kk);   %k_round     = 1               = 0
 TA0 = f1 + k_round*(2*pi);
 
%...Store initial orbital elements (from above) in the vector coe0:
coe0 = [h0 e0 RA0 i0 w0 TA0];

%...Use ODE45 to integrate the Gauss variational equations (Equations
% 12.84) from t0 to tf:
t0 = 0;
tf = 100*days;
nout = 5000; %Number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout);
options = odeset(...
'reltol', 1.e-10, ...
'abstol', 1.e-10, ...
'initialstep', T0/1000);
y0 = coe0';
[t,y] = ode45(@rates, tspan, y0, options);
%...Assign the time histories mnemonic variable names:
h  = y(:,1);
e  = y(:,2);
RA = y(:,3);
i  = y(:,4);
w  = y(:,5);
TA = y(:,6); %True anomaly (theta_dot)

% calculate the rate of change of following variables from the output of ode45 variables (pg.507)
%%% Semi-major axis 
a = (h.^2 ./ mu) .* (1./(1 - e.^2)); 
%%% distance at perigee
rp = a.*(1-e) - RE;

coe = [h e RA i w TA];
for k = 1:length(y)
    [ro(k,:), ~] = sv_from_coe(coe(k,:),mu);%Curtis function to convert from r&v to elements
end
%Position vector components:
x = ro(:,1);
y = ro(:,2);
z = ro(:,3);

figure(1); % Creating Figure
% Plot grid details
hold on
title('Two-Body Trajectory', 'Interpreter', 'Latex')
xlabel('x (km)', 'Interpreter', 'Latex')
ylabel('y (km)', 'Interpreter', 'Latex')
zlabel('z (km)', 'Interpreter', 'Latex')
axis equal
grid minor

% Plotting Earth
opts_units.Units = 'km';
planet3D('earth',opts_units);
light('Position',[10,10,10]);
grid on;
view(-240,30)

% Plotting Trajectory
plot3(x,y,z, 'r','LineWidth',1)
hold off



%...Plot the time histories of the osculatinig elements:
figure(2)
subplot(2,1,1)
plot(t/days,(RA)/deg,'b-','LineWidth',1)
title('Right Ascension (degrees)')
xlabel('days')
grid on
grid minor
axis tight
subplot(2,1,2)
plot(t/days,(w)/deg,'b-','LineWidth',1)
title('Argument of Perigee (degrees)')
xlabel('days')
grid on
grid minor
axis tight
figure(3)
subplot(3,1,1)
plot(t/days,h,'b-','LineWidth',1)
title('Angular Momentum (km^2/s)')
xlabel('days')
grid on
grid minor
axis tight
subplot(3,1,2)
plot(t/days,e,'b-','LineWidth',1)
title('Eccentricity')
xlabel('days')
grid on
grid minor
axis tight
subplot(3,1,3)
plot(t/days,(i)/deg,'b-','LineWidth',1)
title('Inclination (degrees)')
xlabel('days')
grid on
grid minor
axis tight
figure(4)
subplot(3,1,1)
plot(t/days,a,'b-','LineWidth',1)
title('Semi-major axis vs time')
xlabel('days')
ylabel('a (km)')
grid on
grid minor
axis tight
subplot(3,1,2)
plot(t/days,rp,'b-','LineWidth',1)
title('(Distance at perigee - Earth Radius) vs time')
xlabel('days')
ylabel('km')
grid on
grid minor
axis tight
subplot(3,1,3)
plot(t/days,TA/deg/360,'b-','LineWidth',1)
title('True anomaly vs time')
xlabel('days')
ylabel('Degrees')
grid on
grid minor
axis tight
%...Subfunction:
% 
function dfdt = rates(t,f)
% 
%
% This function calculates the time rates of the orbital elements
% from Gauss’s variational equations (Equations 10.84).
% –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
mu = 398600;
radius_earth = 6370;
J2 = 0.00108;

    h     = f(1); % (km^2/s)
    e     = f(2);
    omega = f(3); % [rad]
    i     = f(4); % [rad]
    w     = f(5); % [rad]
    theta = f(6); % [rad]
    u = w + theta; %(rad)
    
    r = h^2/mu/(1 + e*cos(theta)); %The position magnitude

% pg.510 perturbation components in noninertial rsw frame
pr = -3/2*(J2*mu*radius_earth^2)/r^4 * (1 - 3*(sin(i))^2 * (sin(u))^2);
ps = -3/2*(J2*mu*radius_earth^2)/r^4 * ((sin(i))^2 * sin(2*u));
pw = -3/2*(J2*mu*radius_earth^2)/r^4 * ((sin(2*i)) * sin(u));

% Gauss planetary equations pg.506
h_dot = r*ps;
%%%
e_dot = (h/mu) * sin(theta)*pr + 1/(mu*h) *...
    (((h)^2 + mu*r)*cos(theta) + mu*e*r)*ps;
%%%
theta_dot = h/(r^2) + 1/(e * h) *((h^2 / mu * cos(theta) * pr) - ((r + h^2/mu)*sin(theta)*ps));
%%%
omega_dot = r/(h*sin(i)) * sin(u)*pw;
%%%
i_dot = r/h * cos(u)*pw;
%%%
w_dot = -1/(e*h)*((h^2 / mu) *cos(theta)*pr - (r + h^2/mu)*sin(theta)*ps) -...
    ((r*sin(u))/(h*tan(i)))*pw;
    
    dfdt = [h_dot;e_dot;omega_dot;i_dot;w_dot;theta_dot]; 
    %       1-h   2-e   3-omega  4-i    5-w   6-theta
end %rates
% 

% credit: Tamas Kis (2022). 3D Earth and Celestial Bodies (planet3D) 
% https://www.mathworks.com/matlabcentral/fileexchange/86483-3d-earth-and-celestial-bodies-planet3d
% https://github.com/tamaskis/planet3D-MATLAB/releases/tag/v5.2.0
function planet_surface = planet3D(planet,opts)
    
    % ----------------------------
    % Conversion factors and data.
    % ----------------------------
    
    % conversion factors
    factors = {'AU'   1/149597870000;
               'ft'   100/30.48;
               'km'   0.001;
               'm'    1;
               'mi'   100/160934.4;
               'nmi'  1/1852};
    
            % planet/body          radius,       flattening,    obliquity,
            %                      R [m]         f [-]          obl [deg]
    data = {'Sun'                  696000e3      0.000009       0;
            'Moon'                 1738.0e3      0.0012         6.68;
            'Mercury'              2439.0e3      0.0000         0.0;
            'Venus'                6052.0e3      0.000          177.3;
            'Earth'                6370e3        0.0033528131   23.45;
            'Earth Coastlines'     6378.1363e3   0.0033528131   23.45;
            'Earth Cloudy'         6378.1363e3   0.0033528131   23.45;
            'Earth Night'          6378.1363e3   0.0033528131   23.45;
            'Earth Night Cloudy'   6378.1363e3   0.0033528131   23.45;
            'Mars'                 3397.2e3      0.00647630     25.19;
            'Jupiter'              71492.0e3     0.0648744      3.12;
            'Saturn'               60268.0e3     0.0979624      26.73;
            'Uranus'               25559.0e3     0.0229273      97.86;
            'Neptune'              24764.0e3     0.0171         29.56;
            'Pluto'                1151.0e3      0.0            118.0};
    
    % ------------------------------------
    % Sets (or defaults) plotting options.
    % ------------------------------------
    
    % defaults "planet" to 'Earth Cloudy' if not input
    if (nargin == 0) || isempty(planet)
        planet = 'Earth Cloudy';
    end
    
    % sets position of planet's geometric center (defaults to origin)
    if (nargin < 2) || ~isfield(opts,'Position')
        position = [0;0;0];
    else
        position = opts.Position;
    end
    
    % sets rotation angle (defaults to 0)
    if (nargin < 2) || ~isfield(opts,'RotAngle')
        theta = 0;
    else
        theta = opts.RotAngle;
    end
    
    % sets conversion factor (defaults to 1, assuming units of m)
    if (nargin < 2) || ~isfield(opts,'Units')
        units = 'm';
    else
        units = opts.Units;
    end
    
    % sets reference plane (defaults to equatorial plane)
    if (nargin < 2) || ~isfield(opts,'RefPlane')
        reference_plane = 'equatorial';
    else
        reference_plane = opts.RefPlane;
    end
    
    % sets transparency (defaults to 1 so celestial body is solid)
    if (nargin < 2) || ~isfield(opts,'FaceAlpha')
        FaceAlpha = 1;
    else
        FaceAlpha = opts.FaceAlpha;
    end
    
    % determines obliquity
    if strcmpi(reference_plane,'ecliptic')
        obl = data{strcmpi(data(:,1),planet),4};
    else
        obl = 0;
    end
    
    % sets clipping (defaults to 'off')
    if (nargin < 2) || ~isfield(opts,'Clipping')
        Clipping = 'off';
    else
        Clipping = opts.Clipping;
    end
    
    % sets line color (defaults to default MATLAB color)
    if (nargin < 2) || ~isfield(opts,'Color')
        Color = [0,0.4470,0.7410];
    else
        Color = opts.Color;
    end
    
    % sets line style (defaults to solid line)
    if (nargin < 2) || ~isfield(opts,'LineStyle')
        LineStyle = '-';
    else
        LineStyle = opts.LineStyle;
    end
    
    % sets line width (defaults to 0.5)
    if (nargin < 2) || ~isfield(opts,'LineWidth')
        LineWidth = 0.5;
    else
        LineWidth = opts.LineWidth;
    end
    
    % -------------------------------
    % Geometry of the celestial body.
    % -------------------------------
    
    % determines mean equatorial radius and flattening
    R = data{strcmpi(data(:,1),planet),2};
    f = data{strcmpi(data(:,1),planet),3};
    
    % conversion factor to use
    conversion_factor = factors{strcmpi(factors(:,1),units),2};
    
    % determines semi-major and semi-minor axes of body
    a = conversion_factor*R;
    b = a*(1-f);
    
    % coordinates of ellipsoid (uses 400 panels)
    [x,y,z] = ellipsoid(position(1),position(2),position(3),a,a,b,400);
    
    % ------------------------------------------------------------
    % Defining surfaces/coordinates needed to draw celestial body.
    % ------------------------------------------------------------
    
    % not drawing Earth coastlines 
    if ~strcmpi(planet,'Earth Coastlines')
        
        % loads image data
        if strcmpi(planet,'Earth Cloudy')
            cdata = imread('earth.png')+imread('clouds.png');
        elseif strcmpi(planet,'Earth Night Cloudy')
            cdata = imread('earthnight.png')+0.1*imread('clouds.png');
        else
            cdata = imread(strcat('',lower(planet),'.png'));
        end
        
        % draws planet
        planet_surface = surface(x,y,z,'FaceColor','texture',...
            'EdgeColor','none','CData',flipud(cdata),'DiffuseStrength',...
            1,'SpecularStrength',0,'FaceAlpha',FaceAlpha);
        
    end
    
    % drawing Earth coastlines
    if strcmpi(planet,'Earth Coastlines')
        
        % white surface (lines will be plotted on top of this)
        planet_surface = surface(x,y,z,'FaceColor','w','EdgeColor',...
            'none','DiffuseStrength',1,'SpecularStrength',0,...
            'FaceAlpha',FaceAlpha);
        
        % loads coastline data
        coastlines_data = struct2cell(load('coastlines_data'));
        coastlines_data = [coastlines_data{:}];
        
        % extracts ECEF coordinates of coastlines
        x_coast = coastlines_data.X;
        y_coast = coastlines_data.Y;
        z_coast = coastlines_data.Z;
        
    end
    
    % -------------------
    % Performs rotations.
    % -------------------
    
    % transformation matrix for rotation
    R3 = [ cosd(theta)   sind(theta)   0;
          -sind(theta)   cosd(theta)   0;
           0             0             1];
    
    % transformation matrix for tilt
    R1 = [1   0            0;
          0   cosd(obl)   -sind(obl);
          0   sind(obl)    cosd(obl)];
    
    % axes for rotations (must be row vectors)
    alpha1 = [1,0,0];
    alpha2 = (R1*[0;0;1])';
    
    % tilts celestial body if referenced to ecliptic plane
    rotate(planet_surface,alpha1,obl);
    
    % rotates celestial body about its 3rd axis
    rotate(planet_surface,alpha2,theta);
    
    % rotates coordinates of coastlines
    if strcmpi(planet,'Earth Coastlines')
        new_coordinates = R3*R1*[x_coast';y_coast';z_coast'];
        x_coast = new_coordinates(1,:);
        y_coast = new_coordinates(2,:);
        z_coast = new_coordinates(3,:);
    end
    
    % --------------------------------------------------------------
    % Drawing additional lines (i.e. coastlines or rings of Saturn).
    % --------------------------------------------------------------
    
    % draws coastlines
    if strcmpi(planet,'Earth Coastlines')
        hold on;
        plot3(x_coast,y_coast,z_coast,'LineWidth',LineWidth,'LineStyle',...
            LineStyle,'Color',Color);
        hold off;
    end
    
    % draws rings of Saturn
    if strcmpi(planet,'Saturn')
        
        % reads in image
        cdata_rings = imread('saturnrings.png');
        
        % determines number of different colors in ring (if you look at the
        % image, the way it is formatted just looks like horizontal bands
        % of colors)
        n = size(cdata_rings,2);
        
        % preallocates array to store colors
        colors = zeros(n,3);
        
        % extracts rgb values from image data
        for i = 1:n
            colors(i,:) = cdata_rings(1,i,:);
        end
        
        % shrinks the data set of colors down to only 200 (this is for
        % speed - we plot the bands of Saturns rings as individual lines
        % and don't want to plot thousands of lines) - this shrinking
        % process is condensed from the reduced_data_points function (see 
        % https://www.mathworks.com/matlabcentral/fileexchange/86218-reduce
        % -number-of-data-points-reduce_data_points)
        n_new = 200;
        colors = colors(1:round(n/n_new):n,:);
        
        % scales colors to between 0 and 1 (currently between 0 and 255)
        colors = colors/255;
        
        % plots the rings
        theta = 0:0.001:2*pi;
        hold on;
        for i = 1:n_new
            
            % this comes from the fact that Saturns rings extend from 7000
            % to 80000 km (7,000,000 to 8,000,000 km) from the surface of
            % the planet s(https://en.wikipedia.org/wiki/Rings_of_Saturn)
            r = conversion_factor*(R+7000000+((80000000-7000000)/n_new)*i);
            
            % x, y, and z coordinates of Saturns rings in equatorial plane
            x_ring = position(1)+r*cos(theta);
            y_ring = position(2)+r*sin(theta);
            z_ring = position(3)*ones(size(theta));
            
            % rotates rings to equatorial plane (uses same rotation matrix
            % as tilting the planet earlier in code)
            new_coordinates = R1*[x_ring;y_ring;z_ring];
            x_ring = new_coordinates(1,:);
            y_ring = new_coordinates(2,:);
            z_ring = new_coordinates(3,:);
            
            % plots the jth ring
            plot3(x_ring,y_ring,z_ring,'Color',colors(i,:));
            
        end
        hold off;
        
    end
    
    % ----------------------
    % Basic plot formatting.
    % ----------------------
    
    % set axis clipping
    ax = gca;
    ax.Clipping = Clipping;
    
    % equal data unit lengths along each axis
    axis equal;
    
    % 3D view
    view(3);
    
end
