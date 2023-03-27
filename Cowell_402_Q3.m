
clear all; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HW2 Question (2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%...Conversion factors:
hours = 3600; %Hours to seconds
days = 24*hours; %Days to seconds
deg = pi/180; %Degrees to radians

% Constants
mu = 398600; % mu in km^3/s^2
radius_earth = 6370; % in km
period_earth = 86400; %(seconds) for 1 day

%Given Initial Conditions
a = 26600; % semi major axis in km
e_magnitude = 0.74; % eccentricity
i = 1.10654; % radians
%i_degrees = ; % inclination
omega_degrees = 90; %Longitude of ascending node Ω (omega)
w_degrees = 5; %Argument of periapse (w)
M_degrees = 10; % Mean Anomaly


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Convrersion from degrees to radians
%i = deg2rad( i_degrees );
omega = deg2rad( omega_degrees);
w = deg2rad(w_degrees);
M = deg2rad(M_degrees);


I_unit = [1;0;0]; %I unit vector
J_unit = [0;1;0]; %J unit vector
K_unit = [0;0;1]; %K unit vector

h_magnitude = sqrt(mu*a*(1 - e_magnitude^2)); %angular momentum magnitude

E = M;
g = 1;
itr = 0;
while abs(g) > 1e-15
     
g = E-e_magnitude*sin(E) - M; % rad/s Also this is the modified Kepler's equation
dgdE = 1-e_magnitude*cos(E);
E_new = E - g/dgdE;
% Update
E = E_new;
itr = itr + 1;

end


 f1 = (2*atan(sqrt((1 + e_magnitude) / (1 - e_magnitude)) * tan(E/2))); % True Anomaly
    
 kk = E/(2*pi);              %Case 2 = 0.9860   Case 3 = 0.4753
 k_round = round(kk);        %       = 1               = 0
 f = f1 + k_round*(2*pi);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding r and v magnitudes
r0_magnitude = (h_magnitude^2 / mu) / (1 + e_magnitude*cos(f)); 

v0_radial_magnitude = (mu*e_magnitude*sin(f)) / (h_magnitude);
v0_tangential_magnitude = h_magnitude/r0_magnitude;

v0_magnitude = sqrt( (v0_radial_magnitude)^2 + (v0_tangential_magnitude)^2 );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finding r and v vectors
theta = w + f;

r0_vector = ...
   r0_magnitude*(...
      (cos(theta)*cos(omega) - cos(i)*sin(omega)*sin(theta))*I_unit +...
      (cos(theta)*sin(omega) + cos(i)*cos(omega)*sin(theta))*J_unit +...
      (sin(i)*sin(theta))*K_unit);

v0_vector = ...
    - mu/h_magnitude*(cos(omega)*(sin(theta) + e_magnitude*sin(w)) +...
                     sin(omega)*(cos(theta) + e_magnitude*cos(w))*cos(i)...
                     )*I_unit ...
    - mu/h_magnitude*(sin(omega)*(sin(theta) + e_magnitude*sin(w)) -...
                     cos(omega)*(cos(theta) + e_magnitude*cos(w))*cos(i)...
                     )*J_unit ...
    + mu/h_magnitude*(sin(i)*(cos(theta) + e_magnitude*cos(w)))*K_unit;



e_vector = ((v0_magnitude)^2 / mu - 1/r0_magnitude)*r0_vector - (1/mu)*(dot(r0_vector,v0_vector))*v0_vector;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% True anomaly (f)
if dot( r0_vector , v0_vector) < 0
    f_correct1 = 2*pi - acos(dot(r0_vector , e_vector) / (r0_magnitude * e_magnitude) );
else
    f_correct1 = acos(dot(r0_vector , e_vector) / (r0_magnitude * e_magnitude) );
end

%f_degree = 180/pi * f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given Initial Conditions
Y0 = [r0_vector;v0_vector];     % combine r0 and v0 into a single column vector
% or Y0 = [-8903.833; 1208.356; 213.066; -0.971; -6.065; -1.069]; % [x; y; z; vx; vy; vz] [km, km/s]

r0_magnitude = norm(r0_vector); %initial position magnitude
v0_magnitude = norm(v0_vector); %initial velocity magnitude

t0 = 0; % initial time

period = 2*pi*sqrt(a^3 / mu); % Orbit period of sattelite in seconds
tf = 100*period_earth; %required 100 days 
nout = 5000; %Number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout);
options = odeset(...
'reltol', 1.e-10, ...
'abstol', 1.e-10, ...
'initialstep', period);
%options = odeset('RelTol', 1e-10, 'abstol', 1e-10); % Setting tolerance
% Numerical Integration using ode45
[t, Y] = ode45(@ODE2BP, tspan, Y0, options,mu); %(name of function, time span, initial conditions, options)

%Pulling Position and Velocity Data from Output
x = Y(:, 1); % [km]
y = Y(:, 2); % [km]
z = Y(:, 3); % [km]
vx = Y(:, 4); % [km/s]
vy = Y(:, 5); % [km/s]
vz = Y(:, 6); % [km/s]  
   
% J2 perturbation acceleration vector at t = 0
 J2 = 0.00108;
 p_vector = (3/2)*((J2*mu*(radius_earth)^2)/(r0_magnitude)^4 *...
    ((r0_vector(1)/r0_magnitude)*(5*(r0_vector(3)/r0_magnitude)^2 - 1)*I_unit +...
    ((r0_vector(2)/r0_magnitude)*(5*(r0_vector(3)/r0_magnitude)^2 - 1)*J_unit +...
    ((r0_vector(3)/r0_magnitude)*(5*(r0_vector(3)/r0_magnitude)^2 - 3)*K_unit)))); % J2 perturbation acceleration vector



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
plot3(x, y, z, 'color', "#D95319"','LineWidth',1)
hold off



for j= 1:length(Y)
    
r_current = [x(j);y(j);z(j)];    %position
v_current = [vx(j);vy(j);vz(j)]; %velocity

r_current_magnitude = norm(r_current); %initial position magnitude
v_current_magnitude = norm(v_current); %initial velocity magnitude

a_current = 1/(2/r_current_magnitude - (v_current_magnitude)^2 / mu); %semi major axis

e_current_vector = ((v_current_magnitude)^2 / mu - 1/r_current_magnitude)*r_current - (1/mu)*(dot(r_current,v_current))*v_current;
e_current_magnitude = norm(e_current_vector);

h_current_vector = cross(r_current,v_current); %angular momentum
h_current_magnitude = norm(h_current_vector);

I_unit = [1;0;0]; %I unit vector
J_unit = [0;1;0]; %J unit vector
K_unit = [0;0;1]; %K unit vector

node_current_vector = cross(K_unit , h_current_vector/h_current_magnitude); % node vector
node_current_magnitude = norm(node_current_vector);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Inclination
i_current = acos(dot(h_current_vector/h_current_magnitude , K_unit)); %inclination
i_current_degree = 180/pi * i_current;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Longitude of ascending node Ω (omega)
if dot( node_current_vector , J_unit) < 0 
    omega_current = 2*pi - acos((dot(node_current_vector , I_unit))/node_current_magnitude); 
else 
    omega_current = acos((dot(node_current_vector , I_unit))/node_current_magnitude); 
end 
omega_current_degree = 180/pi * omega_current;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Argument of periapse (w)
if dot( e_current_vector , K_unit) < 0
    w_current = 2*pi - acos(dot(node_current_vector , e_current_vector) / (node_current_magnitude * e_current_magnitude) );
else
    w_current = acos(dot(node_current_vector , e_current_vector) / (node_current_magnitude * e_current_magnitude) );
end
w_current_degree = 180/pi * w_current;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% True anomaly (f)
if dot( r_current , v_current) < 0
    f_current = 2*pi - acos(dot(r_current , e_current_vector) / (r_current_magnitude * e_current_magnitude) );
else
    f_current = acos(dot(r_current , e_current_vector) / (r_current_magnitude * e_current_magnitude) );
end
f_current_degree = 180/pi * f_current;

elements(j,:) = [omega_current_degree,w_current_degree,a_current,e_current_magnitude,i_current_degree];


end
%%%%%% output values of final positioin and velocity vectors
fprintf('position vector r: [');
fprintf('%.4f ',r_current);
fprintf('] km\n');

fprintf('velocity vector r: [');
fprintf('%.4f ',v_current);
fprintf('] km/s\n');

fprintf('acceleration a\t = %.4f km/s^2\n',a_current);
fprintf('Eccentricity e \t = %.4f \n',e_current_magnitude);
fprintf('Inclinaion i \t = %.4f deg\n',i_current_degree);
fprintf('Right Ascention \x03A9\t = %.4f deg\n',omega_current_degree);
fprintf('Argument of perigee w\t = %.4f deg\n',w_current_degree);

hold on
%...Plot the time histories of the osculatinig elements:
figure(2)
subplot(2,1,1)
plot(t/days,elements(:,1),'b-','LineWidth',1)
title('Right Ascension (Ω) vs time')
xlabel('days')
ylabel('Ω (degrees)')
grid on
grid minor
axis tight

subplot(2,1,2)
plot(t/days,elements(:,2),'b-','LineWidth',1)
title('Argument of Perigee (w) vs time')
xlabel('days')
ylabel('w (degrees)')
grid on
grid minor
axis tight

figure(3)
subplot(3,1,1)
plot(t/days,elements(:,3),'b-','LineWidth',1)
title('Semi-major axis vs time')
xlabel('days')
ylabel('a (km)')
grid on
grid minor
axis tight

subplot(3,1,2)
plot(t/days,elements(:,4),'b-','LineWidth',1)
title('Eccentricity vs time')
xlabel('days')
ylabel('e')
grid on
grid minor
axis tight

subplot(3,1,3)
plot(t/days,elements(:,5),'b-','LineWidth',1)
title('Inclination vs time')
ylabel(' i (degrees)')
xlabel('days')
grid on
grid minor
axis tight
hold off
%%

% User-Defined ODE Function, ODE 2BodyProblem
function dYdt = ODE2BP(t, Y,mu) %   t and Y are output variables. (t) is numerical integration time step and (Y) is the state vector at each time step.
   radius_earth = 6370;
   J2 = 0.00108;
   
    x = Y(1); % [km]
    y = Y(2); % [km]
    z = Y(3); % [km]
    r = norm(Y(1:3)); %r_magnitude
    
    p_vector = [(-3/2*J2*(mu/r^2)*(radius_earth/r)^2)*((1 - 5*(z/r)^2)*x/r),...
                (-3/2*J2*(mu/r^2)*(radius_earth/r)^2)*((1 - 5*(z/r)^2)*y/r),...
                (-3/2*J2*(mu/r^2)*(radius_earth/r)^2)*((3 - 5*(z/r)^2)*z/r)]; % J2 perturbation acceleration vector
    
    vx = Y(4); % [km/s]
    vy = Y(5); % [km/s]
    vz = Y(6); % [km/s]
    ax = - mu * x / (r^3) + p_vector(1); % [km/s^2]
    ay = - mu * y / (r^3) + p_vector(2); % [km/s^2]
    az = - mu * z / (r^3) + p_vector(3); % [km/s^2]    
    
    dYdt = [vx;vy;vz;ax;ay;az]; % Y' or derivative of Y

end

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



