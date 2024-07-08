function [] = plot_singleRod(center,angle,l,h,c_rods)

% Parameters for the rod-shaped bacterium
length = l; % l: Length of the bacterium
width = h;   % h: Height of the bacterium
x_center = center(1); % X-coordinate of the center of the bacterium
y_center = center(2); % Y-coordinate of the center of the bacterium
%angle: Rotation angle in degrees

% Calculate half-dimensions
half_length = length / 2;
half_width = width / 2;
% Discretization to draw the circles
th_disc=(linspace(-pi/2,pi/2,6))';
th_disc2=(linspace(3*pi/2,pi/2,6))';

corners = [half_length+half_width*cos(th_disc), half_width*sin(th_disc);...
           -half_length+half_width*cos(th_disc2),half_width*sin(-th_disc2)];

% Rotation matrix
theta = angle; 
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];

% Apply rotation and translate to the center position
rotated_corners = (R * corners')';
rotated_corners(:,1) = rotated_corners(:,1) + x_center;
rotated_corners(:,2) = rotated_corners(:,2) + y_center;

% Plot the rod-shaped bacterium
fill(rotated_corners(:,1), rotated_corners(:,2), c_rods,'LineStyle','none'); % Draw the filled rotated rectangle


end