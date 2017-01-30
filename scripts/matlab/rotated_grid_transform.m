function [grid_out] = rotated_grid_transform(grid_in, option, SP_coor)

lon = grid_in(:,1);
lat = grid_in(:,2);

lon = (lon*pi)/180; % Convert degrees to radians
lat = (lat*pi)/180;

SP_lon = SP_coor(1);
SP_lat = SP_coor(2);

%theta = 90+SP_lat; % Rotation around y-axis
theta = SP_lat-90; % Rotation around y-axis
phi = SP_lon; % Rotation around z-axis

phi = (phi*pi)/180; % Convert degrees to radians
theta = (theta*pi)/180;

x = cos(lon).*cos(lat); % Convert from spherical to cartesian coordinates
y = sin(lon).*cos(lat);
z = sin(lat);

if option == 1 % Regular -> Rotated

    x_new = cos(theta).*cos(phi).*x + cos(theta).*sin(phi).*y + sin(theta).*z;
    y_new = -sin(phi).*x + cos(phi).*y;
    z_new = -sin(theta).*cos(phi).*x - sin(theta).*sin(phi).*y + cos(theta).*z;

elseif option == 2 % Rotated -> Regular

    phi = -phi;
    theta = -theta;

    x_new = cos(theta).*cos(phi).*x + sin(phi).*y + sin(theta).*cos(phi).*z;
    y_new = -cos(theta).*sin(phi).*x + cos(phi).*y - sin(theta).*sin(phi).*z;
    z_new = -sin(theta).*x + cos(theta).*z;

end

lon_new = atan2(y_new,x_new); % Convert cartesian back to spherical coordinates
lat_new = asin(z_new);

lon_new = (lon_new*180)/pi; % Convert radians back to degrees
lat_new = (lat_new*180)/pi;

grid_out = [lon_new lat_new];