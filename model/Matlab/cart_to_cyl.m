function [u_r,u_theta] = cart_to_cyl(u,v,theta)
% this function transforms a cartesian vector u,v
% to cylindrical coordinates u_r, u_theta, with angle theta (rad)

u_r = u .* cos(theta) + v .* sin(theta);
u_theta = -u .* sin(theta) + v .* cos(theta);

end