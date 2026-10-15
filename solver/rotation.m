function [d]=rotation(t,theta,d_hat)
d=t*dot(t,d_hat)+cos(theta)*(d_hat-t*dot(t,d_hat))+sin(theta)*cross(t,d_hat);
end