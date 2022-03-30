function [theta]=scaling(theta)
% The purpose of this function is to rescale theta values if chosen to be
% greater than pi or -pi, so the rest of the transformations for
% ObscurinKnotting.m code makes sense
%
% Callie J Miller 3/29/22

if theta > pi
    theta = -(2*pi-theta);
elseif theta < -pi
    theta = 2*pi-abs(theta);
else
    theta=theta;
end