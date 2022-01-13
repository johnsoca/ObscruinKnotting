function [count]=crossovers3D(ax, ay, az, bx, by, bz, cx, cy, cz, dx, dy, dz)
% Callie J Miller
% 01/13/22
% The purpose of this function is to determine if two line segments in 3D
% intersect. Based on the formula referenced in https://math.stackexchange.com/questions/3114665/how-to-find-if-two-line-segments-intersect-in-3d
% Using the equation of a line between points A and B (tA+(1-t)B) where
% 0 <=t<=1, and the equation of a line between points C and D
% (sC+(1-s)D) for 0<=s<=1, in 3D these are 3 equations and 2 unknowns.
% Solve the first two equations (linear algebra with an inverse
% matrix), then check the solutions for t and s to see if the third
% equation is true. If the third equation is true, then the line
% segments intersect.
%
% Inputs: 4 points in 3D that represent 2 line segments from (ax,ay,az) to
% (bz, by, bz) and the line segment from (cx, cy, cz) to (dx, dy, dz).
% Output: count = 0 if the 2 line segments do not intersect, or count = 1
% if the 2 line segments do intersect

count = 0;

A = [ax-bx, dx-cx; ay-by, dy-cy];
y = [dx-bx; dy-by];
x = inv(A)*y;

check_eqn_left=x(1)*(az-bz)+x(2)*(dz-cz);
check_eqn_right=dz-bz;
        
if check_eqn_left==check_eqn_right
    count=count+1;
end