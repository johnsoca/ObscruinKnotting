function [x,y]=multipleEqnSolver(n, l, theta, x1, x2, y1, y2, z1, z2)
% the purpose of this function is to solve 2 equations and 2 unknowns from
% previously identified vector magnitude and vector dot product. See notes
% for details about how the following variables were declared and used to
% solve.

A = (n*l*cos(theta) -(z2-z1)*(z3-z2)+y2*(y2-y1))/(x2-x1)+x2;
B=-(y2-y1)/(x2-x1);
C=l^2-(x3-x2)^2;
D= -(C-A^2+2*A*x2-x2^2-y2^2);
E=2*A*B-2*x2*B-2*y2;
F=B^2-1;

%check if the determinant is imaginary
if E^2<4*F*D
    disp('imaginary determinant');
    return
end

y_pos = (-E+sqrt(E^2-4*F*D))/(2*F);
y_neg = (-E-sqrt(E^2-4*F*D))/(2*F);

if rand()>0.5
    y=y_pos;
else
    y=y_neg;
end

x=A+B*y;