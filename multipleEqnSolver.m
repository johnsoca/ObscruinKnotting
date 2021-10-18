function [x,y]=multipleEqnSolver(n, l, theta, x1, x2, y1, y2, z1, z2, z3)
% the purpose of this function is to solve 2 equations and 2 unknowns from
% previously identified vector magnitude and vector dot product. See notes
% for details about how the following variables were declared and used to
% solve.

A = n*l*cos(theta) -(z3-z2)*(z1-z2)+y2*(y1-y2)+x2*(x1-x2);
B=x1-x2;
C=y1-y2;
D=l^2-(z3-z2)^2;
E=(A/B)-x2;
F=-C/B;
G=-2*y2;
H=E^2-D;
J=F^2+1;
K=G+2*E*F;

%check if the determinant is imaginary
if K^2<4*H*J
    disp('imaginary determinant');
    return
end

y_pos = (-K+sqrt(K^2-4*J*H))/(2*J);
y_neg = (-K-sqrt(K^2-4*J*H))/(2*J);

if rand()>0.5
    y=y_pos;
else
    y=y_neg;
end

x=(A-C*y)/B;