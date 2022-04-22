%Testing.m
A0=[0; 0; 1];
% theta=pi/6;
% gamma=pi/3;
% Ae=A0+30*[cos(theta); sin(theta)*cos(gamma); sin(theta)*sin(gamma)];
% A=Ae-A0;
A=A0;
norm(A);

alpha=acos(A(1)/norm(A));
beta=acos(A(2)/norm(A));
gamma=acos(A(3)/norm(A));
R=[cos(alpha)*cos(beta) cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma)   cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma);
   sin(alpha)*cos(beta) sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma)   sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma);
   -sin(beta)           cos(beta)*sin(gamma)                                    cos(beta)*cos(gamma)];
AP=R*A;

phi=pi/8;
tau=pi/5;

BeP=AP+13*[cos(phi); sin(phi)*cos(tau); sin(phi)*sin(tau)];

Be=inv(R)*BeP;

norm(Be-A)

phi_check=acos(dot(A,Be-A)/(norm(A)*norm(Be-A)))
