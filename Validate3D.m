A=[0; 2; 0];
e_x=[1; 0; 0];

R=rotateAlign(A/norm(A),e_x);


% alpha=acos(A(1)/norm(A));
% beta=acos(A(2)/norm(A));
% gamma=acos(A(3)/norm(A));
% R=[ cos(beta)*cos(gamma)    sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma)   cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma);
%     sin(gamma)*cos(beta)    sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma)   cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma);
%     -sin(beta)              cos(beta)*sin(alpha)                                    cos(beta)*cos(alpha)];

Ap=R*A;
phi=pi/9;
tau=pi/3;
Bp=Ap+13*[cos(phi); cos(tau)*sin(phi); sin(tau)*sin(phi)];

B=inv(R)*Bp;

phi_test=acos(dot(A,B-A)/(norm(A)*norm(B-A)));

function [R]=rotateAlign(v1,v2)
% Create rotation matrix to rotate a given vector v1 to be aligned with
% another vector v2.
%  https://gist.github.com/kevinmoran/b45980723e53edeb8a5a43c49f134724

ax=cross(v1,v2);
cosA = dot(v1, v2);
k = 1/(1+cosA);

R = [ax(1)^2*k+cosA         ax(2)*ax(1)*k-ax(3)     ax(3)*ax(1)*k+ax(2);
     ax(1)*ax(2)*k+ax(3)    ax(2)^2*k+cosA          ax(3)*ax(2)*k-ax(1);
     ax(1)*ax(3)*k-ax(2)    ax(2)*ax(3)*k + ax(1)   ax(3)^2*k+cosA];
end
