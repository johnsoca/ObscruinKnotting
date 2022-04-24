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
