function [dP]=dist3D_Segment_to_Segment(Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Dx,Dy,Dz)
% MATLAB code modified from C++ Code from Dan Sunday
% Copyright 2001, 2012, 2021 Dan Sunday
% This code may be freely used and modified for any purpose providing that
% this copyright notice is included with it. There is no warranty for this
% code, and the author of it cannot be held liable for any real or imagined
% damage from its use. Users of this code must verify correctness for their
% application.

% For my purposes, I have the coordinates of line segments separated into
% x, y, and z components, so I'll define the appropriate vectors here.
u=[Bx-Ax; By-Ay; Bz-Az];
v=[Dx-Cx; Dy-Cy; Dz-Cz];
w=[Cx-Ax; Cy-Ay; Cz-Az];

a=dot(u,u); %this should always be >=0
b=dot(u,v);
c=dot(v,v); %always >=0
d=dot(u,w);
e=dot(v,w);
D=a*c-b*b; %always >=0
sD=D;
tD=D;

sNum=0.00000001;

if D<sNum
    sN=0;
    sD=1;
    tN=e;
    tD=c;
else
    sN=(b*e-c*d);
    tN=(a*e-b*d);
    if sN <0 %end case
        sN=0;
        tN=e;
        tD=c;
    elseif sN>sD
        sN=sD;
        tN=e+b;
        tD=c;
    end
end

if tN <0 %end case
    tN=0;
    if -d<0
        sN=0;
    elseif (-d>a)
        sN=sD;
    else
        sN=-d;
        sD=a;
    end
elseif tN>tD %end case
    tN=tD;
    if (-d+b)<0
        sN=0;
    elseif ((-d+b)>a)
        sN=sD;
    else
        sN=(-d+b);
        sD=a;
    end
end

if abs(sN)<sNum
    sc=0;
else
    sc=sN/sD;
end

if abs(tN)<sNum
    tc=0;
else
    tc=tN/tD;
end

% Difference of the 2 closest points
dP=w+(sc*u)-(tc*v);

dP=norm(dP);