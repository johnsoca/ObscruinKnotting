% ObscurinKnotting3D.m
% Callie J Miller
% 10/01/21
% The purpose of this code is to simulate the node-linker morphologies
% based on data from YASARA for obscurin.
% Definitions:
%   domain- Ig domain
%   linker- Linker between Ig domains
%
% Includes the function(s): plus_minus.m

% Variables:
% phi- angle between the domain and the following linker: randomly generated according to normal
%       distribution mu=126, sigma=18.5 degrees
% x - x position of the tip of the domain
% y - y position of the tip of the domain node
% z - z position of the tip of the domain node
% x_e, y_e, z_e - x, y and z coordinates of the end of the domain node, also the start of
%            the linker- domain is between 32-34 Angstroms in length
% xt, yt, zt - x, y and z coordinates of the end of the linker- linker is 12.8 angstroms
%          in length
% theta - randomly generated according to a bi-modal normal distribution
%           with mu=93.3 and sigma=13.2, and mu=58.4 and sigma=9.44 (degrees) as the
%           angle between the linker to the next domain
% N - number of nodes which include a domain and a linker
% l_l - length of linker, 12.8 Angstroms
% l_d - length of domain, between 32-34 Angstroms, *** choose to be 33 ***
% tau - random with each call, determines the degree of "twist" as defined
%   from the angle between the vector and the z-axis.

%--MONTE CARLO--
%sims=1000000;
%---------------

% Angle distribution values
mu_d2l = deg2rad(126); %mu of node domain to linker
s_d2l = deg2rad(18.5); %sigma of node domain to linker
mu_l2d = deg2rad([93.3,58.4]); %mu of linker to domain bimodal distribution
s_l2d = deg2rad([13.2,9.44]); %sigma of linker to domain bimodal distribution

% Initialize other variables
l_l = 12.8;
l_d = 33;
N = 2; % 18; %number of IG domain-linker pairs, i.e. one node


%-- MONTE CARLO--
%for sim = 1:sims
    % Initialize variables
    phi = zeros(1,N);
    x = zeros(1,N);
    y = zeros(1,N);
    z = zeros(1,N);
    x_e = zeros(1,N);
    y_e = zeros(1,N);
    z_e = zeros(1,N);
    xt = zeros(1,N);
    yt = zeros(1,N);
    zt = zeros(1,N);
    theta = zeros(1,N);
    phi = zeros(1,N);
    tau=zeros(1,2*N);


    % Initialize the first position and angle of the first domain of node 1
    x(1) = rand();
    y(1) = rand();
    z(1) = rand();
    bimodal = round(rand()+1); % 50% probability of choosing bimodal distribution index 1 as 2
    theta(1) = plus_minus*normrnd(mu_l2d(bimodal),s_l2d(bimodal));
    %use spherical coordinate system where tau is the inclination and theta
    %is the azimuth, and the magnitude is the radius
    tau(1)=plus_minus*rand()*pi; % random each time called
    x_e(1)=x(1)+l_d*cos(theta(1))*sin(tau(1));
    y_e(1)=y(1)+l_d*sin(theta(1))*sin(tau(1));
    z_e(1)=z(1)+l_d*cos(tau(1));
     
%###############CHECK#####################
if abs((sqrt((x_e(1)-x(1))^2+(y_e(1)-y(1))^2+(z_e(1)-z(1))^2))-l_d)>0.01 %defined tolerance of "equal"
    disp("first node is not within tolerance of the length of the Ig doman (l_d)");
end
%#########################################

    phi(1) = plus_minus * normrnd(mu_d2l,s_d2l);
    tau(2)=plus_minus*rand()*pi; % random each time called
    xt(1)=x_e(1)+l_l*cos(phi(1))*sin(tau(2));
    yt(1)=y_e(1)+l_l*sin(phi(1))*sin(tau(2));
    zt(1)=z_e(1)+l_l*cos(tau(2));
    
%###############CHECK#####################
if abs((sqrt((xt(1)-x_e(1))^2+(yt(1)-y_e(1))^2+(zt(1)-z_e(1))^2))-l_l)>0.01 %defined tolerance of "equal"
    disp("first linker is not within tolerance of the length of the linker doman (l_d)");
end
%#########################################  

    %update the beginning of the next domain such that it's the same as the
    %termination of the linker (xt,yt,zt)
% maybe the beginning of the loop?
    i=1; %******
    if i ~=N
        x(i+1) = xt(i);
        y(i+1) = yt(i);
        z(i+1) = zt(i);
        bimodal = round(rand()+1); % 50% probability of choosing bimodal distribution index 1 as 2
        theta(i+1) = plus_minus*normrnd(mu_l2d(bimodal),s_l2d(bimodal)); 
    end
    tau(3)=plus_minus*rand()*pi;
    x_e(2)=x(2)+l_d*cos(theta(2))*sin(tau(3));
    y_e(2)=y(2)+l_d*sin(theta(2))*sin(tau(3));
    z_e(2)=z(2)+l_d*cos(tau(3));
   
%###############CHECK#####################
if abs((sqrt((x_e(2)-x(2))^2+(y_e(2)-y(2))^2+(z_e(2)-z(2))^2))-l_d)>0.01 %defined tolerance of "equal"
    disp("second node is not within tolerance of the length of the Ig doman (l_d)");
end
%#########################################    
    phi(2) = plus_minus * normrnd(mu_d2l,s_d2l);
    tau(4)=plus_minus*rand()*pi; % random each time called
    xt(2)=x_e(2)+l_l*cos(phi(2))*sin(tau(4));
    yt(2)=y_e(2)+l_l*sin(phi(2))*sin(tau(4));
    zt(2)=z_e(2)+l_l*cos(tau(4));
%###############CHECK#####################
if abs((sqrt((xt(2)-x_e(2))^2+(yt(2)-y_e(2))^2+(zt(2)-z_e(2))^2))-l_l)>0.01 %defined tolerance of "equal"
    disp("second linker is not within tolerance of the length of the linker doman (l_d)");
end
%######################################### 
% Troubleshooting 
for i=1:N    
    plot3([x(i),x_e(i)],[y(i),y_e(i)],[z(i),z_e(i)],'k');
    hold on
    plot3([x_e(i),xt(i)],[y_e(i),yt(i)],[z_e(i),zt(i)],'r+:');
end
    
    
    
    
    
    
    
    
    
    