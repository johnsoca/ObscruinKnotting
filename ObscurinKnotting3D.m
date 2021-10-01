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
% theta - angle between the linker and the following domain: randomly
%       generated bimodal distribution
% x_e, y_e, z_e - x, y and z coordinates of the end of the domain node, also the start of
%            the linker- domain is between 32-34 Angstroms in length
% xt, yt, zt - x, y and z coordinates of the end of the linker- linker is 12.8 angstroms
%          in length
% alpha - randomly generated according to a bi-modal normal distribution
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
    % **** could define this in terms of theta(1), but for now leaving alone
    % because it doesn't work in my head ************************
    alpha = plus_minus * rand()*pi; % for 3D space, need to randomly declare alpha as angle between the first domain and x-axis
    beta = plus_minus * rand()*pi; % for 3D space, need to randomly declare beta as angle between the first domain and y-axis
    gamma = plus_minus * rand()*pi; % for 3D space, need to randomly declare gamma as angle between the first domain and z-axis
    
    x_e(1)=x(1)+l_d*cos(alpha);
    y_e(1)=y(1)+l_d*cos(beta);
    z_e(1)=z(1)+l_d*cos(gamma);
    phi(1) = plus_minus * normrnd(mu_d2l,s_d2l);
     
    tau(1)=plus_minus*rand()*pi; % random each time called
    zt(1)=l_d*cos(tau(1))+z_e(1);
    [xt(1),yt(1)]=multipleEqnSolver(l_d,l_l,phi(1), x(1), x_e(1), y(1), y_e(1), z(1), z_e(1),zt(1));
    
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
    tau(2)=plus_minus *rand()*pi;
    z_e(2)=l_d*cos(tau(2))+z(2);
    [x_e(2),y_e(2)]=multipleEqnSolver(l_l,l_d,theta(2),x_e(1),x(2), y_e(1), y(2), z_e(1),z(2),z_e(2)); %NOTE switch l_l and l_d because of the vector in question the magnitude changes
    
    phi(2) = plus_minus * normrnd(mu_d2l,s_d2l);
    tau(3)=plus_minus*rand()*pi; % random each time called
    zt(2)=l_d*cos(tau(3))+z_e(1);
    [xt(2),yt(2)]=multipleEqnSolver(l_d,l_l,phi(2), x(2), x_e(2), y(2), y_e(2), z(2), z_e(2),zt(2));

% Troubleshooting 
for i=1:N    
    plot3([x(i),x_e(i)],[y(i),y_e(i)],[z(i),z_e(i)],'k');
    hold on
    plot3([x_e(i),xt(i)],[y_e(i),yt(i)],[z_e(i),zt(i)],'r+:');
end
    
    
    
    
    
    
    
    
    
    