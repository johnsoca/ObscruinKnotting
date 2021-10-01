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
% theta - 
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

%sims=1000000;
%x_range = zeros(1,sims);
%y_range = zeros(1,sims);
%crossings = zeros(1,sims);

% Angle distribution values
mu_d2l = deg2rad(126); %mu of node domain to linker
s_d2l = deg2rad(18.5); %sigma of node domain to linker
mu_l2d = deg2rad([93.3,58.4]); %mu of linker to domain bimodal distribution
s_l2d = deg2rad([13.2,9.44]); %sigma of linker to domain bimodal distribution

% Initialize other variables
l_l = 12.8;
l_d = 33;


%tic
%m=1;
%for sim = 1:sims
    % Initialize variables
    N = 18;
    phi = zeros(1,N);
    x = zeros(1,N);
    y = zeros(1,N);
    z = zeros(1,N);
    theta_h = zeros(1,N);
    x_e = zeros(1,N);
    y_e = zeros(1,N);
    z_e = zeros(1,N);
    xt = zeros(1,N);
    yt = zeros(1,N);
    zt = zeros(1,N);
%     X = zeros(1,37);
%     Y = zeros(1,37);

    % Initialize the first position and angle of the first domain of node 1
    x(1) = rand();
    y(1) = rand();
    theta_h(1) = plus_minus * rand()*pi; % initial orienation of domain doens't matter 

%    %********** COMMENT IF RUNNING MONTE CARLO
%    figure()
%    hold on
%    %*********
    for i = 1: N
        x_e(i) = x(i) + 30*cos(theta_h(i));
        y_e(i) = y(i) + 30*sin(theta_h(i));

        phi(i) = plus_minus * normrnd(mu_d2l,s_d2l) * acos(8/18); 
        xt(i) = x_e(i) + 18*cos(theta_h(i)+phi(i));
        yt(i) = y_e(i) + 18*sin(theta_h(i)+phi(i));
        
%         % ********************
%         plot([x(i),x_e(i)],[y(i),y_e(i)],'k');
%         plot([x_e(i),xt(i)],[y_e(i),yt(i)],'r+:');
%         % ********************
        
        if i ~=N
            x(i+1) = xt(i);
            y(i+1) = yt(i);
            bimodal = round(rand()+1); % 50% probability of choosing bimodal distribution index 1 as 2
            alpha = plus_minus*normrnd(mu_l2d(bimodal),s_l2d(bimodal)); 
            beta = acos((x(i+1)-x_e(i))/18); % angle of the tail linker wrt the x axis
            theta_h(i+1) = pi-(beta+alpha);
        end
    end
    x_range(sim) = max(x)-min(x);
    y_range(sim) = max(y)-min(y);
    

    % -------------------
    % Code below is used to test if a "cross-over" has occurred,
    % and count the instances of cross-overs
    %--------------------
    s = 1;
    for i = 1 : N
        X(s) = x(i);
        X(s+1) = x_e(i);
        Y(s) = y(i);
        Y(s+1) = y_e(i);
        s=s+2;
    end
    X(37) = xt(18);
    Y(37) = yt(18);
    
    count=0;
    for s = 1:length(X)-3
        x1 = X(s);
        y1 = Y(s);
        x2 = X(s+1);
        y2 = Y(s+1);
        for n = s+2: length(X)-1
            x3 = X(n);
            y3 = Y(n);
            x4 = X(n+1);
            y4 = Y(n+1);
            c = crossovers(x1,y1,x2,y2,x3,y3,x4,y4);
            count = count + c;
        end
    end
    
    crossings(sim)=count;
%     if crossings(sim)>100
%         dX(m,:)=X;
%         dY(m,:)=Y;
%         dtheta(m,:)=theta_h;
%         dphi(m,:)=phi;
%         CROSSINGCHECK(m) = crossings(sim);
%         m=m+1;
%     end
end

%------------------------------------
% Plot the histograms of the "spread" in the nodes in the x and y
% directions to determine the likelihood that an 18-domain/linker node is
% stretched or crumpled on itself
%------------------------------------
%figure()
%hist(x_range);
%figure()
%hist(y_range);
figure()
hist(crossings);
figure()
hist(crossings/N);
toc
