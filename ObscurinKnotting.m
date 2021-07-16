% ObscurinKnotting.m
% Callie J Miller
% 06/18/21
% Includes the function(s): plus_minus.m

% Variables:
% phi- angle between the domain and the linker of a single node: randomly generated according to normal
%       distribution mu=126, sigma=18.5 degrees
% x - x position of the tip of the domain
% y - y position of the tip of the domain
% theta_h - angle of the domain
% x_e, y_e - x and y coordinates of the end of the domain, also the start of
%            the linker- domain is 30 Angstroms in length
% xt, yt - x and y coordinates of the end of the linker- linker is 18 angstroms
%          in length
% alpha - randomly generated according to a bi-modal normal distribution
%           with mu=93.3 and sigma=13.2, and mu=58.4 and sigma=9.44 (degrees) as the
%           angle between the linker to the next domain
% N - number of nodes which include a domain and a linker

sims=1000000;
x_range = zeros(1,sims);
y_range = zeros(1,sims);

% Angle distribution values
mu_d2l = deg2rad(126); %mu of domain to linker
s_d2l = deg2rad(18.5); %sigma of domain to linker
mu_l2d = deg2rad([93.3,58.4]); %mu of linker to domain bimodal distribution
s_l2d = deg2rad([13.2,9.44]); %sigma of linker to domain bimodal distribution

tic
for sim = 1:sims
    % Initialize variables
    N = 18;
    phi = zeros(1,18);
    x = zeros(1,18);
    y = zeros(1,18);
    theta_h = zeros(1,18);
    x_e = zeros(1,18);
    y_e = zeros(1,18);
    xt = zeros(1,18);
    yt = zeros(1,18);

    % Initialize the first position and angle of the first domain of node 1
    x(1) = rand();
    y(1) = rand();
    theta_h(1) = plus_minus * rand()*pi; % initial orienation of domain doens't matter 

%    figure()
%    hold on

    for i = 1: N
        x_e(i) = x(i) + 30*cos(theta_h(i));
        y_e(i) = y(i) + 30*sin(theta_h(i));

        phi(i) = plus_minus * normrnd(mu_d2l,s_d2l) * acos(8/18); 
        xt(i) = x_e(i) + 18*cos(theta_h(i)+phi(i));
        yt(i) = y_e(i) + 18*sin(theta_h(i)+phi(i));

%         plot([x(i),x_e(i)],[y(i),y_e(i)],'k');
%         plot([x_e(i),xt(i)],[y_e(i),yt(i)],'r+:');
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
end
figure()
hist(x_range);
figure()
hist(y_range);
toc
