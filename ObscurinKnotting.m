% ObscurinKnotting.m
% Callie J Miller
% 06/18/21
% Includes the function(s): plus_minus.m

% Variables:
% phi- angle between the head and the tail of a single node
% x - x position of the tip of the head
% y - y position of the tip of the head
% theta_h - angle of the head
% x_e, y_e - x and y coordinates of the end of the head, also the start of
%            the tail linker- head is 30 Angstroms in length
% xt, yt - x and y coordinates of the end of the tail- tail is 18 angstroms
%          in length
% alpha - randomly generated to be the angle between 75 and 180 between the
%         tail of one linker and the head of the next linker
% N - number of nodes which include a head and a tail

sims=1;
x_range = zeros(1,sims);
y_range = zeros(1,sims);

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

    % Initialize the first position and angle of the first head of node 1
    x(1) = rand();
    y(1) = rand();
    theta_h(1) = plus_minus * rand()*pi; 

    figure()
    hold on

    for i = 1: N
        x_e(i) = x(i) + 30*cos(theta_h(i));
        y_e(i) = y(i) + 30*sin(theta_h(i));

        phi(i) = plus_minus * rand() * acos(8/18); 
        xt(i) = x_e(i) + 18*cos(theta_h(i)+phi(i));
        yt(i) = y_e(i) + 18*sin(theta_h(i)+phi(i));

        plot([x(i),x_e(i)],[y(i),y_e(i)],'k');
        plot([x_e(i),xt(i)],[y_e(i),yt(i)],'r+:');
        if i ~=N
            x(i+1) = xt(i);
            y(i+1) = yt(i);
            alpha = deg2rad(plus_minus * (105*rand()+75)); % randomly generate an angle (in radians) between 75 and 180 degrees for the angle between the tail of a node and the head of the next node
            beta = acos((x(i+1)-x_e(i))/18); % angle of the tail linker wrt the x axis
            theta_h(i+1) = pi-(beta+alpha);
        end
    end
    x_range(sim) = max(x)-min(x);
    y_range(sim) = max(y)-min(y);
end
