% ObscurinKnotting.m
% Callie J Miller
% 06/18/21
% Includes the function(s): plus_minus.m

% Variables:
% phi- angle between the domain and the linker of a single node: randomly generated according to normal
%       distribution mu=126, sigma=18.5 degrees
% X - x & y position of the tip of the domain
% E - x and y coordinates of the end of the domain, also the start of
%            the linker- domain is 30 Angstroms in length
% alpha - randomly generated according to a bi-modal normal distribution
%           with mu=93.3 and sigma=13.2, and mu=58.4 and sigma=9.44 (degrees) as the
%           angle between the linker to the next domain
% N - number of nodes which include a domain and a linker

sims=1;
x_range = zeros(1,sims);
y_range = zeros(1,sims);
crossings = zeros(1,sims);

% Angle distribution values
mu_d2l = deg2rad(126); %mu of domain to linker
s_d2l = deg2rad(18.5); %sigma of domain to linker
mu_l2d = deg2rad([93.3,58.4]); %mu of linker to domain bimodal distribution
s_l2d = deg2rad([13.2,9.44]); %sigma of linker to domain bimodal distribution

tic
m=1;
for sim = 1:sims
    % Initialize variables
    N = 18;

    X = zeros(2,N+1);
    E = zeros(2,N);
    
    % Initialize the first position and angle of the first domain of node 1
    X(:,1) = [rand(); rand()];
    alpha = plus_minus * rand()*pi; % initial orienation of head domain wrt x axis doens't matter 

    E(:,1)=X(:,1)+30*[cos(alpha); sin(alpha)]; %coordinates of the end of the head node
    
%----------------------------------------------------
    figure() %comment out if running monte carlo sims
    plot([X(1,1),E(1,1)],[X(2,1),E(2,1)],'b');
    hold on
%----------------------------------------------------

    % Rotate the coordinate system with respect to alpha using
    % transformation/rotation matrix
    eP=[cos(alpha), sin(alpha); -sin(alpha), cos(alpha)]*E(:,1);
        
    phi = plus_minus*normrnd(mu_d2l,s_d2l);
    tP=eP+13*[cos(phi); sin(phi)];
        
    X(:,2)=[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*tP; %coordinates of the end of the tail node, same as start of next head node
%----------------------------------------------------
% CHECK: magnitude of X to E is 30, magnitude of E to X(i+1) is 13, angle 
% between EX and EX(i+1)is phi
norm(E(:,1)-X(:,1)) %should be close to 30
norm(X(:,2)-E(:,1)) % should be close to 13
A=X(:,1)-E(:,1);
B=X(:,2)-E(:,1);
pi-acos(dot(A,B)/(norm(A)*norm(B))) %should be close to pi-phi
phi
%----------------------------------------------------

%----------------------------------------------------
    plot([E(1,1),X(1,2)],[E(2,1),X(2,2)],'r');
%----------------------------------------------------

    for i=2:N
        bimodal = round(rand()+1); % 50% probability of choosing bimodal distribution index 1 as 2
        alpha=normrnd(mu_l2d(bimodal),s_l2d(bimodal));
        E(:,i)=X(:,i)+30*[cos(alpha); sin(alpha)]; %coordinates of the end of the head node

        % Rotate the coordinate system with respect to alpha using
        % transformation/rotation matrix
        eP=[cos(alpha), sin(alpha); -sin(alpha), cos(alpha)]*E(:,i);
        
        phi = normrnd(mu_d2l,s_d2l);
        tP=eP+13*[cos(phi); sin(phi)];
        
        X(:,i+1)=[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*tP; %coordinates of the end of the tail node, same as start of next head node
%----------------------------------------------------
    plot([X(1,i),E(1,i)],[X(2,i),E(2,i)],'b');
    plot([E(1,i),X(1,i+1)],[E(2,i),X(2,i+1)],'r');
%----------------------------------------------------
%----------------------------------------------------
% CHECK: magnitude of X to E is 30, magnitude of E to X(i+1) is 13, angle 
% between EX and EX(i+1)is phi, angle between XE(i-1) and XE is alpha
A=X(:,i)-E(:,i);
B=X(:,i+1)-E(:,i);
norm(A) %should be close to 30
norm(B) % should be close to 13
pi-acos(dot(A,B)/(norm(A)*norm(B))) %should be close to phi
phi

C=E(:,i-1)-X(:,i);
D=E(:,i)-X(:,i);
norm(C) % should be close to 13
norm(D) % should be close to 30
pi-acos(dot(C,D)/(norm(C)*norm(D))) %should be close to alpha
alpha
%----------------------------------------------------
    end
    x_range(sim) = max(X(1,:))-min(X(1,:));
    y_range(sim) = max(X(2,:))-min(X(2,:));
    

% % % %     % -------------------
% % % %     % Code below is used to test if a "cross-over" has occurred,
% % % %     % and count the instances of cross-overs
% % % %     %--------------------
% % % %     s = 1;
% % % %     for i = 1 : N
% % % %         X(s) = x(i);
% % % %         X(s+1) = x_e(i);
% % % %         Y(s) = y(i);
% % % %         Y(s+1) = y_e(i);
% % % %         s=s+2;
% % % %     end
% % % %     X(37) = xt(18);
% % % %     Y(37) = yt(18);
% % % %     
% % % %     count=0;
% % % %     for s = 1:length(X)-3
% % % %         x1 = X(s);
% % % %         y1 = Y(s);
% % % %         x2 = X(s+1);
% % % %         y2 = Y(s+1);
% % % %         for n = s+2: length(X)-1
% % % %             x3 = X(n);
% % % %             y3 = Y(n);
% % % %             x4 = X(n+1);
% % % %             y4 = Y(n+1);
% % % %             c = crossovers(x1,y1,x2,y2,x3,y3,x4,y4);
% % % %             count = count + c;
% % % %         end
% % % %     end
% % % %     
% % % %     crossings(sim)=count;
% % % % %     if crossings(sim)>100
% % % % %         dX(m,:)=X;
% % % % %         dY(m,:)=Y;
% % % % %         dtheta(m,:)=theta_h;
% % % % %         dphi(m,:)=phi;
% % % % %         CROSSINGCHECK(m) = crossings(sim);
% % % % %         m=m+1;
% % % % %     end
end

%------------------------------------
% Plot the histograms of the "spread" in the nodes in the x and y
% directions to determine the likelihood that an 18-domain/linker node is
% stretched or crumpled on itself
%------------------------------------
% figure()
% hist(x_range);
% figure()
% hist(y_range);
% % % % figure()
% % % % hist(crossings);
% % % % figure()
% % % % hist(crossings/N);
% % % % toc
