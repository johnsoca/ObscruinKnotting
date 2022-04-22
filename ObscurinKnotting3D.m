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
% X - x & y position of the tip of the domain
% E - x and y coordinates of the end of the domain, also the start of
%            the linker- domain is 30 Angstroms in length
% alpha - randomly generated according to a bi-modal normal distribution
%           with mu=93.3 and sigma=13.2, and mu=58.4 and sigma=9.44 (degrees) as the
%           angle between the linker to the next domain
% N - number of nodes which include a domain and a linker


%--MONTE CARLO--
sims=1;
%---------------

% Angle distribution values
mu_d2l = deg2rad(126); %mu of node domain to linker
s_d2l = deg2rad(18.5); %sigma of node domain to linker
mu_l2d = deg2rad([93.3,58.4]); %mu of linker to domain bimodal distribution
s_l2d = deg2rad([13.2,9.44]); %sigma of linker to domain bimodal distribution

% Initialize other variables
l_l = 13;
l_d = 33;
N = 18; %number of IG domain-linker pairs, i.e. one node
thresDist = l_d/2; %threshold distance for crossings, the radius of the Ig Domain if assumed to be a spherical shape

R=l_d/2; % set radius of circle to search for clusters
pt_num4=4;
pt_num5=5;
m=1;

tic
%for plotting the range of x, y, and z values
x_range=zeros(1,sims);
y_range=zeros(1,sims);
z_range=zeros(1,sims);
crossings=zeros(1,sims);


%-- MONTE CARLO--
for sim = 1:sims
    % Initialize variables
    X = zeros(3, N+1);
    E = zeros(3,N);
    
    count=0;


    % Initialize the first position and angle of the first domain of node 1
    X(:,1) = [rand(); rand(); rand()];
    theta = plus_minus * rand()*pi; % random theta wrt the x axis
    gamma = plus_minus * rand()*pi; % random gamma wrt to the y-axis
    
    E(:,1)=X(:,1) + l_d * [cos(theta); sin(theta)*cos(gamma); sin(theta)*sin(gamma)];
    A=E(:,1)-X(:,1);
norm(A) %should be l_d
    
    alpha=acos(A(1)/norm(A)); % angle of vector A wrt the x-axis for the rotation.
    eP=[cos(alpha) sin(alpha) 0; -sin(alpha) cos(alpha) 0; 0 0 1]*E(:,1);
    
    phi = plus_minus * normrnd(mu_d2l,s_d2l); %angle wrt the rotated x-axis
    phi = scaling(phi);
    tau=plus_minus*rand()*pi; % random angle wrt the y-axis
    
    xP=eP+l_l*[cos(phi); sin(phi)*cos(tau); sin(phi)*sin(tau)];
    
    X(:,2)=[cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1]*xP;
    
    B=X(:,2)-E(:,1);
    norm(B) % should be l_l
    phi
    acos(dot(A,B)/(norm(A)*norm(B)))
    
    xt(1)=x_e(1)+l_l*cos(phi(1))*sin(tau);
    yt(1)=y_e(1)+l_l*sin(phi(1))*sin(tau);
    zt(1)=z_e(1)+l_l*cos(tau);

    %###############CHECK#####################
    if abs((sqrt((xt(1)-x_e(1))^2+(yt(1)-y_e(1))^2+(zt(1)-z_e(1))^2))-l_l)>0.01 %defined tolerance of "equal"
        disp("first linker is not within tolerance of the length of the linker doman (l_d)");
    end
    %#########################################  

    %update the beginning of the next domain such that it's the same as the
    %termination of the linker (xt,yt,zt)
    for i=2:N
        x(i) = xt(i-1);
        y(i) = yt(i-1);
        z(i) = zt(i-1);
        bimodal = round(rand()+1); % 50% probability of choosing bimodal distribution index 1 as 2
        theta(i) = plus_minus*normrnd(mu_l2d(bimodal),s_l2d(bimodal)); 
        tau=plus_minus*rand()*pi;
        x_e(i)=x(i)+l_d*cos(theta(i))*sin(tau);
        y_e(i)=y(i)+l_d*sin(theta(i))*sin(tau);
        z_e(i)=z(i)+l_d*cos(tau);

    %###############CHECK#####################
    if abs((sqrt((x_e(i)-x(i))^2+(y_e(i)-y(i))^2+(z_e(i)-z(i))^2))-l_d)>0.01 %defined tolerance of "equal"
        disp("second node is not within tolerance of the length of the Ig doman (l_d)",i);
    end
    %#########################################    
        phi(i) = plus_minus * normrnd(mu_d2l,s_d2l);
        tau=plus_minus*rand()*pi; % random each time called
        xt(i)=x_e(i)+l_l*cos(phi(i))*sin(tau);
        yt(i)=y_e(i)+l_l*sin(phi(i))*sin(tau);
        zt(i)=z_e(i)+l_l*cos(tau);
    %###############CHECK#####################
    if abs((sqrt((xt(i)-x_e(i))^2+(yt(i)-y_e(i))^2+(zt(i)-z_e(i))^2))-l_l)>0.01 %defined tolerance of "equal"
        disp("second linker is not within tolerance of the length of the linker doman (l_d)",i);
    end
    %######################################### 

    end
    
    % for purposes of figuring out the range, create an additional x, y, z
    % coordinate
    x(N+1)=xt(N);
    y(N+1)=yt(N);
    z(N+1)=zt(N);
    
    x_range(sim)=max(x)-min(x);
    y_range(sim)=max(y)-min(y);
    z_range(sim)=max(z)-min(z);

    % COMMENT OUT FOR MONTE CARLO SIMS
%     % Plot the result of 18 linkers 
%     for i=1:N    
%         plot3([x(i),x_e(i)],[y(i),y_e(i)],[z(i),z_e(i)],'k');
%         hold on
%         plot3([x_e(i),xt(i)],[y_e(i),yt(i)],[z_e(i),zt(i)],'r+:');
%     end

    % Add code to count the number of cross overs in 3D
    % 
    count=0;
    % I'm going to iterate through the linkers such that linker 1 is
    % compared to linkers 3-18 for determining cross overs (there's no way
    % linkers next to each other with indices can cross over). To help with
    % indexing in the nested loop, I'm going to create a new variable to
    % put the poitns of the linkers end to end.
    L=zeros(3,N*2);
    ind=1;
    for i=1:N
        L(:,ind)=[x(i); y(i); z(i)];
        L(:,ind+1)=[x_e(i); y_e(i); z_e(i)];
        ind=ind+2;
    end
    
    for i=1:(N*2)-3
        for j=(i+2):(N*2)-1
            [dist vShort]=DistBetween2Segment(L(:,i), L(:,i+1),L(:,j),L(:,j+1));
            %cross = crossovers3D(L(1,i), L(2,i), L(3,i), L(1,i+1), L(2,i+1), L(3,i+1), L(1,j), L(2,j), L(3,j), L(1,j+1),L(2,j+1),L(3,j+1));
            if dist <= thresDist % if linkers and domains are within a specified threshold distance, it counts as a crossing
                count=count+1;
            end
        end
    end
    crossings(sim)=count;
    
end
toc
figure()
hist(x_range)
figure()
hist(y_range)
figure()
hist(z_range)
figure()
hist(crossings)

% Save data into txt file
M=zeros(sims,4);
M(:,1)=x_range;
M(:,2)=y_range;
M(:,3)=z_range;
M(:,4)=crossings;
csvwrite('MillionSim3D.txt',M); %NOTE! This will overwrite exisiting file

    
    
    
    
    
    
    
    