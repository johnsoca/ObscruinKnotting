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
pt_num4=4;
pt_num5=5;
m=1;
e_x=[1; 0; 0]; % The unit vector for the x-axis

tic
%for plotting the range of x, y, and z values
x_range=zeros(1,sims);
y_range=zeros(1,sims);
z_range=zeros(1,sims);
crossings=zeros(1,sims);

%-- MONTE CARLO--
for sim = 1:sims
    % Initialize variables
    X = zeros(3,N+1);
    E = zeros(3,N);
    count=0;

    % Initialize the first position and angle of the first domain of node 1
    X(:,1)=[rand(); rand(); rand()];
        alpha = plus_minus * rand()*pi; % initial orientaiton random
        tau = plus_minus * rand()*pi; % initial twist random
    E(:,1)=X(:,1)+l_d*[cos(alpha); cos(tau)*sin(alpha); sin(tau)*sin(alpha)];
    
% %     %---------------------------------------------------
% %     figure() % comment out if running Monte Carlo sims
% %     plot3([X(1,1),E(1,1)],[X(2,1), E(2,1)],[X(3,1),E(3,1)],'b');
% %     hold on
% %     %---------------------------------------------------
    
    % Rotate the coordinate system to make the x-axis align with the vector
    % created by X(:,1) and E(:,1) pointing from X to E
    A=E(:,1)-X(:,1);
        R=rotateAlign(A/norm(A),e_x);
        eP=R*E(:,1);
        phi = plus_minus * normrnd(mu_d2l, s_d2l);
        phi = scaling(phi);
        tau = plus_minus * pi* rand();
        tP = eP+l_l*[cos(phi); cos(tau)*sin(phi); sin(tau)*sin(phi)];
    X(:,2)=inv(R)*tP;
    
% %     %--------------------------------------------------
% %     plot3([E(1,1),X(1,2)],[E(2,1),X(2,2)],[E(3,1),X(3,2)],'r');
% %     %--------------------------------------------------
    
    for i=2:N
        % Rotate the coordinate system wrt the prior segment
        A=X(:,i)-E(:,i-1);
            R=rotateAlign(A/norm(A),e_x);
            xP=R*X(:,i);
            bimodal = round(rand()+1); % 50% probability of choosing bimodal distribution index 1 as 2
            alpha=scaling(plus_minus*normrnd(mu_l2d(bimodal),s_l2d(bimodal)));
            tau = plus_minus * pi* rand();
            eP = xP+l_d*[cos(alpha); cos(tau)*sin(alpha); sin(tau)*sin(alpha)];
        E(:,i)=inv(R)*eP;
        
        A=E(:,i)-X(:,i);
            R=rotateAlign(A/norm(A),e_x);
            eP=R*E(:,i);
            phi = plus_minus * normrnd(mu_d2l, s_d2l);
            phi = scaling(phi);
            tau = plus_minus * pi* rand();
            tP = eP+l_l*[cos(phi); cos(tau)*sin(phi); sin(tau)*sin(phi)];
        X(:,i+1)=inv(R)*tP;
        
% %         %--------------------------------------------------
% %         plot3([X(1,i),E(1,i)],[X(2,i),E(2,i)],[X(3,i),E(3,i)],'b');
% %         plot3([E(1,i),X(1,i+1)],[E(2,i),X(2,i+1)],[E(3,i),X(3,i+1)],'r');
% %         %--------------------------------------------------


    end
% figure() %comment out if running monte carlo sims
% for i=1:N
%     plot3([X(1,i),E(1,i)],[X(2,i),E(2,i)],[X(3,i),E(3,i)],'k'); %plot one simulation for the final figure where the IG domain is green
%     hold on
%     plot3([E(1,i),X(1,i+1)],[E(2,i),X(2,i+1)],[E(3,i),X(3,i+1)],'color',[0,1,1]);% linker is cyan
% end    
    
    x_range(sim)=max(X(1,:))-min(X(1,:));
    y_range(sim)=max(X(2,:))-min(X(2,:));
    z_range(sim)=max(X(3,:))-min(X(3,:));

    % -------------------
    % Code below is used to test if a "cross-over" has occurred,
    % and count the instances of cross-overs
    %--------------------
    L=zeros(3,N*2);
    ind=1;
    for i=1:N
        L(:,ind)=X(:,i);
        L(:,ind+1)=E(:,i);
        ind=ind+2;
    end
    
    % COUNT the number of crossings
    for i=1:(N*2)-3
        for j=(i+2):(N*2)-1
            [dist vShort]=DistBetween2Segment(L(:,i), L(:,i+1),L(:,j),L(:,j+1));
            if dist <= thresDist % if linkers and domains are within a specified distance, it "crosses"
                count=count+1;
            end
        end
    end
    crossings(sim)=count;
    
    % COUNT the number of clusters, i.e. tangles
    D=zeros(1,N*2);
    P=[1:1:N*2];
    num4Cluster(sim)=0;
    
%     C={'k','g','m','b','r','c',[0.5 0.6 0.7],[0.8 0.2 0.6],'k','g','m','b','r','c',[0.5 0.6 0.7], [0.8 0.2 0.6]};
    for i=1:N*2
        if P(i) ~=0
            D=sqrt((L(1,P(i))-L(1,:)).^2+(L(2,P(i))-L(2,:)).^2+(L(3,P(i))-L(3,:)).^2);
            repeat=logical(P~=0);
            threshold=logical(D<=thresDist);
            cluster=logical(repeat+threshold==2);
            clusterSize=nnz(cluster);
            if clusterSize >= pt_num4
                num4Cluster(sim)=num4Cluster(sim)+1;
                ind=find(cluster); % indices of points within the cluster
%                 for j=1:length(ind)
%                     plot3(L(1,P(ind(j))),L(2,P(ind(j))),L(3,P(ind(j))),'color',C{num4Cluster},'marker','*');
%                 end
                P(ind)=0;
            end

        end
    end
    D=zeros(1,N*2);
    P=[1:1:N*2];
    num5Cluster(sim)=0;
% figure() %comment out if running monte carlo sims
% for i=1:N
%     plot3([X(1,i),E(1,i)],[X(2,i),E(2,i)],[X(3,i),E(3,i)],[60,180,75]); %plot one simulation for the final figure where the IG domain is green
%     hold on
%     plot3([E(1,i),X(1,i+1)],[E(2,i),X(2,i+1)],[E(3,i),X(3,i+1)],[70,240,240]);% linker is cyan
% end
    for i=1:N*2
        if P(i) ~=0
            D=sqrt((L(1,P(i))-L(1,:)).^2+(L(2,P(i))-L(2,:)).^2+(L(3,P(i))-L(3,:)).^2);
            repeat=logical(P~=0);
            threshold=logical(D<=thresDist);
            cluster=logical(repeat+threshold==2);
            clusterSize=nnz(cluster);
            if clusterSize >= pt_num5
                num5Cluster(sim)=num5Cluster(sim)+1;
                ind=find(cluster); % indices of points within the cluster
%                 for j=1:length(ind)
%                     plot3(L(1,P(ind(j))),L(2,P(ind(j))),L(3,P(ind(j))),'color',C{num5Cluster},'marker','*');
%                 end
                P(ind)=0;
            end
        end
    end
    
    
end
toc
% figure()
% hist(x_range)
% figure()
% hist(y_range)
% figure()
% hist(z_range)
% figure()
% hist(crossings)
% figure()
% hist(num4Cluster)
% figure()
% hist(num5Cluster)
% 
% % Save data into txt file
% M=zeros(sims,6);
% M(:,1)=x_range;
% M(:,2)=y_range;
% M(:,3)=z_range;
% M(:,4)=crossings;
% M(:,5)=num4Cluster;
% M(:,6)=num5Cluster;
% csvwrite('MillionSim3D.txt',M); %NOTE! This will overwrite exisiting file

    
    
    
    
    
    
    
    