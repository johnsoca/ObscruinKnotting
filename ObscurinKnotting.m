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

sims=1000000;
x_range = zeros(1,sims);
y_range = zeros(1,sims);
crossings = zeros(1,sims);

% Angle distribution values
mu_d2l = deg2rad(126); %mu of domain to linker
s_d2l = deg2rad(18.5); %sigma of domain to linker
mu_l2d = deg2rad([93.3,58.4]); %mu of linker to domain bimodal distribution
s_l2d = deg2rad([13.2,9.44]); %sigma of linker to domain bimodal distribution

R=33/2; % set radius of circle to search for clusters as 33 Angstroms.
pt_num4=4; % number of points within a circle of radius R to qualify as a cluster
pt_num5=5;
tic
m=1;
for sim = 1:sims
    % Initialize variables
    N = 18;

    X = zeros(3,N+1);
    E = zeros(3,N);
    
    count=0;
    thresDist=33/2; % threshold distance for crossing
    
    % Initialize the first position and angle of the first domain of node 1
    X(1:2,1) = [rand(); rand()];
    alpha = plus_minus * rand()*pi; % initial orienation of head domain wrt x axis doens't matter 

    E(1:2,1)=X(1:2,1)+30*[cos(alpha); sin(alpha)]; %coordinates of the end of the head node
    
% %----------------------------------------------------
%     figure() %comment out if running monte carlo sims
%     plot([X(1,1),E(1,1)],[X(2,1),E(2,1)],'b');
%     hold on
% %----------------------------------------------------

    % Rotate the coordinate system with respect to alpha using
    % transformation/rotation matrix
    eP=[cos(alpha), sin(alpha); -sin(alpha), cos(alpha)]*E(1:2,1);
        
    phi = plus_minus * normrnd(mu_d2l,s_d2l);
    phi = scaling(phi);
    tP=eP+13*[cos(phi); sin(phi)];
        
    X(1:2,2)=[cos(alpha), -sin(alpha); sin(alpha), cos(alpha)]*tP; %coordinates of the end of the tail node, same as start of next head node
% %----------------------------------------------------
% % CHECK: magnitude of X to E is 30, magnitude of E to X(i+1) is 13, angle 
% % between EX and EX(i+1)is phi
% norm(E(:,1)-X(:,1)) %should be close to 30
% norm(X(:,2)-E(:,1)) % should be close to 13
% A=X(:,1)-E(:,1);
% B=X(:,2)-E(:,1);
% pi-acos(dot(A,B)/(norm(A)*norm(B))) %should be close to pi-phi
% phi
% %----------------------------------------------------

% %----------------------------------------------------
%     plot([E(1,1),X(1,2)],[E(2,1),X(2,2)],'r');
% %----------------------------------------------------

    for i=2:N
        % Rotate the coordinate system with respect to the prior segment using
        % transformation/rotation matrix
        theta=atan((X(2,i)-E(2,i-1))/(X(1,i)-E(1,i-1)));
        xP = [cos(theta), sin(theta); -sin(theta), cos(theta)]*X(1:2,i);
        bimodal = round(rand()+1); % 50% probability of choosing bimodal distribution index 1 as 2
        alpha=scaling(plus_minus*normrnd(mu_l2d(bimodal),s_l2d(bimodal)));
        eP=xP+30*[cos(alpha); sin(alpha)]; %coordinates of the end of the head node
        E(1:2,i)=[cos(theta), -sin(theta); sin(theta), cos(theta)]*eP; % coordinates of the point between the domain and linker
        
        % Rotate the coordinate system with respect to the prior segment
        theta=atan((E(2,i)-X(2,i))/(E(1,i)-X(1,i)));
        eP= [cos(theta), sin(theta); -sin(theta), cos(theta)]*E(1:2,i);
        phi = scaling(plus_minus*normrnd(mu_d2l,s_d2l));
        tP=eP+13*[cos(phi); sin(phi)];
        X(1:2,i+1)=[cos(theta), -sin(theta); sin(theta), cos(theta)]*tP; %coordinates of the end of the tail node, same as start of next head node
% %----------------------------------------------------
%     plot([X(1,i),E(1,i)],[X(2,i),E(2,i)],'b');
%     plot([E(1,i),X(1,i+1)],[E(2,i),X(2,i+1)],'r');
% %----------------------------------------------------
% %----------------------------------------------------
% % CHECK: magnitude of X to E is 30, magnitude of E to X(i+1) is 13, angle 
% % between EX and EX(i+1)is phi, angle between XE(i-1) and XE is alpha
% A=X(:,i)-E(:,i);
% B=X(:,i+1)-E(:,i);
% norm(A) %should be close to 30
% norm(B) % should be close to 13
% phi_test = acos(dot(A,B)/(norm(A)*norm(B))) %should be close to phi
% phi
% 
% C=E(:,i-1)-X(:,i);
% D=E(:,i)-X(:,i);
% alpha_test = acos(dot(C,D)/(norm(C)*norm(D))) %should be close to alpha
% alpha
% %----------------------------------------------------
        
    end
    x_range(sim) = max(X(1,:))-min(X(1,:));
    y_range(sim) = max(X(2,:))-min(X(2,:));
    

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
    num5Cluster(sim)=0;
%     C={'k','g','m','b','r','c',[0.5 0.6 0.7],[0.8 0.2 0.6],'k','g','m','b','r','c',[0.5 0.6 0.7], [0.8 0.2 0.6]};
    for i=1:N*2
        if P(i) ~=0
            D=sqrt((L(1,P(i))-L(1,:)).^2+(L(2,P(i))-L(2,:)).^2+(L(3,P(i))-L(3,:)).^2);
            repeat=logical(P~=0);
            threshold=logical(D<=R);
            cluster=logical(repeat+threshold==2);
            clusterSize=nnz(cluster);
            if clusterSize >= pt_num4
                num4Cluster(sim)=num4Cluster(sim)+1;
                ind=find(cluster); % indices of points within the cluster
%                 for j=1:length(ind)
%                     plot(L(1,P(ind(j))),L(2,P(ind(j))),'color',C{numCluster},'marker','*');
%                 end
                P(ind)=0;
            end
            if clusterSize >= pt_num5
                num5Cluster(sim)=num5Cluster(sim)+1;
                ind=find(cluster); % indices of points within the cluster
%                 for j=1:length(ind)
%                     plot(L(1,P(ind(j))),L(2,P(ind(j))),'color',C{numCluster},'marker','*');
%                 end
                P(ind)=0;
            end
        end
    end
end

%------------------------------------
% Plot the histograms of the "spread" in the nodes in the x and y
% directions to determine the likelihood that an 18-domain/linker node is
% stretched or crumpled on itself
%------------------------------------
figure()
hist(x_range);
figure()
hist(y_range);
figure()
hist(crossings);
figure()
hist(num4Cluster);
figure()
hist(num5Cluster);
% % % % % figure()
% % % % % hist(crossings/N);
toc


% Save data into txt file
M=zeros(sims,6);
M(:,1)=x_range;
M(:,2)=y_range;
M(:,3)=crossings;
M(:,4)=crossings/N;
M(:,5)=num4Cluster;
M(:,6)=num5Cluster;
csvwrite('MillionSim2D.txt',M); %NOTE! This will overwrite exisiting file

