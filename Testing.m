figure() %comment out if running monte carlo sims
hold on
for i=1:18
    plot([X(1,i),E(1,i)],[X(2,i),E(2,i)],'b');
    plot([E(1,i),X(1,i+1)],[E(2,i),X(2,i+1)],'r');
end

R=33/2;
pt_num=5;

% COUNT the number of clusters, i.e. tangles
    D=zeros(1,N*2);
    P=[1:1:N*2];
    numCluster(sim)=0;
    C={'k','g','m','b','r','c',[0.5 0.6 0.7],[0.8 0.2 0.6]};
    for i=1:N*2
        if P(i) ~=0
            D=sqrt((L(1,P(i))-L(1,:)).^2+(L(2,P(i))-L(2,:)).^2+(L(3,P(i))-L(3,:)).^2);
            repeat=logical(P~=0);
            threshold=logical(D<=R);
            cluster=logical(repeat+threshold==2);
            clusterSize=nnz(cluster);
            if clusterSize >= pt_num % there's the potential for tangling because there are greater than pt_num points in the area
                %Now, determine if there are a minimum of 2 crossings for the identified points within the circle    
                ind=find(cluster); % indices of points within the cluster
                numCluster(sim)=numCluster(sim)+1;
                for j=1:length(ind)
                    plot(L(1,P(ind(j))),L(2,P(ind(j))),'color', C{numCluster},'marker','*');
                end
                P(ind)=0;
            end
        end
    end