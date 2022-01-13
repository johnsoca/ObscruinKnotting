% CheckCrossing3D.m
% Callie J Miller
% 01/13/22
%
% The purpose of this code is to check that the code to count cross overs
% in 3D is correct. Generate data where cross overs or intersections must
% occur, then count how many intersections exist.
N=18;
count=0;

L=rand(3,N*2+1);

for i=1:(N*2)-3
    for j=(i+2):(N*2)
        cross = crossovers3D(L(1,i), L(2,i), L(3,i), L(1,i+1), L(2,i+1), L(3,i+1), L(1,j), L(2,j), L(3,j), L(1,j+1),L(2,j+1),L(3,j+1));
        count=count+cross;
    end
end

disp(count);
figure()
for i=1:2:N*2    
    plot3([L(1,i),L(1,i+1)],[L(2,i),L(2,i+1)],[L(3,i),L(3,i+1)],'k');
    hold on
    plot3([L(1,i+1),L(1,i+2)],[L(2,i+1),L(2,i+2)],[L(3,i+1),L(3,i+2)],'r+:');
end