figure()
hold on
for its=1:N
    plot3([x(1,its) x_e(1,its)],[y(1,its) y_e(1,its)],[z(1,its) z_e(1,its)],'b-', 'LineWidth',3);
    plot3([x_e(1,its) xt(1,its)],[y_e(1,its) yt(1,its)],[z_e(1,its) zt(1,its)],'g-','LineWidth',3);
end
view(211,19)
xlabel('x-axis');
ylabel('y-axis');
zlabel('z-axis');

% figure()
% hold on
% for its=1:N
%     plot3([x(1,its) x_e(1,its)],[y(1,its) y_e(1,its)],[z(1,its) z_e(1,its)],'b-');
%     plot3([x_e(1,its) xt(1,its)],[y_e(1,its) yt(1,its)],[z_e(1,its) zt(1,its)],'g-');
% end
% % Add code to count the number of cross overs in 3D
%     % 
%     count=0;
%     % I'm going to iterate through the linkers such that linker 1 is
%     % compared to linkers 3-18 for determining cross overs (there's no way
%     % linkers next to each other with indices can cross over). To help with
%     % indexing in the nested loop, I'm going to create a new variable to
%     % put the poitns of the linkers end to end.
%     L=zeros(3,N*2);
%     ind=1;
%     for i=1:N
%         L(:,ind)=[x(i); y(i); z(i)];
%         L(:,ind+1)=[x_e(i); y_e(i); z_e(i)];
%         ind=ind+2;
%     end
%     
%     for i=1:(N*2)-3
%         for j=(i+2):(N*2)-1
%             [dist vShort]=DistBetween2Segment(L(:,i), L(:,i+1),L(:,j),L(:,j+1));
%             %cross = crossovers3D(L(1,i), L(2,i), L(3,i), L(1,i+1), L(2,i+1), L(3,i+1), L(1,j), L(2,j), L(3,j), L(1,j+1),L(2,j+1),L(3,j+1));
%             if dist <= thresDist % if linkers and domains are within a specified threshold distance, it counts as a crossing
%                 count=count+1;
%                 plot3([L(1,i) L(1,i+1)],[L(2,i) L(2,i+1)],[L(3,i) L(3,i+1)],'ro');
%                 plot3([L(1,j) L(1,j+1)],[L(2,j) L(2,j+1)],[L(3,j) L(3,j+1)],'ko');
%                 
%             end
%         end
%     end
%     crossings(sim)=count;