% The purpose of this code is to test how to plot saved data for specific
% instances from ObscurinKnotting.m
% Callie J Miller
% 07/19/21

% d = size(dphi);
% 
% for i=1: 1%d(1)
%     figure();
%     hold on
%     for n=1:2:d(2)-2
%         plot([dX(i,n),dX(i,n+1)],[dY(i,n),dY(i,n+1)],'k');
%         plot([dX(i,n+1),dX(i,n+2)],[dY(i,n+1),dY(i,n+2)],'r+:');
%     end
% end

figure()
hold on
% test counting cross overs
    count=0;
    for s = 1:37-3
        x1 = dX(1,s);
        y1 = dY(1,s);
        x2 = dX(1,s+1);
        y2 = dY(1,s+1);
        plot([x1,x2],[y1,y2],'b');
        for n = s+2: 36
            x3 = dX(1,n);
            y3 = dY(1,n);
            x4 = dX(1,n+1);
            y4 = dY(1,n+1);
            plot([x3,x4],[y3,y4],':');
            c = crossovers(x1,y1,x2,y2,x3,y3,x4,y4);
            count = count + c;
        end
    end

% % Test if the distance between (x,y) and (x_e,y_e) is 30 A and distance
% % between (x_e,y_e) and (x,y) of next domain is 18 Angstroms
% s=1;
% for i=1:17
%     dist30 = sqrt((dX(1,s)-dX(1,s+1))^2+((dY(1,s)-dY(1,s+1))^2))
%     dist18 = sqrt((dX(1,s+1)-dX(1,s+2))^2+(dY(1,s+1)-dY(1,s+2))^2)
%     s=s+2;
% end