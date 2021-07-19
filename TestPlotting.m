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

% test counting cross overs
    count=0;
    for s = 1:37-3
        x1 = dX(s);
        y1 = dY(s);
        x2 = dX(s+1);
        y2 = dY(s+1);
        for n = s+2: 36
            x3 = dX(n);
            y3 = dY(n);
            x4 = dX(n+1);
            y4 = dY(n+1);
            c = crossovers(x1,y1,x2,y2,x3,y3,x4,y4);
            count = count + c;
        end
    end