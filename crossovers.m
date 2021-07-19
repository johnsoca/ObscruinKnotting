function [count]=crossovers(x1,y1,x2,y2,x3,y3,x4,y4)

fx3y3 = (x3-x1)*(y2-y1)-(y3-y1)*(x2-x1);
fx4y4 = (x4-x1)*(y2-y1)-(y4-y1)*(x2-x1);
gx1y1 = (x1-x3)*(y4-y3)-(y1-y3)*(x4-x3);
gx2y2 = (x2-x3)*(y4-y3)-(y2-y3)*(x4-x3);

if (fx3y3*fx4y4) < 0 && (gx1y1*gx2y2) < 0
    count = 1;
else
    count=0;
end
