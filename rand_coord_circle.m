function [x, y]=rand_coord_circle(rc)
%x1,y1,
a=2*pi*rand;% get the angle that is a fraction of 2pi 
r=sqrt(rand);% get the radius vector magnitude (r^2 = x^2 + y^2, so r = sqrt(some number))
x=(rc*r)*cos(a);%scale the x to be a fraction of selected radius
y=(rc*r)*sin(a);%scale the y to be a fraction of selected radius
end