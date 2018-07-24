function [cells1,cells2]=seed_cells_rand(rc,n,cellsz_pxl,celltyperatio_equal)

close all;axis equal;hold on
clear tmp1; clear tmp2;clear cells1, clear cells2
[x,y,~] = cylinder(rc);plot(x(1,:),y(1,:),'--m');
% get the total number of randomly generated points
if celltyperatio_equal ==1
    nn = 2;
else
    nn= 3;
end
for i=1:n %n points inside the circle
[x, y]=rand_coord_circle(rc);%
if mod(i,nn)==0 % even   % 3
tmp1(i,1:2) = [x y];
else           % uneven
tmp2(i,1:2) = [x y];
end
end
% make sure there are no centroids within cellsz_pxl pixels distance from any other
% centroid (no overlapping cells), this effectively gives cell a size
% works with already randomly generated cells, no new data points are
% generated here
indx1 = (1:size(tmp1))';
tmp1 = cat(2,tmp1,indx1);% assign indixes to cells
indx2 = (1:size(tmp2))';
tmp2 = cat(2,tmp2,indx2);% assign indixes to cells
[cells1,cells2] = check_close_cells(cellsz_pxl,tmp1,tmp2);% usethis when initializing the cells
% this function will remove the cells that are too close to each other
indx1 = (1:size(cells1))';
cells1 = cat(2,cells1(:,1:2),indx1);% assign indixes to cells
indx2 = (1:size(cells2))';
cells2 = cat(2,cells2(:,1:2),indx2);% assign indixes to cells

% CLEAN UP plot
close all
[x,y,~]=cylinder(rc);
figure(1),plot(x(1,:),y(1,:),'--m');hold on
axis equal
figure(1),scatter(cells1(:,1),cells1(:,2),'LineWidth',2);hold on
scatter(cells2(:,1),cells2(:,2),'LineWidth',2);
box on
title(['After clean neighbors: type1: ' num2str(size(cells1,1)) 'pts   type2: ' num2str(size(cells2,1)) 'pts'])

end