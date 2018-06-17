% simulation of sorting cells
% generate random points within circle, setup two different groups of
% cells, make sure they are not too close
clf
n = 200;% 
rc=700;
cellsz_pxl =40;% needs to be put in the paramfile for sorting simulations
close all;axis equal;hold on
clear tmp1; clear tmp2;clear cells1, clear cells2
[x,y,~] = cylinder(rc);plot(x(1,:),y(1,:),'--m');
% get the total number of randomly generated points
for i=1:n %n points inside the circle
[x, y]=rand_coord_circle(rc);%
if mod(i,2)==0 % even
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

[cells1,cells2] = check_close_cells(cellsz_pxl,tmp1,tmp2);

indx1 = (1:size(cells1))';
cells1 = cat(2,cells1(:,1:2),indx1);% assign indixes to cells
indx2 = (1:size(cells2))';
cells2 = cat(2,cells2(:,1:2),indx2);% assign indixes to cells
%%
[x,y,~]=cylinder(rc);
figure(2),plot(x(1,:),y(1,:),'--m');hold on
axis equal
figure(2),scatter(cells1(:,1),cells1(:,2),'LineWidth',2);hold on
scatter(cells2(:,1),cells2(:,2),'LineWidth',2);
box on
title(['After clean neighbors: type1: ' num2str(size(cells1,1)) 'pts   type2: ' num2str(size(cells2,1)) 'pts'])
%% next: calculate the quantity that will be changed as cells move randomly (closer or further to each other)
% then generate the increments to each cell coordinate such that they stay within the circle and not move too close to neighboring cells
% calculate the quantity again, iterate. This Q will be smaller if like
% cells are coming together, and larger, if unlike cells are close together
%TODO: make sure the cells don't move outside the circle
cells1_new = [];
cells2_new = [];
scale_rnd = 35;
coord_tmp =[];
cells1_new(:,1:2) = cells1(:,1:2) +scale_rnd*rand(size(cells1(:,1:2)));% increment coordinates of cells1 by randon numbrs
%make sure cells don't go out of the circle 
for jj=1:size(cells1_new,1)
    coord_tmp(jj) = (cells1_new(jj,1).*cells1_new(jj,1) + cells1_new(jj,2).*cells1_new(jj,2));
   % disp(power(coord_tmp(jj),0.5))
    if power(coord_tmp(jj),0.5)> rc
       disp('outside');
       cells1_new(jj,1:2)=cells1(jj,1:2);
       coord_tmp(jj) = (cells1_new(jj,1).*cells1_new(jj,1) + cells1_new(jj,2).*cells1_new(jj,2));
       disp(power(coord_tmp(jj),0.5))
    end
end
indx1 = (1:size(cells1_new))';
cells1_new = cat(2,cells1_new,indx1);% assign indixes to cells
cells2_new(:,1:2) = cells2(:,1:2) + scale_rnd*rand(size(cells2(:,1:2)));% increment coordinates of cells2 by randon numbrs

indx2 = (1:size(cells2_new))';
cells2_new = cat(2,cells2_new,indx2);

% figure(2),hold on,scatter(cells1_new(:,1),cells1_new(:,2),'LineWidth',3);hold on
% scatter(cells2_new(:,1),cells2_new(:,2),'LineWidth',3);
% text(cells2_new(:,1),cells2_new(:,2), num2str(cells2_new(:,3)));
% box on
%%
close all
[cells1_tmp,cells2_tmp] = check_close_cells(cellsz_pxl,cells1_new,cells2_new);
[x,y,~]=cylinder(rc);
figure(2),plot(x(1,:),y(1,:),'--m');hold on
axis equal
figure(2),hold on,scatter(cells1_tmp(:,1),cells1_tmp(:,2),'LineWidth',1.5);hold on
scatter(cells2_tmp(:,1),cells2_tmp(:,2),'LineWidth',1.5);
box on
% now the two cell coordinate cells to use for calculation of the quantity
% to be minimized
% cells1_tmp, cells2_tmp: cell positions after the increment was made
% cells1,cells2:          cell positions before the increment was made
% TODO: get the local neighborhood of each cell, calculate the Ising-like
% Hamiltonion (without external field)

nbrh_pxl_sz = 90;
synthetic_trackID = 2;% will be a loop variable
[like_cells,unlike_cells]=find_neighbors_simulation(synthetic_trackID,cells1_tmp,cells2_tmp,nbrh_pxl_sz);
figure(2), hold on
plot(cells1_tmp(synthetic_trackID,1),cells1_tmp(synthetic_trackID,2),'r.');hold on
text(cells1_tmp(synthetic_trackID,1)-5,cells1_tmp(synthetic_trackID,2)-5,['current, ' num2str(cells1_tmp(synthetic_trackID,3))],'FontSize',9)
if ~isempty(like_cells)
plot(like_cells(:,1),like_cells(:,2),'r.');hold on
text(like_cells(:,1)+5,like_cells(:,2)+5,num2str(like_cells(:,3)),'Color','r','FontSize',9);hold on

end
if ~isempty(unlike_cells)
plot(unlike_cells(:,1),unlike_cells(:,2),'b.');hold on
text(unlike_cells(:,1)+5,unlike_cells(:,2)+5,num2str(unlike_cells(:,3)),'Color','b','FontSize',9);hold on
end

% calculate the energy for each cell 'synthetic_trackID'  based on its neighborhood
% 

