%% initialize simulation
%  simulation of sorting cells
% generate random points within circle, setup two different groups of
% cells, make sure they are not too close
clf
close all
n = 200;% total number of cells
rc=500; % radius of circle in pixels
cellsz_pxl =50;% cell size in pixels 
%synthetic_trackID = 2;% cell ID to watch
nbrh_pxl_sz = 75; % local neighborhood size
scale_rnd = 15; % multiply rand number btw 0 and 1  (randomly generated cell displacements)
J = 1;
[cells1,cells2]=seed_cells_rand(rc,n,cellsz_pxl);% generates a set of cells of two tyoes, intexed with ID and not overlapping
 fig = 1;
% calculate H1 for each cell with its neighborhood for cells1 set
 H1 = zeros(size(cells1,1),1);
 for synthetic_trackID=1:size(cells1,1)
[like_cells,unlike_cells,no_neighbors]=find_neighbors_simulation(synthetic_trackID,cells1,cells2,nbrh_pxl_sz);
%disp([ 'before' num2str(cells1(synthetic_trackID,:))]);
axis equal;hold on
[x,y,~] = cylinder(rc);plot(x(1,:),y(1,:),'--m');hold on
scatter(cells1(:,1),cells1(:,2),'LineWidth',2);hold on
scatter(cells2(:,1),cells2(:,2),'LineWidth',2);hold on
plot_like_unlike(synthetic_trackID,cells1,like_cells,unlike_cells,fig);hold on;box on
if no_neighbors == 0
H1(synthetic_trackID)=calculate_interactionE(like_cells,unlike_cells,J);
end

 end
 initE1=H1 ;
% same for cells2

% calculate H1 for each cell with its neighborhood for cells2 set
 H1 = zeros(size(cells2,1),1);
 for synthetic_trackID=1:size(cells2,1)
[like_cells,unlike_cells,no_neighbors]=find_neighbors_simulation(synthetic_trackID,cells2,cells1,nbrh_pxl_sz);
%disp([ 'before' num2str(cells1(synthetic_trackID,:))]);
axis equal;hold on
[x,y,~] = cylinder(rc);plot(x(1,:),y(1,:),'--m');hold on
scatter(cells1(:,1),cells1(:,2),'LineWidth',2);hold on
scatter(cells2(:,1),cells2(:,2),'LineWidth',2);hold on
plot_like_unlike(synthetic_trackID,cells2,like_cells,unlike_cells,fig);
if no_neighbors == 0
H1(synthetic_trackID)=calculate_interactionE(like_cells,unlike_cells,J);
end

 end
  initE2=H1 ;
  
cells1 = cat(2,cells1,initE1);
cells2 = cat(2,cells2,initE2);

%% loop to update cell's position and check if will accept ot reject thenew coordinate

colormap = jet;
kT =  2.5; %27e-3 meV, kT at 37 C
MCstep = 40;
cellstomove = 10;
dH1 = zeros(MCstep,1);% size as the number of cells
dH2 = zeros(MCstep,1);% size as the number of cells
MC_iter_data = struct;
MC_iter_data2 = struct;
clear acceptance_p2
clear acceptance_p1
for MCiter = 1:MCstep % Monte Carlo iterations
    % cells1 - cell positions before random move
    % cells2 - cell positions before random move
[cells1_new,random_move1]=displace_cells_randomly(cells1,rc,scale_rnd);
[cells2_new,random_move2]=displace_cells_randomly(cells2,rc,scale_rnd);

[cells1_tmp,cells2_tmp] = check_close_cells2(cellsz_pxl,cells1_new,cells2_new,cells1,cells2,random_move1,random_move2);
cells1_tmp = cat(2,cells1_tmp,zeros(size(cells1_tmp,1),1));% init the individual interaction energy column with zeros
cells2_tmp = cat(2,cells2_tmp,zeros(size(cells2_tmp,1),1));
%cellstomove
for synthetic_trackID = 1:size(cells1_tmp,1) % tracked cells is always from from the second argument of find_neighbors_simulation
acceptance_p1 = [];
[like_cells,unlike_cells,no_neighbors]=find_neighbors_simulation(synthetic_trackID,cells1_tmp,cells2_tmp,nbrh_pxl_sz);
if no_neighbors == 1 % always accept the move if cell had no neighbors
    cells1(synthetic_trackID,:)=cells1_tmp(synthetic_trackID,:); % update the new cell coordinate and energy 
else
H2(synthetic_trackID)=calculate_interactionE(like_cells,unlike_cells,J);% 
dH1(synthetic_trackID) = H2(synthetic_trackID)-cells1(synthetic_trackID,4);% interaction energy difference for given cell
% disp(['prev. cell count. type1: ' num2str(size(cells1,1)) 'pts   type2: ' num2str(size(cells2,1)) 'pts']);
if dH1(synthetic_trackID)<=0 % if the local cell interaction energy decreased
     cells1(synthetic_trackID,:)=cells1_tmp(synthetic_trackID,:); % update the new cell coordinate and energy 
else 
     acceptance_p1 = exp(-dH1(synthetic_trackID)/kT);  % acceptance ratio
     %disp(acceptance_p1);
     u = rand;
     if u<= acceptance_p1
         cells1(synthetic_trackID,:)=cells1_tmp(synthetic_trackID,:); %
%      else
%          cells1(synthetic_trackID,:)=cells1(synthetic_trackID,:); %
         
     end
end
end
MC_iter_data(MCiter).dat = cells1;
[x,y,~]=cylinder(rc);
figure(2),plot(x(1,:),y(1,:),'--m');hold on
axis equal
figure(2),hold on,scatter(cells1(:,1),cells1(:,2),'LineWidth',1.5);hold on
scatter(cells2(:,1),cells2(:,2),'LineWidth',1.5);
box on
fig = 2;
plot_like_unlike(synthetic_trackID,cells1,like_cells,unlike_cells,fig)
end
% ---- other cell type 
% clear necessary stuff here
H2 = [];
like_cells = [];
unlike_cells = [];
acceptance_p2=[];
clear synthetic_trackID

for synthetic_trackID = 1:size(cells2_tmp,1) % tracked cells is always from from the second argument of find_neighbors_simulation
[like_cells,unlike_cells,no_neighbors]=find_neighbors_simulation(synthetic_trackID,cells2_tmp,cells1_tmp,nbrh_pxl_sz);
if no_neighbors == 1 % always accept the move if cell had no neighbors
         cells2(synthetic_trackID,:)=cells2_tmp(synthetic_trackID,:); % update the new cell coordinate TODO: this will update all cell coordinates though, regardless of dH for those cells
else
[H2(synthetic_trackID)]=calculate_interactionE(like_cells,unlike_cells,J);% new interaction energy
dH2(synthetic_trackID) = H2(synthetic_trackID)-cells2(synthetic_trackID,4);% interaction energy difference for given cell

if dH2(synthetic_trackID)<=0 % if the local cell interaction energy decreased
     cells2(synthetic_trackID,:)=cells2_tmp(synthetic_trackID,:); % update the new cell coordinate TODO: this will update all cell coordinates though, regardless of dH for those cells
else 
     acceptance_p2= exp(-(dH2(synthetic_trackID))/kT);  % acceptance probability
     %disp(acceptance_p2);
     u = rand;
     if u<=acceptance_p2% accept the new coordinate
         cells2(synthetic_trackID,:)=cells2_tmp(synthetic_trackID,:); %
%      else
%          cells2(synthetic_trackID,:) =  cells2(synthetic_trackID,:);
     end
end
end
MC_iter_data2(MCiter).dat = cells2;
% disp(['current cell count. type1: ' num2str(size(cells1,1)) 'pts   type2: ' num2str(size(cells2,1)) 'pts'])
[x,y,~]=cylinder(rc);
figure(3),plot(x(1,:),y(1,:),'--m');hold on
axis equal
figure(3),hold on,scatter(cells1(:,1),cells1(:,2),'LineWidth',1.5);hold on
fig = 3;
scatter(cells2(:,1),cells2(:,2),'LineWidth',1.5);
box on
plot_like_unlike(synthetic_trackID,cells2,like_cells,unlike_cells,fig)
end
end
%% plot cells1
close all
N = 4;
celltype1_simulatied = struct;
for syn_ID =1:size(cells1_tmp,1)%cellstomove
maketrack = [];
for MCiter=1:MCstep-2
    maketrack(MCiter,1:2) = [MC_iter_data(MCiter).dat(syn_ID,1) MC_iter_data(MCiter).dat(syn_ID,2)];
    maketrack(MCiter,3)= MC_iter_data(MCiter).dat(syn_ID,4); % interaction energy

end
celltype1_simulatied(syn_ID).dat = maketrack;

[x,y,~]=cylinder(rc);
figure(4),plot(x(1,:),y(1,:),'--m');hold on
axis equal
figure(4),hold on,scatter(cells1(:,1),cells1(:,2),'LineWidth',1.5);hold on
scatter(cells2(:,1),cells2(:,2),'LineWidth',1.5);hold on
box on
figure(4), plot(maketrack(1,1),maketrack(1,2),'rp','MarkerFaceColor','r');hold on
figure(4), plot(maketrack(1:end,1),maketrack(1:end,2),'-m.');hold on
text(maketrack(1,1)+3,maketrack(1,2)+3,num2str(syn_ID),'Color','r');
figure(7), plot(maketrack(1,1),maketrack(1,2),'rp','MarkerFaceColor','b','Markersize',15);hold on
figure(8), plot(maketrack(end,1),maketrack(end,2),'rp','MarkerFaceColor','b','Markersize',15);hold on
%xlim([-300 300]);ylim([-300 300]);
end
%title('All cell moves are accepted')
title(['Monte Carlo steps: ' num2str(MCstep) 'Cells moving: ' num2str(cellstomove) ]);
%% plot cells2
%close all
celltype2_simulatied = struct;
for syn_ID =1:size(cells2_tmp,1)%cellstomove
maketrack = [];
for MCiter=1:MCstep-2
    maketrack(MCiter,1:2) = [MC_iter_data2(MCiter).dat(syn_ID,1) MC_iter_data2(MCiter).dat(syn_ID,2)];
    maketrack(MCiter,3)= MC_iter_data2(MCiter).dat(syn_ID,4); % interaction energy
end
celltype2_simulatied(syn_ID).dat = maketrack;

[x,y,~]=cylinder(rc);
figure(4),plot(x(1,:),y(1,:),'--m');hold on
axis equal
figure(4),hold on,scatter(cells1(:,1),cells1(:,2),'LineWidth',1.5);hold on
scatter(cells2(:,1),cells2(:,2),'LineWidth',1.5);hold on
box on
figure(4), plot(maketrack(1,1),maketrack(1,2),'rp','MarkerFaceColor','r');
figure(7), plot(maketrack(1,1),maketrack(1,2),'rp','MarkerFaceColor','r','Markersize',15);hold on
figure(8), plot(maketrack(end,1),maketrack(end,2),'rp','MarkerFaceColor','r','Markersize',15);hold on

figure(4), plot(maketrack(1:end,1),maketrack(1:end,2),'-b.');
text(maketrack(1,1)+4,maketrack(1,2)+4,num2str(syn_ID),'Color','r');
%xlim([-300 300]);ylim([-300 300]);
end
%title('All cell moves are accepted')
title(['Monte Carlo steps: ' num2str(MCstep) 'Cells moving: ' num2str(cellstomove) ]);
title(['Monte Carlo steps: ' num2str(MCstep) 'Cells moving: All' ]);

[x,y,~]=cylinder(rc);
figure(7),plot(x(1,:),y(1,:),'--m');hold on
axis equal
title('before sorting')
[x,y,~]=cylinder(rc);
figure(8),plot(x(1,:),y(1,:),'--m');hold on
axis equal
title('after sorting')
%% analyze motion of simulated cells
%load('simulated_data_allmovesallowed.mat')
paramfile = 'C:\Users\Nastya\Desktop\FromGithub\CellTracker\paramFiles\setUserParamTrackSortingAN_20X.m';
delta_t = 60; % then the time will be in iteration number (not hours)
for syn_ID2 =1:size(celltype2_simulatied,2)
for syn_ID1 = 1:size(celltype1_simulatied,2)
[d_btw]=get_dist_btw_cellcenters(celltype1_simulatied(syn_ID1).dat(:,1:2),celltype2_simulatied(syn_ID2).dat(:,1:2),paramfile,delta_t);
if d_btw(1)<40
figure(5), plot(d_btw(:,2),d_btw(:,1),'-*m');hold on;box on
end
end
end
ylim([0 150])
title('Distance between two simulated cells as they move')
xlabel('time,Monte Carlo Iterations')
ylabel('UnLike cells separation, microns')

%% distance between like cells 

for syn_ID1 = 1:size(celltype1_simulatied,2)-1
[d_btw]=get_dist_btw_cellcenters(celltype1_simulatied(syn_ID1).dat(:,1:2),celltype1_simulatied(syn_ID1+1).dat(:,1:2),paramfile,delta_t);
if (d_btw(1)>30) && (d_btw(1)<150)
figure(6), plot(d_btw(:,2),d_btw(:,1),'-*m');hold on;box on
end
end
%end
ylim([0 150])
title('Distance between two simulated cells as they move')
xlabel('time,Monte Carlo Iterations')
ylabel('Like cells separation, microns')

for syn_ID2 = 1:size(celltype2_simulatied,2)-1
[d_btw]=get_dist_btw_cellcenters(celltype2_simulatied(syn_ID2).dat(:,1:2),celltype2_simulatied(syn_ID2+1).dat(:,1:2),paramfile,delta_t);
if (d_btw(1)>30) && (d_btw(1)<150)
figure(7), plot(d_btw(:,2),d_btw(:,1),'-*m');hold on;box on
end
end
%end
ylim([0 150])
title('Distance between two simulated cells as they move')
xlabel('time,Monte Carlo Iterations')
ylabel('Like cells separation, microns')

%% plot MSD vs time for all simulated tracks
for synID=1:size(celltype1_simulatied,2)
[msd1]=getMSD(celltype1_simulatied,paramfile,delta_t,synID);
figure(8), plot(msd1(synID).dat(:,2),msd1(synID).dat(:,1),'-p');hold on
box on
end
for synID=1:size(celltype2_simulatied,2)
[msd2]=getMSD(celltype2_simulatied,paramfile,delta_t,synID);
figure(8), plot(msd2(synID).dat(:,2),msd2(synID).dat(:,1),'-p');hold on
box on
end
xlabel('time,Monte Carlo Iterations')
ylabel('MSD, microns^2')
title('Simulated cell motion')



