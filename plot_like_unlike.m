function plot_like_unlike(synthetic_trackID,cells1,like_cells,unlike_cells,fig)
% synthetic_trackID - ID of a cell from data cells1 (second arg)
% like cells = same as cells1 (second arg)
figure(fig), hold on
plot(cells1(synthetic_trackID,1),cells1(synthetic_trackID,2),'r.');hold on
text(cells1(synthetic_trackID,1)-5,cells1(synthetic_trackID,2)-5,['current, ' num2str(cells1(synthetic_trackID,3))],'FontSize',9)
if ~isempty(like_cells)
plot(like_cells(:,1),like_cells(:,2),'r.');hold on
text(like_cells(:,1)+5,like_cells(:,2)+5,num2str(like_cells(:,3)),'Color','r','FontSize',9);hold on
end
if ~isempty(unlike_cells)
plot(unlike_cells(:,1),unlike_cells(:,2),'b.');hold on
text(unlike_cells(:,1)+5,unlike_cells(:,2)+5,num2str(unlike_cells(:,3)),'Color','b','FontSize',9);hold on
end
end