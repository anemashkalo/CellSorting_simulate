function  [cells1,cells2] = check_close_cells2(cellsz_pxl,tmp1,tmp2,cells1prev,cells2prev,random_move1,random_move2)
% make sure there are no centroids within N pixel distance from any other
% centroid (no overlapping cells) after the random dispalcement. If they
% overlap, return the coordinate of the cell to its previous position

clear type1; clear type2;clear cells1;clear cells2
new_sz1 = size(nonzeros(tmp1(:,1)),1);
type1(1:new_sz1,1) = nonzeros(tmp1(:,1));
type1(1:new_sz1,2) = nonzeros(tmp1(:,2));
type1(1:new_sz1,3) = tmp1(1:new_sz1,3);

new_sz2 = size(nonzeros(tmp2(:,1)),1);
type2(1:new_sz2,1) = nonzeros(tmp2(:,1));
type2(1:new_sz2,2) = nonzeros(tmp2(:,2));
type2(1:new_sz2,3) = tmp2(1:new_sz2,3);

clear too_close_diff
clear too_close_1
clear too_close_2

too_close_diff= ipdm(type1,type2,'Result','Structure','Subset','Maximum','Limit',cellsz_pxl);% get cells that are closer then realistic cell-cell distance
too_close_1= ipdm(type1,'Result','Structure','Subset','Maximum','Limit',cellsz_pxl);% get cells that are closer then realistic cell-cell distance
too_close_2= ipdm(type2,'Result','Structure','Subset','Maximum','Limit',cellsz_pxl);% get cells that are closer then realistic cell-cell distance

title(['type1: ' num2str(new_sz1) 'pts   type2: ' num2str(new_sz2) 'pts'])
% check which centroids are too close and keep them at their previous position

for jj=1:size(too_close_diff.distance,1)
% figure(1), hold on
if too_close_diff.distance(jj)>0
% make the move much smaller, for cells that end up 
% running into neighboring cell
type2(too_close_diff.columnindex(jj),1:2)=cells2prev(too_close_diff.columnindex(jj),1:2)-random_move2(too_close_diff.columnindex(jj),1:2).*0.2;
%[checked_data]= stay_within_circle(cells_new, rc,cells);
type1(too_close_diff.rowindex(jj),1:2)=cells1prev(too_close_diff.rowindex(jj),1:2)-random_move1(too_close_diff.rowindex(jj),1:2).*0.2;
end
end
% make the move much smaller, for cells that end up 
% running into neighboring cell
[var1,~] = find(too_close_1.distance>0);
for jj=1:size(too_close_1.distance,1)
    if too_close_1.distance(jj)>0        
 type1(too_close_1.rowindex(jj),1:2)=cells1prev(too_close_1.rowindex(jj),1:2)-random_move1(too_close_1.rowindex(jj),1:2).*0.2;
        end
end
% make the move much smaller (previous displacement * 0.2), for cells that end up 
% running into neighboring cell
for jj=1:size(too_close_2.distance,1)
    if too_close_2.distance(jj)>0        
 type2(too_close_2.rowindex(jj),1:2)=cells2prev(too_close_2.rowindex(jj),1:2)-random_move2(too_close_2.rowindex(jj),1:2).*0.2;
   
    end
end

% final simulated cell mix, 50:50 ratio before cleanup
cells1=zeros(size(nonzeros(type1(:,1)),1),2);
cells2=zeros(size(nonzeros(type2(:,1)),1),2);

cells1(:,1)=nonzeros(type1(:,1));cells1(:,2)=nonzeros(type1(:,2));cells1(:,3)=nonzeros(type1(:,3));
cells2(:,1)=nonzeros(type2(:,1));cells2(:,2)=nonzeros(type2(:,2));cells2(:,3)=nonzeros(type2(:,3));

end