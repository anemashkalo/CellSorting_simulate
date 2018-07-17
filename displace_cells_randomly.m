function [checked_data,random_move]=displace_cells_randomly(cells,rc,scale_rnd)
%move cells and make sure the cells don't move outside the circle
% input arguments: cells previous locations
% output = new locations
cells_new = [];
random_move=[];
random_move = scale_rnd*randn(size(cells(:,1:2)));
cells_new(:,1:2) = cells(:,1:2) +random_move;% increment coordinates of cells1 by randon numbrs

%make sure cells don't go out of the circle 
[checked_data]= stay_within_circle(cells_new, rc,cells);
indx = (1:size(checked_data))';
checked_data = cat(2,checked_data,indx);% assign indixes to cells

end