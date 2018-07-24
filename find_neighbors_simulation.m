function [like_cells,unlike_cells,no_neighbors]=find_neighbors_simulation(synthetic_trackID,next_iter_cells1,next_iter_cells2,nbrh_pxl_sz)
clear like_cells; clear unlike_cells
clear totest1
allcells = cat(1,next_iter_cells1,next_iter_cells2);
allcells_type2 = next_iter_cells2;
coordinate_same = struct;
coordinate_other = struct;
coordinate_same = [];
coordinate_other = [];

no_neighbors = 0;
local_neighbors = ipdm(next_iter_cells1(synthetic_trackID,1:2),allcells(:,1:2),'Result','Structure','Subset','Maximum','Limit',nbrh_pxl_sz);% get all the cells, closer than local_sz
[~,c]=find(local_neighbors.distance==min(local_neighbors.distance));
        local_neighbors.columnindex(c)=[];
        % local_neighbors.columnindex - cells within the rad of local_sz
        % now find which of these neighbors are cfp and wich are pluri
        totest1 = allcells(local_neighbors.columnindex',:);
        if isempty(totest1)
           % disp('no local neighbors found; increase local neighborhood R')
            like_cells=[];
            unlike_cells=[];
            no_neighbors = 1;
            return
        end
        totest = totest1;% 
        q1 = 1;
        q2 = 1;        
        for h=1:size(totest,1)
            % if one of the found local neighborhood cells are in the set of
            % celltype2, the nearest neighbor to it will be at zero distance (the cell is closest to itself);
            % counting how many of those, gives the number of neighbors of cell
            % type2, the rest (nearest-counter) is the other cell type
            tmp = ipdm(totest(h,1:2),allcells_type2(:,1:2),'Result','Structure','Subset','NearestNeighbor');%
            if tmp.distance == 0
                %disp('same')
                coordinate_other(q1).xyID = totest(h,1:3);% save the coordinate and index of the nearest like cell
                q1 = q1+1;
            else
                %disp('other')
                coordinate_same(q2).xyID = totest(h,1:3);% save the coordinate and index of the nearest unlike cell
                q2 = q2+1;
            end
        end
        
        if isempty(coordinate_same) && isempty(coordinate_other)
           % disp('no local neighbors found; increase local neighborhood R')
            like_cells=[];
            unlike_cells=[];
            no_neighbors = 1;
            return
        end
        
        if size(coordinate_same,2) == size(totest,1)% if all cells are same type
            coordinate_other = [];
            unlike_cells = [];
        else 
              unlike_cells = cat(1,coordinate_other.xyID);   

        end
                
        if size(coordinate_other,2) == size(totest,1)% if all cells are other type
            coordinate_same = [];
            like_cells = [];
        else
            like_cells =cat(1,coordinate_same.xyID);
        end
end