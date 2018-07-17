function [H]=calculate_interactionE(like_cells,unlike_cells,J)
% calculate the energy for each cell 'synthetic_trackID'  based on its neighborhood
% 1. H(t1) - energy before move
% 2. H(2) - energy after move
% 3. dH - energy gain upon move
% 4. H = J*sum(sisj)   J, interaction strength; sisj=0 if it;s a pair...
% of like cells, and one, if it's a pair of unlike cells
% 5. if dH upon move is <=0, the accept move and update the new cell
% coordinate
% 5a. if dH is >0, accept/reject move with probability = exp(-dH/kT)
clear H

H = J*(size(like_cells,1)+size(unlike_cells,1));
end