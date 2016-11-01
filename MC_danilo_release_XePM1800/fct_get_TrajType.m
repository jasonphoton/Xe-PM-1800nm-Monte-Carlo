function [rescind dirind] = fct_get_TrajType(trind_mat)

trind_mat_ = trind_mat(:,1);
% identify direct trajectories
ind    = find(trind_mat_~=0);
rescind = ind;

ind    = find(trind_mat_==0);
dirind = ind;


end