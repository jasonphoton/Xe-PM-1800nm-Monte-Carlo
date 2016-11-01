function [dirind] = fct_get_dirTraj(trind_mat)

trind_mat_ = trind_mat(:,1);

% identify direct trajectories
ind = find(trind_mat_==0);

end