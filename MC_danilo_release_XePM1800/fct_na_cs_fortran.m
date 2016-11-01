function [ result ] = fct_na_cs_fortran( q_input,theta,energy,cs )
% function takes a q_input return momentum (vector) and gives out a matrix with
% cross sections corresponding to the return momenta.
% result has length(q_input) lines and length(angle_grid) columns.  the
% angle grid was defined in the fortran cs program and MUST be the same as
% the anglegrid in fct_gen_resc2D_fromList_v1_fortran.m. this is achieved
% by passing the anglegrid to fct_gen_resc2D_fromList_v1_fortran.m as a
% variable called winkel_total in fct_gen_resc2D_fromList_v1_fortran.m
% by Danilo Zille

 energy = energy / 27.211; %convert eV to au

result=zeros(length(q_input),length(cs(1,:)));

for ind = 1:length(q_input) % scan over array of return momenta
    
energy_input =0.5*q_input(ind)^2; % return energy for this iteration

try
    inds = find((energy-energy_input)<=0);
    choose = inds(end);    % this is the index on the energy grid, which the returning electron has    
     result(ind,:) = cs(choose,:);  %return the cs at given return momentum as a vector
catch % if incoming energy smaller than smallest fortran energy, take the 1st cs value (usually 1eV)
      result(ind,:) = cs(1,:);  %return the cs at given return momentum as a vector


end


end

