function [vxf vyf vzf] = fct_get_vscatt(A,texit_ind_sample,v0x_sample,v0y_sample,v0z_sample,PlotOpt)

% calculate final momentum direct
vxf = v0x_sample;
vyf = v0y_sample;
vzf = v0z_sample + A(texit_ind_sample);

vfgrid = -3:0.05:3;
hist_vxf = hist(vxf,vfgrid);
hist_vyf = hist(vyf,vfgrid);
hist_vzf = hist(vzf,vfgrid);

if PlotOpt==1
    
end

end