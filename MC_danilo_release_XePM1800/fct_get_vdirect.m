function [vxf vyf vzf] = fct_get_vdirect(A,texit_ind_sample,v0x_sample,v0y_sample,v0z_sample,Up,PlotOpt)

% calculate final momentum direct
vxf = v0x_sample;
vyf = v0y_sample;
vzf = v0z_sample - A(texit_ind_sample)';

if PlotOpt==1
    figure;
    vfgrid = -1.*max(sqrt(vxf.^2+vyf.^2+vzf.^2)):1e-3:max(sqrt(vxf.^2+vyf.^2+vzf.^2));
    hist_vxf = hist(vxf,vfgrid);
    hist_vyf = hist(vyf,vfgrid);
    hist_vzf = hist(vzf,vfgrid);

    % calculate the energy
    efgrid = 0:1e-2:(max(0.5*vzf.^2));
    ezf      = 0.5.*vzf.^2;
    ind      = find(vzf<0);
    ezf_left = ezf(ind);
    ind      = find(vzf>0);
    ezf_right= ezf(ind);
    hist_efz = hist(ezf,efgrid);
    hist_efz_left  = hist(ezf_left,efgrid);
    hist_efz_right = hist(ezf_right,efgrid);
    
    subplot(2,1,1)
    plot(vfgrid,hist_vxf./max(hist_vxf),'b'); hold on
    plot(vfgrid,hist_vyf./max(hist_vyf),'r'); hold on
    plot(vfgrid,hist_vzf./max(hist_vzf),'k'); hold on
    title('velocity direct')
    grid on
    
    subplot(2,1,2)
    semilogy(efgrid./Up,(hist_efz),'b'); hold on
    semilogy(efgrid./Up,(hist_efz_left),'r'); hold on
    semilogy(efgrid./Up,(hist_efz_right),'k'); hold on
    title('energy parallel to pol Up')
    grid on
    xlabel('energy Up'); ylabel('counts');
end

end