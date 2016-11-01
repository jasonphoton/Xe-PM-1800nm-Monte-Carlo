function [hist_e_all hist_e_left hist_e_right] = fct_get_lrespec(v_sample,egrid,PlotOpt)

e_sample = 0.5.*v_sample.^2;

ind = find(v_sample<0);
e_sample_left = e_sample(ind);

ind = find(v_sample>0);
e_sample_right = e_sample(ind);

hist_e_all   = hist(e_sample,egrid);
hist_e_left  = hist(e_sample_left,egrid);
hist_e_right = hist(e_sample_right,egrid);

if PlotOpt==1
    figure;
    semilogy(egrid,(hist_e_all),'b'); hold on
    semilogy(egrid,(hist_e_left),'r'); hold on
    semilogy(egrid,(hist_e_right),'k'); hold on
    title('energy (au)')
    grid on
    xlabel('energy');
    ylabel('counts');
end