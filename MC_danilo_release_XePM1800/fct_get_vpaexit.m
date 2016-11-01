function [vpaxit] = fct_get_vpaexit(vpexit,PlotOpt)

% vpexit, perpendicular velocity at the exit of the tunnel
vpaxit = zeros(1,length(vpexit));

% vpaxit discret
vpaxitgrid = -2*max(vpexit):1e-2:2.*vpexit;

if PlotOpt==1
    
    hist_vpaxit      = hist(vpaxit,vpaxitgrid);
    
    figure;
    subplot(1,1,1)
    plot(vpaxitgrid,hist_vpaxit,'r'); hold on
    xlabel('angle (rand)');
    ylabel('counts');
    grid on
    axis tight
end

end