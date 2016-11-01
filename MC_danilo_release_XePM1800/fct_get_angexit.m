function [chi] = fct_get_angexit(vpexit,PlotOpt)

% vpexit, perpendicular velocity at the exit of the tunnel
chi = rand(1,length(vpexit)).*2.*pi - pi;

% v0x = vpexit.*cos(chi_dir);
% v0y = vpexit.*sin(chi_dir);

% discretization of the angle
anglegrid=-pi:1e-2:pi;

if PlotOpt==1
    
    hist_chi      = hist(chi,anglegrid);
    
    figure;
    subplot(1,1,1)
    plot(anglegrid,hist_chi,'r'); hold on
    xlabel('angle (rand)');
    ylabel('counts');
    grid on
    axis tight
end

end