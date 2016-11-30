function [texit_ind texit] = fct_gen_texit(nsample,t_grid,IonAmp,PlotOpt)

% t_grid, temporal grid
% IonAmp, ionization amplitude on the temporal grid

% generate random indices according to the ionization amplitude
Pnorm    = IonAmp(1:1:length(IonAmp))./sum(IonAmp(1:1:length(IonAmp)));
texit_ind= fct_gen_distr(Pnorm,1,nsample);



% generate random texit by interpolation
texit    = t_grid(texit_ind)+rand(1,nsample)*(t_grid(2)-t_grid(1));

if PlotOpt==1
    
    hist_IonAmp = hist(texit,t_grid);

    figure;
    plot(t_grid,IonAmp./max(IonAmp),'b'); hold on
    plot(t_grid,hist_IonAmp./max(hist_IonAmp),'r'); hold on
    xlabel('time (au)');
    ylabel('rate / counts (norm)');
    grid on
end

end