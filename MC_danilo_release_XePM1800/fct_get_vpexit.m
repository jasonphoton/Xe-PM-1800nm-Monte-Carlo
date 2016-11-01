function [vpexit_ind vpexit vp_grid] = fct_get_vpexit(Eexit_sample,Ip_sample,MeanOpt,PlotOpt)

% Eexit_sample, ionization field strength 
% Ip_sample, related ionization potential 
kappa_sample = sqrt(2.*abs(Ip_sample));
absE_sample  = abs(Eexit_sample);
c_sample     = sqrt(absE_sample./(2*kappa_sample));
vpexit_ind   = zeros(1,length(Eexit_sample));

% discretization of the perp momentum
vp_grid = 0:1e-3:(5*((2*sqrt(2*log(2)))*max(c_sample)));

if MeanOpt==0
    % loop over all E_samples and Ip_sample - this is the "correct" way
    for k=1:length(Eexit_sample)
        kappak = kappa_sample(k);     
        Ek     = absE_sample(k);
        Pk     = vp_grid.*(((4.*pi.*kappak)./(Ek)).*exp(- (kappak.*(vp_grid.^2))./(Ek)));
        Pknorm = Pk./sum(Pk);
        vpexit_ind(k) = fct_gen_distr(Pknorm,1,1);
    end
    vpexit = vp_grid(vpexit_ind)+rand(1,length(vpexit_ind))*(vp_grid(2)-vp_grid(1));
end

kappa_sample_mean = mean(kappa_sample);
absE_sample_mean  = mean(absE_sample);
Pmean = vp_grid.*(((4.*pi.*kappa_sample_mean./(absE_sample_mean)).*exp(- (kappa_sample_mean.*(vp_grid.^2))./(absE_sample_mean))));
    
if MeanOpt==1
    % option two, take an averaged field strength and generate it in one step
    Pmean_norm = Pmean./sum(Pmean);
    vpexit_ind_mean = fct_gen_distr(Pmean_norm,1,length(Eexit_sample));
    vpexit_mean = vp_grid(vpexit_ind_mean)+rand(1,length(vpexit_ind_mean))*(vp_grid(2)-vp_grid(1));
    vpexit_ind      = vpexit_ind_mean;
    vpexit          = vpexit_mean;
end

if PlotOpt==1
    
    hist_vpexit      = hist(vpexit,vp_grid);
%    hist_vpexit_mean = hist(vpexit_mean,vp_grid);
     
    figure;
    subplot(1,1,1)
    plot(vp_grid,Pmean./max(Pmean),'r'); hold on
    plot(vp_grid,hist_vpexit./max(hist_vpexit),'k'); hold on
%    plot(vp_grid,hist_vpexit_mean./max(hist_vpexit_mean),'b'); hold on
    xlabel('vp (au)');
    ylabel('counts');
    grid on
end

end