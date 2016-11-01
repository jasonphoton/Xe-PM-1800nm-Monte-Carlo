
function [IonAmp Yield] = fct_TolRate_GenAtom_TBIcor(E,t_grid,Ip,Z,Cnl,l,m,beta,SaturationSwitch,TbiSwitch,PlotOpt)

%% ionisation rate as shown in
% PHYSICAL REVIEW A 84, 053423(R) (2011)
% Theory of tunneling ionization of molecules: Weak-field asymptotics including dipole effects
% Oleg I. Tolstikhin, Toru Morishita, and Lars Bjoer Madsen

% define kappa
kappa = sqrt(2*abs(Ip));

%% structure factor
Qlm = (-1.^((abs(m)-m)*0.5)).*sqrt( ((2*l+1).*factorial(l+abs(m)))./(2.*factorial(l-abs(m))) );
a = ((((-1).^l).*Cnl.*(2.^(1-Z./kappa)).*Qlm))./(sqrt( (kappa.^(abs(m)+1)).*factorial(abs(m))) );      % equ. (73)
Gnm = (abs(a).^2);
%% field factor
W00   = 0.5.*kappa.*(((4.*kappa.^2)./abs(E)).^( (2*Z./kappa)-abs(m)-1)).*exp(-(2.*kappa.^3)./(3.*abs(E)));  % equ. (60), mu=0, nksi=0 

%% obi factor
WTBI  = exp(-beta.*((2*Z.^2)./kappa.^2).*(abs(E)./kappa.^3));

%% ionization amplitudes of all starting times
dt    = t_grid(2)-t_grid(1);

if SaturationSwitch ==0 && TbiSwitch==0
    IonAmp = Gnm.*W00;
    Yield  = sum(IonAmp.*dt);
end

if SaturationSwitch ==1 && TbiSwitch==0
    IonAmp = Gnm.*W00;
    IonAmp = IonAmp.*exp(-1.*cumsum(IonAmp).*dt);       % saturation
    Yield  = sum(IonAmp.*dt);
end

if SaturationSwitch ==0 && TbiSwitch==1
    IonAmp = Gnm.*W00.*WTBI;
    Yield  = sum(IonAmp.*dt);
end

if SaturationSwitch ==1 && TbiSwitch==1
    IonAmp = Gnm.*W00.*WTBI;
    IonAmp = IonAmp.*exp(-1.*cumsum(IonAmp).*dt);       % saturation
    Yield  = sum(IonAmp.*dt);
end


% plot the fields and so on
if PlotOpt==1
   figure;
   
   subplot(3,1,1)
   plot(t_grid,E,'.')
   xlabel('time (au)');
   ylabel('E (au)');
   title('E field');
   grid on

   subplot(3,1,2)
   plot(t_grid,IonAmp./max(IonAmp),'b.'); hold on
   %plot(t_grid,log10(IonAmp./max(IonAmp)),'r.'); hold on
   xlabel('time (au)');
   ylabel('rate (au)');
   title('ionization rate (au)');
   grid on
   
   subplot(3,1,3)
   plot(t_grid,IonAmp./max(IonAmp),'b-'); hold on
   plot(t_grid,E./max(E),'r-'); hold on
   xlabel('time (au)');
   ylabel('rate / field (norm)');
   title('ionization rate & field (au)');
   grid on
   
end


%figure;
% %% save
% savestr = ['sim_tot_yield_',num2str(IWcm,'%10.1e\n'),'_dt_',num2str(IWcm,4),'.mat'];
% save(savestr,'IWcm','n','CEP','Gleft','Gright');


