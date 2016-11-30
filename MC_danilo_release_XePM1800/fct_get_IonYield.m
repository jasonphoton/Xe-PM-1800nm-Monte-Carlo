
function [IonYield] = fct_get_IonYield(atom,IWcm,CEP,pulse_)
%% parameters for the laser pulse 
c0       = 299792458./(2.1876912633e6);          % au, speed of light 
n        = pulse_(1);                            % propto pulse duration
N        = pulse_(2);
wvlm     = pulse_(3);                            % 
wvl      = wvlm/(5.2917720859e-11);              % [au]
omega    = (2*pi*c0)/wvl;                        % [au], angular frequency
% IWcm     = 1e14;               %
CEP      = 0;                                    % au
dt       = 1e-1;                                 % au
Up       = (0.09337.*IWcm.*(wvlm.^2))./27.211;   % au

%% define the field
[Env E I A ALPHA BETA v r tgrid E_fh] = fct_get_EnvEIAAlphBetvr_sin2_nN(dt,wvlm,IWcm,CEP,n,N);

%% function handle for bsi corrections
Fc = @(kappa_,Z_,m_) (kappa_.^4)./(8.*(2.*Z_-kappa_.*(m_+1)));

%%  All atoms, pure rate, hydrogen rate, gen atom, wo_BI, wo_me
SaturationSwitch = 1;
PlotOpt          = 0;
tbiSwitch        = 1;
AtomString       = atom;
[Ip kappa Z Cnl l beta alphaN alphaI] = fct_get_Atom(AtomString);

Ip    = Ip.*ones(1,length(E));
Z     = Z.*ones(1,length(E));
Cnl   = Cnl.*ones(1,length(E));
l     = l*ones(1,length(E));
beta  = beta*ones(1,length(E));
m     = 0*ones(1,length(E));
alphaN= alphaN.*ones(1,length(E));
alphaI= alphaI.*ones(1,length(E));
%Ip    = Ip-0.5.*(alphaN-alphaI).*E.^2;
Fc_GenA_H1s_wo_tbi_wo_me  = Fc(kappa,Z,m);
[IonAmp_GenA_H1s_wo_tbi_wo_me Yield_GenA_H1s_wo_tbi_wo_me] = fct_TolRate_GenAtom_TBIcor(E,tgrid,Ip,Z,Cnl,l,m,beta,SaturationSwitch,tbiSwitch,PlotOpt);
IonYield = Yield_GenA_H1s_wo_tbi_wo_me;
end
