function Main_loop_save_v7_spineffects
% close all
 clear all

% pulse characteristics
n        = 5;                                   % propto pulse duration, (cos^2)
N        = 0;                                   % propto flat part
wvlm     = 9000e-9;                            %
atom     = 'He1s2s-triplet'; % %atom     = 'He1s2s-singlet';    %'He1s2s-triplet'
s        = 3/2;
pulse    = [n N wvlm];
dt       = 0.1;               % temporaresolution
kmax     = 200;              % number of cep/int settings for vol avg
nsample  = 1e8;             % number of particles per cep/int setting


% focal characteristics
I0    = 1e12;             % I0, [W/cm^2]
w0    = 0.087;            % w0, [mm]
Zr    = 3.3;              % Zr, [mm]

% get hist of intensity distribution in focus
% IWcm_grid_Vavg       = 5e11:0.1e10:I0;
IWcm_grid_Vavg       = 1e11:1e11:I0;   % for single intensity!!!
NI_Vavg              = length(IWcm_grid_Vavg);
PlotOpt              =0;
[Iaxis_Vavg,hist_Isample_Vavg] = fct_get_IhistVavg_zR(w0,I0,Zr,IWcm_grid_Vavg,PlotOpt);

% ion yield in igrid from intensity distribution
IonYield = zeros(1,length(hist_Isample_Vavg));
for l=1:1:NI_Vavg
    [Yield]  = fct_get_IonYield(atom,Iaxis_Vavg(l),0,pulse);        % CEP = 0 ist a good estimate
    IonYield(l) = Yield;  
end
hist_Isample_Vavg = hist_Isample_Vavg./max(hist_Isample_Vavg);
IonYield          = IonYield./max(IonYield);
Y_Vavg_Iyield     = (hist_Isample_Vavg.*IonYield)./sum(hist_Isample_Vavg.*IonYield);

% prep vol averaged ion yield for big monte carlo
Pnorm       = Y_Vavg_Iyield;
intens_ind  = fct_gen_distr(Pnorm,1,1e3);
Isample     = Iaxis_Vavg(intens_ind); 
hist_Isample= hist(Isample,Iaxis_Vavg);
%   plot(Iaxis_Vavg,hist_Isample)
 
% prepare name1 for this run
savestr =  ['n',num2str(n),'_N',num2str(N),'_',atom,'_',num2str(wvlm./1e-9),'_Ip',num2str(I0,'%10.1e\n')];

% prepare histograms for final i averaged and cep averaged result
Up       = (0.09337.*1e12.*(wvlm.^2))./27.211;   % au


resol = 0.01; %0.002 ist ok, aber ehr zu fein, 0.0857 == 0.1eV
lim_correction = 1.4;
if wvlm == 1800e-9
y_lim = sqrt(2*10*Up)-sqrt(2*3.17*Up);
vpagrid = -0.5:resol:0.5;
vpegrid = -lim_correction*y_lim:resol:lim_correction*y_lim;
else
y_lim = sqrt(2*10*Up)-sqrt(2*3.17*Up);
vpagrid_01 = -2.7:resol:2.7;
vpegrid_01 = -lim_correction*y_lim:resol:lim_correction*y_lim;
end

resol = 0.020; %0.002 ist ok, aber ehr zu fein, 0.0857 == 0.1eV
lim_correction = 1.4;
if wvlm == 1800e-9
y_lim = sqrt(2*10*Up)-sqrt(2*3.17*Up);
vpagrid = -0.5:resol:0.5;
vpegrid = -lim_correction*y_lim:resol:lim_correction*y_lim;
else
y_lim = sqrt(2*10*Up)-sqrt(2*3.17*Up);
vpagrid_025 = -2.7:resol:2.7;
vpegrid_025 = -lim_correction*y_lim:resol:lim_correction*y_lim;
end
% % % % 
% % % % resol = 0.05; %0.002 ist ok, aber ehr zu fein, 0.0857 == 0.1eV
% % % % lim_correction = 1.4;
% % % % if wvlm == 1800e-9
% % % % y_lim = sqrt(2*10*Up)-sqrt(2*3.17*Up);
% % % % vpagrid = -0.5:resol:0.5;
% % % % vpegrid = -lim_correction*y_lim:resol:lim_correction*y_lim;
% % % % else
% % % % y_lim = sqrt(2*10*Up)-sqrt(2*3.17*Up);
% % % % vpagrid_05 = -2.7:resol:2.7;
% % % % vpegrid_05 = -lim_correction*y_lim:resol:lim_correction*y_lim;
% % % % end
% % % % 
% % % % resol = 0.1; %0.002 ist ok, aber ehr zu fein, 0.0857 == 0.1eV
% % % % lim_correction = 1.4;
% % % % if wvlm == 1800e-9
% % % % y_lim = sqrt(2*10*Up)-sqrt(2*3.17*Up);
% % % % vpagrid = -0.5:resol:0.5;
% % % % vpegrid = -lim_correction*y_lim:resol:lim_correction*y_lim;
% % % % else
% % % % y_lim = sqrt(2*10*Up)-sqrt(2*3.17*Up);
% % % % vpagrid_1 = -2.7:resol:2.7;
% % % % vpegrid_1 = -lim_correction*y_lim:resol:lim_correction*y_lim;
% % % % end

hist_vyvz_resc_loop_01 = zeros(length(vpegrid_01),length(vpagrid_01));
hist_vyvz_resc_loop_025 = zeros(length(vpegrid_025),length(vpagrid_025));
% % % % hist_vyvz_resc_loop_05 = zeros(length(vpegrid_05),length(vpagrid_05));
% % % % hist_vyvz_resc_loop_1 = zeros(length(vpegrid_1),length(vpagrid_1));

savename2    = atom;

for k=1:kmax
    
    tic
    Pnorm          = Y_Vavg_Iyield;
    intens_ind     = fct_gen_distr(Pnorm,1,1);
    Isample(k)     = 5.5e11;%Iaxis_Vavg(intens_ind);  %!!!!!!!!!!!!!!!! PUT SINGLE INT HERE !!!!!!!!!!
    CEPsample(k)   = 0;%rand(1)*2.*pi;
    savename1      = [num2str(k),'_',savestr,'.mat'];
    
    % get the momenta for one set of CEP and intensity
    
    %nsample = 1e5;
    [vxf_all vyf_all vzf_all vxf_dir vyf_dir vzf_dir vxf_resc vyf_resc vzf_resc T] = singleICEP_v0(Isample(k),CEPsample(k),pulse,atom,nsample,savename1,dt,s); 
   
    % histogram for current i & cep averaged results
    PlotOpt = 0;
    [hist_vxvy_resc hist_vxvz_resc hist_vyvz_resc_01] = fct_illstr_Projections(vpagrid_01,vpegrid_01,vxf_resc,vyf_resc,vzf_resc,PlotOpt);
    [hist_vxvy_resc hist_vxvz_resc hist_vyvz_resc_025] = fct_illstr_Projections(vpagrid_025,vpegrid_025,vxf_resc,vyf_resc,vzf_resc,PlotOpt);
% % % % %     [hist_vxvy_resc hist_vxvz_resc hist_vyvz_resc_05] = fct_illstr_Projections(vpagrid_05,vpegrid_05,vxf_resc,vyf_resc,vzf_resc,PlotOpt);
% % % % %     [hist_vxvy_resc hist_vxvz_resc hist_vyvz_resc_1] = fct_illstr_Projections(vpagrid_1,vpegrid_1,vxf_resc,vyf_resc,vzf_resc,PlotOpt);

    hist_vyvz_resc_loop_01 = hist_vyvz_resc_loop_01 + hist_vyvz_resc_01;
    hist_vyvz_resc_loop_025 = hist_vyvz_resc_loop_025 + hist_vyvz_resc_025;
% % % % %     hist_vyvz_resc_loop_05 = hist_vyvz_resc_loop_05 + hist_vyvz_resc_05;
% % % % %     hist_vyvz_resc_loop_1 = hist_vyvz_resc_loop_1 + hist_vyvz_resc_1;
% % % % %     
if mod(k,1)==0 || k==kmax  
display(strcat('k loop ',num2str(k),' done!'));

   % %     savename2    = [num2str(k),'_all',savestr,'.mat'];
    save(savename2,...
         'vpagrid_01','vpegrid_01','vpagrid_025','vpegrid_025','hist_vyvz_resc_loop_01','hist_vyvz_resc_loop_025');
    
end
toc
end
foldername = ['new_s',num2str(s),'_n',num2str(n),'N',num2str(N),'wvl',num2str(wvlm./1e-9),'I',num2str(I0,'%10.1e\n'),...
             '_dt',num2str(dt),'kmax',num2str(kmax),'nsample',num2str(nsample,'%10.0e\n')];
save('params.mat');       
mkdir(foldername);
movefile('*.mat',foldername)  
end

function [vxf_all vyf_all vzf_all vxf_dir vyf_dir vzf_dir vxf_resc vyf_resc vzf_resc T] = singleICEP_v0(IWcm,CEP,pulse_,atom,nsample,savename,dt,s) 


%% parameters for the laser pulse 
c0       = 299792458./(2.1876912633e6);          % au, speed of light 
n        = pulse_(1);                            % propto pulse duration
N        = pulse_(2);
wvlm     = pulse_(3);                            % 
wvl      = wvlm/(5.2917720859e-11);              % [au]
omega    = (2*pi*c0)/wvl;                        % [au], angular frequency
T        = 2.*pi./omega;

% IWcm     = 1e14;               %
% CEP      = pi/2;                               % au
% dt       = 1;                                  % au
Up       = (0.09337.*IWcm.*(wvlm.^2))./27.211;   % au

%% define the field
[Env E I A ALPHA BETA v r tgrid E_fh] = fct_get_EnvEIAAlphBetvr_sin2_nN(dt,wvlm,IWcm,CEP,n,N);
% ind = find(I>0.5.*max(I));
% tp  = tgrid(max(ind))-tgrid(min(ind));
%disp(['n= ',num2str(n),'i.e. tp = ',num2str(tp.*24.2e-18/1e-15),' fs at a wvl of ',num2str(wvlm./1e-6),' mum'])
% plot(tgrid,Env,tgrid,E,'LineWidth',4,'Color','r')
% set(gcf,'color','w'); % set figure background to white
%% function handle for bsi
Fc = @(kappa_,Z_,m_) (kappa_.^4)./(8.*(2.*Z_-kappa_.*(m_+1)));

%%  All atoms, pure rate, hydrogen rate, gen atom, wo_BI, wo_me
SaturationSwitch = 0;
PlotOpt          = 0;
tbiSwitch        = 0;
AtomString       = atom;
[Ip kappa Z Cnl l beta alphaN alphaI] = fct_get_Atom(AtomString);
gamma            = sqrt(abs(Ip)/(2*Up));

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
IonAmp = IonAmp_GenA_H1s_wo_tbi_wo_me;
IonAmp=IonAmp/max(IonAmp);

%% settings for the monte carlo
% nsample = 1e3;

%% generate starting times
PlotOpt = 0;
[texit_ind texit] = fct_gen_texit(nsample,tgrid,IonAmp,PlotOpt);
hist_texit =  hist(texit,tgrid);
%texit_ind = 1:1:(length(tgrid)-10);

%% all trajectories - set tunneling conditions for all trajectories
Eexit_sample  = E(texit_ind);                                    % field at tunneling
Ip_sample     = Ip(texit_ind);                                   % Ip at tunneling
vpexit        = zeros(1,length(texit));                          % perp velocity
[chi_sample]  = fct_get_angexit(vpexit,PlotOpt);                 % emission angle
vpaxit        = zeros(1,length(texit));                          % para velocity
v0x_sample    = vpexit.*cos(chi_sample);
v0y_sample    = vpexit.*sin(chi_sample);
v0z_sample    = vpaxit;
rpaxit        = zeros(1,length(texit));                          % parallel tunnel exit
rpexit        = zeros(1,length(texit));                          % perp tunnel exit
r0x_sample    = rpexit.*cos(chi_sample);
r0y_sample    = rpexit.*sin(chi_sample);
r0z_sample    = rpaxit;

% all trajectories - condense initial conditions in matrix
[texit v0mat r0mat] = fct_get_initMat(texit,v0x_sample,v0y_sample,v0z_sample,r0x_sample,r0y_sample,r0z_sample);

% all trajectories -  illustrate tunneling 
PlotOpt= 0;
v0grid = -0.5:0.005:0.5;
r0grid = -15:0.005:15;
fct_illstr_Tunnel(v0grid,r0grid,v0mat,r0mat,PlotOpt)

% dice an actual return according to the scattering prob and so on
PlotOpt = 0;
%s       = 3/2;
sigwp   = @(ts,tr)   (tr-ts).^(-s); %sigwp   = @(ts,tr)   (tr-ts).^(-s);
% sigtg   = @(q,theta) abs(q);
% sigtg   = @(q,theta) abs(q); %volumeelement for filled circle

%% cross-section for noble gas potential?
scXe = [54.00 3.604];
scAr = [18.00 2.722];
scNe = [10.00 2.038];


if strcmp(atom,'Xe')
    sc   = scXe;
    rc   = 2.4566; % covalent radius, https://de.wikipedia.org/wiki/Kovalenter_Radius 130e-12 m converted to atomic units;
end
if strcmp(atom,'Ne')
    sc   = scNe;
    rc   = 1.304;  % 69e-12 m
end
if strcmp(atom,'Ar')
    sc   = scAr;
    rc   = 1.833;   % 97e-12 m
end

Z=2;
mu      = 0.06;
mu_prime = mu/2;


if strcmp(atom,'He1s2s-singlet')==1
     sigtg_para_mu   = @(q,theta)   4*((Z-1)./(4*q.^2.*sin(theta/2).^2+mu^2) + (2*Z^2+q.^2.*sin(theta/2).^2)./4./(Z^2+q.^2.*sin(theta/2).^2).^2 - Z.^4./((q.^2+mu_prime^2).*(Z^2+q.^2.*sin(theta/2).^2).^2) ).^2;

    sigtg =@(q,theta) sigtg_para_mu(q,theta);                       
end 

if strcmp(atom,'He1s2s-triplet')==1
       sigtg_ortho_mu =  @(q,theta)   4*((Z-1)./(4*q.^2.*sin(theta/2).^2+mu^2) + (2*Z^2+q.^2.*sin(theta/2).^2)./4./(Z^2+q.^2.*sin(theta/2).^2).^2 + Z.^4./((q.^2+mu_prime^2).*(Z^2+q.^2.*sin(theta/2).^2).^2) ).^2;

      sigtg =@(q,theta) sigtg_ortho_mu(q,theta); 
end

   rc   = 0.5871;   % 31e-12 m
    rcov = rc;

%% cross-section without angle dependence 
%sigtg=@(q,theta) (abs(q));      % abs(q) is equal to homogenously filled cylinder 


%% %% find return list, vx0=vy0=vz0=x0=y0=z0=0
PlotOpt = 0;
NrMax   = 5;
texit_ind = 1:1:length(tgrid);%2.65e4:1:2.93e4;%1:1:length(tgrid);
 

load_table = 0;
if exist('dt0.01_1ion.mat', 'file') == 2
load_table = 1;
end

if load_table==1
    table_name ='dt0.01_1ion.mat';
    load(table_name);
else
[trind_mat tr_mat tt_mat rescind dirind] = fct_get_return1Traj1D_v2(tgrid,texit_ind,zeros(1,length(tgrid)),zeros(1,length(tgrid)),r,rcov,NrMax,PlotOpt);
toc
save('dt0.01_1ion.mat','trind_mat','tr_mat','tt_mat','rescind','dirind');
end
disp('found return list')
%[trind_mat tr_mat tt_mat rescind dirind] = fct_get_return1Traj1D_v1(tgrid,texit_ind,zeros(1,length(tgrid)),zeros(1,length(tgrid)),r,v,NrMax,PlotOpt);
%% scattering & nonscattering
Pne_all        = IonAmp./sum(IonAmp);
te_ind_all     = fct_gen_distr(Pne_all,1,nsample);
te_sample_all  = tgrid(te_ind_all);                         % list of all starting times
Ne_all         = hist(te_sample_all,tgrid);

nsample_resc = sum(Ne_all(texit_ind(rescind)));                        % number of resc electrons
nsample_dir  = sum(Ne_all(texit_ind(dirind)));                         % number of direct electrons

Pne_dir    = zeros(1,length(tgrid));
Pne_dir(texit_ind(dirind))    = Pne_all(texit_ind(dirind));
te_ind_dir = fct_gen_distr(Pne_dir./sum(Pne_dir),1,nsample_dir);
te_sample_dir = tgrid(te_ind_dir);        

%% distribute rescattering
PlotOpt =0;
[vx_resc_list vy_resc_list vz_resc_list te_resc_list tr_resc_list tt_resc_list vr_resc_list th_resc_list trind_resc_list] = fct_gen_resc2D_fromList_v1(A,IonAmp,tgrid,rescind,texit_ind(rescind),trind_mat,nsample_resc,sigwp,sigtg,PlotOpt);

% calculate final velocity of rescattered
PlotOpt  = 0;                    
[vxf_resc vyf_resc vzf_resc] = fct_get_vfinal(A,trind_resc_list',vx_resc_list,vy_resc_list,vz_resc_list,Up,PlotOpt); 

%
texit_resc = te_resc_list;
tr_resc    = tr_resc_list;
angle_resc = th_resc_list;
vr_resc    = vr_resc_list;
vxsc       = vx_resc_list;
vysc       = vy_resc_list;
vzsc       = vz_resc_list;
vxf_resc   = vxf_resc';  
vyf_resc   = vyf_resc';
vzf_resc   = vzf_resc';

% vpegrid = -3:0.025:3;
% vpagrid = -3:0.025:3;
% [hist_vxvy_all hist_vxvz_all hist_vyvz_all]    = fct_illstr_Projections(vpagrid,vpegrid,vxf_resc',vyf_resc',vzf_resc',1);

%% direct trajectories
% here comes an option, either you take all dirind only or you take
% boths... dirind & rescind for this

% % % % texit_ind_dir= te_ind_dir';        % index on tgrid of the starting time
% % % % texit_dir    = te_sample_dir;      % starting time of direct 
% % % % 
% % % % % dice perp velocities
% % % % MeanOpt = 1;
% % % % PlotOpt = 0;
% % % % [vpexit_ind_dir vpexit_dir vpgrid] = fct_get_vpexit(Eexit_sample(texit_ind_dir),Ip_sample(texit_ind_dir),MeanOpt,PlotOpt);
% % % %  
% % % % % set/ calculate para velocity 
% % % % PlotOpt = 0; 
% % % % [vpaxit_dir] = fct_get_vpaexit(vpexit_dir,PlotOpt);
% % % % 
% % % % % % dice the emission angle 
% % % % PlotOpt = 0;
% % % % [chi_dir] = fct_get_angexit(vpexit_dir,PlotOpt);
% % % % 
% % % % % set the initial velocities
% % % % vsx_dir = (vpexit_dir.*cos(chi_dir))';                 
% % % % vsy_dir = (vpexit_dir.*sin(chi_dir))';                 
% % % % vsz_dir = (vpaxit_dir)';
% % % % 
% % % % % calculate the final velocity of direct electrons
% % % % PlotOpt = 0;
% % % % [vxf_dir vyf_dir vzf_dir] = fct_get_vfinal(A,texit_ind_dir,vsx_dir',vsy_dir',vsz_dir',Up,PlotOpt); 

%% results direct
% % % % texit_dir;  % starting times of direct
% % % % vxf_dir; 
% % % % vyf_dir; 
% % % % vzf_dir;

%% all thogether
vxf_dir = 0;%vxf_dir';
vyf_dir = 0;% vyf_dir';
vzf_dir = 0;%vzf_dir';

vxf_all = 0;%[vxf_dir;vxf_resc];
vyf_all = 0;%[vyf_dir;vyf_resc];
vzf_all = 0;%[vzf_dir;vzf_resc];

%% save everything for potential use later... 
% % % % save(savename,...
% % % %                 'tgrid','E','I','wvlm',...                                                     % the field somehow
% % % %                 'texit_dir','vxf_dir','vyf_dir','vzf_dir',...                           % directs
% % % %                 'texit_resc','tr_resc','angle_resc','vr_resc','vxsc','vysc','vzsc',...  % rescattering
% % % %                 'vxf_resc','vyf_resc','vzf_resc');                                      % final resc
% % % % 


end
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
SaturationSwitch = 0;
PlotOpt          = 0;
tbiSwitch        = 0;
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

