function [vx_resc_list vy_resc_list vz_resc_list te_resc_list tr_resc_list tt_resc_list vr_resc_list th_resc_list trind_resc_list] = fct_gen_resc2D_fromList_v1(A,IonAmp,tgrid,original_rescind,rescind,trind_mat,nsample_resc,sigwp,sigtg,PlotOpt)
% sigwp cross-section wavepacket dispersion
% sigtg cross-section wavepacket dispersion

%% distribute rescattering
tsind_resc     = rescind;
trind_mat_resc = trind_mat(original_rescind,:);           % index of return

% % % angle_limit = 0.1;
% % % 
% % % anglegrid1  =-pi:1e-3:-angle_limit;
% % % anglegrid2  =angle_limit:1e-3:pi;
% % % anglegrid  = [anglegrid1,anglegrid2];

anglegrid  =-pi:1e-3:pi;
killerangle       = 3;
inglgeidkillerind = find(abs(anglegrid)<deg2rad(killerangle));
killer_grid = anglegrid(inglgeidkillerind);

NrMax      = length(trind_mat(1,:));
te_ind_mat = zeros(length(tsind_resc),NrMax);
tr_ind_mat = zeros(length(tsind_resc),NrMax);
te_mat     = zeros(length(tsind_resc),NrMax);
tr_mat     = zeros(length(tsind_resc),NrMax);
tt_mat     = zeros(length(tsind_resc),NrMax);
vr_mat     = zeros(length(tsind_resc),NrMax);
Wte_mat    = zeros(length(tsind_resc),NrMax);   
Wtt_mat    = zeros(length(tsind_resc),NrMax);   
Wth_mat    = zeros(length(tsind_resc),NrMax);

% make matrices
if PlotOpt==1
figure;
end
for m=1:NrMax
    trind_mat_resc_m = trind_mat_resc(:,m);
    ind              = find(trind_mat_resc_m>0);
    
    % indizes
    te_indm = tsind_resc(ind);
    tr_indm = trind_mat_resc(ind,m);
    
    te_ind_mat(ind,m) = te_indm;
    tr_ind_mat(ind,m) = tr_indm;
    
    % times
    tem = tgrid(tsind_resc(ind));
    trm = tgrid(trind_mat_resc(ind,m));
    ttm = tgrid(trind_mat_resc(ind,m))-tgrid(tsind_resc(ind));
    
    te_mat(ind,m) = tem;
    tr_mat(ind,m) = trm;
    tt_mat(ind,m) = ttm;
    
    
    % velocity
    vrm = A(tr_indm)-A(te_indm);
    vr_mat(ind,m) = vrm;

    % Wte
    Wtem = IonAmp(te_indm);
    Wte_mat(ind,m) = Wtem;
    
    % Wtt
    Wttm  = sigwp(tem,trm);
    Wtt_mat(ind,m) = Wttm;
    
    % Wfs sc
    Wthm = zeros(1,length(vrm));
    for l=1:length(vrm)
        Wthm(l) = sum(sigtg(vrm(l).*ones(1,length(anglegrid)),anglegrid));
    end
    Wth_mat(ind,m) = Wthm;
    

    if PlotOpt==1
    
    subplot(3,3,1)
    plot(tem,trm,'.'); hold on
    xlabel('texit')
    ylabel('treturn')
    
    subplot(3,3,4)
    plot(tem,ttm,'.'); hold on
    xlabel('texit')
    ylabel('traveltime')
    
    subplot(3,3,7)
    plot(tem,vrm,'.'); hold on
    xlabel('texit')
    ylabel('vr')
    
    subplot(3,3,2)
    plot(tem,Wtem,'.'); hold on
    xlabel('texit')
    ylabel('Wtexit')
    
    subplot(3,3,5)
    semilogy(tem,Wttm,'.'); hold on
    xlabel('texit')
    ylabel('Wtt')
    
    subplot(3,3,8)
    plot(tem,Wthm,'.'); hold on
    xlabel('texit')
    ylabel('Wtheta')
    
    subplot(3,3,3)
    semilogy(tem,Wtem,'.'); hold on
    xlabel('texit')
    ylabel('Wtt')
    
    subplot(3,3,6)
    semilogy(tem,Wtem.*Wttm,'.'); hold on
    xlabel('texit')
    ylabel('Wte.*Wtt')

    subplot(3,3,9)
    semilogy(tem,Wtem.*Wttm.*Wthm,'.'); hold on
    xlabel('texit')
    ylabel('Wte.*Wtt.*Wtheta')
    
    end

end
%% introduce a condition to kill low energy returns
indkiller = find(abs(vr_mat)<0.15);  
Wte_mat(indkiller) = 0;   
Wtt_mat(indkiller) = 0;   
Wth_mat(indkiller) = 0;


%% distribute particles 
Pne_allreturns = Wte_mat(:,1)./sum(Wte_mat(:,1));
te_ind_sample  = fct_gen_distr(Pne_allreturns',1,nsample_resc);
te_sample_resc = te_mat(te_ind_sample,1);                   % list ot rescattering starting times

% illustration on histogram on tgrid
hist_te_sample_resc = hist(te_sample_resc,tgrid);

if PlotOpt ==1
plot(tgrid,hist_te_sample_resc./max(hist_te_sample_resc)); hold on
plot(te_mat(:,1),Wte_mat(:,1)./max(Wte_mat(:,1)),'.');
end
N_te_sample_resc    = hist(te_sample_resc,te_mat(:,1));
if PlotOpt ==1
plot(te_mat(:,1),Wte_mat(:,1)./max(Wte_mat(:,1)),'-'); hold on
plot(te_mat(:,1),N_te_sample_resc./max(N_te_sample_resc),'.'); 
end;

N_te_tr = zeros(length(tsind_resc),NrMax);
nr_ind  = 1:1:NrMax;
for k=1:length(tsind_resc)
    nek   = N_te_sample_resc(k);
    
    if nek>0
        wttk     = Wtt_mat(k,:);
        wthk     = Wth_mat(k,:);
        vrk      = vr_mat(k,:);
        Pttth    = (wttk.*wthk)/sum(wttk.*wthk);
        ttth_ind = fct_gen_distr(Pttth,1,nek);
        N_te_tr(k,:) = hist(nr_ind(ttth_ind),nr_ind);
    else
        N_te_tr(k,:) = zeros(1,NrMax);
    end
end

if PlotOpt==1
figure;
subplot(2,1,1)
pcolor(te_mat(:,1),nr_ind,N_te_tr'); shading flat
xlabel('te (au)')
ylabel('nr (#)')
colorbar
title('N(te,nr)')

subplot(2,1,2)
plot(te_mat(:,1),N_te_tr,'.'); shading flat
xlabel('te (au)')
ylabel('nr (#)')
title('N(te,nr)')
legend(num2str(nr_ind))

end

te_resc_list = zeros(1,nsample_resc);
trind_resc_list = zeros(1,nsample_resc);
tr_resc_list = zeros(1,nsample_resc);
tt_resc_list = zeros(1,nsample_resc);
vr_resc_list = zeros(1,nsample_resc);
th_resc_list = zeros(1,nsample_resc);

n_=0;
for l=1:length(tsind_resc)
    for g=1:NrMax
        if N_te_tr(l,g)>0
            n_te_tr = N_te_tr(l,g);
            te_lg   = te_mat(l,g);
            tr_lg   = tr_mat(l,g);
            tt_lg   = tt_mat(l,g);
            trind_lg= tr_ind_mat(l,g);
            vr_lg   = vr_mat(l,g);
            Pth         = sigtg(vr_lg.*ones(1,length(anglegrid)),anglegrid);
            Pth         = Pth./sum(Pth);       
            Pth(inglgeidkillerind) = 0;             % kill forward scttering to avoid divergence

            angle_ind   = fct_gen_distr(Pth,1,n_te_tr);
            
            p=length(angle_ind);
            trind_resc_list((n_+1):1:(n_+p)) =  trind_lg.*ones(1,p);
            te_resc_list((n_+1):1:(n_+p)) = te_lg.*ones(1,p);
            tr_resc_list((n_+1):1:(n_+p)) = tr_lg.*ones(1,p);
            tt_resc_list((n_+1):1:(n_+p)) = tt_lg.*ones(1,p);
            vr_resc_list((n_+1):1:(n_+p)) = vr_lg.*ones(1,p);
            th_resc_list((n_+1):1:(n_+p)) = anglegrid(angle_ind);
            n_      = n_+n_te_tr;   
            
        end
    end
end
vx_resc_list = zeros(1,nsample_resc);
vy_resc_list = sin(th_resc_list).*(vr_resc_list);
vz_resc_list = cos(th_resc_list).*(vr_resc_list);