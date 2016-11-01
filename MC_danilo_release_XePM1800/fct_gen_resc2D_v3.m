function [tresc_ind vxsc vysc vzsc angle_resc v_resc] = fct_gen_resc2D_v3(A,tgrid,texit_ind,trind_mat,sigwp,sigtg,PlotOpt)
% sigwp cross-section wavepacket dispersion
% sigtg cross-section wavepacket dispersion

Pwpdis  = zeros(length(trind_mat(:,1)),length(trind_mat(1,:)));
vreturn = zeros(length(trind_mat(:,1)),length(trind_mat(1,:)));
% angleret= zeros(length(trind_mat(:,1)),length(trind_mat(1,:)));

% get return velocity & calc prob due to wave packet dispersion 
for m=1:length(trind_mat(1,:))
    trind_mat_m    = trind_mat(:,m);
    indm           = find(trind_mat_m>0);
    texit_ind_m    = texit_ind(indm);
    Pwpdis(indm,m) = sigwp(tgrid(texit_ind(indm)),tgrid(trind_mat_m(indm)));
    vreturn(indm,m)= A(trind_mat_m(indm))-A(texit_ind(indm));                   % A(tr)-A(ts);
end

anglegrid  =-pi:1e-3:pi;
%Pscatt     = zeros(length(anglegrid),length(trind_mat(1,:)));
%Pnorm_     = Pscatt(:);
[mG,aG]    = meshgrid(1:1:length(trind_mat(1,:)),anglegrid);
angle_resc = zeros(length(trind_mat(:,1)),1);
v_resc     = zeros(length(trind_mat(:,1)),1);
tresc_ind  = zeros(length(trind_mat(:,1)),1);
for l=1:length(trind_mat(:,1))
    indl        = find(abs(vreturn(l,:))>0);
    Pwpdis_l    = Pwpdis(l,indl);
    vreturn_l   = vreturn(l,indl);
    trind_mat_l = trind_mat(l,indl);
    Pscatt      = zeros(length(anglegrid),length(vreturn_l));
    sPscatt_    = zeros(1,length(vreturn_l));
    ssPscatt    = zeros(1,length(vreturn_l));
    for m=1:length(vreturn_l)
        sPscatt_(m)= sum(sigtg(vreturn_l(m),anglegrid)); 
        Pscatt(:,m)= Pwpdis_l(m).*(sigtg(vreturn_l(m),anglegrid)./sPscatt_(m));
        ssPscatt(m)= sum(Pscatt(:,m));
    end
    Pnorm         = Pscatt(:)./sum(sum(Pscatt(:)));
    asc_ind       = fct_gen_distr(Pnorm',1,1);
    angle_resc(l) = aG(asc_ind);
    v_resc(l)     = vreturn_l(mG(asc_ind));
    tresc_ind(l)  = trind_mat_l(mG(asc_ind)); 
end

vxsc = zeros(1,length(trind_mat(:,1)))';
vysc = sin(angle_resc).*(v_resc);
vzsc = cos(angle_resc).*(v_resc);

if PlotOpt==1
    vpegrid = -3:0.05:3;
    vpagrid = -3:0.05:3;
    [hist_v0xv0y hist_v0xv0z hist_v0yv0z] = fct_illstr_VMI(vpagrid,vpegrid,vxsc,vysc,vzsc);
end