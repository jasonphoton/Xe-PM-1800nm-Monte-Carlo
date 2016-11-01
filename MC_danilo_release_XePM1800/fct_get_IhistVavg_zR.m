
function [Iaxis,hist_Isample] = fct_get_IhistVavg_zR(w0,I0,Zr,Itheogrid,PlotOpt)

if I0<Itheogrid(1)
    'think about what you do, the max I smaller than the min Itheogrid'
end;

%% define int. function
F_w0I0Zr_Izr = @(w0,Zr,I0,Z,R) I0./(1+(Z/Zr).^2) .* exp(-2*R.^2./(w0^2*(1+(Z/Zr).^2)));


%% plot the intensity function for reference on a mm scale

% define int. function in mm 
% w0 = 0.5.*(22.44e-6+25.16e-6)./1e-3;        % [mm],    e.g. H2+ 1693
% I0 = 4e14;                                  % [W/cm2], e.g. H2+ 1693
% Zr = 4.21;                                  % [mm],    e.g. H2+ 1693
F_Izr = @(Z,R) F_w0I0Zr_Izr(w0,Zr,I0,Z,R);

%% do a monte carlo
N    = 1000000;
Zmin = -2*Zr;
Zmax = 2*Zr;
Rmax = 4.*w0;

% do the monte carlo

Xtemp = Rmax*rand([1,2*N]);
Ytemp = Rmax*rand([1,2*N]); 
Rtemp = sqrt(Xtemp.^2+Ytemp.^2);
indx    = find(Rtemp<Rmax+eps);
Xsample = Xtemp(indx(1:N));
Ysample = Ytemp(indx(1:N));
Rsample = Rtemp(indx(1:N));
Zsample = Zmin+(Zmax-Zmin)*rand([1,N]);
Isample = F_Izr(Zsample,Rsample);
Istep   = Itheogrid(2)-Itheogrid(1);
Iaxis   = (Itheogrid(1)-Istep):Istep:(Itheogrid(end)+Istep);          % linspace(1e13,I0,100);
hist_Isample = hist(Isample,Iaxis);
hist_Isample = hist_Isample./max(hist_Isample);
hist_Isample = hist_Isample(2:1:(end-1));
Iaxis        = Iaxis(2:1:(end-1));
if PlotOpt ==1

% plot the intensity in this region
figure;
subplot(2,3,1);
Z=Zmin:.01:Zmax;
R=0*Z;
plot(Z,F_Izr(Z,R));
xlabel('z [mm]');
ylabel('Intensity');
title('I(z) W/cm^2');
grid on
axis tight
clear 'Z' 'R';

subplot(2,3,2);
R=-Rmax:1e-3:Rmax;
Z=0*R;
plot(R,(F_Izr(Z,R)));
xlabel('\rho [mm]');
ylabel('I(\rho)  W/cm^2');
grid on
title('I(\rho)');
clear 'Z' 'R';

subplot(2,3,4);
z=Zmin:.01:Zmax;
r=-Rmax:1e-3:Rmax;
[Z,R] = meshgrid(z,r);
imagesc(z,r,F_Izr(Z,R));
xlabel('z [mm]');
ylabel('\rho [mm]');
colorbar('location','southoutside');
title('I [W/cm^2]');

subplot(2,3,5);
imagesc(z,r,log10(F_Izr(Z,R)));
xlabel('z [mm]');
ylabel('\rho [mm]');
colorbar('location','southoutside');
title('log_{10}(I(z,\rho))')
clear 'z' 'r' 'Z' 'R';

subplot(2,3,3)
plot(Iaxis,hist_Isample,'.');
xlabel('I bin [W/cm2]');
ylabel('Counts');
grid on

subplot(2,3,6)
semilogy(Iaxis,hist_Isample,'.');
xlabel('I bin [W/cm2]');
ylabel('log(Counts)');
grid on

end

% %% plot the intensity function for reference
% 
% subplot(2,2,1);
% Z=-10:.01:10;
% R=0*Z;
% plot(Z,F_w0I0Zr_Izr(1,1,1,Z,R));
% xlabel('z');
% ylabel('Intensity');
% title('I(z)')
% clear 'Z' 'R';
% 
% subplot(2,2,2);
% R=-10:.01:10;
% Z=0*R;
% plot(R,F_w0I0Zr_Izr(1,1,1,Z,R));
% xlabel('\rho');
% ylabel('Intensity');
% title('I(\rho)')
% clear 'Z' 'R';
% 
% subplot(2,2,3);
% z= -10:.01:10;
% r=-10:.01:10;
% [Z,R] = meshgrid(z,r);
% imagesc(z,r,F_w0I0Zr_Izr(1,1,1,Z,R));
% xlabel('z');
% ylabel('\rho');
% title('I(z,\rho)')
% colorbar('location','southoutside');
% 
% subplot(2,2,4);
% imagesc(z,r,log10(F_w0I0Zr_Izr(1,1,1,Z,R)));
% xlabel('z');
% ylabel('\rho');
% colorbar('location','southoutside');
% title('log_{10}(I(z,\rho))')
% clear 'z' 'r' 'Z' 'R';

end

