function [hist_v0xv0y hist_v0xv0z hist_v0yv0z] = fct_illstr_Projections(vpagrid,vpegrid,vxf,vyf,vzf,PlotOpt)

hist_v0xv0y = hist3([vxf vyf], {vpegrid' vpegrid'});
hist_v0xv0z = hist3([vxf vzf], {vpegrid' vpagrid'});
hist_v0yv0z = hist3([vyf vzf], {vpegrid' vpagrid'});

if PlotOpt==1
figure;
subplot(2,3,1)
pcolor(vpegrid,vpegrid,hist_v0xv0y'); shading flat
xlabel('vx'); ylabel('vy');
axis equal
colorbar

subplot(2,3,2)
pcolor(vpegrid,vpagrid,hist_v0xv0z'); shading flat
xlabel('vx'); ylabel('vz');
axis equal
colorbar

subplot(2,3,3)
pcolor(vpegrid,vpagrid,hist_v0yv0z'); shading flat
xlabel('vy'); ylabel('vz');
axis equal
colorbar

subplot(2,3,4)
pcolor(vpegrid,vpegrid,log10(hist_v0xv0y')); shading flat
xlabel('vx'); ylabel('vy');
axis equal
colorbar

subplot(2,3,5)
pcolor(vpegrid,vpagrid,log10(hist_v0xv0z')); shading flat
xlabel('vx'); ylabel('vz');
axis equal
colorbar

subplot(2,3,6)
pcolor(vpegrid,vpagrid,log10(hist_v0yv0z')); shading flat
xlabel('vy'); ylabel('vz');
axis equal
colorbar
end
% subplot(2,3,4)
% imagesc(r0grid,r0grid,hist_r0xr0y)
% xlabel('r0x'); ylabel('r0y');
% 
% subplot(2,3,5)
% imagesc(r0grid,r0grid,hist_r0xr0z)
% xlabel('r0x'); ylabel('r0z');
% 
% subplot(2,3,6)
% imagesc(r0grid,r0grid,hist_r0yr0z)
% xlabel('r0y'); ylabel('r0z');

end