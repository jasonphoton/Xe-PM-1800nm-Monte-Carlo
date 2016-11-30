function fct_illstr_Tunnel(v0grid,r0grid,v0mat,r0mat,PlotOpt)

v0x_sample = v0mat(:,1);
v0y_sample = v0mat(:,2);
v0z_sample = v0mat(:,3);
r0x_sample = r0mat(:,1);
r0y_sample = r0mat(:,2);
r0z_sample = r0mat(:,3);


hist_v0xv0y = hist3([v0x_sample v0y_sample], {v0grid' v0grid'});
hist_v0xv0z = hist3([v0x_sample v0z_sample], {v0grid' v0grid'});
hist_v0yv0z = hist3([v0y_sample v0z_sample], {v0grid' v0grid'});

r0grid = -15:0.2:15;
hist_r0xr0y = hist3([r0x_sample r0y_sample], {r0grid' r0grid'});
hist_r0xr0z = hist3([r0x_sample r0z_sample], {r0grid' r0grid'});
hist_r0yr0z = hist3([r0y_sample r0z_sample], {r0grid' r0grid'});

if PlotOpt ==1
figure;
subplot(2,3,1)
imagesc(v0grid,v0grid,hist_v0xv0y)
xlabel('v0x'); ylabel('v0y');

subplot(2,3,2)
imagesc(v0grid,v0grid,hist_v0xv0z)
xlabel('v0x'); ylabel('v0z');

subplot(2,3,3)
imagesc(v0grid,v0grid,hist_v0yv0z)
xlabel('v0y'); ylabel('v0z');

subplot(2,3,4)
imagesc(r0grid,r0grid,hist_r0xr0y)
xlabel('r0x'); ylabel('r0y');

subplot(2,3,5)
imagesc(r0grid,r0grid,hist_r0xr0z)
xlabel('r0x'); ylabel('r0z');

subplot(2,3,6)
imagesc(r0grid,r0grid,hist_r0yr0z)
xlabel('r0y'); ylabel('r0z');
end
end