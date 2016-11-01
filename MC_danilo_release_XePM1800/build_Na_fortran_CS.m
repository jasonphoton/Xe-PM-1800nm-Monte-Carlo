clear all
load('na_5_5ev_cs')

ha = @(xq) interp1(Winkel,cs,xq,'pchip');

semilogy(Winkel,cs,'.r')
hold on
ang = 0:0.1:180;
semilogy(ang,ha(ang),'.b')
hold off

%% programm

%% programm

clear all
nr_energies = 80;
%load the cross sections
for i=1:nr_energies    
filename_used = strcat('elastic_',num2str(i),'.out');
delimiterIn = ' ';
headerlinesIn = 0;
A = importdata(filename_used,delimiterIn,headerlinesIn);

if i==1
Winkel = A(:,1);
% set boundaries right
Winkel(1) = 0;
Winkel(end) = 180;
end
cs(i,:) = A(:,2);
sprintf(strcat('cs ',num2str(i),' is done!'))
end
% load the energies used
for i=1:nr_energies    
    filename_used = strcat('nrg_used_',num2str(i),'.out');
delimiterIn = ' ';
headerlinesIn = 0;
A = importdata(filename_used,delimiterIn,headerlinesIn);
energy(i) = A; %in eV
sprintf(strcat('energy',num2str(i),' is done!'))
end

%
save('na_cs.mat','energy','Winkel','cs')