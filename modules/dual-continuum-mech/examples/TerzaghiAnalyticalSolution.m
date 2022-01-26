% Terzaghi: simple analytical solution (1D problem)
% author: Ruslan Rin (iskhakov@stanford.edu)
% 5/13/2016
close all;clc;clear all;
bar = 1e5;        %bar to Pa
md = 1e-15;       %md to m2
cp = 1e-3;        %cp to Pa*s
b = 1;            %Biot's coefficient
poro = 0.375;     %porosity
nu = 0.25;        %Poisson ratio
E = 1e4*bar;      %Young's modulus, Pa
cf = 4.4e-5/bar;  %compressibility, 1/Pa
cbulk = 0;        %grain compressibility, 1/Pa
mu = 0.0981*cp;   %viscosity, Pa*s
alpha = 1;
k = 1*md;         %permeability, m2
pL = 100*bar;     %applied pressure, Pa
L = 100;          %length of the domain
Nb = 10;         %number of blocks
h = L/Nb;         %discretized, m
K = E/3/(1-2*nu); %Bulk modulus
G = E/2/(1+nu);   %shear modulus
M = 1/(poro*cf + (1-poro)*cbulk); %Biot's modulues
K1D = K+4/3*G;    %uniaxial drained bulk modulus
c = k/mu*M * K1D/(K1D+alpha^2*M); %consolidation coefficient
tc = L^2/c;       %characteristic time
Nt = 100;         %number of timesteps
N = 500;          %number of elements in series
m = [0:1:N];      %fourier series sequence
p0 = alpha*M*pL/(K1D+alpha^2*M);  % undrained pressure, Pa
u0 = L*pL/(K1D+alpha^2*M);        %max initial dispacement, m
t = [1:1:Nt]*tc/Nt;
z = [h/2:h:L];    % center of our blocks
z2 = [0:h:L];
if (Nb < 20)
  %we need to shift since numerical schemes uses 9 blocks to mimic
  %boundaries capabilities
  z(Nb) = z(Nb-1);    
  z2(1) = -h/8; 
end;
%analytical solution
for i = 1:Nt  
  tmp_exp = exp(-pi^2*c*t(i)*(2*m+1).^2/(4*L^2));
  for iz = 1:Nb
    xtmp = sum(1./(2*m+1).*tmp_exp.*sin((2*m+1)*pi.*z(iz)/2/L));
    p(iz,i) = 4/pi*p0*xtmp; %Pa
    xtmp = sum(1./(2*m+1).^2.*tmp_exp.*cos((2*m+1)*pi.*z2(iz)/2/L));
    u(iz,i) = alpha*p0/K1D*(L-z2(iz)-8*L/pi^2*xtmp) + u0;
  end;
end;
p = p/bar;
tt = t/tc;

% pressure and dispacement 
pmax = 85;
umax = 0.8;

figure; set(gcf,'color','white'); set(gcf, 'Position', [400 200 600 500])
[hax,h1,h2] = plotyy(tt,p(Nb,:), tt,u(1,:));
set(h1,'color','red')
set(h2,'color','blue')
ylim(hax(1),[0 pmax]);
ylabel(hax(2), 'dispacement, m');
set(hax(2),'YTick',0:0.2:0.8); 
ylim(hax(2),[0 umax]);
xlabel('time, t/t_c');
ylabel(hax(1), 'pressure, bar');
set(hax(1),'YTick',0:10:85); 
hold on;

[hax2,h3,h4] = plotyy(tt,p_vtk(Nb,:),tt,u_vtk(1,:));
set(hax2(2),'ytick',[])
set(hax2(2),'yticklabel',[])
set(h3,'marker', '.', 'color','red')
set(h4,'marker', '.', 'color','blue')
ylim(hax2(1),[0 pmax]);
ylim(hax2(2),[0 umax]);

[hax3,h5,h6] = plotyy(tt,p_vtk_seq(Nb,:),tt,u_vtk_seq(1,:));
set(hax3(2),'ytick',[])
set(hax3(2),'yticklabel',[])
set(h5,'marker', 'd', 'color','red')
set(h6,'marker', 'd', 'color','blue')
ylim(hax3(1),[0 pmax]);
ylim(hax3(2),[0 umax]);

v = [h1 h3 h2 h4 h5 h6]';
legend(v, 'p_{bot,analytical}', 'p_{bot,numerical}','p_{bot,sequential}',...
          'u_{top,analytical}','u_{top,numerical}','u_{top,sequential}',...
          'Location', 'SouthWest');
Image = getframe(gcf); imwrite(Image.cdata, 'terzagi10_comparison.png');

% show numerical solution
%pressure images
i1 = 1; i2 = 10; i3 = 100;
figure; set(gcf,'color','white'); set(gcf, 'Position', [400 200 400 400]);
ax(1) = subplot(1,3,1); imagesc(p_vtk(:,i1)); caxis([0 pmax]); title('t = 0.01 t_c');  set(ax(1),'XTickLabel','');
ax(2) = subplot(1,3,2); imagesc(p_vtk(:,i2));  caxis([0 pmax]); title('t = 0.1 t_c');   axis off;
ax(3) = subplot(1,3,3); imagesc(p_vtk(:,i3));   caxis([0 pmax]); title('t = t_c');       axis off;

x0 = get(ax(1), 'Position'); dx = 0.6 * x0(3); dshift = dx/2;
for i=1:3 
  xleft = x0(1) + (dx+dshift)*(i-1);
  set(ax(i), 'Position', [xleft x0(2) dx x0(4)]);
end;
xleft = x0(1) + (dx+dshift)*3;
h=colorbar; set(h, 'Position', [xleft+dshift/2 x0(2) dx*.6 x0(4)]);
Image = getframe(gcf); imwrite(Image.cdata, 'terzagi10_pressure.png');

%displacement images
figure; set(gcf,'color','white'); set(gcf, 'Position', [400 200 400 400]);
ax(1) = subplot(1,3,1); imagesc(u_vtk(:,i1)); caxis([0 umax]); title('t = 0.01 t_c');  set(ax(1),'XTickLabel','');
ax(2) = subplot(1,3,2); imagesc(u_vtk(:,i2));  caxis([0 umax]); title('t = 0.1 t_c');   axis off;
ax(3) = subplot(1,3,3); imagesc(u_vtk(:,i3));   caxis([0 umax]); title('t = t_c');       axis off;

x0 = get(ax(1), 'Position'); dx = 0.6 * x0(3); dshift = dx/2;
for i=1:3 
  xleft = x0(1) + (dx+dshift)*(i-1);
  set(ax(i), 'Position', [xleft x0(2) dx x0(4)]);
end;
xleft = x0(1) + (dx+dshift)*3;
h=colorbar; set(h, 'Position', [xleft+dshift/2 x0(2) dx*.6 x0(4)]);
Image = getframe(gcf); imwrite(Image.cdata, 'terzagi10_displacement.png');



