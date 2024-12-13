clear; close all; clc
set(0,'defaultTextInterpreter','latex')

nsteps = 2000; dt = 0.0005;
[r,wx,scl,usq] = importfrdata('fr.dat',nsteps+1);

figure(1)
subplot(1,3,1)
p1 = plot(r(:,1),usq(:,1),'b','LineWidth',2);
xlabel('$r$','FontSize',24)
ylabel('$u_{\theta}$','FontSize',24)
ylim([0 max(usq(:,1))])

subplot(1,3,2)
p2 = plot(r(:,1),sqrt(wx(:,1).^2),'r','LineWidth',2);
xlabel('$r$','FontSize',24)
ylabel('$|\omega_x|$','FontSize',24)
ylim([0 max(abs(wx(:,1)))])

subplot(1,3,3)
p3 = plot(r(:,1),scl(:,1),'color',[0 0.5 0],'LineWidth',2);
xlabel('$r$','FontSize',24)
ylabel('$\gamma$','FontSize',24)
ylim([0 max(scl(:,1))])

for ii = 1:nsteps
    p1.YData = usq(:,ii+1);
    p2.YData = abs(wx(:,ii+1));
    p3.YData = scl(:,ii+1);
    drawnow
end

%%
clc
nu = 0.01; sigma0 = 0.1; lambda = 0.1; % Example values
tvec = (0:dt:dt*nsteps)/lambda;
rvec = r(:,1)/sigma0;

% Vortex core
vcore = sqrt(nu*tvec + (sigma0)^2);
[rm,tm] = meshgrid(rvec',tvec');
figure(2)
clf
contourf(tm,rm,scl','LineStyle','none')
hold on
plot(tvec,vcore/sigma0,'r','LineWidth',2)
set(gca,'FontWeight','bold','FontSize',18)
cb = colorbar(gca,"eastoutside");
%cb.Label = '$\overline{\gamma}$';
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 18;
ylim([0 5])
xlabel('$t/\lambda$','fontsize',24)
ylabel('$r/\sigma_0$','FontSize',24)
title('$\overline{\gamma}$','FontSize',24)