% Use this script to plot 2 results from driver_MultiMeal.m
clear all;

fname1 = './MultiMealSim/26-Jan-2024_driverMultiMeal_doinsulin-1_doFF-1_doMKX-0_notes-try1.mat';
fname2 = './MultiMealSim/26-Jan-2024_driverMultiMeal_doinsulin-1_doFF-1_doMKX-0_notes-higherKamt.mat';
fname3 = './MultiMealSim/26-Jan-2024_driverMultiMeal_doinsulin-1_doFF-1_doMKX-0_notes-highK2.mat';

dat1 = load(fname1);
dat2 = load(fname2);
dat3 = load(fname3);

lab1 = 'normal Kamt';
lab2 = 'high K 1';
lab3 = 'high K 2';


%----------------------
% plot results
%---------------------

fprintf('plotting results \n')
figure(1);
clf;
nr = 2; nc = 2;
lw = 3; lwgray = 2; lsgray = '--';
cmap = parula(7);
c1 = cmap(1,:); c2 = cmap(4,:); c3 = cmap(6,:);
cgraymap = gray(5);
cgray = cgraymap(3,:);
subplot(nr,nc,1)
hold on
plot(dat1.T, dat1.Y(:,1), 'linewidth', lw, color = c1)
plot(dat2.T, dat2.Y(:,1), 'linewidth', lw, color = c2)
plot(dat3.T, dat3.Y(:,1), 'linewidth', lw, color = c3)
xlabel('Time (hrs)')
ylabel('Gut amount')
title('Gut amount')
set(gca,'fontsize',18)
xlim([0,24])
grid on

subplot(nr,nc,2)
hold on
plot(dat1.T,dat1.Y(:,2)/dat1.pars.V_plasma, 'linewidth',lw,'color',c1)
plot(dat2.T,dat2.Y(:,2)/dat2.pars.V_plasma, 'linewidth',lw,'color',c2)
plot(dat3.T,dat3.Y(:,2)/dat3.pars.V_plasma, 'linewidth',lw,'color',c3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (hrs)')
ylabel('Plasma [K^+]')
title('Plasma [K^+]')
set(gca,'fontsize',18)
xlim([0,24])
grid on

subplot(nr,nc,3)
hold on
plot(dat1.T,dat1.Y(:,3)/dat1.pars.V_interstitial,'linewidth',lw,'color',c1)
plot(dat2.T,dat2.Y(:,3)/dat2.pars.V_interstitial,'linewidth',lw,'color',c2)
plot(dat3.T,dat3.Y(:,3)/dat3.pars.V_interstitial,'linewidth',lw,'color',c3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (hrs)')
ylabel('Interstitial [K^+]')
title('Interstitial [K^+]')
set(gca,'fontsize',18)
xlim([0,24])
grid on

subplot(nr,nc,4)
hold on
plot(dat1.T,dat1.Y(:,4)/dat1.pars.V_muscle,'linewidth',lw,'color',c1)
plot(dat2.T,dat2.Y(:,4)/dat2.pars.V_muscle,'linewidth',lw,'color',c2)
plot(dat3.T,dat3.Y(:,4)/dat3.pars.V_muscle,'linewidth',lw,'color',c3)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (hrs)')
ylabel('Intracellular [K^+]')
title('Intracellular [K^+]')
set(gca,'fontsize',18)
xlim([0,24])
grid on

legend(lab1, lab2, lab3)