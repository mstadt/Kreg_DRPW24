% Plots 3 outputs from driverMeal.m with compute_vars output
% Choose MealOnly, KCl Only and Meal + KCl to show the 3 Preston
% experiments
% This is used in manuscript

clear all;
%% load data
% Meal Only simulation results
f_MealOnly = './MealSim/18-Dec-2023_driverMeal_insulin-1_Kin-0_notes-InsulinOnly.mat';
f_KClOnly = './MealSim/18-Dec-2023_driverMeal_insulin-0_Kin-35_notes-KClOnly.mat';
f_MealKCl = './MealSim/18-Dec-2023_driverMeal_insulin-1_Kin-35_notes-MealKCl.mat';

dat1 = load(f_MealOnly);
dat2 = load(f_KClOnly);
dat3 = load(f_MealKCl);

% Preston et al data file name
f_PrestonDat = './PrestonData/20-Jun-2023_PrestonData.mat';
PrestonDat = load(f_PrestonDat);

lab1 = 'Meal Only';
lab2 = 'KCl Only';
lab3 = 'Meal + KCl';

%% make figures
fprintf('making figures \n')
% figure specs
lw = 3;
f.xlab = 20; f.ylab = 20; f.title = 22;
f.leg = 16; f.gca = 16; f.figlab = 18;
cmap = parula(7);
c1 = cmap(1,:);
c2 = cmap(3,:);
c3 = cmap(5,:);
ls1 = '-'; ls2 = ':'; ls3 = '-.';
cgraymap = gray(5);
cgray = cgraymap(3,:);
lwgray = 2; lsgray = '--';
leglabs = {lab1, lab2, lab3};
xlims = [-6, 8];
xtickvals = -6:2:8;

% time in hours
meal_start = 6;
t1_hrs = dat1.t/60 - meal_start;
t2_hrs = dat2.t/60 - meal_start;
t3_hrs = dat3.t/60 - meal_start;

%% concentrations
figure(1)
clf
nrows = 1; ncols  = 3;
subplot(nrows,ncols,1)
hold on
plot(t1_hrs,dat1.y(:,2)/dat1.pars.V_plasma,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(t2_hrs,dat2.y(:,2)/dat2.pars.V_plasma,'linewidth',lw,'color',c2, 'linestyle',ls2)
plot(t3_hrs,dat3.y(:,2)/dat3.pars.V_plasma,'linewidth',lw,'color',c3, 'linestyle',ls3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize',f.gca)
ylabel('[K^+]_{plasma}', 'fontsize', f.ylab)
xlabel('time (hours)', 'fontsize', f.xlab)
xlim(xlims)
xticks(xtickvals)
yticks(3.5:0.25:5.0)
title('Plasma [K^+]', 'fontsize', f.title)
grid on
legend(leglabs, 'fontsize', f.leg)

subplot(nrows,ncols,2)
hold on
plot(t1_hrs,dat1.y(:,3)/dat1.pars.V_interstitial,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(t2_hrs,dat2.y(:,3)/dat2.pars.V_interstitial,'linewidth',lw,'color',c2, 'linestyle',ls2)
plot(t3_hrs,dat3.y(:,3)/dat3.pars.V_interstitial,'linewidth',lw,'color',c3, 'linestyle',ls3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{inter}', 'fontsize', f.ylab)
xlabel('time (hours)', 'fontsize', f.xlab)
set(gca, 'fontsize',f.gca)
xlim(xlims)
xticks(xtickvals)
yticks(3.5:0.25:5.0)
title('Interstitial [K^+]', 'fontsize', f.title)
grid on
legend(leglabs, 'fontsize', f.leg)

subplot(nrows,ncols,3)
hold on
plot(t1_hrs,dat1.y(:,4)/dat1.pars.V_muscle,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(t2_hrs,dat2.y(:,4)/dat2.pars.V_muscle,'linewidth',lw,'color',c2, 'linestyle',ls2)
plot(t3_hrs,dat3.y(:,4)/dat3.pars.V_muscle,'linewidth',lw,'color',c3, 'linestyle',ls3)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
set(gca, 'fontsize',f.gca)
ylabel('[K^+]_{intracellular}', 'fontsize', f.ylab)
xlabel('time (hours)', 'fontsize', f.xlab)
xlim(xlims)
xticks(xtickvals)
title('Intracellular [K^+]', 'fontsize', f.title)
grid on
legend(leglabs, 'fontsize', f.leg)

AddLetters2Plots(figure(1), {'(a)', '(b)','(c)'},...
                'HShift', -0.05, 'VShift', -0.06, ...
                'fontsize', f.figlab)

