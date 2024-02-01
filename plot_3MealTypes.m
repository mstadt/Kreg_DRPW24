% Plots 3 outputs from driverMeal.m with compute_vars output
% Choose MealOnly, KCl Only and Meal + KCl to show the 3 Preston
% experiments

clear all;
%---------------------
% Begin user input
%---------------------
% File names where simulation results are stored
% Meal Only simulation results
f_MealOnly = './MealSim/01-Feb-2024_driverMeal_insulin-1_Kin-0_notes-insonly.mat';
lab1 = 'Meal Only'; % label for figures


% KCl Only simulation results
f_KClOnly = './MealSim/01-Feb-2024_driverMeal_insulin-0_Kin-35_notes-kclonly.mat';
lab2 = 'KCl Only'; % label for figures

% Meal + KCl simulation results
f_MealKCl = './MealSim/01-Feb-2024_driverMeal_insulin-1_Kin-35_notes-mealkcl.mat';
lab3 = 'Meal + KCl'; % label for figures
%----------------------
% End user input
%----------------------


%% Load Data
dat1 = load(f_MealOnly);
dat2 = load(f_KClOnly);
dat3 = load(f_MealKCl);

% Preston et al data file name
f_PrestonDat = './PrestonData/20-Jun-2023_PrestonData.mat';
PrestonDat = load(f_PrestonDat);


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


%% other figures
%% concentrations
figure(11)
clf
nrows = 1; ncols  = 3;
subplot(nrows,ncols,1)
hold on
plot(dat1.t,dat1.y(:,2)/dat1.pars.V_plasma,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,2)/dat2.pars.V_plasma,'linewidth',lw,'color',c2, 'linestyle',ls2)
plot(dat3.t,dat3.y(:,2)/dat3.pars.V_plasma,'linewidth',lw,'color',c3, 'linestyle',ls3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{plasma}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Plasma K concentration', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
hold on
plot(dat1.t,dat1.y(:,3)/dat1.pars.V_interstitial,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,3)/dat2.pars.V_interstitial,'linewidth',lw,'color',c2, 'linestyle',ls2)
plot(dat3.t,dat3.y(:,3)/dat3.pars.V_interstitial,'linewidth',lw,'color',c3, 'linestyle',ls3)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{inter}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Interstitial K concentration', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
hold on
plot(dat1.t,dat1.y(:,4)/dat1.pars.V_muscle,'linewidth',lw,'color',c1, 'linestyle',ls1)
plot(dat2.t,dat2.y(:,4)/dat2.pars.V_muscle,'linewidth',lw,'color',c2, 'linestyle',ls2)
plot(dat3.t,dat3.y(:,4)/dat3.pars.V_muscle,'linewidth',lw,'color',c3, 'linestyle',ls3)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{muscle}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Muscle K concentration', 'fontsize', f.title)
grid on
legend(leglabs, 'fontsize', f.leg)

%% feedforward and feedback response
figure(15)
clf
nrows = 2; ncols = 2;
subplot(nrows,ncols,1)
hold on
rhoins1 = [dat1.vals1.rho_insulin; dat1.vals2.rho_insulin; dat1.vals3.rho_insulin];
rhoins2 = [dat2.vals1.rho_insulin; dat2.vals2.rho_insulin; dat2.vals3.rho_insulin];
rhoins3 = [dat3.vals1.rho_insulin; dat3.vals2.rho_insulin; dat3.vals3.rho_insulin];
plot(dat1.t, rhoins1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, rhoins2, 'linewidth',lw,'color',c2,'linestyle',ls2)
plot(dat3.t, rhoins3, 'linewidth',lw,'color',c3,'linestyle',ls3)
ylabel('\rho_{ins}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('[insulin] effect on \Phi_{ECtoIC} (\rho_{insulin})', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
hold on
rhoal1 = [dat1.vals1.rho_al; dat1.vals2.rho_al; dat1.vals3.rho_al];
rhoal2 = [dat2.vals1.rho_al; dat2.vals2.rho_al; dat2.vals3.rho_al];
rhoal3 = [dat3.vals1.rho_al; dat3.vals2.rho_al; dat3.vals3.rho_al];
plot(dat1.t, rhoal1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, rhoal2,'linewidth',lw,'color',c2,'linestyle',ls2)
plot(dat3.t, rhoal3,'linewidth',lw,'color',c3,'linestyle',ls3)
ylabel('\rho_{al}', 'fontsize', f.ylab)
xlabel('t','fontsize',f.xlab)
title('[ALD] effect on \Phi_{ECtoIC} (\rho_{al})', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
hold on
gamkin1 = [dat1.vals1.gamma_Kin; dat1.vals2.gamma_Kin; dat1.vals3.gamma_Kin];
gamkin2 = [dat2.vals1.gamma_Kin; dat2.vals2.gamma_Kin; dat2.vals3.gamma_Kin];
gamkin3 = [dat3.vals1.gamma_Kin; dat3.vals2.gamma_Kin; dat3.vals3.gamma_Kin];
plot(dat1.t, gamkin1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, gamkin2,'linewidth',lw,'color',c2,'linestyle',ls2)
plot(dat3.t, gamkin3, 'linewidth',lw,'color',c3,'linestyle',ls3)
ylabel('\gamma_{Kin}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('GI feedforward effect (\gamma_{Kin})', 'fontsize', f.title)
grid on

subplot(nrows,ncols,4)
hold on
gamal1 = [dat1.vals1.gamma_al; dat1.vals2.gamma_al; dat1.vals3.gamma_al];
gamal2 = [dat2.vals1.gamma_al; dat2.vals2.gamma_al; dat2.vals3.gamma_al];
gamal3 = [dat3.vals1.gamma_al; dat3.vals2.gamma_al; dat3.vals3.gamma_al];
plot(dat1.t, gamal1, 'linewidth', lw,'color',c1,'linestyle',ls1)
plot(dat2.t, gamal2, 'linewidth',lw,'color',c2,'linestyle',ls2)
plot(dat3.t, gamal3, 'linewidth',lw,'color',c3,'linestyle',ls3)
ylabel('gamma_{al}', 'fontsize', f.ylab)
xlabel('t','fontsize',f.xlab)
title('[ALD] effect on \Phi_{dt-Ksec} (\gamma_{al})', 'fontsize', f.title)
grid on

legend(leglabs, 'fontsize', f.leg)

%% kidney
figure(21)
clf
nrows = 2; ncols = 2;
subplot(nrows,ncols,1)
hold on
filK1 = [dat1.vals1.filK; dat1.vals2.filK; dat1.vals3.filK];
filK2 = [dat2.vals1.filK; dat2.vals2.filK; dat2.vals3.filK];
filK3 = [dat3.vals1.filK; dat3.vals2.filK; dat3.vals3.filK];
plot(dat1.t, filK1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, filK2, 'linewidth',lw,'color',c2,'linestyle',ls2)
plot(dat3.t, filK3, 'linewidth',lw,'color',c3,'linestyle',ls3)
ylabel('fil_K', 'fontsize', f.ylab)
xlabel('t','fontsize',f.xlab)
title('Filtered K^+ load', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
hold on
dtKsec1 = [dat1.vals1.dtKsec; dat1.vals2.dtKsec; dat1.vals3.dtKsec];
dtKsec2 = [dat2.vals1.dtKsec; dat2.vals2.dtKsec; dat2.vals3.dtKsec];
dtKsec3 = [dat3.vals1.dtKsec; dat3.vals2.dtKsec; dat3.vals3.dtKsec];
plot(dat1.t, dtKsec1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, dtKsec2, 'linewidth',lw,'color',c2,'linestyle',ls2)
plot(dat3.t, dtKsec3, 'linewidth',lw,'color',c3,'linestyle',ls3)
ylabel('\Phi_{dt-Ksec}', 'fontsize', f.ylab)
xlabel('t','fontsize',f.xlab)
title('DT K^+ secretion', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
hold on
cdKtrans1 = [dat1.vals1.cdKsec; dat1.vals2.cdKsec; dat1.vals3.cdKsec] ...
    - [dat1.vals1.cdKreab; dat1.vals2.cdKreab; dat1.vals3.cdKreab];
cdKtrans2 = [dat2.vals1.cdKsec; dat2.vals2.cdKsec; dat2.vals3.cdKsec] ...
    - [dat2.vals1.cdKreab; dat2.vals2.cdKreab; dat2.vals3.cdKreab];
cdKtrans3 = [dat3.vals1.cdKsec; dat3.vals2.cdKsec; dat3.vals3.cdKsec]...
    - [dat3.vals1.cdKreab; dat3.vals2.cdKreab; dat3.vals3.cdKreab];
plot(dat1.t, cdKtrans1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, cdKtrans2, 'linewidth',lw,'color',c2,'linestyle',ls2)
plot(dat3.t, cdKtrans3, 'linewidth',lw,'color',c3,'linestyle',ls3)
ylabel('\Phi_{cd-Ksec} - \Phi_{cd-Kreab}', 'fontsize', f.ylab)
xlabel('t','fontsize',f.xlab)
title('CD K^+ transport', 'fontsize', f.title)
grid on

subplot(nrows,ncols,4)
hold on
uK1 = [dat1.vals1.UrineK; dat1.vals2.UrineK; dat1.vals3.UrineK];
uK2 = [dat2.vals1.UrineK; dat2.vals2.UrineK; dat2.vals3.UrineK];
uK3 = [dat3.vals1.UrineK; dat3.vals2.UrineK; dat3.vals3.UrineK];
plot(dat1.t, uK1, 'linewidth',lw,'color',c1,'linestyle',ls1)
plot(dat2.t, uK2, 'linewidth',lw,'color',c2,'linestyle',ls2)
plot(dat3.t, uK3, 'linewidth',lw,'color',c3,'linestyle',ls3)
ylabel('\Phi_{uK}', 'fontsize', f.ylab)
xlabel('t','fontsize',f.xlab)
title('Urine K^+ excretion', 'fontsize', f.title)
grid on

