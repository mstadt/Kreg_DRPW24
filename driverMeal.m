% This script can be used to simulate a single meal.
% The output is saved to ./MealSim/
% Save output to ./MealSim/
clear all;
%--------------------
% User input
%-------------------
MealInsulin = 0; % set to 0 for no insulin
Kamt = 35; % amount of K in meal

% muscle-kidney cross talk options
MKX = 0; MKXslope = 0; 
%-------------------
% End of user input
%-------------------

%% initialize parameter values
fprintf('loading params \n')
pars = set_params();
[params, parnames] = pars2vector(pars,0);

%% set initial conditions
temp = load('./SS/SS_4vars.mat');
SS = temp.SS;
[IC, ~, ~] = getSS(SS, params, 'do_figs', 0); % start at SS

%% Fasting state
% ODE options
t0 = 0;
tf = 6*60; % 6 hours of fasting
tspan = [t0, tf]; % time of simulation
options = odeset('RelTol',1.0e-6,'AbsTol',1e-9); % ode solver settings

% simulation settings
do_insulin = 0; % set to 1 when doing meal
do_FF      = 1; % set to 1 unless no FF effect
alt_sim    = 0; % only if have other options

fprintf('solving ODEs \n')
[t1,y1] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                            'alt_sim', alt_sim,...
                            'do_MKX', [MKX, MKXslope],...
                            'do_insulin', do_insulin,...
                            'do_FF', do_FF,...
                            'Kintake', 0), ...
                            tspan, IC, options);

vals1 = compute_vars(t1,y1,params,...
                        'do_insulin', do_insulin,...
                        'do_FF', do_FF,...
                        'do_MKX', [MKX, MKXslope]);

%% add a meal with glucose
% NOTE: need to restart t from 0 to get insulin dynamics, shift later
t0 = 0; 
tf = t0 + 30;
tspan = [t0, tf];
IC = y1(end,:);
Kintake = Kamt/(tf - t0);
do_insulin = MealInsulin; % add glucose componenent
[t2,y2] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                            'alt_sim', alt_sim,...
                            'do_MKX', [MKX, MKXslope],...
                            'do_insulin', do_insulin,...
                            'do_FF', do_FF, ...
                            'Kintake', Kintake), ...
                            tspan, IC, options);

vals2 = compute_vars(t2,y2,params,...
                        'do_insulin', do_insulin,...
                        'do_FF', do_FF,...
                        'do_MKX', [MKX, MKXslope]);
% done K intake
IC = y2(end,:);
t0 = t2(end);
tf = t0 + 1000 - 30;
tspan = [t0, tf];
Kintake = 0;
[t3,y3] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                            'alt_sim', alt_sim,...
                            'do_MKX', [MKX, MKXslope],...
                            'do_insulin', do_insulin,...
                            'do_FF', do_FF, ...
                            'Kintake', Kintake), ...
                            tspan, IC, options);

t2 = t2 + t1(end); % shift by t1
vals3 = compute_vars(t3,y3,params,...
                        'do_insulin', do_insulin,...
                        'do_FF', do_FF,...
                        'do_MKX', [MKX, MKXslope]);


t3 = t3 + t1(end); % shift by t1
t = [t1;t2;t3];
y = [y1;y2;y3];


%% make figures
% load old data
fprintf('making figures \n')
% figure specs
lw = 3;
f.xlab = 16; f.ylab = 16; f.title = 18;
f.leg = 16;
cmap = spring(5);
c1 = cmap(1,:);
c2 = cmap(3,:);
ls1 = '-'; ls2 = '-';
cgraymap = gray(5);
cgray = cgraymap(3,:);
lwgray = 2; lsgray = '--';
leglabs = {'old data', 'new sims'};

% variables
figure(1)
clf
nrows = 2; ncols  = 3;
subplot(nrows,ncols,1)
hold on
plot(t,y(:,1),'linewidth',lw,'color',c2, 'linestyle',ls2)
ylabel('M_{Kgut}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Gut K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
hold on
plot(t,y(:,2),'linewidth',lw,'color',c2, 'linestyle',ls2)
ylabel('M_{Kplas}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Plasma K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
hold on
plot(t,y(:,3),'linewidth',lw,'color',c2, 'linestyle',ls2)
ylabel('M_{Kinter}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Interstitial K', 'fontsize', f.title)
grid on

subplot(nrows,ncols,4)
hold on
plot(t,y(:,4),'linewidth',lw,'color',c2, 'linestyle',ls2)
ylabel('M_{Kmuscle}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Muscle K', 'fontsize', f.title)
grid on

% subplot(nrows,ncols,5)
% hold on
% plot(t,y(:,5),'linewidth',lw,'color',c2, 'linestyle',ls2)
% ylabel('N_{al}', 'fontsize', f.ylab)
% xlabel('t', 'fontsize', f.xlab)
% title('Normalized ALD', 'fontsize', f.title)
% grid on

% concentrations
figure(2)
clf
nrows = 1; ncols  = 3;
subplot(nrows,ncols,1)
hold on
plot(t,y(:,2)/pars.V_plasma,'linewidth',lw,'color',c2, 'linestyle',ls2)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{plasma}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Plasma K concentration', 'fontsize', f.title)
grid on

subplot(nrows,ncols,2)
hold on
plot(t,y(:,3)/pars.V_interstitial,'linewidth',lw,'color',c2, 'linestyle',ls2)
yline(3.5,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(5.0,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{inter}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Interstitial K concentration', 'fontsize', f.title)
grid on

subplot(nrows,ncols,3)
hold on
plot(t,y(:,4)/pars.V_muscle,'linewidth',lw,'color',c2, 'linestyle',ls2)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
ylabel('[K^+]_{muscle}', 'fontsize', f.ylab)
xlabel('t', 'fontsize', f.xlab)
title('Muscle K concentration', 'fontsize', f.title)
grid on



%% save simulations options
save_sim = input('Do you want to save the simulation? (0 - no/1 - yes) ');
if save_sim
    notes = input('notes: ');
    fname = strcat('./MealSim/', date, '_driverMeal', ...
                    '_insulin-',num2str(MealInsulin),...
                    '_Kin-', num2str(Kamt), ...
                    '_notes-', notes, ...
                    '.mat');
    save(fname, 't', 'y', ...
                'pars', 'params', 'parnames',...
                'Kamt', 'MealInsulin', 'MKX', 'MKXslope',...
                'IC',...
                'vals1', 'vals2', 'vals3',...
                'options')
    fprintf('results saved to: \n %s \n', fname)
end
