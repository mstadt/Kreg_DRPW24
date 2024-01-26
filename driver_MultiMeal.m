% This script can be used to simulate multiple meals throughout a day
% The output is saved to ./MultiMealSim/

% clear
clear all;

%--------------------
% User input
%-------------------
Kamt_meal = 78/3; % amount of K in each meal
len_meal = 30;

doFF = 1; % do GI feedforward effect?
doinsulin = 1; % set to 0 for no insulin effect

% muscle-kidney cross talk options
doMKX = 0; % MKX option -- 0: no MKX; 1: dtKsec; 2: cdKsec, 3: cdKreab
%-------------------
% End of user input
%-------------------

%% initialize parameter values
fprintf('loading params \n')
pars = set_params();
[params, parnames] = pars2vector(pars,0);

%% simulation settings
% Meal times
MealTimes = [6, 12, 18] * 60; % meal times 
Meal_Kamts = [Kamt_meal, Kamt_meal, Kamt_meal];%Kamt_meal*ones(size(MealTimes)); % amounts of K per meal
opts.do_insulin = doinsulin;
opts.do_FF = doFF;
if doMKX > 0
    if doMKX == 1
        % slope tries (0.005, 0.01, 0.025, 0.05, 0.075, 0.1)
        MKXslope = 0.1; % dtKsec slope
    elseif doMKX == 2
        MKXslope = 0.1; % cdKsec slope
    elseif doMKX == 3
        MKXslope = -0.1; % cdKreab slope
    end
else
    MKXslope = -1;
end
opts.do_MKX = [doMKX, MKXslope];

%% set initial conditions
temp = load('./SS/SS_4vars.mat');
SS = temp.SS;
[IC, ~, ~] = getSS(SS, params, ...
    'do_insulin', opts.do_insulin,...
    'do_FF',opts.do_FF,...
    'do_MKX', opts.do_MKX,...
    'do_figs', 0); % start at SS

%IC0 = IC;

%% Run simulations of intake for multiple meals
fprintf('run 3 meals \n')
% Run 1 day of meals
% Fasting simulation to start
[tf, yf, ~] = fast_sim(IC, ...
                    [0,MealTimes(1)],...
                    -60*60, ... % insulin is low
                    params,...
                    opts);

% store full simulation in T and Y
T = tf; 
Y = yf;

% Do meal simulations for the day
for ii = 1:length(MealTimes)
    % Run meal simulation
    [tm,ym,~] = meal_sim(Y(end,:), MealTimes(ii), len_meal, Meal_Kamts(ii),...
                                params, opts);

    % set end point for next fasting step
    if ii<length(MealTimes)
        t_end = MealTimes(ii+1);
    else
        t_end = 24 * 60; % end of day time
    end
    % Run fasting simulation
    [tf,yf, ~] = fast_sim(ym(end,:),[tm(end),t_end],tm(1),...
                        params,opts);
    % append meal and fasting simulation to T and Y
    T = [T;tm;tf]; Y = [Y;ym;yf];
end
fprintf('simulation done \n')

%----------------------
% plot results
%---------------------
T = T./60; % change time to hours
fprintf('plotting results \n')
figure(1);
clf;
nr = 2; nc = 2;
lw = 3; lwgray = 2; lsgray = '--';
cmap = parula(6);
c1 = cmap(1,:); c2 = cmap(2,:); c3 = cmap(3,:); c4 = cmap(4,:);
cgraymap = gray(5);
cgray = cgraymap(3,:);
subplot(nr,nc,1)
plot(T, Y(:,1), 'linewidth', lw, color = c1)
xlabel('Time (hrs)')
ylabel('Gut amount')
title('Gut amount')
set(gca,'fontsize',18)
xlim([0,24])
grid on

subplot(nr,nc,2)
hold on
plot(T,Y(:,2)/pars.V_plasma, 'linewidth',lw,'color',c2)
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
plot(T,Y(:,3)/pars.V_interstitial,'linewidth',lw,'color',c3)
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
plot(T,Y(:,4)/pars.V_muscle,'linewidth',lw,'color',c4)
yline(120,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
yline(140,'color',cgray,'linestyle',lsgray, 'linewidth', lwgray)
xlabel('Time (hrs)')
ylabel('Intracellular [K^+]')
title('Intracellular [K^+]')
set(gca,'fontsize',18)
xlim([0,24])
grid on

%% save simulations options
save_sim = input('Do you want to save the simulation? (0 - no/1 - yes) ');
if save_sim
    notes = input('notes: ');
    fname = strcat('./MultiMealSim/', date, '_driverMultiMeal', ...
                    '_doinsulin-',num2str(doinsulin),...
                    '_doFF-', num2str(doFF),...
                    '_doMKX-',num2str(doMKX),...
                    '_notes-', notes, ...
                    '.mat');
    save(fname, 'T', 'Y', ...
                'pars', 'params', 'parnames',...
                'MealTimes', 'Meal_Kamts', ... % meal settings
                'opts') % simulation settings
    fprintf('results saved to: \n %s \n', fname)
end


%% Functions used
%---------------------
% Functions
%---------------------
% These functions are used within this script
% Meal simulation
function [t, y, vals] = meal_sim(IC, t0, len_meal, Kamt, params, opts)
    % Meal simulation
    % Inputs:
    %    IC -- initial condition
    %    t0 -- time to start simulation
    %    len_meal -- how long meal is (minutes)
    %    Kamt --- how much K in the meal
    %    params -- parameters
    %    opts -- simulation settings

    % Output:
    %   t -- timepoints
    %   y -- variable outpus
    %   vals -- extra variable outputs
    tf = t0 + len_meal; % meal length
    tspan = [t0, tf];
    Kintake = Kamt / (tf - t0);
    options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9); % ode solver settings
    [t, y] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                                'do_insulin', opts.do_insulin,...
                                'do_FF', opts.do_FF,...
                                'meal_time', t0,... % t0 is start of meal
                                'Kintake', Kintake,...
                                'do_MKX', opts.do_MKX),...
                                tspan, IC, options);
    vals = compute_vars(t,y,params,...
                            'do_insulin', opts.do_insulin,...
                            'do_FF', opts.do_FF, ...
                            'meal_time', t0,...
                            'do_MKX', opts.do_MKX);
end

% Fasting simulation
function [t, y, vals] = fast_sim(IC, tspan, last_meal, params, opts)
    % Fasting simulation
    % IC -- initial condition
    % tspan -- [start_time, end_time]
    % last_meal -- start time of the last meal
    % params -- parameters
    % opts -- simulation settings

    %  t -- timesteps
    %  y -- variable outputs
    %  vals -- extra variable outputs
    options = odeset('RelTol', 1.0e-6, 'AbsTol', 1e-9); % ode solver settings
    [t, y] = ode15s(@(t,y) kreg_eqns(t,y,params,...
                                'do_insulin', opts.do_insulin,...
                                'do_FF', opts.do_FF,...
                                'meal_time', last_meal,...
                                'Kintake', 0,... % fasting state
                                'do_MKX', opts.do_MKX),... 
                                tspan, IC, options);
    vals = compute_vars(t,y,params,...
                            'do_insulin', opts.do_insulin,...
                            'do_FF', opts.do_FF,...
                            'meal_time', last_meal,...
                            'Kintake', 0,...
                            'do_MKX', opts.do_MKX); 
end