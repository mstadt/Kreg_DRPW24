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

%% Run simulations
% Run 1 day of meals
% if MealTimes(1) > 0
%     % start with a fasting simulation if not starting at 0
%     [tf, yf, ~] = fast_sim(IC0, [0,MealTimes(1)], -60*6, params, opts);
%     T = tf;
%     Y = yf;
% else
%     error('MealTimes(1) = 0') % should not have meal at 0
% end
% % Simulate meals
% for ii = 1:length(MealTimes)
%     [tm, ym, ~] = meal_sim(Y(end,:), MealTimes(ii), len_meal, Meal_Kamts(ii), ...
%                                     params, opts);
%     if ii < length(MealTimes)
%             t_end = MealTimes(ii+1);
%         else
%             t_end = 24 * 60;
%     end
%     [tf, yf, ~] = fast_sim(ym(end,:), [tm(end), t_end], tm(1), ...
%                                     params, opts);
%     T = [T; tm; tf]; Y = [Y; ym; yf];
% end

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