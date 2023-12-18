function v = compute_vars(tvals, yvals, params, varargin)
% computes algebraic equations given output of ODE solver
% Input: tvals, yvals -- output from ODE solver, params -- parameter vector

% set variable names
MKgut_vals    = yvals(:,1); % amount of K in gut
MKplas_vals   = yvals(:,2); % amount of K in plasma
MKinter_vals  = yvals(:,3); % amount of K in interstitial space
MKmuscle_vals = yvals(:,4); % amount of K in muscle

% fprintf('WARNING>>> MAY NEED TO UPDATE EQUATIONS!')

% set parameter names
Phi_Kin_ss = params(1);
t_insulin_ss = params(2);
fecal_excretion = params(3);
kgut = params(4);
MKgutSS = params(5);
V_plasma = params(6);
V_interstitial = params(7);
V_muscle = params(8);
Kecf_total = params(9);
P_ECF = params(10);
Kmuscle_baseline = params(11);
Vmax = params(12);
Km = params(13);
P_muscle = params(14);
GFR_base = params(15);
eta_ptKreab_base = params(16);
eta_LoHKreab = params(17);
dtKsec_eq = params(18);
A_dtKsec = params(19);
B_dtKsec = params(20);
cdKsec_eq = params(21);
A_cdKsec = params(22);
B_cdKsec = params(23);
alpha_TGF = params(24);
A_cdKreab = params(25);
ALD_eq = params(26);
m_K_ALDO = params(27);
FF = params(28);
A_insulin = params(29);
B_insulin = params(30);

%% Get variable inputs
% default settings, varargin is used to change settings
SS = false; % compute SS solution
alt_sim = false; % use alternate equations
do_insulin = true;
do_FF = true;
MKX = 0;
Kintake = 0;
meal_start = 0;
highK_eff = 0;
TGF_eff = 0;
for i = 1:2:length(varargin)
    temp = varargin{i+1};
    if strcmp(varargin{i}, 'SS')
        SS = temp;
    elseif strcmp(varargin{i}, 'alt_sim')
        alt_sim = temp;
    elseif strcmp(varargin{i}, 'do_MKX')
        MKX = temp(1);
        MKslope = temp(2);
    elseif strcmp(varargin{i}, 'do_insulin')
        do_insulin = temp(1);
    elseif strcmp(varargin{i}, 'do_FF')
        do_FF = temp(1);
    elseif strcmp(varargin{i}, 'Kintake')
        Kintake = temp(1);
    elseif strcmp(varargin{i}, 'meal_time')
        meal_start = temp(1);
%     elseif strcmp(varargin{i}, 'highK_eff')
%         highK_eff = temp(1);
    elseif strcmp(varargin{i}, 'TGF_eff')
        TGF_eff = temp(1);
        alpha_TGF = temp(2);
        eta_ptKreab = temp(3);
%         if TGF_eff
%             fprintf('doing TGF_eff \n')
%         end
    else
        disp('WRONG VARARGIN INPUT')
        fprintf('What is this varargin input? %s \n', varargin{i})
        error('wrong varargin input')
    end % if
end %for

% set insulin level
if do_insulin
    if SS
        t_insulin = t_insulin_ss.*ones(size(tvals));
    else
        t_insulin = tvals - meal_start;
    end
    v.C_insulin = zeros(size(t_insulin));
    for ii = 1:length(v.C_insulin)
        v.C_insulin(ii) = get_Cinsulin(t_insulin(ii));
    end
else
    v.C_insulin = ones(size(tvals)).*22.6/1000; % steady state insulin
end

% concentrations
v.K_plas = MKplas_vals./V_plasma;
v.K_inter = MKinter_vals./V_interstitial;
v.K_muscle = MKmuscle_vals./V_muscle;
v.K_ECFtot = (MKplas_vals + MKinter_vals)./(V_plasma + V_interstitial);

% ALD
%v.xi_ksod = max(0,((v.K_ECFtot./Csod)./(Kecf_total/144/(xi_par+1))-xi_par));
%v.N_als = v.xi_ksod;

% Gut K
v.Gut2plasma = kgut.*MKgut_vals;

% ALD impact
v.Nal_vals = exp(m_K_ALDO .* (v.K_ECFtot - Kecf_total));
v.C_al = v.Nal_vals .* ALD_eq;
v.gamma_al = A_dtKsec .* v.C_al .^B_dtKsec;
v.lambda_al = A_cdKsec .* v.C_al.^B_cdKsec;

if do_FF
    v.gamma_Kin = max(1, FF.*(MKgut_vals - MKgutSS));
else
    v.gamma_Kin = ones(size(MKgut_vals));
end

if MKX > 0
    v.omegaKic = max(0, (MKslope.*(v.K_muscle - Kmuscle_baseline) + 1));
else
    v.omegaKic = ones(size(v.K_muscle));
end

% renal K handling
% if highK_eff
%     if highK_eff == 1
%         GFR = (1 - 0.29) * 0.125;
%         etapsKreab = 0.36 + 0.25; % PT + TAL part
%     elseif highK_eff == 2
%         etapsKreab = 0.36 + 0.25; % PT + TAL part
%     elseif highK_eff == 3
%         GFR = (1 - 0.29) * 0.125;
%     end
% end

eta_psKreab_base = eta_ptKreab_base + eta_LoHKreab;
if TGF_eff == 1
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    v.GFR = GFR_base + alpha_TGF * (eta_psKreab - eta_psKreab_base);
elseif TGF_eff == 2 % GFR only
    eta_ptKreab = eta_ptKreab_base; % PT K reab is baseline value
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    v.GFR = GFR_base + alpha_TGF * (eta_psKreab - eta_psKreab_base);
elseif TGF_eff == 3 % PT only
    eta_psKreab = eta_ptKreab + eta_LOHKreab;
    v.GFR = GFR_base;
else
    eta_ptKreab = eta_ptKreab_base; % PT K reab is baseline value
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    v.GFR = GFR_base; 
end

v.filK = v.GFR .* v.K_plas;

v.psKreab = eta_psKreab * v.filK;


% v.filK = GFR .* v.K_plas;
% v.psKreab = etapsKreab .* v.filK;

% distal tubule
if MKX == 1
    v.eta_dtKsec = v.gamma_al .* v.gamma_Kin .* v.omegaKic;
else
    v.eta_dtKsec = v.gamma_al .* v.gamma_Kin;
end
v.dtKsec = dtKsec_eq .* v.eta_dtKsec;

% collecting duct
if MKX == 2
    v.eta_cdKsec = v.lambda_al .* v.omegaKic;
else
    v.eta_cdKsec = v.lambda_al;
end
v.cdKsec = cdKsec_eq .* v.eta_cdKsec;

if MKX == 3
    v.eta_cdKreab = v.omegaKic;
else
    v.eta_cdKreab = ones(size(v.omegaKic));
end

dtK = v.filK - v.psKreab + v.dtKsec;
v.cdKreab = dtK .* A_cdKreab.*v.eta_cdKreab;

v.UrineK = dtK + v.cdKsec - v.cdKreab;

% interstitial K
v.rho_al = (66.4 + 0.273.* v.C_al)./89.6050;
% insulin
L = 100*ones(size(v.C_insulin)); x0 = 0.5381 * ones(size(v.C_insulin)); k = 1.069;
ins_A = A_insulin * ones(size(v.C_insulin)); ins_B = 100*B_insulin * ones(size(v.C_insulin));
temp = (ins_A.*(L./(1+exp(-k.*(log10(v.C_insulin)-log10(x0)))))+ ins_B)./100;
if do_insulin
    v.rho_insulin = max(1.0, temp);
else
    v.rho_insulin = ones(size(v.C_insulin));
end
v.eta_NKA = v.rho_insulin .* v.rho_al;

v.Inter2Muscle = v.eta_NKA .*((Vmax * v.K_inter)./(Km + v.K_inter));
v.Muscle2Inter = P_muscle.*(v.K_muscle - v.K_inter);







end % compute_vars