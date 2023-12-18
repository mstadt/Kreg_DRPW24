function pars = set_params()
% this file sets the current parameter values 
%% K intake at SS
pars.Phi_Kin_ss        = 70/1440; %mEq/min, steady state for Phi_Kin (Preston 2015)
pars.t_insulin_ss      = 270; % ss t_insulin value

%% gut parameters
pars.fecal_excretion = 0.1;
pars.kgut = 0.01;
pars.MKgutSS = (0.9*pars.Phi_Kin_ss)/pars.kgut;
%% volumes
pars.V_plasma          = 4.5; %plasma fluid volume (L)
pars.V_interstitial    = 10; % interstitial ECF volume (L)
pars.V_muscle          = 24; %19.0; % intracellular fluid volume (L)

%% baseline concentrations
Kecf_baseline      = 4.3;% baseline ECF K concentration (total) mEq/L
pars.Kecf_total    = 4.2;
pars.P_ECF             = 0.3; % this parameter will have to be fit I think

pars.Kmuscle_baseline       = 130; % baseline muscle concentration mEq/L
%% NKA activity values
pars.Vmax              = 130;% mmol/min Cheng 2013
pars.Km                = 1.4; % mmol/L (Cheng 2013 gives between 0.8 and 1.5)

%% compute permeability values
NKA_baseline = pars.Vmax*Kecf_baseline/(pars.Km + Kecf_baseline);
pars.P_muscle = (NKA_baseline)/(pars.Kmuscle_baseline - Kecf_baseline);


%% Kidney
pars.GFR_base   = 0.125; %baseline GFR L/min

pars.eta_ptKreab_base = 0.67; % fractional PT K reabsorption baseline

pars.eta_LoHKreab = 0.25; % fractional LoH K reabsorption, fixed

pars.dtKsec_eq = 0.041;
pars.A_dtKsec = 0.3475;
pars.B_dtKsec = 0.23792;

pars.cdKsec_eq = 0.0022;
pars.A_cdKsec = 0.161275;
pars.B_cdKsec = 0.410711;

%% TGF response

GFR0 = pars.GFR_base; % baseline GFR
Tong_HighKeff = 0.29; % Reduction of GFR with high K from Tong
HighKeff_etaPT = 0.36; % fractional PT K reabsorption at high K (from Tong)
Base_etaPT = 0.67; % baseline fractional PT K reabsorption

HighKeff_etaPS = HighKeff_etaPT + pars.eta_LoHKreab; % proximal segment reabsorption (high K)
Base_etaPS = Base_etaPT + pars.eta_LoHKreab; % proximal segment reabsorption (base)
% This is baseline alpha_TGF
pars.alpha_TGF = ((1-Tong_HighKeff)*GFR0 - GFR0) / (HighKeff_etaPS - Base_etaPS);



%% parameters A and B are divided by 1000 and 100 respectively in k_reg_mod
% because otherwise, when fitting the parameters, the steps would be too small. 
pars.A_cdKreab = 0.499994223625298; 

%% ALD
pars.ALD_eq = 85; % ng/L
% pars.T_al = 60; % ALD half life (min)
% pars.Csod = 144; % sodium concentration mEq/L
% pars.xi_par = 2; %lower xi_pars makes C_al less sensitive

pars.m_K_ALDO = 0.5; %951.2/1000; % adjust for mL to L

%% effects
pars.FF = 0.250274; 

pars.A_insulin = 0.999789;
pars.B_insulin =  0.6645;
end %set_params