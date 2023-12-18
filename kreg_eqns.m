function dydt = kreg_eqns(t,y,params,varargin)
% K regulation model equations

%% variable names
M_Kgut    = y(1); % amount of K in gut
M_Kplas   = y(2); % amount of K in plasma
M_Kinter  = y(3); % amount of K in interstitial space
M_Kmuscle = y(4); % amount of K in muscle


%% set parameter names
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
        t_insulin = t_insulin_ss;
    else
        t_insulin = t - meal_start;
    end
    C_insulin = get_Cinsulin(t_insulin);
else
    C_insulin = 22.6/1000; % steady state insulin
end


%% model equations
dydt = zeros(length(y),1);

% concentrations
K_plas    = M_Kplas/V_plasma; % plasma K concentration
K_inter   = M_Kinter/V_interstitial; % interstitial K concentration
K_muscle  = M_Kmuscle/V_muscle; % intracellular K concentration
K_ECFtot  = (M_Kplas + M_Kinter)/(V_plasma + V_interstitial); % total ECF concentration

%% ALD (N_al)
% equation from Maddah & Hallow 2022
N_al = exp(m_K_ALDO * (K_ECFtot - Kecf_total));
C_al = N_al*ALD_eq;

%% Gut K (M_Kgut)
if SS
    Phi_Kin = Phi_Kin_ss;
else
    Phi_Kin = Kintake;
end
K_intake   = (1-fecal_excretion)*Phi_Kin;
Gut2plasma = kgut*M_Kgut;
dydt(1) = K_intake - Gut2plasma;


%% Plasma K (M_Kplas)
Plas2ECF = P_ECF*(K_plas - K_inter);

% ALD impact
gamma_al = A_dtKsec * C_al.^B_dtKsec;
lambda_al = A_cdKsec * C_al.^B_cdKsec;

% GI feedforward effect
if do_FF
    gamma_Kin = max(1, FF*(M_Kgut - MKgutSS));
else
    gamma_Kin = 1;
end

% muscle-kidney crosstalk
if MKX > 0
    omegaKic = max(0, (MKslope*(K_muscle - Kmuscle_baseline) + 1));
else
    omegaKic = 1;
end

% renal K handling

% if highK_eff > 0
%     if highK_eff == 1
%         GFR = (1 - 0.29) * 0.125;
%         eta_ptKreab = 0.36; % lower fractional PT reab
%     elseif highK_eff == 2
%         eta_ptKreab = 0.36; % lower fractional PT reab only
%     elseif highK_eff == 3
%         GFR = (1 - 0.29) * 0.125; % GFR change only
%     else
%         fprintf('What is this highK_eff? %i', highK_eff)
%     end
% end

% NOTE: should be able to make highK_eff happen by changing
% TGF instead

eta_psKreab_base = eta_ptKreab_base + eta_LoHKreab;
if TGF_eff == 1
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    GFR = GFR_base + alpha_TGF * (eta_psKreab - eta_psKreab_base);
elseif TGF_eff == 2 % GFR only
    eta_ptKreab = eta_ptKreab_base; % PT K reab is baseline value
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    GFR = GFR_base + alpha_TGF * (eta_psKreab - eta_psKreab_base);
elseif TGF_eff == 3 % PT only
    eta_psKreab = eta_ptKreab + eta_LOHKreab;
    GFR = GFR_base;
else
    eta_ptKreab = eta_ptKreab_base; % PT K reab is baseline value
    eta_psKreab = eta_ptKreab + eta_LoHKreab;
    GFR = GFR_base; 
end

filK = GFR*K_plas;

psKreab = eta_psKreab * filK;

% distal tubule
if MKX == 1
    eta_dtKsec = gamma_al * gamma_Kin * omegaKic;
else
    eta_dtKsec = gamma_al * gamma_Kin;
end
dtKsec = dtKsec_eq * eta_dtKsec;

% collecting duct
if MKX == 2
    eta_cdKsec = lambda_al * omegaKic;
else
    eta_cdKsec = lambda_al;
end
cdKsec = cdKsec_eq * eta_cdKsec;

if MKX == 3
    eta_cdKreab = omegaKic;
else
    eta_cdKreab = 1;
end
dtK= filK - psKreab + dtKsec; % flow from dt
cdKreab = dtK*A_cdKreab*eta_cdKreab;

UrineK   = dtK + cdKsec - cdKreab;

dydt(2) = Gut2plasma - Plas2ECF - UrineK;

%% Interstitial K (M_Kinter)
rho_al = (66.4 + 0.273*C_al)./89.6050;
% insulin
L = 100; x0 = 0.5381; k = 1.069;
ins_A = A_insulin; ins_B = 100*B_insulin;
temp = (ins_A.*(L./(1+exp(-k.*(log10(C_insulin)-log10(x0)))))+ ins_B)./100;
if do_insulin
    rho_insulin = max(1.0,temp);
    %disp(C_insulin)
    %disp(temp)
    %disp(rho_insulin)
else
    rho_insulin = 1;
end
eta_NKA = rho_insulin * rho_al;

Inter2Muscle = eta_NKA* ((Vmax * K_inter)/(Km + K_inter));
Muscle2Inter = P_muscle*(K_muscle - K_inter);

dydt(3) = Plas2ECF - Inter2Muscle + Muscle2Inter;

%% Intracellular K (M_Kmuscle)
dydt(4) = Inter2Muscle - Muscle2Inter;

end