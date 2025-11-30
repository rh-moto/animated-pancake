%-------------------------------------------------------------------------------------------------%
% File Name --- model_desc.mod
%
% README (written by Ryuichiro Hashimoto, Oct. 2021)
% This code section declares the model.
% The structure of this section is:
%   0) External functions
%   1) Exogenous variables
%   2) Parameters
%   3) Parameter values
%   4) Steady-State values
%   5) Endogenous variables
%   6) Model equations
%   7) Initial values
%
% This section uses the following files:
%   1) data/params_original.mat     - contains parameter values
%-------------------------------------------------------------------------------------------------%


%-------------------------------------------------------------------------------------------------%
                                        % Exogenous Variables
%-------------------------------------------------------------------------------------------------%
varexo

// shock variables
nu_aa                   % TFP level
nu_ad                   % Inv-specific technology level
nu_ea                   % TFP growth
nu_ed                   % Inv-specific technology growth
nu_nb                   % banking sector BS
nu_nc                   % firm sector BS
nu_r                    % monetary policy rule
nu_sb                   % bank-level productivity
nu_sc                   % firm-level productivity
nu_gc                   % government spending (ration w.r.t. private consumption)
nu_kc                   % investment adj. costs
nu_pc                   % price markup
nu_wc                   % wage markup
nu_ut                   % capital utilization
nu_d                    % discount rate
nu_pi                   % taylor-rule inflation
// shocks added
nu_fdr                  % flood shocks
nu_edr                  % earthquake shocks
nu_fdr1                 % anticipated flood shocks (+1)
//nu_fdr_p                % flood shocks (to public infrastructure)
//nu_edr_p                % earthquake shocks (to public infrastructure)
//nu_fdr_p1               % anticipated flood shocks (to public infrastructure)
;

%-------------------------------------------------------------------------------------------------%
                                        % Parameters
%-------------------------------------------------------------------------------------------------%
parameters

// Params in model equations
Share_Labor             % labor share (sum of HH + E + FI)
kappa_wc                % adj. costs in wage
Fc                      % Fixed costs in production
kappaC                  % adj. costs in capital stock
gammab                  % non-exit probability in banking sector
mu_ec                   % monitoring costs in FE contract
mu_b                    % monitoring costs in IF contract
gama                    % share of intermediate goods input
gammaec                 % non-exit probability in firm sector
vv                      % elasticity w.r.t. labor suply
ksi                     % weight of HH labor disutility
theta_wc                % elasticity w.r.t. labor input at steady-state
phi_c                   % coefficient on marginal costs (?)
theta_c                 % elasticity w.r.t. intermediary input at steady-state
kappa_pc                % adj. costs for price
rpi                     % taylor rule (inflation)
delta                   % depreciation rate (quarterly)
alpha                   % labor share (HH)
alpha_E                 % labor share (Entrepreneur)
alpha_FI                % labor share (Financial Intermediary)
beta                    % discount rate
A1                      % the power of TFP growth (Okazaki-Sudo p.49)
A2                      % the power of Inv-specific tech growth (see above)
kappa_uc                % adj. costs in capital utilization
phi_u                   % inverse of adj. costs in capital utilization
FcCg                    % ratio of production fixed costs to final goods production
rho_NSR                 % interest rate smoothing
habit                   % consumption habitat
theta_fdr               % elasticity of TFP w.r.t. flood
theta_edr               % elasticity of TFP w.r.t. earthquake
//theta_fdr_p               % elasticity of TFP w.r.t. flood (public)
//theta_edr_p               % elasticity of TFP w.r.t. earthquake (public)

// AR(1) coefficient in shock process
rho_aa                  % TFP level
rho_ad                  % Inv-specific technology level
rho_ea                  % TFP growth
rho_ed                  % Inv-specific technology growth
rho_nc                  % banking sector BS
rho_nb                  % firm sector BS
rho_r                   % monetary policy rule
rho_sb                  % bank-level productivity
rho_sc                  % firm-level productivity
rho_gc                  % government spending (ratio w.r.t. private consumption)
rho_kc                  % investment adj. costs
rho_pc                  % price markup
rho_wc                  % wage markup
rho_ut                  % capital utilization
rho_d                   % discount rate
rho_pi                  % taylor-rule inflation
rho_fdr                 % disaster shocks to TFP level
rho_edr                 % disaster shocks to TFP level
//rho_fdr_p                 % disaster shocks to TFP level (public)
//rho_edr_p                 % disaster shocks to TFP level (public)

// Std error in shock process
// Used in main.mod if you choose to give shock
sigma_aa                % TFP level
sigma_ad                % Inv-specific technology level
sigma_ea                % TFP growth
sigma_ed                % Inv-specific technology growth
sigma_nc                % banking sector BS
sigma_nb                % firm sector BS
sigma_r                 % monetary policy rule
sigma_sb                % bank-level productivity
sigma_sc                % firm-level productivity
sigma_gc                % government spending (ratio w.r.t. private consumption)
sigma_kc                % investment adj. costs
sigma_pc                % price markup
sigma_wc                % wage markup
sigma_ut                % capital utilization
sigma_d                 % discount rate
sigma_pi                % taylor-rule inflation

// Variables at Steady-State
gec_SS                  % ratio of government spending w.r.t. consumption
NSR_SS                  % nominal short rate
RSR_SS                  % real short-term rate
Pc_SS                   % price
eta_d_SS                % Inv-specific technology growth rate
eta_a_SS                % TFP growth rate
NEC_SS                  % firm net worth
NB_SS                   % bank net worth
C_SS                    % consumption
REC_SS                  % gross nominal rate of capital stock (utilization rate adjusted)
GRKC_SS                 % scaling param for adj. costs in capital utilization
QC_SS                   % Tobin Q
IC_SS                   % investment
KC_SS                   % capital stock
LEC_SS                  % lending to entrepreneur
omegabarUB_SS           % cut-off value of idiosyncratic shock in IF contracts
omegabarUEC_SS          % cut-off value of idiosyncratic shock in FE contracts
gcdfUB_SS               % G in IF contracts (see BGG, CMR, and/or Muto-Sudo-Yoneyama)
gcdfUEC_SS              % G in FE contracts
gammacdfUB_SS           % gamma in IF contracts
gammacdfUEC_SS          % gamma in FE contracts
d_gcdfUB_SS             % derivative of G w.r.t. omegabar in IF contracts
d_gcdfUEC_SS            % derivative of G w.r.t. omegabar in FE contracts
d_gammacdfUB_SS         % derivative of gamma w.r.t. omegabar in IF contracts
d_gammacdfUEC_SS        % derivative of gamma w.r.t. omegabar in FE contracts
sigmaUB_SS              % uncertainty in bank-level productivity shock
sigmaUEC_SS             % uncertainty in firm-level productivity shock
NECQKC_SS               % ratio of firm net worth to tobin Q*capital stock
NBQK_SS                 % ratio of bank net worth to tobin Q*capital stock
Cg_SS                   % final goods production
MCC_SS                  % nominal marginal cost
WC_SS                   % nominal wage for HH
WEC_SS                  % nominal wage for Entrepreneur
WFIC_SS                 % nominal wage for Financial Intermediary
HC_SS                   % nominal labor input
C_C_SS                  % intermediary input
Output_SS               % nominal output
Solow_SS                % nominal Solow residual
lambda_b_SS             % nominal HH Lagrange multiplier
RB_SS                   % return on financial intermediation
pi_c_SS                 % inflation at steady-state
UC_SS                   % capital utilization

// params for flood damage rate
sigma_fdr1              % sigma for anticipated flood damage rate (+1)
sigma_fdr               % sigma for damage rate
sigma_edr               % sigma for damage rate
//sigma_fdr_p1            % sigma for anticipated flood damage rate (to public infrastructure)
//sigma_fdr_p             % sigma for damage rate (to public infrastructure)
//sigma_edr_p             % sigma for damage rate (to public infrastructure)
sqrtfdr_avg             % average value of (fdr)^(1/2)
sqrtedr_avg             % average value of (edr)^(1/2)
;

%-------------------------------------------------------------------------------------------------%
                                        % Parameterization (Fixed, calibrated)
%-------------------------------------------------------------------------------------------------%
load(strcat(basedir, "/13_Data/params_original.mat"));

// Calibrated params
alpha           = 0.6000;
alpha_E         = 0.0200;
alpha_FI        = 0.0200;
gama            = 0.58333333;
theta_c         = 7;
theta_wc        = 7;
GRKC_SS         = parametersI(9); //the model somehow breaks when GRKC_SS is explicitly set at 0.05
ksi             = 0.2000;
beta            = 0.9980;
delta           = 0.0280;
gec_SS          = 0.2400;
omegabarUB_SS   = 0.73722330;
omegabarUEC_SS  = 0.40008466;
NECQKC_SS       = 0.6000;
NBQK_SS         = 0.1000;

if loadOriginal == 1;
    // Estimated params (either from params_original or previous result)
    vv              = 1.0;
    kappaC          = 0.5;
    kappa_pc        = 25;
    kappa_wc        = 25;
    rpi             = 2.75;
    phi_u           = 5;
    rho_NSR         = 0.5;
    sigmaUEC_SS     = XC(4);
    sigmaUB_SS      = XC(5);
    mu_ec           = XC(6);
    mu_b            = XC(7);
    gammaec         = XC(8);
    gammab          = XC(9);
    eta_d_SS        = 1.0;
    eta_a_SS        = 1.0;
    pi_c_SS         = 1;
    habit           = parametersI(17);

    // steady-state values (either from params_original or previous result)
    UC_SS           = 1;
    QC_SS           = 1;
    HC_SS           = var_ss1(1);
    C_SS            = var_ss1(2);
    C_C_SS          = var_ss1(3);
    KC_SS           = var_ss1(4);
    MCC_SS          = var_ss1(5);
    Pc_SS           = var_ss1(6);
    WC_SS           = var_ss1(7);
    WFIC_SS         = var_ss1(8);
    WEC_SS          = var_ss1(9);
    Cg_SS           = var_ss1(10);
end;

if loadPrevResult == 1;
    load(strcat(basedir, "/13_Data/est_params.csv"));
    load(strcat(basedir, "/13_Data/est_ss.csv"));

    // Estimated params (either from params_original or previous result)
    vv              = est_params(1);
    kappaC          = est_params(2);
    kappa_pc        = est_params(3);
    kappa_wc        = est_params(4);
    rpi             = est_params(5);
    phi_u           = est_params(6);
    rho_NSR         = est_params(7);
    sigmaUEC_SS     = est_params(8);
    sigmaUB_SS      = est_params(9);
    mu_ec           = est_params(10);
    mu_b            = est_params(11);
    gammaec         = est_params(12);
    gammab          = est_params(13);
    eta_a_SS        = est_params(14);
    eta_d_SS        = est_params(15);
    pi_c_SS         = est_params(16);
    habit           = est_params(17);

    // steady-state values (either from params_original or previous result)
    UC_SS           = 1;
    QC_SS           = 1;
    Pc_SS           = 1;
    HC_SS           = est_ss(1);
    C_SS            = est_ss(2);
    C_C_SS          = est_ss(3);
    Cg_SS           = est_ss(4);
    KC_SS           = est_ss(5);
    MCC_SS          = est_ss(6);
    WC_SS           = est_ss(7);
    WFIC_SS         = est_ss(8);
    WEC_SS          = est_ss(9);
end;

%-------------------------------------------------------------------------------------------------%
                                        % Steady-State Values
%-------------------------------------------------------------------------------------------------%
// steady state value of credit contracts
// auxiliary variables (to use standard normal distribution)
// log(omega) ~ N(-sigma^2/2, sigma^2) a la BGG [1999]
zUB_SS              = (log(omegabarUB_SS) + sigmaUB_SS^2/2)/sigmaUB_SS;
gcdfUB_SS           = normcdf(zUB_SS - sigmaUB_SS);
d_gcdfUB_SS         = normpdf(zUB_SS - sigmaUB_SS)/(omegabarUB_SS*sigmaUB_SS);
gammacdfUB_SS       = normcdf(zUB_SS - sigmaUB_SS) + omegabarUB_SS*(1-normcdf(zUB_SS));
d_gammacdfUB_SS     = d_gcdfUB_SS + (1-normcdf(zUB_SS)) - normpdf(zUB_SS)/sigmaUB_SS;

zUEC_SS             = (log(omegabarUEC_SS) + sigmaUEC_SS^2/2)/sigmaUEC_SS;
gcdfUEC_SS          = normcdf(zUEC_SS - sigmaUEC_SS);
d_gcdfUEC_SS        = normpdf(zUEC_SS - sigmaUEC_SS)/(omegabarUEC_SS*sigmaUEC_SS);
gammacdfUEC_SS      = normcdf(zUEC_SS - sigmaUEC_SS) + omegabarUEC_SS*(1-normcdf(zUEC_SS));
d_gammacdfUEC_SS    = d_gcdfUEC_SS + (1-normcdf(zUEC_SS)) - normpdf(zUEC_SS)/sigmaUEC_SS;

// constant
phi_c = (gama^-gama*(alpha*(1-gama))^(-alpha*(1-gama)))*
            (alpha_E*(1-gama))^(-alpha_E*(1-gama))*
            (alpha_FI*(1-gama))^(-alpha_FI*(1-gama))*
            ((1-alpha-alpha_E-alpha_FI)*(1-gama))^(-(1-alpha-alpha_E-alpha_FI)*(1-gama));

A1 = 1/(1-gama)/(alpha+alpha_E+alpha_FI);
A2 = (1-alpha-alpha_E-alpha_FI)/(alpha+alpha_E+alpha_FI);

// other variables
IC_SS           = delta*KC_SS;
Output_SS       = (1+gec_SS)*C_SS+IC_SS;
Share_Labor     = alpha+alpha_E+alpha_FI;
Solow_SS        = Output_SS/HC_SS^(Share_Labor)/(KC_SS)^(1-Share_Labor);
FcCg            = theta_c/(theta_c-1)-1;
Fc              = FcCg*Cg_SS;
kappa_uc        = GRKC_SS/Pc_SS/QC_SS;
NSR_SS          = (pi_c_SS*eta_a_SS^A1*eta_d_SS^A2)/beta;
RSR_SS          = (eta_a_SS^A1*eta_d_SS^A2)/beta;
REC_SS          = GRKC_SS/QC_SS/eta_d_SS+(1-delta)/eta_d_SS;
RB_SS           = 1/(gammacdfUB_SS-mu_b*gcdfUB_SS)*(1-NECQKC_SS-NBQK_SS)/(1-NECQKC_SS)*RSR_SS;
NEC_SS          = NECQKC_SS*(QC_SS*KC_SS);
NB_SS           = NBQK_SS*(QC_SS*KC_SS);
LEC_SS          = QC_SS*KC_SS - NEC_SS;
lambda_b_SS     = 1/Pc_SS/C_SS*(1-beta*habit)/(1-habit);

//shock process (AR1)
rho_aa = 0.9237;
rho_ad = 0.9161;
rho_ea = 0.1136;
rho_ed = 0.5928;
rho_nb = 0.3231;
rho_nc = 0.793;
rho_gc = 0.7723;
rho_kc = 0.3014;
rho_pc = 0.8551;
rho_wc = 0.9621;
rho_d = 0.6692;
rho_pi = 0.4287;
rho_ut = 0.9651;
rho_r           = 0;
rho_sc          = 0;
rho_sb          = 0;

//shock process (SD)
sigma_aa = 0.006;
sigma_ad = 0.0061;
sigma_ea = 0.0024;
sigma_ed = 0.0021;
sigma_r = 0.0017;
sigma_nb = 0.0029;
sigma_nc = 0.0039;
sigma_gc = 0.0096;
sigma_kc = 0.0294;
sigma_pc = 0.0278;
sigma_wc = 0.152;
sigma_d = 0.0042;
sigma_pi = 0.0036;
sigma_ut = 0.0045;
sigma_sc        = 0;
sigma_sb        = 0;

// flood damage rate
sigma_fdr       = 0.02818;
sigma_fdr1      = 0.02818;
//sigma_fdr_p     = 0.0007;
//sigma_fdr_p1    = 0.0007;
sigma_edr       = 0.05;
//sigma_edr_p     = 0.0010;
//sigma_edr_p1    = 0.0010;
rho_fdr = 0.4636;
rho_edr = 0.5562;
//rho_fdr_p = 0.3712;
//rho_edr_p = 0.2248;
theta_fdr = 0.5887;
//theta_fdr_p = 0.2416;
theta_edr = 0.405;
//theta_edr_p = 0.5017;
sqrtfdr_avg     = 0.00540;
sqrtedr_avg     = 0.00177;


%-------------------------------------------------------------------------------------------------%
                                        % Endogenous Variables
%-------------------------------------------------------------------------------------------------%
var
//
pi_c                % inflation
eta_a               % TFP growth rate
eta_d               % Inv-specific technology growth rate
HC                  % nominal labor input
RSR                 % real short-term rate
NSR                 % nominal short-term rate
RLR                 % real long-term rate
NLR                 % nominal long-term rate
REC                 % gross nominal rate of capital stock (utilization rate adjusted)
EX_REC              % expected gross nominal rate of capital stock (utilization rate adjusted)
NECQKC              % ratio of firm net worth to tobin Q*capital stock
NBQK                % ratio of bank net worth to tobin Q*capital stock
pi_d                % detrended inflation
gec                 % ratio of government spending w.r.t. consumption
omegabarUB          % cut-off value of idiosyncratic shock in IF contracts
omegabarUEC         % cut-off value of idiosyncratic shock in FE contracts
gcdfUB              % G in IF contracts (see BGG, CMR, and/or Muto-Sudo-Yoneyama)
gcdfUEC             % G in FE contracts
gammacdfUB          % gamma in IF contracts
gammacdfUEC         % gamma in FE contracts
d_gcdfUB            % derivative of G w.r.t. omegabar in IF contracts
d_gcdfUEC           % derivative of G w.r.t. omegabar in FE contracts
d_gammacdfUB        % derivative of gamma w.r.t. omegabar in IF contracts
d_gammacdfUEC       % derivative of gamma w.r.t. omegabar in FE contracts
sigmaUB             % uncertainty in bank-level productivity shock in IF contracts
sigmaUEC            % uncertainty in firm-level productivity shock in FE contracts
ZB                  % loan rate in IF contracts
ZEC                 % loan rate in FE contracts
RB                  % return on financial intermediation
UC                  % capital utilization
SP_B                % interest rate spread in IF contracts
SP_EC               % interest rate spread in FE contracts
SP_RE               % spread between REC and real short-term interest rate
SP_EX_RE            % spread between EX_REC and real short-term interest rate
RE_P                % spread between pi_c and pi_d
pi_star             % policy inflation target
LEC                 % lending to entrepreneur

// shock process
e_aa                % TFP level
e_ad                % Inv-specific technology level
e_nc                % banking sector BS
e_nb                % firm sector BS
e_r                 % monetary policy rule
e_kc                % investment adj. costs
e_pc                % price markup
e_wc                % wage markup
e_ut                % capital utilization
e_d                 % discount rate
e_fdr               % flood shocks on TFP level
e_edr               % earthquake shocks on TFP level
//e_fdr_p             % flood shocks on TFP level (public)
//e_edr_p             % earthquake shocks on TFP level (public)

// First log difference
l_C                 % consumption
l_C_C               % interemediary goods consumption
l_Cg                % final goods consumption
l_NB                % bank net worth
l_NEC               % firm net worth
l_IC                % investment
l_KC                % capital stock
l_GRKC              % scaling param for adj. costs in capital utilization
l_MCC               % nominal marginal cost
l_WC                % nominal wage
l_WEC               % nominal wage for HH
l_WFIC              % nominal wage for Entrepreneur
l_HC                % nominal wage for Financial Intermediary
l_QC                % Tobin Q
l_lambda_b          % nominal HH Lagrange multiplier
l_Output            % nominal output
l_Solow             % nominal Solow residualw
l_pi_c              % inflation
l_pi_d              % detrended inflation
l_RSR               % real short-term rate
l_NSR               % nominal short-term rate
l_RLR               % real long-term rate
l_NLR               % nominal long-term rate
l_UC                % capital utilization
l_REC               % gross nominal rate of capital stock (utilization rate adjusted)


// Climate-related variables (% of capital depreciated)
fdr                 % flood damage ratio
edr                 % earthquake damage ratio
//fdr_p               % flood damage ratio (public)
//edr_p               % earthquake damage ratio (public)
;

trend_var(growth_factor = pi_c) Pc;	    % detrended price
trend_var(growth_factor = eta_a) Za ;	% detrended TFP
trend_var(growth_factor = eta_d) Zd ;	% detrended Inv-specific

% deterented consumption, output, net worth
var(deflator = Za^A1*Zd^A2) C C_C Cg Output NB NEC;

% detrended Investment & Capital
var(deflator = Za^A1*Zd^(A2+1)) IC KC KKC;

%
var(deflator = Pc*Zd^(-1)) Pd GRKC;

%
var(deflator = Pc*Za^A1*Zd^A2) WC WEC WFIC;

%
var(deflator = Zd^(-1)) QC;

%
var(deflator = (Pc*Za^A1*Zd^A2)^(-1)) lambda_b;

%
var(deflator = Pc) MCC;

%
var(deflator = Za^(A1*Share_Labor)*Zd^(A2*Share_Labor-(1-Share_Labor))) Solow;


%-------------------------------------------------------------------------------------------------%
                                        % Model Equations
%-------------------------------------------------------------------------------------------------%
model;

// Households
    [name = 'Euler equation']
    1 = lambda_b(+1)/lambda_b*NSR*beta;

    [name = 'HH shadow price']
    lambda_b = exp(e_d)/(C-habit*C(-1))/Pc-beta*habit*exp(e_d(+1))/(C(+1)-habit*C)/Pc;


    [name = 'Nominal wage in intermediate goods sector']
    exp(e_d)*theta_wc*exp(e_wc)*ksi*(HC)^vv =
        - lambda_b*Pc*WC/Pc*(1-theta_wc*exp(e_wc))
        + lambda_b*Pc*WC/Pc*kappa_wc*(WC/WC(-1)-pi_c_SS*(eta_a_SS^A1*eta_d_SS^A2))*(WC/WC(-1))
        - beta*lambda_b(+1)*Pc(+1)*WC(+1)/Pc(+1)*kappa_wc*
            (WC(+1)/WC-pi_c_SS*(eta_a_SS^A1*eta_d_SS^A2))*HC(+1)/HC*WC(+1)/WC;

// Intermediate goods sector
    @#if withIns
        [name = 'Net capital return for intermediate goods sector']
        // followed the util adj. cost by Sugo & Ueda 2008JJIE P498(explanation on CEE2005)
        REC = (UC*GRKC/Pc-kappa_uc*((exp(e_ut)*UC)^(1+phi_u)-1)/(1+phi_u)
                + QC*(1-delta))/QC(-1);
    @#else
        [name = 'Net capital return for intermediate goods sector']
        // followed the util adj. cost by Sugo & Ueda 2008JJIE P498(explanation on CEE2005)
        REC = (UC*GRKC/Pc-kappa_uc*((exp(e_ut)*UC)^(1+phi_u)-1)/(1+phi_u)
                + QC*(1-delta))*(1-fdr-edr)/QC(-1);
    @#endif

    // this adj. cost variation is following Sugo & Ueda 2008JJIE P480.
    [name = 'FOC w.r.t. investment']
    QC*(1-kappaC*(IC*exp(e_kc)/IC(-1)-(eta_a_SS^A1*eta_d_SS^(A2+1)))^2/2-
            (IC*exp(e_kc)/IC(-1))*kappaC*(IC*exp(e_kc)/IC(-1)-(eta_a_SS^A1*eta_d_SS^(A2+1))))
        - 1/exp(e_ad)/Zd =
            - beta*lambda_b(+1)*Pc(+1)/lambda_b/Pc(+1)*QC(+1)*(IC(+1)/IC)^2*
                exp(e_kc(+1))*kappaC*(IC(+1)*exp(e_kc(+1))/IC-(eta_a_SS^A1*eta_d_SS^(A2+1)));

    [name = 'Capital accumulation']
    KC = (1-kappaC*(IC*exp(e_kc)/IC(-1)-(eta_a_SS^A1*eta_d_SS^(A2+1)))^2/2)*IC
            + (1-delta)*(1-fdr-edr)*KC(-1);

    [name = 'Entrepreneur profit maximization problem for intermediate goods']
    GRKC/Pc=kappa_uc*(exp(e_ut))^(phi_u+1)*(UC)^phi_u;

    [name = 'Intermediate goods producer demand function (intermeidate input)']
    gama*MCC*(Cg+Fc)=Pc*C_C;

    [name = 'Intermediate goods producer demand function (capital)']
    (1-gama)*(1-alpha-alpha_E-alpha_FI)*MCC*(Cg+Fc)=(UC)*GRKC*(1-fdr-edr)*KC(-1);

    [name = 'Intermediate goods producer demand function (household labor)']
    (1-gama)*alpha*MCC*(Cg+Fc)=WC*HC;

    [name = 'Intermediate goods producer demand function (entrepreneur labor)']
    (1-gama)*alpha_E*MCC*(Cg+Fc)=WEC;

    [name = 'Intermediate goods producer demand function (FI labor)']
    (1-gama)*alpha_FI*MCC*(Cg+Fc)=WFIC;

    @#if noTFP
        [name = 'Nominal marginal cost in intermediate goods sector']
        MCC = (exp(0*e_edr)*exp(0*e_fdr))/
                (exp(e_aa))/Za*phi_c*Pc^gama*
                (WEC^alpha_E*WFIC^alpha_FI*WC^alpha*GRKC^(1-alpha-alpha_E-alpha_FI))^(1-gama);
    @#else
        [name = 'Nominal marginal cost in intermediate goods sector']
        MCC = (exp(theta_edr*e_edr)*exp(theta_fdr*e_fdr))/
                (exp(e_aa))/Za*phi_c*Pc^gama*
                (WEC^alpha_E*WFIC^alpha_FI*WC^alpha*GRKC^(1-alpha-alpha_E-alpha_FI))^(1-gama);
    @#endif

    [name = 'Price dynamics for intermediate goods goods sector']
    (1-theta_c*exp(e_pc)) =
        - theta_c*exp(e_pc)*MCC/Pc
        + kappa_pc*(Pc/Pc(-1)-pi_c_SS)*Pc/Pc(-1)
        - lambda_b(+1)*Pc(+1)/lambda_b/Pc*beta*kappa_pc*(Pc(+1)/Pc-pi_c_SS)*Pc(+1)/Pc*Cg(+1)/Cg;


// Chained credit contracts
    // G distribution
    // auxiliary variables (to use standard normal distribution)
    // log(omega) ~ N(-sigma^2/2, sigma^2) a la BGG [1999]
    # zUB   = (log(omegabarUB) + sigmaUB(-1)^2/2)/sigmaUB(-1);
    # zUEC  = (log(omegabarUEC) + sigmaUEC(-1)^2/2)/sigmaUEC(-1);

    [name = 'G in IF contracts']
    gcdfUB = normcdf(zUB - sigmaUB(-1));

    [name = 'G in FE contracts']
    gcdfUEC = normcdf(zUEC - sigmaUEC(-1));

    [name = 'Derivative of G w.r.t. omegabarUB in IF contracts']
    d_gcdfUB = normpdf(zUB - sigmaUB(-1))/(omegabarUB*sigmaUB(-1));

    [name = 'Derivative of G w.r.t. omegabarUEC in FE contracts']
    d_gcdfUEC = normpdf(zUEC - sigmaUEC(-1))/(omegabarUEC*sigmaUEC(-1));

    [name = 'Gamma in IF contracts']
    gammacdfUB = gcdfUB +omegabarUB*(1-normcdf(zUB));

    [name = 'Gamma in FE contracts']
    gammacdfUEC = gcdfUEC +omegabarUEC*(1-normcdf(zUEC));

    [name = 'Derivative of Gamma w.r.t. omegabarUB in IF contracts']
    d_gammacdfUB = d_gcdfUB + (1-normcdf(zUB)) - normpdf(zUB)/sigmaUB(-1);

    [name = 'Derivative of Gamma w.r.t. omegabarUEC in FE contracts']
    d_gammacdfUEC = d_gcdfUEC + (1-normcdf(zUEC)) - normpdf(zUEC)/sigmaUEC(-1);

    [name = 'Zero profit condition for IF contract']
    (gammacdfUB -mu_b*gcdfUB)*(gammacdfUEC -mu_ec*gcdfUEC)*REC = RSR(-1)*(1-NECQKC(-1)-NBQK(-1));

    [name = 'Entrepreners return lower limit in FE contracts']
    1 - gammacdfUEC = NECQKC(-1);

    [name = 'Bank profit maximization in FE contracts']
    0 = (1-gammacdfUB(+1))*(gammacdfUEC(+1)-mu_ec*gcdfUEC(+1))*REC(+1)
        + d_gammacdfUB(+1)/(d_gammacdfUB(+1)-mu_b*d_gcdfUB(+1))*
            ((gammacdfUB(+1)-mu_b*gcdfUB(+1))*(gammacdfUEC(+1)-mu_ec*gcdfUEC(+1))*REC(+1))
        - d_gammacdfUB(+1)/(d_gammacdfUB(+1)-mu_b*d_gcdfUB(+1))*RSR
        + (1-gammacdfUB(+1))*(d_gammacdfUEC(+1)-mu_ec*d_gcdfUEC(+1))*
            (1-gammacdfUEC(+1))/d_gammacdfUEC(+1)*REC(+1)
        + d_gammacdfUB(+1)*
            (gammacdfUB(+1)-mu_b*gcdfUB(+1))*(d_gammacdfUEC(+1)-mu_ec*d_gcdfUEC(+1))*
            (1-gammacdfUEC(+1))*REC(+1)/(d_gammacdfUB(+1)-mu_b*d_gcdfUB(+1))/d_gammacdfUEC(+1);

    [name = 'Bank net worth']
    NBQK*QC*KC/(QC(-1)*KC(-1)) =
        gammab*(1-gammacdfUB)*(gammacdfUEC -mu_ec*gcdfUEC)*(REC)
        + (e_nb*exp(Output_SS))/(QC(-1)*KC(-1))
        + alpha_FI/(1-alpha-alpha_E-alpha_FI)*(UC)*GRKC/QC(-1)/Pc;

    [name = 'Entrepreners net worth']
    NECQKC*(QC*KC/(QC(-1)*KC(-1))) =
        gammaec*(1-gammacdfUEC)*(REC)
        + (e_nc*exp(Output_SS))/(QC(-1)*KC(-1))
        + alpha_E/(1-alpha-alpha_E-alpha_FI)*(UC)*GRKC/QC(-1)/Pc;

    [name = 'Lending to Entrepreneurs']
    LEC = QC*KC - NEC;


// Market clearing, monetary policy rule, etc.
    [name = 'Fisher equation']
    RSR = NSR/pi_c(+1);

    [name = 'Market clearing condition']
    Cg = (1+gec)*C
            + C_C
            + IC/exp(e_ad)/Zd
            + kappa_uc*((exp(e_ut)*UC)^(1+phi_u)-1)/(1+phi_u)*(1-fdr-edr)*KC(-1)
            + mu_ec*gcdfUEC*(REC)*QC(-1)*KC(-1)
	        + mu_b*gcdfUB*(gammacdfUEC -mu_ec*gcdfUEC)*(REC)*QC(-1)*KC(-1)
            + (1-gammaec)*(1-gammacdfUEC)*(REC)*QC(-1)*KC(-1)
            + (1-gammab)*(1-gammacdfUB)*(gammacdfUEC -mu_ec*gcdfUEC)*(REC)*QC(-1)*KC(-1);

    [name = 'Output']
    Output = (1+gec)*C+IC/exp(e_ad)/Zd;

    [name = 'Solow residuals']
    //Solow = Output/(((1-fdr-edr)*KC(-1)*UC)^(1-Share_Labor)*HC^(Share_Labor));
    Solow = Output/((1-fdr-edr)*KC(-1)*UC)^(1-Share_Labor)/(HC^Share_Labor);

    [name = 'Monetary policy rule']
    NSR = NSR(-1)^rho_NSR*((pi_c/pi_star)^(rpi*(1-rho_NSR)))*
            ((pi_c_SS*eta_a_SS^A1*eta_d_SS^A2)/beta)^(1-rho_NSR)*exp(e_r);

    [name = 'Shock on Taylor rule']
    log(pi_star/pi_c_SS) = rho_pi*log(pi_star(-1)/pi_c_SS) + nu_pi;

    [name = 'Inflation']
    pi_d=Pd/Pd(-1);

    [name = 'Detrended price']
    Pd=Pc/exp(e_ad)/Zd;


// Interest & loan rate
    [name = 'Long term real interest rate']
    l_RLR = (1/40)*(
    @#for i in 0:39
        +l_RSR(@{i})
    @#endfor
    );

    [name = 'Long term real interest rate (exp)']
    RLR = exp(l_RLR);

    [name = 'Long term nominal interest rate']
    l_NLR = (1/40)*(
    @#for i in 0:39
        +l_NSR(@{i})
    @#endfor
    );

    [name = 'Long term nominal interest rate (exp)']
    NLR = exp(l_NLR);

    [name = 'Loan rate in FE contracts']
    ZEC = 1/(1-NECQKC)*omegabarUEC(+1)*(REC(+1));

    [name = 'Loan rate in IF contracts']
    ZB = (1-NECQKC)/(1-NECQKC-NBQK)*omegabarUB(+1)*(RB(+1));

    [name = 'Return on financial intermediation']
    RB/RSR(-1) = 1/(gammacdfUB-mu_b*gcdfUB)*(1-NBQK(-1)-NECQKC(-1))/(1-NECQKC(-1));


// Shock process
    [name = 'Shock process (TFP level)']
    e_aa = rho_aa*e_aa(-1) + nu_aa;

    [name = 'Shock process (Inv-specific technology level)']
    e_ad = rho_ad*e_ad(-1) + nu_ad;

    @#if bankLoss
        [name = 'Shock process (FI sector BS)']
        e_nb = rho_nb*e_nb(-1)+ nu_nb - QC_SS*KC_SS*(fdr+edr)*0.25;
    @#else
        [name = 'Shock process (FI sector BS)']
        e_nb = rho_nb*e_nb(-1)+ nu_nb;
    @#endif

    [name = 'Shock process (Firm sector BS)']
    e_nc = rho_nc*e_nc(-1)+ nu_nc;

    [name = 'Shock process (Investment adj. costs)']
    e_kc = rho_kc*e_kc(-1)+ nu_kc;

    [name = 'Shock process (Price markup)']
    e_pc = rho_pc*e_pc(-1)+ nu_pc;

    [name = 'Shock process (Wage markup)']
    e_wc = rho_wc*e_wc(-1)+ nu_wc;

    [name = 'Shock process (Capital utilization)']
    e_ut = rho_ut*e_ut(-1)+ nu_ut;

    [name = 'Shock process (Discount rate)']
    e_d  = rho_d*e_d(-1)+ nu_d;

    [name = 'Shock process (Monetary policy rule)']
    e_r  = rho_r*e_r(-1) + nu_r;

    [name = 'Shock process (Sigma in FE contracts)']
    log(sigmaUEC/sigmaUEC_SS)   = rho_sc*log(sigmaUEC(-1)/sigmaUEC_SS) + nu_sc;

    [name = 'Shock process (Sigma in IF contracts)']
    log(sigmaUB/sigmaUB_SS)     = rho_sb*log(sigmaUB(-1)/sigmaUB_SS) + nu_sb;

    [name = 'Shock process (Inv-specific techology growth rate)']
    log(eta_d/eta_d_SS)         = rho_ed*log(eta_d(-1)/eta_d_SS) + nu_ed;

    [name = 'Shock process (TFP growth rate)']
    log(eta_a/eta_a_SS)         = rho_ea*log(eta_a(-1)/eta_a_SS) + nu_ea;

    [name = 'Shock process (Govt spending share)']
    log(gec/gec_SS)             = rho_gc*log(gec(-1)/gec_SS) + nu_gc;


// Auxiliary Variables
    [name = 'Auxiliary variables (Capital pre-period)']
    KKC = (1-delta)*(1-fdr-edr)*KC(-1);

    [name = 'Auxiliary variables (NB)']
    NB          = NBQK*(QC*KC);

    [name = 'Auxiliary variables (NEC)']
    NEC         = NECQKC*(QC*KC);

    [name = 'Auxiliary variables (C growth)']
    l_C         = log(C/C(-1));

    [name = 'Auxiliary variables (C_C growth)']
    l_C_C       = log(C_C/C_C(-1));

    [name = 'Auxiliary variables (Cg growth)']
    l_Cg        = log(Cg/Cg(-1));

    [name = 'Auxiliary variables (NB growth)']
    l_NB        = log(NB/NB(-1));

    [name = 'Auxiliary variables (NEC growth)']
    l_NEC       = log(NEC/NEC(-1));

    [name = 'Auxiliary variables (REC)']
    l_REC       = log(REC);

    [name = 'Auxiliary variables (IC growth)']
    l_IC        = log(IC/IC(-1));

    [name = 'Auxiliary variables (KC growth)']
    l_KC        = log((1-fdr(+1)-edr(+1))*KC/((1-fdr-edr)*KC(-1)));

    [name = 'Auxiliary variables (GRKC growth)']
    l_GRKC      = log(GRKC/GRKC(-1));

    [name = 'Auxiliary variables (NSR)']
    l_NSR       = log(NSR);

    [name = 'Auxiliary variables (RSR)']
    l_RSR       = log(RSR);

    [name = 'Auxiliary variables (MCC growth)']
    l_MCC       = log(MCC/MCC(-1));

    [name = 'Auxiliary variables (Household wage growth)']
    l_WC        = log(WC/WC(-1));

    [name = 'Auxiliary variables (Entrepreneur wage growth)']
    l_WEC       = log(WEC/WEC(-1));

    [name = 'Auxiliary variables (FI wage growth)']
    l_WFIC      = log(WFIC/WFIC(-1));

    [name = 'Auxiliary variables (Labor input growth)']
    l_HC        = log(HC/HC(-1));

    [name = 'Auxiliary variables (QC growth)']
    l_QC        = log(QC/QC(-1));

    [name = 'Auxiliary variables (Inflation)']
    l_pi_c      = log(pi_c);

    [name = 'Auxiliary variables (Detrended inflation)']
    l_pi_d      = log(pi_d);

    [name = 'Auxiliary variables (Capital utilization growth)']
    l_UC        = log(UC/UC(-1));

    [name = 'Auxiliary variables (Household shadow price growth)']
    l_lambda_b  = log(lambda_b/lambda_b(-1));

    [name = 'Auxiliary variables (Output growth)']
    l_Output    = log(Output/Output(-1));

    [name = 'Auxiliary variables (Solow residual growth)']
    l_Solow     = log(Solow/Solow(-1));


// Macro variables (growth, level, and spread)
    [name = 'Interest rate spread in IF contracts']
    SP_B        = ZB - RSR;

    [name = 'Interest rate spread in FE contracts']
    SP_EC       = ZEC - ZB;

    [name = 'Interest rate spread between capital return and risk-free rate']
    SP_RE       = REC - RSR;

    [name = 'Inflation rate spread between detrended and non-adjusted']
    RE_P        = l_pi_d - l_pi_c;

    [name = 'Expected return on capital']
    EX_REC      = REC(+1);

    [name = 'Expected interest rate spread between capital return and risk-free rate']
    SP_EX_RE    = EX_REC - RSR;


// Climate-related components
// flood damage ratio is set by quarterly shocks
    [name = 'Private capital depreciation rate due to flood']
    fdr         = (sqrtfdr_avg + nu_fdr + nu_fdr1(-1))^2;

    [name = 'Private capital depreciation rate due to earthquake']
    edr         = (sqrtedr_avg + nu_edr)^2;

    //[name = 'Public capital depreciation rate due to flood']
    //fdr_p       = nu_fdr_p + nu_fdr_p1(-1);

    //[name = 'Public capital depreciation rate due to earthquake']
    //edr_p       = nu_edr_p;

    [name = 'Flood shock in TFP level']
    e_fdr       = rho_fdr*e_fdr(-1) + fdr;

    //[name = 'Flood shock in TFP level (public)']
    //e_fdr_p     = rho_fdr_p*e_fdr_p(-1) + fdr_p;

    [name = 'Earthquake shock in TFP level']
    e_edr       = rho_edr*e_edr(-1) + edr;

    //[name = 'Earthquake shock in TFP level (public)']
    //e_edr_p     = rho_edr_p*e_edr_p(-1) + edr_p;

end;


%-------------------------------------------------------------------------------------------------%
                                        % Initial Values
%-------------------------------------------------------------------------------------------------%
initval;
    eta_a           = eta_a_SS;
    eta_d           = eta_d_SS;
    C               = C_SS;
    RSR             = RSR_SS;
    NSR             = NSR_SS;
    RLR             = RSR_SS;
    NLR             = NSR_SS;
    Pd              = 1;
    GRKC            = GRKC_SS;
    REC             = REC_SS;
    QC              = QC_SS;
    IC              = IC_SS;
    KC              = KC_SS;
    KKC             = (1-delta)*KC_SS;
    gcdfUB          = gcdfUB_SS;
    gcdfUEC         = gcdfUEC_SS;
    omegabarUB      = omegabarUB_SS;
    omegabarUEC     = omegabarUEC_SS;
    sigmaUB         = sigmaUB_SS;
    sigmaUEC        = sigmaUEC_SS;
    d_gcdfUB        = d_gcdfUB_SS;
    d_gcdfUEC       = d_gcdfUEC_SS;
    gammacdfUB      = gammacdfUB_SS;
    gammacdfUEC     = gammacdfUEC_SS;
    d_gammacdfUB    = d_gammacdfUB_SS;
    d_gammacdfUEC   = d_gammacdfUEC_SS;
    NECQKC          = NECQKC_SS;
    ZEC             = 1/(1-NECQKC_SS)*omegabarUEC_SS*REC_SS;
    ZB              = (1-NECQKC_SS)/(1-NECQKC_SS-NBQK_SS)*omegabarUB_SS*(RB_SS);
    RB              = 1/(gammacdfUB_SS-mu_b*gcdfUB_SS)*(1-NECQKC_SS-NBQK_SS)/(1-NECQKC_SS)*RSR_SS;
    NBQK            = NBQK_SS;
    Cg              = Cg_SS;
    MCC             = MCC_SS;
    NB              = NB_SS;
    NEC             = NEC_SS;
    LEC             = LEC_SS;
    gec             = gec_SS;
    WC              = WC_SS;
    WEC             = WEC_SS;
    WFIC            = WFIC_SS;
    HC              = HC_SS;
    UC              = UC_SS;
    C_C             = C_C_SS;
    Output          = Output_SS;
    Solow           = Solow_SS;
    lambda_b        = lambda_b_SS;
    pi_c            = pi_c_SS;
    pi_d            = pi_c_SS/eta_d_SS;
    pi_star         = pi_c_SS;
    SP_EC           = 1/(1-NECQKC_SS)*omegabarUEC_SS*REC_SS
                        - (1-NECQKC_SS)/(1-NECQKC_SS-NBQK_SS)*omegabarUB_SS*(RB_SS);
    SP_B            = (1-NECQKC_SS)/(1-NECQKC_SS-NBQK_SS)*omegabarUB_SS*(RB_SS) - RSR_SS;
    EX_REC          = REC_SS;
    SP_RE           = REC_SS - RSR_SS;
    SP_EX_RE        = REC_SS - RSR_SS;
    RE_P            = log(pi_c_SS/eta_d_SS) - log(pi_c_SS);

    l_pi_c          = log(pi_c_SS);
    l_MCC           = log(pi_c_SS);
    l_WC            = log(pi_c_SS*eta_a_SS^A1*eta_d_SS^A2);
    l_WEC           = log(pi_c_SS*eta_a_SS^A1*eta_d_SS^A2);
    l_WFIC          = log(pi_c_SS*eta_a_SS^A1*eta_d_SS^A2);
    l_RLR           = log(RSR_SS);
    l_NLR           = log(NSR_SS);
    l_pi_d          = log(pi_c_SS/eta_d_SS);
    l_NSR           = log(NSR_SS);
    l_RSR           = log(RSR_SS);
    l_REC           = log(REC_SS);
    l_UC            = log(UC_SS);
    l_C             = log(eta_a_SS^A1*eta_d_SS^A2);
    l_C_C           = log(eta_a_SS^A1*eta_d_SS^A2);
    l_Cg            = log(eta_a_SS^A1*eta_d_SS^A2);
    l_NB            = log(eta_a_SS^A1*eta_d_SS^A2);
    l_NEC           = log(eta_a_SS^A1*eta_d_SS^A2);
    l_IC            = log(eta_a_SS^A1*eta_d_SS^(A2+1));
    l_KC            = log(eta_a_SS^A1*eta_d_SS^(A2+1));
    l_GRKC          = log(pi_c_SS*eta_d_SS^-1);
    l_QC            = log(1/eta_d_SS);
    l_lambda_b      = log((pi_c_SS*eta_a_SS^A1*eta_d_SS^A2)^-1);
    l_Output        = log(eta_a_SS^A1*eta_d_SS^A2);
    l_Solow         = log(eta_a_SS^(A1*Share_Labor)*eta_d_SS^(A2*Share_Labor-(1-Share_Labor)));

    e_aa            = 0;
    e_ad            = 0;
    e_r             = 0;
    e_nb            = 0;
    e_nc            = 0;
    e_kc            = 0;
    e_ut            = 0;
    e_pc            = 0;
    e_wc            = 0;
    e_d             = 0;
    e_fdr           = 0;
    e_edr           = 0;
    //e_fdr_p         = 0;
    //e_edr_p         = 0;
    e_d             = 0;

    nu_aa           = 0;
    nu_ad           = 0;
    nu_r            = 0;
    nu_sb           = 0;
    nu_sc           = 0;
    nu_nb           = 0;
    nu_nc           = 0;
    nu_ea           = 0;
    nu_ed           = 0;
    nu_gc           = 0;
    nu_d            = 0;
    nu_kc           = 0;
    nu_pc           = 0;
    nu_wc           = 0;
end;
