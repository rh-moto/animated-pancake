//2017/07/03 Okazaki Changed Price,Wage,Monetary Policy�̎���Pc_SS��pi_c_SS�֕ύX
//2017/07/03re Okazaki Changed �ڕW�C���t���V���b�N��ǉ�

//���܂��Ȃ� // 2016/09/21 OKazaki Added
//opengl('save','software')

//2017/09/12 Okazaki Changed
// Exogenous variables
//varexo nu_aa, nu_ad, nu_r,  nu_sb,  nu_sc,  nu_nb,  nu_nc,  nu_ea, nu_ed, nu_gc, nu_kc, nu_pc, nu_wc, nu_ut, nu_d, nu_l, nu_pop, nu_pi;
//varexo nu_aa, nu_ad, nu_r,  nu_sb,  nu_sc,  nu_nb,  nu_nc,  nu_ea, nu_ed, nu_gc, nu_kc, nu_pc, nu_wc, nu_ut, nu_d, nu_l, nu_pop, nu_pi, nu_1, nu_2, nu_3, nu_4;
//varexo nu_aa, nu_ad, nu_r,  nu_sb,  nu_sc,  nu_nb,  nu_nc,  nu_ea, nu_ed, nu_gc, nu_kc, nu_pc, nu_wc, nu_ut, nu_d, nu_pop, nu_pi, nu_1, nu_2, nu_3, nu_4;
varexo nu_aa, nu_ad, nu_r,  nu_sb,  nu_sc,  nu_nb,  nu_nc,  nu_ea, nu_ed, nu_gc, nu_kc, nu_pc, nu_wc, nu_ut, nu_d, nu_pop, nu_pi, nu_1, nu_2, nu_3, nu_4, nu_5, nu_6, nu_7, nu_8, nu_9, nu_10, nu_11, nu_12;


// parameters
parameters Share_Labor kappa_wc NSR_SS Pc_SS sigmaUB_SS sigmaUEC_SS Fc kappaC gammab mu_ec mu_b gama gammaec vv ksi theta_wc phi_c theta_c kappa_pc rpi delta alpha alpha_E alpha_FI beta A1 A2  NEC_SS kappa_uc phi_u Fc_Cg gec_SS eta_d_SS eta_a_SS;

parameters rho_aa rho_ad rho_ea rho_ed rho_nc rho_nb rho_gc rho_kc rho_sc rho_sb rho_r rho_pc rho_wc rho_ut rho_NSR;

parameters sigma_aa sigma_nc  sigma_nb  sigma_sc  sigma_sb  sigma_r  sigma_ea  sigma_ed sigma_kc sigma_gc  sigma_pc  sigma_wc  sigma_ut;

parameters C_SS REC_SS GRKC_SS QC_SS IC_SS KC_SS gcdfUB_SS gcdfUEC_SS omegabarUB_SS omegabarUEC_SS d_gcdfUB_SS d_gcdfUEC_SS gammacdfUB_SS gammacdfUEC_SS d_gammacdfUB_SS d_gammacdfUEC_SS RSR_SS NECQKC_SS NBQK_SS CG_SS MCC_SS WC_SS WEC_SS WFIC_SS HC_SS C_C_SS Output_SS Solow_SS lambda_b_SS;

parameters C_FL_SS  REC_FL_SS  GRKC_FL_SS QC_FL_SS IC_FL_SS KC_FL_SS gcdfUB_FL_SS gcdfUEC_FL_SS omegabarUB_FL_SS omegabarUEC_FL_SS d_gcdfUB_FL_SS d_gcdfUEC_FL_SS gammacdfUB_FL_SS gammacdfUEC_FL_SS d_gammacdfUB_FL_SS d_gammacdfUEC_FL_SS RSR_FL_SS NECQKC_FL_SS NBQK_FL_SS CG_FL_SS MCC_FL_SS WC_FL_SS WEC_FL_SS  WFIC_FL_SS HC_FL_SS C_C_FL_SS Output_FL_SS  Solow_FL_SS lambda_b_FL_SS NEC_FL_SS;

parameters MNC_SS MNC_FL_SS;
//2017/06/05 Okazaki Changed
parameters rho_d sigma_d;
//parameters sigma_d;
parameters habit;
//parameters RLQ_SS RLQ_FL_SS NLQ_SS rho_l sigma_l phi_l_SS;
parameters eta_wc eta_pc;
parameters rho_pop POP_SS;
parameters NB_SS NB_FL_SS RB_SS RB_FL_SS;

//2017/07/03 Okazaki Added
parameters pi_c_SS;

//2017/07/03re Okazaki Added
parameters rho_pi sigma_pi;

//2017/08/09 Okazaki Added
parameters UC_SS UC_FL_SS;

//2017/09/12 Okazaki Added
parameters sigma_1 sigma_2 sigma_3 sigma_4;

//2017/11/08 Okazaki Added
parameters sigma_5 sigma_6 sigma_7 sigma_8 sigma_9 sigma_10 sigma_11 sigma_12;

//2017/09/21 Okazaki Added
parameters RLR_SS RLR_FL_SS;

load data_for_dynare0330;

theta_c  = parametersI(1);
CPI_SS   = parametersI(2);
delta    = parametersI(3);
alpha    = parametersI(4);
alpha_E  = parametersI(5);
alpha_FI = parametersI(6);
gama     = parametersI(7);
theta_wc = parametersI(8);
GRKC_SS  = parametersI(9);
beta     = parametersI(10);
ksi      = parametersI(11);
//2017/07/03 Okazaki Changed
//vv       = 1.5;
vv = 1.0;
gec_SS   = parametersI(13);
habit    = parametersI(17);

//phi_l_SS = 0;
kappaC   = 0.5;
kappa_pc = 25;
kappa_wc = 25;
phi_u    = 5; // a bit big

eta_wc   = 0;
eta_pc   = 0;

phi_c=gama^-gama*(alpha*(1-gama))^(-alpha*(1-gama));
phi_c=phi_c*(alpha_E*(1-gama))^(-alpha_E*(1-gama));
phi_c=phi_c*(alpha_FI*(1-gama))^(-alpha_FI*(1-gama));
phi_c=phi_c*((1-alpha-alpha_E-alpha_FI)*(1-gama))^(-(1-alpha-alpha_E-alpha_FI)*(1-gama));

A1 = 1/(1-gama)/(alpha+alpha_E+alpha_FI);
A2 = (1-alpha-alpha_E-alpha_FI)/(alpha+alpha_E+alpha_FI);


//calculated ss

eta_d_SS=1.0;
eta_a_SS=1.0;

HC_SS=var_ss1(1);
C_SS=var_ss1(2);
C_C_SS=var_ss1(3);
KC_SS=var_ss1(4);
MCC_SS=var_ss1(5);
Pc_SS=var_ss1(6);
WC_SS=var_ss1(7);
WFIC_SS=var_ss1(8);
WEC_SS=var_ss1(9);
CG_SS=var_ss1(10);
IC_SS=delta*KC_SS;
Output_SS=(1+gec_SS)*C_SS+IC_SS; // deleted fic on 2.23.15
Share_Labor=alpha+alpha_E+alpha_FI;
Solow_SS=Output_SS/HC_SS^(Share_Labor)/KC_SS^(1-Share_Labor);

lambda_b_SS=1/Pc_SS/C_SS*(1-beta*habit)/(1-habit);

NECQKC_SS=parametersFI(5);
NBQK_SS=parametersFI(6);

omegabarUEC_SS = XC(2);
omegabarUB_SS = XC(3);
sigmaUEC_SS = XC(4);
sigmaUB_SS = XC(5);
mu_ec = XC(6);
mu_b = XC(7);
gammaec = XC(8);
gammab = XC(9);

//steady state value of other varaibles
gcdfUB_SS = normcdf(((log(omegabarUB_SS)-sigmaUB_SS^2/2)/sigmaUB_SS));
gcdfUEC_SS = normcdf(((log(omegabarUEC_SS)-sigmaUEC_SS^2/2)/sigmaUEC_SS));

d_gcdfUB_SS = 1/sigmaUB_SS/omegabarUB_SS*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUB_SS)-sigmaUB_SS^2/2)/sigmaUB_SS)^2);

d_gcdfUEC_SS = 1/sigmaUEC_SS/omegabarUEC_SS*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUEC_SS)-sigmaUEC_SS^2/2)/sigmaUEC_SS)^2);

gammacdfUB_SS = normcdf(((log(omegabarUB_SS)-sigmaUB_SS^2/2)/sigmaUB_SS))+omegabarUB_SS*(1-normcdf(((log(omegabarUB_SS)+sigmaUB_SS^2/2)/sigmaUB_SS)));

gammacdfUEC_SS = normcdf(((log(omegabarUEC_SS)-sigmaUEC_SS^2/2)/sigmaUEC_SS))+omegabarUEC_SS*(1-normcdf(((log(omegabarUEC_SS)+sigmaUEC_SS^2/2)/sigmaUEC_SS)));

d_gammacdfUB_SS = 1/sigmaUB_SS/omegabarUB_SS*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUB_SS)-sigmaUB_SS^2/2)/sigmaUB_SS)^2)+ (1-normcdf(((log(omegabarUB_SS)+sigmaUB_SS^2/2)/sigmaUB_SS)))- 1/sigmaUB_SS*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUB_SS)+ sigmaUB_SS^2/2)/sigmaUB_SS)^2);

d_gammacdfUEC_SS = 1/sigmaUEC_SS/omegabarUEC_SS*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUEC_SS)-sigmaUEC_SS^2/2)/sigmaUEC_SS)^2) + (1-normcdf(((log(omegabarUEC_SS)+sigmaUEC_SS^2/2)/sigmaUEC_SS)))- 1/sigmaUEC_SS*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUEC_SS)+ sigmaUEC_SS^2/2)/sigmaUEC_SS)^2);

//2017/08/09 Okazaki Added
UC_SS = 1;
UC_FL_SS = UC_SS;

QC_SS=1;
Fc_Cg = theta_c/(theta_c-1)-1;
Fc=Fc_Cg*CG_SS;

POP_SS = 1.0;

//2017/07/03 Okazaki Added
//2017/09/06 Okazaki Changed
//2017/09/07 Okazaki Changed
//pi_c_SS = 1.01^0.25; //�N���łP��
pi_c_SS = 1.0; //�N���O��
//pi_c_SS = 1.0015683218; //�`1998�̃T���v���ł̐��v�l

//NLQ_SS=0;
//RLQ_SS=NLQ_SS;

kappa_uc = GRKC_SS/Pc_SS/QC_SS;

// monetary policy
//2017/09/06 Okazaki Changed
//rho_NSR=0.8087     ;
//rpi=2.1240     ;
rho_NSR=0.5;
rpi=2.75;
//rho_NSR=0.6441; //�`1998�̃T���v���ł̐��v�l
//rpi=2.4913; //�`1998�̃T���v���ł̐��v�l


//2017/07/03 Okazaki Changed
//2017/08/08 Okaaki Changed
//NSR_SS=pi_c_SS/beta - NLQ_SS;
//RSR_SS=1/beta - RLQ_SS;
//NSR_SS=(pi_c_SS*eta_a_SS^A1*eta_d_SS^A2)/beta - NLQ_SS;
//RSR_SS=(eta_a_SS^A1*eta_d_SS^A2)/beta - RLQ_SS;
NSR_SS=(pi_c_SS*eta_a_SS^A1*eta_d_SS^A2)/beta;
RSR_SS=(eta_a_SS^A1*eta_d_SS^A2)/beta;

//2017/09/21 Okazaki Added
//RLR_SS = 1/beta - RLQ_SS;
//RLR_FL_SS = 1/beta - RLQ_SS;
RLR_SS = RSR_SS;
RLR_FL_SS = RLR_SS;

//2017/08/08 Okazaki Changed
//REC_SS=(GRKC_SS+(1-delta)*QC_SS)/QC_SS;
REC_SS=GRKC_SS/QC_SS/eta_d_SS+(1-delta)/eta_d_SS;

//RB_SS = 1/(gammacdfUB_SS-mu_b*gcdfUB_SS)*(1-NECQKC_SS-NBQK_SS)/(1-NECQKC_SS)*(RSR_SS + RLQ_SS);
RB_SS = 1/(gammacdfUB_SS-mu_b*gcdfUB_SS)*(1-NECQKC_SS-NBQK_SS)/(1-NECQKC_SS)*RSR_SS;

//shock process
rho_aa = 0;
rho_ad = 0;
rho_r  = 0;

//rho_sc = 0.8;
//rho_sb = 0.8;
//rho_nb = 0.25;
//rho_nc = 0.34;
//rho_ea = 0.14;
//rho_ed = 0.23;
//rho_gc = 0.80;
//rho_kc = 0.92;
//rho_pc = 0.47;
//rho_wc = 0.43;
rho_sc = 0;
rho_sb = 0;
rho_nb = 0;
rho_nc = 0;
rho_ea = 0;
rho_ed = 0;
rho_gc = 0;
rho_kc = 0;
rho_pc = 0;
rho_wc = 0;
rho_ut = 0;
//2017/06/05 Okazaki Deleted
rho_d  = 0;
//rho_l  = 0.5;
//rho_l  = 0;
rho_pop = 0;
rho_pi = 0;

sigma_aa=0.79;
sigma_ad=0.79;
sigma_nc=0.07;
sigma_nb=0.07;
//sigma_sc=0.001;
//sigma_sb=.001;
sigma_sc=0;
sigma_sb=0;
sigma_r=.01;
sigma_ea=0.06;
sigma_ed=0.06;
sigma_gc=0.14;
sigma_kc=0.54;
sigma_pc=0.11;
sigma_wc=0.39;
sigma_ut=0;
sigma_d = 0.05;
//sigma_d = 0.01;
//sigma_d = 0;
//sigma_l = 0.05;
sigma_pi = 0.01;

//2017/09/12 OKazaki Added
sigma_1 = 0.01;
sigma_2 = 0.01;
sigma_3 = 0.01;
sigma_4 = 0.01;

//2017/11/08 Okazaki Added
sigma_5 = 0.01;
sigma_6 = 0.01;
sigma_7 = 0.01;
sigma_8 = 0.01;
sigma_9 = 0.01;
sigma_10 = 0.01;
sigma_11 = 0.01;
sigma_12 = 0.01;


NEC_SS=parametersFI(5)*(QC_SS*KC_SS);
NB_SS=parametersFI(6)*(QC_SS*KC_SS);

MNC_SS=mu_ec*gcdfUEC_SS*(REC_SS)*QC_SS*KC_SS+mu_b*gcdfUB_SS*((gammacdfUEC_SS -mu_ec*gcdfUEC_SS)*(REC_SS)*QC_SS*KC_SS);


C_FL_SS=C_SS;
REC_FL_SS=REC_SS;
GRKC_FL_SS=GRKC_SS;
QC_FL_SS=QC_SS;
IC_FL_SS=IC_SS;
KC_FL_SS=KC_SS;
gcdfUB_FL_SS=gcdfUB_SS;
gcdfUEC_FL_SS=gcdfUEC_SS;
omegabarUB_FL_SS=omegabarUB_SS;
omegabarUEC_FL_SS=omegabarUEC_SS;
d_gcdfUB_FL_SS=d_gcdfUB_SS;
d_gcdfUEC_FL_SS=d_gcdfUEC_SS;
gammacdfUB_FL_SS=gammacdfUB_SS;
gammacdfUEC_FL_SS=gammacdfUEC_SS;
d_gammacdfUB_FL_SS=d_gammacdfUB_SS;
d_gammacdfUEC_FL_SS=d_gammacdfUEC_SS;
RSR_FL_SS=RSR_SS;
NECQKC_FL_SS=NECQKC_SS;
NBQK_FL_SS=NBQK_SS;
CG_FL_SS=CG_SS;
MCC_FL_SS=MCC_SS;
WC_FL_SS=WC_SS;
WEC_FL_SS=WEC_SS;
WFIC_FL_SS=WFIC_SS;
HC_FL_SS=HC_SS;
C_C_FL_SS=C_C_SS;
Output_FL_SS=Output_SS;
Solow_FL_SS=Solow_SS;
lambda_b_FL_SS=lambda_b_SS;
NEC_FL_SS=NEC_SS;
MNC_FL_SS=MNC_SS;
NB_FL_SS=NB_SS;
RB_FL_SS=RB_SS;
//RLQ_FL_SS=NLQ_SS;


// Endogenous variab1les
var pi_c eta_a eta_d;
trend_var(growth_factor = pi_c) Pc;
trend_var(growth_factor = eta_a) Za ;
trend_var(growth_factor = eta_d) Zd ;

var POP l_POP;

var(deflator = Za^A1*Zd^A2) C C_C Cg Output NB NEC C_FL C_C_FL Cg_FL Output_FL NB_FL NEC_FL;

var(deflator = Za^A1*Zd^(A2+1)) IC KC IC_FL KC_FL;

var(deflator = Pc*Zd^(-1)) Pd GRKC;
var(deflator = Zd^(-1)) GRKC_FL;

var(deflator = Pc*Za^A1*Zd^A2) WC WEC WFIC;
var(deflator = Za^A1*Zd^A2) WC_FL WEC_FL WFIC_FL;

var(deflator = Zd^(-1)) QC QC_FL;

var(deflator = (Pc*Za^A1*Zd^A2)^(-1)) lambda_b;
var(deflator = (Za^A1*Zd^A2)^(-1)) lambda_b_FL;

var(deflator = Pc) MCC;
var MCC_FL;

var(deflator = Za^(A1*Share_Labor)*Zd^(A2*Share_Labor-(1-Share_Labor))) Solow Solow_FL;

var NSR RSR HC REC NECQKC NBQK pi_d gec RSR_FL HC_FL REC_FL NECQKC_FL NBQK_FL;

var gcdfUB gcdfUEC omegabarUB omegabarUEC sigmaUB sigmaUEC d_gcdfUB d_gcdfUEC gammacdfUB gammacdfUEC d_gammacdfUB d_gammacdfUEC gcdfUB_FL gcdfUEC_FL omegabarUB_FL omegabarUEC_FL d_gcdfUB_FL d_gcdfUEC_FL gammacdfUB_FL gammacdfUEC_FL d_gammacdfUB_FL d_gammacdfUEC_FL;

var e_aa e_ad e_r e_nb e_nc e_kc e_pc e_wc e_ut e_d;

var l_C l_C_C l_Cg l_NB l_NEC l_IC l_KC l_GRKC l_MCC l_WC l_WEC l_WFIC l_HC l_QC l_lambda_b l_Output l_Solow l_pi_c l_pi_d l_NSR l_RSR l_UC l_C_FL l_C_C_FL l_Cg_FL l_NB_FL l_NEC_FL l_IC_FL l_KC_FL l_GRKC_FL l_MCC_FL l_WC_FL l_WEC_FL l_WFIC_FL l_HC_FL l_QC_FL l_lambda_b_FL l_Output_FL l_Solow_FL l_RSR_FL l_UC_FL;

var ZEC ZB RB ZEC_FL ZB_FL RB_FL;
var mnc_cg l_REC mnc_cg_FL l_REC_FL;
var UC UC_FL;

//var NLQ RLQ RLQ_FL;
//var(deflator = (Za^A1*Zd^A2)^(-1)) phi_l;


var L_Output L_C L_IC L_HC L_KC L_NB L_NEC;
var L_Output_FL L_C_FL L_IC_FL L_HC_FL L_KC_FL L_NB_FL L_NEC_FL;

var SP_B SP_EC SP_B_FL SP_EC_FL;
var EX_REC EX_REC_FL SP_RE SP_RE_FL SP_EX_RE SP_EX_RE_FL RE_P;

var GAP RGAP;

//2017/07/03re Okazaki Added
var pi_star;

//2017/09/12 Okazaki Added
var e_1 e_2 e_3 e_4;

//2017/09/15 Okazaki Added
var EXP_R1 EXP_R2 EXP_R3 EXP_R4;

//2017/11/08 Okazaki Added
var EXP_R5 EXP_R6 EXP_R7 EXP_R8 EXP_R9 EXP_R10 EXP_R11 EXP_R12 e_5 e_6 e_7 e_8 e_9 e_10 e_11 e_12;

//2017/09/21 Okazaki Added
var RLR RLR_FL l_RLR l_RLR_FL;

//2017/10/19 Okazaki Added
var NLR l_NLR;

//2017/11/09 Okazaki Added
var RSR_MA RSR_FL_MA l_RSR_MA l_RSR_FL_MA;
var N1R N5R l_N1R l_N5R;
var R1R R5R R1R_FL R5R_FL l_R1R l_R1R_FL l_R5R l_R5R_FL;

model;

////////////////////////////////////////////////////////////
////////////Sticky price - wage economy
////////////////////////////////////////////////////////////


// Eular equation
//1 = lambda_b(+1)/lambda_b*(NSR+NLQ)*beta; //(#1)
1 = lambda_b(+1)/lambda_b*NSR*beta; //(#1)

lambda_b = exp(e_d)/(C-habit*C(-1))/Pc-beta*habit*exp(e_d(+1))*POP(+1)/(C(+1)-habit*C)/Pc; //(#2)


// net capital return for non-durable sector
REC = (UC*GRKC/Pc-kappa_uc*((exp(e_ut)*UC)^(1+phi_u)-1)/(1+phi_u)+QC*(1-delta))/QC(-1); //(#3)
// followed the util adj. cost by Sugo & Ueda 2008JJIE P498(explanation on CEE2005)


//investment for non-durable sector
//first order condition of capital producer that provides capital goods used for non-durable sector
//2017/08/08 Okazaki Changed
//QC*(1-kappaC*(IC*exp(e_kc)*POP/IC(-1)-1)^2/2-(IC*exp(e_kc)*POP/IC(-1))*kappaC*(IC*exp(e_kc)*POP/IC(-1)-1))-1/exp(e_ad)/Zd = -beta*lambda_b(+1)*Pc(+1)/lambda_b/Pc*QC(+1)*(IC(+1)*POP(+1)/IC)^2*exp(e_kc(+1))*kappaC*(IC(+1)*exp(e_kc(+1))*POP(+1)/IC-1); //(#)
QC*(1-kappaC*(IC*exp(e_kc)*POP/IC(-1)-(eta_a_SS^A1*eta_d_SS^(A2+1)))^2/2-(IC*exp(e_kc)*POP/IC(-1))*kappaC*(IC*exp(e_kc)*POP/IC(-1)-(eta_a_SS^A1*eta_d_SS^(A2+1))))-1/exp(e_ad)/Zd = -beta*lambda_b(+1)*Pc(+1)/lambda_b/Pc*QC(+1)*(IC(+1)*POP(+1)/IC)^2*exp(e_kc(+1))*kappaC*(IC(+1)*exp(e_kc(+1))*POP(+1)/IC-(eta_a_SS^A1*eta_d_SS^(A2+1))); //(#4)

// this adj. cost variation is following Sugo & Ueda 2008JJIE P480.


//capital accumulation in non-durable sector
//2017/08/08 Okaaki Changed
//KC = (1-kappaC*(IC*exp(e_kc)*POP/IC(-1)-1)^2/2)*IC + (1-delta)*KC(-1)/POP; //(#)
KC = (1-kappaC*(IC*exp(e_kc)*POP/IC(-1)-(eta_a_SS^A1*eta_d_SS^(A2+1)))^2/2)*IC + (1-delta)*KC(-1)/POP; //(#5)

// G distribution
gcdfUB = normcdf(((log(omegabarUB)-sigmaUB(-1)^2/2)/sigmaUB(-1))); //(#6)

gcdfUEC = normcdf(((log(omegabarUEC)-sigmaUEC(-1)^2/2)/sigmaUEC(-1))); //(#7)

d_gcdfUB = 1/sigmaUB(-1)/omegabarUB*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUB)-sigmaUB(-1)^2/2)/sigmaUB(-1))^2); //(#8)

d_gcdfUEC = 1/sigmaUEC(-1)/omegabarUEC*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUEC)-sigmaUEC(-1)^2/2)/sigmaUEC(-1))^2); //(#9)


// gamma distribution
gammacdfUB = gcdfUB +omegabarUB*(1-normcdf(((log(omegabarUB)+sigmaUB(-1)^2/2)/sigmaUB(-1)))); //(#10)

gammacdfUEC = gcdfUEC +omegabarUEC*(1-normcdf(((log(omegabarUEC)+sigmaUEC(-1)^2/2)/sigmaUEC(-1)))); //(#11)

d_gammacdfUB = 1/sigmaUB(-1)/omegabarUB*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUB)-sigmaUB(-1)^2/2)/sigmaUB(-1))^2)+(1-normcdf(((log(omegabarUB)+sigmaUB(-1)^2/2)/sigmaUB(-1)))) - 1/sigmaUB(-1)*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUB)+sigmaUB(-1)^2/2)/sigmaUB(-1))^2); //(#12)

d_gammacdfUEC = 1/sigmaUEC(-1)/omegabarUEC*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUEC)-sigmaUEC(-1)^2/2)/sigmaUEC(-1))^2) + (1-normcdf(((log(omegabarUEC)+sigmaUEC(-1)^2/2)/sigmaUEC(-1)))) - 1/sigmaUEC(-1)*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUEC)+sigmaUEC(-1)^2/2)/sigmaUEC(-1))^2); //(#13)


//zero profit condition (state by state) for FI contract
//(gammacdfUB -mu_b*gcdfUB)*(gammacdfUEC -mu_ec*gcdfUEC)*REC = (RSR(-1)+RLQ(-1))*(1-NECQKC(-1)-NBQK(-1)); //(#14)
(gammacdfUB -mu_b*gcdfUB)*(gammacdfUEC -mu_ec*gcdfUEC)*REC = RSR(-1)*(1-NECQKC(-1)-NBQK(-1)); //(#14)

//entrepreners' return lower limit in non-durable sector
1 - gammacdfUEC = NECQKC(-1); //(#15)

//bank's profit maximization1 in non-durable sector
//0=(1-gammacdfUB(+1))*(gammacdfUEC(+1)-mu_ec*gcdfUEC(+1))*REC(+1) + d_gammacdfUB(+1)/(d_gammacdfUB(+1)-mu_b*d_gcdfUB(+1))*((gammacdfUB(+1)-mu_b*gcdfUB(+1))*(gammacdfUEC(+1)-mu_ec*gcdfUEC(+1))*REC(+1)) - d_gammacdfUB(+1)/(d_gammacdfUB(+1)-mu_b*d_gcdfUB(+1))*(RSR+RLQ)+(1-gammacdfUB(+1))*(d_gammacdfUEC(+1)-mu_ec*d_gcdfUEC(+1))*(1-gammacdfUEC(+1))/d_gammacdfUEC(+1)*REC(+1)+d_gammacdfUB(+1)*(gammacdfUB(+1)-mu_b*gcdfUB(+1))*(d_gammacdfUEC(+1)-mu_ec*d_gcdfUEC(+1))*(1-gammacdfUEC(+1))*REC(+1)/(d_gammacdfUB(+1)-mu_b*d_gcdfUB(+1))/d_gammacdfUEC(+1); //(#16)
0=(1-gammacdfUB(+1))*(gammacdfUEC(+1)-mu_ec*gcdfUEC(+1))*REC(+1) + d_gammacdfUB(+1)/(d_gammacdfUB(+1)-mu_b*d_gcdfUB(+1))*((gammacdfUB(+1)-mu_b*gcdfUB(+1))*(gammacdfUEC(+1)-mu_ec*gcdfUEC(+1))*REC(+1)) - d_gammacdfUB(+1)/(d_gammacdfUB(+1)-mu_b*d_gcdfUB(+1))*RSR+(1-gammacdfUB(+1))*(d_gammacdfUEC(+1)-mu_ec*d_gcdfUEC(+1))*(1-gammacdfUEC(+1))/d_gammacdfUEC(+1)*REC(+1)+d_gammacdfUB(+1)*(gammacdfUB(+1)-mu_b*gcdfUB(+1))*(d_gammacdfUEC(+1)-mu_ec*d_gcdfUEC(+1))*(1-gammacdfUEC(+1))*REC(+1)/(d_gammacdfUB(+1)-mu_b*d_gcdfUB(+1))/d_gammacdfUEC(+1); //(#16)

//net worth of bank
NBQK*QC*KC/(QC(-1)*KC(-1))*POP = gammab*(1-gammacdfUB)*(gammacdfUEC -mu_ec*gcdfUEC)*(REC) + (e_nb*exp(Output_SS))/(QC(-1)*KC(-1))+alpha_FI/(1-alpha-alpha_E-alpha_FI)*(UC)*GRKC/QC(-1)/Pc; //(#17)

//net worth of entrepreneur in non-durable sector
NECQKC*(QC*KC/QC(-1)/KC(-1))*POP = gammaec*(1-gammacdfUEC)*(REC)+(e_nc*exp(Output_SS))/QC(-1)/KC(-1)+alpha_E/(1-alpha-alpha_E-alpha_FI)*(UC)*GRKC/QC(-1)/Pc; //(#18)

// Fisher equation
RSR = NSR/pi_c(+1); //(#19)


//market clearing condition for non-durables
Cg=(1+gec)*C+C_C+IC/exp(e_ad)/Zd+mu_ec*gcdfUEC*(REC)*QC(-1)*KC(-1)/POP+mu_b*gcdfUB*(gammacdfUEC -mu_ec*gcdfUEC)*(REC)*QC(-1)*KC(-1)/POP+kappa_uc*((exp(e_ut)*UC)^(1+phi_u)-1)/(1+phi_u)*KC(-1)/POP+(1-gammaec)*(1-gammacdfUEC)*(REC)*QC(-1)*KC(-1)/POP+(1-gammab)*(1-gammacdfUB)*(gammacdfUEC -mu_ec*gcdfUEC)*(REC)*QC(-1)*KC(-1)/POP; //(#20)

//Aggregate Variables (output, inflation and Solow residual)
Output = (1+gec)*C+IC/exp(e_ad)/Zd; //(#21)


pi_d=Pd/Pd(-1); //(#22)

Solow=Output/(KC(-1)^(1-Share_Labor)*HC^(Share_Labor))*POP^(1-Share_Labor); //(#23)


//production function for non-durable sector
gama*MCC*(Cg+Fc)=Pc*C_C; //(#24)
(1-gama)*(1-alpha-alpha_E-alpha_FI)*MCC*(Cg+Fc)=(UC)*GRKC*KC(-1)/POP; //(#25)
(1-gama)*alpha*MCC*(Cg+Fc)=WC*HC; //(#26)
(1-gama)*alpha_E*MCC*(Cg+Fc)=WEC; //(#27)
(1-gama)*alpha_FI*MCC*(Cg+Fc)=WFIC; //(#28)


//nominal wage in non-durable sector
//2017/07/03 Okazaki Changed
//2017/08/09 Okazaki Changed
//exp(e_d)*theta_wc*exp(e_wc)*ksi*(HC)^vv=-lambda_b*Pc*WC/Pc*(1-theta_wc*exp(e_wc))+lambda_b*Pc*WC/Pc*kappa_wc*(WC/WC(-1)-eta_wc*WC(-1)/WC(-2)-(1-eta_wc)*Pc_SS)*(WC/WC(-1))-beta*POP(+1)*lambda_b(+1)*Pc(+1)*WC(+1)/Pc(+1)*kappa_wc*(WC(+1)/WC-eta_wc*WC/WC(-1)-(1-eta_wc)*Pc_SS)*HC(+1)/HC*WC(+1)/WC; //(#)
//exp(e_d)*theta_wc*exp(e_wc)*ksi*(HC)^vv=-lambda_b*Pc*WC/Pc*(1-theta_wc*exp(e_wc))+lambda_b*Pc*WC/Pc*kappa_wc*(WC/WC(-1)-eta_wc*WC(-1)/WC(-2)-(1-eta_wc)*pi_c_SS)*(WC/WC(-1))-beta*POP(+1)*lambda_b(+1)*Pc(+1)*WC(+1)/Pc(+1)*kappa_wc*(WC(+1)/WC-eta_wc*WC/WC(-1)-(1-eta_wc)*pi_c_SS)*HC(+1)/HC*WC(+1)/WC; //(#)
exp(e_d)*theta_wc*exp(e_wc)*ksi*(HC)^vv=-lambda_b*Pc*WC/Pc*(1-theta_wc*exp(e_wc))+lambda_b*Pc*WC/Pc*kappa_wc*(WC/WC(-1)-eta_wc*WC(-1)/WC(-2)-(1-eta_wc)*pi_c_SS*(eta_a_SS^A1*eta_d_SS^A2))*(WC/WC(-1))-beta*POP(+1)*lambda_b(+1)*Pc(+1)*WC(+1)/Pc(+1)*kappa_wc*(WC(+1)/WC-eta_wc*WC/WC(-1)-(1-eta_wc)*pi_c_SS*(eta_a_SS^A1*eta_d_SS^A2))*HC(+1)/HC*WC(+1)/WC; //(#29)

//nominal marginal cost in non-durable sector
MCC=1/exp(e_aa)/Za*phi_c*Pc^gama*(WEC^alpha_E*WFIC^alpha_FI*WC^alpha*GRKC^(1-alpha-alpha_E-alpha_FI))^(1-gama); //(#30)


// Price Dynamics for non-durable goods sector
//2017/07/03 Okazaki Changed
//(1-theta_c*exp(e_pc))=-theta_c*exp(e_pc)*MCC/Pc+kappa_pc*(Pc/Pc(-1)-eta_pc*Pc(-1)/Pc(-2)-(1-eta_pc)*Pc_SS)*Pc/Pc(-1)-lambda_b(+1)*Pc(+1)/lambda_b/Pc*beta*POP(+1)*kappa_pc*(Pc(+1)/Pc-eta_pc*Pc/Pc(-1)-(1-eta_pc)*Pc_SS)*Pc(+1)/Pc*Cg(+1)/Cg; //(#)
(1-theta_c*exp(e_pc))=-theta_c*exp(e_pc)*MCC/Pc+kappa_pc*(Pc/Pc(-1)-eta_pc*Pc(-1)/Pc(-2)-(1-eta_pc)*pi_c_SS)*Pc/Pc(-1)-lambda_b(+1)*Pc(+1)/lambda_b/Pc*beta*POP(+1)*kappa_pc*(Pc(+1)/Pc-eta_pc*Pc/Pc(-1)-(1-eta_pc)*pi_c_SS)*Pc(+1)/Pc*Cg(+1)/Cg; //(#31)

//Monetary Policy Rule
//2017/07/03 Okazaki Changed
//NSR = NSR(-1)^rho_NSR*((pi_c/Pc_SS)^(rpi*(1-rho_NSR)))*NSR_SS^(1-rho_NSR)*exp(e_r); //(#)

//2017/07/03re Okazaki Changed & Added
//NSR = NSR(-1)^rho_NSR*((pi_c/pi_c_SS)^(rpi*(1-rho_NSR)))*NSR_SS^(1-rho_NSR)*exp(e_r); //(#)

//2017/08/09 Okazaki Changed
//NSR = NSR(-1)^rho_NSR*((pi_c/pi_star)^(rpi*(1-rho_NSR)))*NSR_SS^(1-rho_NSR)*exp(e_r); //(#)
//NSR = NSR(-1)^rho_NSR*((pi_c/pi_star)^(rpi*(1-rho_NSR)))*((pi_c_SS*eta_a_SS^A1*eta_d_SS^A2)/beta - NLQ_SS)^(1-rho_NSR)*exp(e_r); //(#32)
NSR = NSR(-1)^rho_NSR*((pi_c/pi_star)^(rpi*(1-rho_NSR)))*((pi_c_SS*eta_a_SS^A1*eta_d_SS^A2)/beta)^(1-rho_NSR)*exp(e_r); //(#32)

log(pi_star/pi_c_SS) = rho_pi*log(pi_star(-1)/pi_c_SS) + nu_pi; //(#33)


//2017/09/21 OKazaki Added
l_RLR = (l_RSR + l_RSR(+1) + l_RSR(+2)+ l_RSR(+3)+ l_RSR(+4)+ l_RSR(+5)+ l_RSR(+6)+ l_RSR(+7)+ l_RSR(+8)+ l_RSR(+9)+ l_RSR(+10)+ l_RSR(+11)+ l_RSR(+12)+ l_RSR(+13)+ l_RSR(+14)+ l_RSR(+15)+ l_RSR(+16)+ l_RSR(+17)+ l_RSR(+18)+ l_RSR(+19)+ l_RSR(+20)+ l_RSR(+21)+ l_RSR(+22)+ l_RSR(+23)+ l_RSR(+24)+ l_RSR(+25)+ l_RSR(+26)+ l_RSR(+27)+ l_RSR(+28)+ l_RSR(+29)+ l_RSR(+30)+ l_RSR(+31)+ l_RSR(+32)+ l_RSR(+33)+ l_RSR(+34)+ l_RSR(+35)+ l_RSR(+36)+ l_RSR(+37)+ l_RSR(+38)+ l_RSR(+39))/40;
RLR = exp(l_RLR);


//2017/10/19 Okazaki Added
l_NLR =  (l_NSR + l_NSR(+1) + l_NSR(+2)+ l_NSR(+3)+ l_NSR(+4)+ l_NSR(+5)+ l_NSR(+6)+ l_NSR(+7)+ l_NSR(+8)+ l_NSR(+9)+ l_NSR(+10)+ l_NSR(+11)+ l_NSR(+12)+ l_NSR(+13)+ l_NSR(+14)+ l_NSR(+15)+ l_NSR(+16)+ l_NSR(+17)+ l_NSR(+18)+ l_NSR(+19)+ l_NSR(+20)+ l_NSR(+21)+ l_NSR(+22)+ l_NSR(+23)+ l_NSR(+24)+ l_NSR(+25)+ l_NSR(+26)+ l_NSR(+27)+ l_NSR(+28)+ l_NSR(+29)+ l_NSR(+30)+ l_NSR(+31)+ l_NSR(+32)+ l_NSR(+33)+ l_NSR(+34)+ l_NSR(+35)+ l_NSR(+36)+ l_NSR(+37)+ l_NSR(+38)+ l_NSR(+39))/40;
NLR = exp(l_NLR);

//2017/11/09 Okazaki Added
l_RSR_MA = (l_RSR(-2) + l_RSR(-1) + l_RSR + l_RSR(+1) + l_RSR(+2))/5;
RSR_MA = exp(l_RSR_MA);
l_RSR_FL_MA = (l_RSR_FL(-2) + l_RSR_FL(-1) + l_RSR_FL + l_RSR_FL(+1) + l_RSR_FL(+2))/5;
RSR_FL_MA = exp(l_RSR_FL_MA);

l_N1R =  (l_NSR + l_NSR(+1) + l_NSR(+2)+ l_NSR(+3))/4;
N1R = exp(l_N1R);
l_N5R =  (l_NSR + l_NSR(+1) + l_NSR(+2)+ l_NSR(+3)+ l_NSR(+4)+ l_NSR(+5)+ l_NSR(+6)+ l_NSR(+7)+ l_NSR(+8)+ l_NSR(+9)+ l_NSR(+10)+ l_NSR(+11)+ l_NSR(+12)+ l_NSR(+13)+ l_NSR(+14)+ l_NSR(+15)+ l_NSR(+16)+ l_NSR(+17)+ l_NSR(+18)+ l_NSR(+19))/20;
N5R = exp(l_N5R);

l_R1R =  (l_RSR + l_RSR(+1) + l_RSR(+2)+ l_RSR(+3))/4;
R1R = exp(l_R1R);
l_R5R =  (l_RSR + l_RSR(+1) + l_RSR(+2)+ l_RSR(+3)+ l_RSR(+4)+ l_RSR(+5)+ l_RSR(+6)+ l_RSR(+7)+ l_RSR(+8)+ l_RSR(+9)+ l_RSR(+10)+ l_RSR(+11)+ l_RSR(+12)+ l_RSR(+13)+ l_RSR(+14)+ l_RSR(+15)+ l_RSR(+16)+ l_RSR(+17)+ l_RSR(+18)+ l_RSR(+19))/20;
R5R = exp(l_R5R);

l_R1R_FL =  (l_RSR_FL + l_RSR_FL(+1) + l_RSR_FL(+2)+ l_RSR_FL(+3))/4;
R1R_FL = exp(l_R1R_FL);
l_R5R_FL =  (l_RSR_FL + l_RSR_FL(+1) + l_RSR_FL(+2)+ l_RSR_FL(+3)+ l_RSR_FL(+4)+ l_RSR_FL(+5)+ l_RSR_FL(+6)+ l_RSR_FL(+7)+ l_RSR_FL(+8)+ l_RSR_FL(+9)+ l_RSR_FL(+10)+ l_RSR_FL(+11)+ l_RSR_FL(+12)+ l_RSR_FL(+13)+ l_RSR_FL(+14)+ l_RSR_FL(+15)+ l_RSR_FL(+16)+ l_RSR_FL(+17)+ l_RSR_FL(+18)+ l_RSR_FL(+19))/20;
R5R_FL = exp(l_R5R_FL);



ZEC = 1/(1-NECQKC)*omegabarUEC(+1)*(REC(+1)); //(#34)
ZB = (1-NECQKC)/(1-NECQKC-NBQK)*omegabarUB(+1)*(RB(+1)); //(#35)
//RB/(RSR(-1)+RLQ(-1)) = 1/(gammacdfUB-mu_b*gcdfUB)*(1-NBQK(-1)-NECQKC(-1))/(1-NECQKC(-1)); //(#36)
RB/RSR(-1) = 1/(gammacdfUB-mu_b*gcdfUB)*(1-NBQK(-1)-NECQKC(-1))/(1-NECQKC(-1)); //(#36)
Pd=Pc/exp(e_ad)/Zd; //(#37)


//shocks
e_aa = rho_aa*e_aa(-1) + nu_aa; //(#38)
e_ad = rho_ad*e_ad(-1) + nu_ad; //(#39)

//2017/09/12 Okazaki Changed(ZLB)
//e_r = rho_r*e_r(-1) + nu_r; //(#)
e_r = rho_r*e_r(-1) + nu_r + e_1(-1); //(#40)
e_1 = e_2(-1) + nu_1;
e_2 = e_3(-1) + nu_2;
e_3 = e_4(-1) + nu_3;

//2017/11/08 Okazaki Changed
//e_4 = nu_4;
e_4 = e_5(-1) + nu_4;

//2017/09/15 Okazaki Added
EXP_R1 = l_NSR(+1);
EXP_R2 = l_NSR(+2);
EXP_R3 = l_NSR(+3);
EXP_R4 = l_NSR(+4);

//2017/11/08 Okazaki Added
e_5 = e_6(-1) + nu_5;
e_6 = e_7(-1) + nu_6;
e_7 = e_8(-1) + nu_7;
e_8 = e_9(-1) + nu_8;
e_9 = e_10(-1) + nu_9;
e_10 = e_11(-1) + nu_10;
e_11 = e_12(-1) + nu_11;
e_12 = nu_12;

EXP_R5 = l_NSR(+5);
EXP_R6 = l_NSR(+6);
EXP_R7 = l_NSR(+7);
EXP_R8 = l_NSR(+8);
EXP_R9 = l_NSR(+9);
EXP_R10 = l_NSR(+10);
EXP_R11 = l_NSR(+11);
EXP_R12 = l_NSR(+12);

log(sigmaUEC/sigmaUEC_SS) = rho_sc*log(sigmaUEC(-1)/sigmaUEC_SS) + nu_sc; //(#41)
log(sigmaUB/sigmaUB_SS) = rho_sb*log(sigmaUB(-1)/sigmaUB_SS) + nu_sb; //(#42)
e_nb = rho_nb*e_nb(-1)+ nu_nb; //(#43)
e_nc = rho_nc*e_nc(-1)+ nu_nc; //(#44)
log(eta_d/eta_d_SS) = rho_ed*log(eta_d(-1)/eta_d_SS) + nu_ed; //(#45)
log(eta_a/eta_a_SS) = rho_ea*log(eta_a(-1)/eta_a_SS) + nu_ea; //(#46)
log(gec/gec_SS) = rho_gc*log(gec(-1)/gec_SS) + nu_gc; //(#47)
e_kc = rho_kc*e_kc(-1)+ nu_kc; //(#48)
e_pc = rho_pc*e_pc(-1)+ nu_pc; //(#49)
e_wc = rho_wc*e_wc(-1)+ nu_wc; //(#50)
e_ut = rho_ut*e_ut(-1)+ nu_ut; //(#51)
//2017/06/05 Okazaki Changed
e_d = rho_d*e_d(-1)+ nu_d; //(#52)
//e_d = nu_d; //(#)


// Burasagari Variables
NB  = NBQK*(QC*KC); //(#53)
NEC = NECQKC*(QC*KC); //(#54)
l_C = log(C/C(-1)); //(#55)
l_C_C = log(C_C/C_C(-1)); //(#56)
l_Cg = log(Cg/Cg(-1)); //(#57)
l_NB = log(NB/NB(-1)); //(#58)
l_NEC = log(NEC/NEC(-1)); //(#59)
l_IC = log(IC/IC(-1)); //(#60)
l_KC = log(KC/KC(-1)); //(#61)
l_GRKC = log(GRKC/GRKC(-1)); //(#62)
l_NSR = log(NSR); //(#63)
l_RSR = log(RSR); //(#64)
l_MCC = log(MCC/MCC(-1)); //(#65)
l_WC = log(WC/WC(-1)); //(#66)
l_WEC = log(WEC/WEC(-1)); //(#67)
l_WFIC = log(WFIC/WFIC(-1)); //(#68)
l_HC = log(HC/HC(-1)); //(#69)
l_QC = log(QC/QC(-1)); //(#70)
l_lambda_b = log(lambda_b/lambda_b(-1)); //(#71)
l_Output = log(Output/Output(-1)); //(#72)
l_Solow = log(Solow/Solow(-1)); //(#73)
l_pi_c = log(pi_c); //(#74)
l_pi_d = log(pi_d); //(#75)

mnc_cg =(mu_ec*gcdfUEC*(REC)*QC(-1)*KC(-1)/POP+mu_b*gcdfUB*(gammacdfUEC -mu_ec*gcdfUEC)*(REC)*QC(-1)*KC(-1)/POP)/Cg; //(#76)

l_REC=log(REC); //(#77)


// Entrepreners' profit maximization problem for those in non-durable sector
GRKC/Pc=kappa_uc*(exp(e_ut))^(phi_u+1)*(UC)^phi_u; //(#78)

l_UC = log(UC/UC(-1)); //(#79)


//F.O.C. of B_L
//1 = lambda_b(+1)/lambda_b*NSR*beta+phi_l*exp(e_d)/(lambda_b*Pc); //(#80)

//RLQ = NLQ/pi_c(+1); //(#81)

//phi_l/(lambda_b*Pc) = rho_l*(phi_l(-1)/(lambda_b(-1)*Pc(-1))) + nu_l; //(#82)

l_POP = nu_pop; //(#83)
POP = exp(l_POP); //(#84)


L_Output = l_Output + l_POP; //(#85)
L_C = l_C + l_POP; //(#86)
L_IC = l_IC + l_POP; //(#87)
L_HC = l_HC + l_POP; //(#88)
L_KC = l_KC + l_POP; //(#89)
L_NB = l_NB + l_POP; //(#90)
L_NEC = l_NEC + l_POP; //(#91)

SP_B = ZB - RSR; //(#92)
SP_EC = ZEC - ZB; //(#93)

EX_REC = REC(+1); //(#94)
SP_RE = REC - RSR; //(#95)
SP_EX_RE = EX_REC - RSR; //(#96)
RE_P = l_pi_d - l_pi_c; //(#97)

//2017/09/07 Okazaki Changed
//GAP = log(Output/Output_FL); //(#)
GAP = Output/Output_FL; //(#98)
RGAP = RSR - RSR_FL; //(#99)

////////////////////////////////////////////////////////////
////////////Flexible economy
////////////////////////////////////////////////////////////

//1 = lambda_b_FL(+1)/lambda_b_FL*(RSR_FL+RLQ_FL)*beta; //(#100)
1 = lambda_b_FL(+1)/lambda_b_FL*RSR_FL*beta; //(#100)

lambda_b_FL = exp(e_d)/(C_FL-habit*C_FL(-1))-beta*habit*exp(e_d(+1))*POP(+1)/(C_FL(+1)-habit*C_FL); //(#101)


REC_FL = (UC_FL*GRKC_FL-kappa_uc*((exp(e_ut)*UC_FL)^(1+phi_u)-1)/(1+phi_u)+QC_FL*(1-delta))/QC_FL(-1); //(#102)

//2017/08/09 Okazaki Changed
//QC_FL*(1-kappaC*(IC_FL*exp(e_kc)*POP/IC_FL(-1)-1)^2/2-(IC_FL*exp(e_kc)*POP/IC_FL(-1))*kappaC*(IC_FL*exp(e_kc)*POP/IC_FL(-1)-1))-1/exp(e_ad)/Zd = -beta*lambda_b_FL(+1)/lambda_b_FL*QC_FL(+1)*(IC_FL(+1)*POP(+1)/IC_FL)^2*exp(e_kc(+1))*kappaC*(IC_FL(+1)*exp(e_kc(+1))*POP(+1)/IC_FL-1); //(#)
QC_FL*(1-kappaC*(IC_FL*exp(e_kc)*POP/IC_FL(-1)-(eta_a_SS^A1*eta_d_SS^(A2+1)))^2/2-(IC_FL*exp(e_kc)*POP/IC_FL(-1))*kappaC*(IC_FL*exp(e_kc)*POP/IC_FL(-1)-(eta_a_SS^A1*eta_d_SS^(A2+1))))-1/exp(e_ad)/Zd = -beta*lambda_b_FL(+1)/lambda_b_FL*QC_FL(+1)*(IC_FL(+1)*POP(+1)/IC_FL)^2*exp(e_kc(+1))*kappaC*(IC_FL(+1)*exp(e_kc(+1))*POP(+1)/IC_FL-(eta_a_SS^A1*eta_d_SS^(A2+1))); //(#103)

//2017/08/09 Okazaki Changed
//KC_FL = (1-kappaC*(IC_FL*exp(e_kc)*POP/IC_FL(-1)-1)^2/2)*IC_FL + (1-delta)*KC_FL(-1)/POP; //(#)
KC_FL = (1-kappaC*(IC_FL*exp(e_kc)*POP/IC_FL(-1)-(eta_a_SS^A1*eta_d_SS^(A2+1)))^2/2)*IC_FL + (1-delta)*KC_FL(-1)/POP; //(#104)

gcdfUB_FL = normcdf(((log(omegabarUB_FL)-sigmaUB(-1)^2/2)/sigmaUB(-1))); //(#105)
gcdfUEC_FL = normcdf(((log(omegabarUEC_FL)-sigmaUEC(-1)^2/2)/sigmaUEC(-1))); //(#106)

d_gcdfUB_FL = 1/sigmaUB(-1)/omegabarUB_FL*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUB_FL)-sigmaUB(-1)^2/2)/sigmaUB(-1))^2); //(#107)

d_gcdfUEC_FL = 1/sigmaUEC(-1)/omegabarUEC_FL*(1/sqrt(2*3.14159265358979))*exp(-.5*((log(omegabarUEC_FL)-sigmaUEC(-1)^2/2)/sigmaUEC(-1))^2); //(#108)


gammacdfUB_FL = gcdfUB_FL +omegabarUB_FL*(1-normcdf(((log(omegabarUB_FL)+sigmaUB(-1)^2/2)/sigmaUB(-1)))); //(#109)

gammacdfUEC_FL = gcdfUEC_FL +omegabarUEC_FL*(1-normcdf(((log(omegabarUEC_FL)+sigmaUEC(-1)^2/2)/sigmaUEC(-1)))); //(#110)

d_gammacdfUB_FL = 1/sigmaUB(-1)/omegabarUB_FL *(1/sqrt(2*3.14159265358979))  *exp(-.5*((log(omegabarUB_FL)-sigmaUB(-1)^2/2)/sigmaUB(-1))^2) + (1-normcdf(((log(omegabarUB_FL)+sigmaUB(-1)^2/2)/sigmaUB(-1)))) - 1/sigmaUB(-1) *(1/sqrt(2*3.14159265358979)) *exp(-.5*((log(omegabarUB_FL)+sigmaUB(-1)^2/2)/sigmaUB(-1))^2); //(#111)

d_gammacdfUEC_FL = 1/sigmaUEC(-1)/omegabarUEC_FL *(1/sqrt(2*3.14159265358979)) *exp(-.5*((log(omegabarUEC_FL)-sigmaUEC(-1)^2/2)/sigmaUEC(-1))^2) + (1-normcdf(((log(omegabarUEC_FL)+sigmaUEC(-1)^2/2)/sigmaUEC(-1)))) - 1/sigmaUEC(-1) *(1/sqrt(2*3.14159265358979)) *exp(-.5*((log(omegabarUEC_FL)+sigmaUEC(-1)^2/2)/sigmaUEC(-1))^2); //(#112)


//(gammacdfUB_FL -mu_b*gcdfUB_FL)*(gammacdfUEC_FL -mu_ec*gcdfUEC_FL)*REC_FL = (RSR_FL(-1)+RLQ_FL(-1))*(1-NECQKC_FL(-1)-NBQK_FL(-1)); //(#113)
(gammacdfUB_FL -mu_b*gcdfUB_FL)*(gammacdfUEC_FL -mu_ec*gcdfUEC_FL)*REC_FL = RSR_FL(-1)*(1-NECQKC_FL(-1)-NBQK_FL(-1)); //(#113)

1 - gammacdfUEC_FL = NECQKC_FL(-1); //(#114)


//0=(1-gammacdfUB_FL(+1))*(gammacdfUEC_FL(+1)-mu_ec*gcdfUEC_FL(+1))*REC_FL(+1) + d_gammacdfUB_FL(+1)/(d_gammacdfUB_FL(+1)-mu_b*d_gcdfUB_FL(+1))*((gammacdfUB_FL(+1)-mu_b*gcdfUB_FL(+1))*(gammacdfUEC_FL(+1)-mu_ec*gcdfUEC_FL(+1))*REC_FL(+1)) - d_gammacdfUB_FL(+1)/(d_gammacdfUB_FL(+1)-mu_b*d_gcdfUB_FL(+1))*(RSR_FL+RLQ_FL)+(1-gammacdfUB_FL(+1))*(d_gammacdfUEC_FL(+1)-mu_ec*d_gcdfUEC_FL(+1))*(1-gammacdfUEC_FL(+1))/d_gammacdfUEC_FL(+1)*REC_FL(+1)+d_gammacdfUB_FL(+1)*(gammacdfUB_FL(+1)-mu_b*gcdfUB_FL(+1))*(d_gammacdfUEC_FL(+1)-mu_ec*d_gcdfUEC_FL(+1))*(1-gammacdfUEC_FL(+1))*REC_FL(+1)/(d_gammacdfUB_FL(+1)-mu_b*d_gcdfUB_FL(+1))/d_gammacdfUEC_FL(+1); //(#115)
0=(1-gammacdfUB_FL(+1))*(gammacdfUEC_FL(+1)-mu_ec*gcdfUEC_FL(+1))*REC_FL(+1) + d_gammacdfUB_FL(+1)/(d_gammacdfUB_FL(+1)-mu_b*d_gcdfUB_FL(+1))*((gammacdfUB_FL(+1)-mu_b*gcdfUB_FL(+1))*(gammacdfUEC_FL(+1)-mu_ec*gcdfUEC_FL(+1))*REC_FL(+1)) - d_gammacdfUB_FL(+1)/(d_gammacdfUB_FL(+1)-mu_b*d_gcdfUB_FL(+1))*RSR_FL+(1-gammacdfUB_FL(+1))*(d_gammacdfUEC_FL(+1)-mu_ec*d_gcdfUEC_FL(+1))*(1-gammacdfUEC_FL(+1))/d_gammacdfUEC_FL(+1)*REC_FL(+1)+d_gammacdfUB_FL(+1)*(gammacdfUB_FL(+1)-mu_b*gcdfUB_FL(+1))*(d_gammacdfUEC_FL(+1)-mu_ec*d_gcdfUEC_FL(+1))*(1-gammacdfUEC_FL(+1))*REC_FL(+1)/(d_gammacdfUB_FL(+1)-mu_b*d_gcdfUB_FL(+1))/d_gammacdfUEC_FL(+1); //(#115)

NBQK_FL*QC_FL*KC_FL/(QC_FL(-1)*KC_FL(-1))*POP = gammab*(1-gammacdfUB_FL)*(gammacdfUEC_FL -mu_ec*gcdfUEC_FL)*(REC_FL) + (e_nb*exp(Output_FL_SS))/(QC_FL(-1)*KC_FL(-1))+alpha_FI/(1-alpha-alpha_E-alpha_FI)*(UC_FL)*GRKC_FL/QC_FL(-1); //(#116)


NECQKC_FL*(QC_FL*KC_FL/QC_FL(-1)/KC_FL(-1))*POP = gammaec*(1-gammacdfUEC_FL)*(REC_FL)+(e_nc*exp(Output_FL_SS))/QC_FL(-1)/KC_FL(-1)+alpha_E/(1-alpha-alpha_E-alpha_FI)*(UC_FL)*GRKC_FL/QC_FL(-1); //(#117)


Cg_FL=(1+gec)*C_FL+C_C_FL+IC_FL/exp(e_ad)/Zd+mu_ec*gcdfUEC_FL*(REC_FL)*QC_FL(-1)*KC_FL(-1)/POP+mu_b*gcdfUB_FL*(gammacdfUEC_FL-mu_ec*gcdfUEC_FL)*(REC_FL)*QC_FL(-1)*KC_FL(-1)/POP+kappa_uc*((exp(e_ut)*UC_FL)^(1+phi_u)-1)/(1+phi_u)*KC_FL(-1)/POP+(1-gammaec)*(1-gammacdfUEC_FL)*(REC_FL)*QC_FL(-1)*KC_FL(-1)/POP+(1-gammab)*(1-gammacdfUB_FL)*(gammacdfUEC_FL -mu_ec*gcdfUEC_FL)*(REC_FL)*QC_FL(-1)*KC_FL(-1)/POP; //(#118)


Output_FL = (1+gec)*C_FL+IC_FL/exp(e_ad)/Zd; //(#119)

Solow_FL=Output_FL/(KC_FL(-1)^(1-Share_Labor)*HC_FL^(Share_Labor))*POP^(1-Share_Labor); //(#120)

gama*MCC_FL*(Cg_FL+Fc)=C_C_FL; //(#121)
(1-gama)*(1-alpha-alpha_E-alpha_FI)*MCC_FL*(Cg_FL+Fc)=(UC_FL)*GRKC_FL*KC_FL(-1)/POP; //(#122)
(1-gama)*alpha*MCC_FL*(Cg_FL+Fc)=WC_FL*HC_FL; //(#123)
(1-gama)*alpha_E*MCC_FL*(Cg_FL+Fc)=WEC_FL; //(#124)
(1-gama)*alpha_FI*MCC_FL*(Cg_FL+Fc)=WFIC_FL; //(#125)


//exp(e_d)*theta_wc*exp(e_wc)*ksi*(HC_FL)^vv=-lambda_b_FL*WC_FL*(1-theta_wc*exp(e_wc)); //(#126)
exp(e_d)*theta_wc*exp(0)*ksi*(HC_FL)^vv=-lambda_b_FL*WC_FL*(1-theta_wc*exp(0)); //(#126)


MCC_FL=1/exp(e_aa)/Za*phi_c*(WEC_FL^alpha_E*WFIC_FL^alpha_FI*WC_FL^alpha*GRKC_FL^(1-alpha-alpha_E-alpha_FI))^(1-gama); //(#127)


//(1-theta_c*exp(e_pc))=-theta_c*exp(e_pc)*MCC_FL; //(#128)
(1-theta_c*exp(0))=-theta_c*exp(0)*MCC_FL; //(#128)


//2017/09/21 Okazaki Added
l_RLR_FL = (l_RSR_FL + l_RSR_FL(+1) + l_RSR_FL(+2)+ l_RSR_FL(+3)+ l_RSR_FL(+4)+ l_RSR_FL(+5)+ l_RSR_FL(+6)+ l_RSR_FL(+7)+ l_RSR_FL(+8)+ l_RSR_FL(+9)+ l_RSR_FL(+10)+ l_RSR_FL(+11)+ l_RSR_FL(+12)+ l_RSR_FL(+13)+ l_RSR_FL(+14)+ l_RSR_FL(+15)+ l_RSR_FL(+16)+ l_RSR_FL(+17)+ l_RSR_FL(+18)+ l_RSR_FL(+19)+ l_RSR_FL(+20)+ l_RSR_FL(+21)+ l_RSR_FL(+22)+ l_RSR_FL(+23)+ l_RSR_FL(+24)+ l_RSR_FL(+25)+ l_RSR_FL(+26)+ l_RSR_FL(+27)+ l_RSR_FL(+28)+ l_RSR_FL(+29)+ l_RSR_FL(+30)+ l_RSR_FL(+31)+ l_RSR_FL(+32)+ l_RSR_FL(+33)+ l_RSR_FL(+34)+ l_RSR_FL(+35)+ l_RSR_FL(+36)+ l_RSR_FL(+37)+ l_RSR_FL(+38)+ l_RSR_FL(+39))/40;

RLR_FL = exp(l_RLR_FL);

ZEC_FL = 1/(1-NECQKC_FL)*omegabarUEC_FL(+1)*(REC_FL(+1)); //(#129)
ZB_FL = (1-NECQKC_FL)/(1-NECQKC_FL-NBQK_FL)*omegabarUB_FL(+1)*(RB_FL(+1)); //(#130)
//RB_FL/(RSR_FL(-1)+RLQ_FL(-1)) = 1/(gammacdfUB_FL-mu_b*gcdfUB_FL)*(1-NBQK_FL(-1)-NECQKC_FL(-1))/(1-NECQKC_FL(-1)); //(#131)
RB_FL/RSR_FL(-1) = 1/(gammacdfUB_FL-mu_b*gcdfUB_FL)*(1-NBQK_FL(-1)-NECQKC_FL(-1))/(1-NECQKC_FL(-1)); //(#131)

NB_FL  = NBQK_FL*(QC_FL*KC_FL); //(#132)
NEC_FL = NECQKC_FL*(QC_FL*KC_FL); //(#133)
l_C_FL = log(C_FL/C_FL(-1)); //(#134)
l_C_C_FL = log(C_C_FL/C_C_FL(-1)); //(#135)
l_Cg_FL = log(Cg_FL/Cg_FL(-1)); //(#136)
l_NB_FL = log(NB_FL/NB_FL(-1)); //(#137)
l_NEC_FL = log(NEC_FL/NEC_FL(-1)); //(#138)
l_IC_FL = log(IC_FL/IC_FL(-1)); //(#139)
l_KC_FL = log(KC_FL/KC_FL(-1)); //(#140)
l_GRKC_FL = log(GRKC_FL/GRKC_FL(-1)); //(#141)
//l_NSR_FL = log(NSR_FL); //(#)
l_RSR_FL = log(RSR_FL); //(#142)
l_MCC_FL = log(MCC_FL/MCC_FL(-1)); //(#143)
l_WC_FL = log(WC_FL/WC_FL(-1)); //(#144)
l_WEC_FL = log(WEC_FL/WEC_FL(-1)); //(#145)
l_WFIC_FL = log(WFIC_FL/WFIC_FL(-1)); //(#146)
l_HC_FL = log(HC_FL/HC_FL(-1)); //(#147)
l_QC_FL = log(QC_FL/QC_FL(-1)); //(#148)
l_lambda_b_FL = log(lambda_b_FL/lambda_b_FL(-1)); //(#149)
l_Output_FL = log(Output_FL/Output_FL(-1)); //(#150)
l_Solow_FL = log(Solow_FL/Solow_FL(-1)); //(#151)

mnc_cg_FL =(mu_ec*gcdfUEC_FL*(REC_FL)*QC_FL(-1)*KC_FL(-1)/POP+mu_b*gcdfUB_FL*(gammacdfUEC_FL -mu_ec*gcdfUEC_FL)*(REC_FL)*QC_FL(-1)*KC_FL(-1)/POP)/Cg_FL; //(#152)

l_REC_FL=log(REC_FL); //(#153)

GRKC_FL=kappa_uc*(exp(e_ut))^(phi_u+1)*(UC_FL)^phi_u; //(#154)

l_UC_FL = log(UC_FL/UC_FL(-1)); //(#155)

//1 = lambda_b_FL(+1)/lambda_b_FL*RSR_FL*beta+phi_l*exp(e_d)/lambda_b_FL; //(#156)

L_Output_FL = l_Output_FL + l_POP; //(#157)
L_C_FL = l_C_FL + l_POP; //(#158)
L_IC_FL = l_IC_FL + l_POP; //(#159)
L_HC_FL = l_HC_FL + l_POP; //(#160)
L_KC_FL = l_KC_FL + l_POP; //(#161)
L_NB_FL = l_NB_FL + l_POP; //(#162)
L_NEC_FL = l_NEC_FL + l_POP; //(#163)

SP_B_FL = ZB_FL - RSR_FL; //(#164)
SP_EC_FL = ZEC_FL - ZB_FL; //(#165)

EX_REC_FL = REC_FL(+1); //(#166)
SP_RE_FL = REC_FL - RSR_FL; //(#167)
SP_EX_RE_FL = EX_REC_FL - RSR_FL; //(#168)

end;

initval;
C=C_SS;
NSR=NSR_SS;
Pd=1;
GRKC=GRKC_SS;
REC=REC_SS;
QC=QC_SS;
IC=IC_SS;
KC=KC_SS;
gcdfUB=gcdfUB_SS;
gcdfUEC=gcdfUEC_SS;
omegabarUB=omegabarUB_SS;
omegabarUEC=omegabarUEC_SS;
sigmaUB=sigmaUB_SS;
sigmaUEC=sigmaUEC_SS;
d_gcdfUB=d_gcdfUB_SS;
d_gcdfUEC=d_gcdfUEC_SS;
gammacdfUB=gammacdfUB_SS;
gammacdfUEC=gammacdfUEC_SS;
d_gammacdfUB=d_gammacdfUB_SS;
d_gammacdfUEC=d_gammacdfUEC_SS;
NECQKC=NECQKC_SS;
ZEC = 1/(1-NECQKC_SS)*omegabarUEC_SS*REC_SS;
ZB = (1-NECQKC_SS)/(1-NECQKC_SS-NBQK_SS)*omegabarUB_SS*(RB_SS);
//RB = 1/(gammacdfUB_SS-mu_b*gcdfUB_SS)*(1-NECQKC_SS-NBQK_SS)/(1-NECQKC_SS)*(RSR_SS+RLQ_SS);
RB = 1/(gammacdfUB_SS-mu_b*gcdfUB_SS)*(1-NECQKC_SS-NBQK_SS)/(1-NECQKC_SS)*RSR_SS;
NBQK=NBQK_SS;
RSR=RSR_SS;

Cg=CG_SS;
MCC=MCC_SS;
WC=WC_SS;
WEC=WEC_SS;
WFIC=WFIC_SS;
HC=HC_SS;
C_C=C_C_SS;
Output=Output_SS;
Solow=Solow_SS;
lambda_b=lambda_b_SS;
//2017/07/03 Okazaki Changed
//pi_c=1;
pi_c=pi_c_SS;
l_pi_c = log(pi_c_SS);
l_MCC = log(pi_c_SS);

//2017/08/09 OKazaki Changed
//l_WC = log(pi_c_SS);
l_WC = log(pi_c_SS*eta_a_SS^A1*eta_d_SS^A2);
//l_WEC = log(pi_c_SS);
l_WEC = log(pi_c_SS*eta_a_SS^A1*eta_d_SS^A2);
//l_WFIC = log(pi_c_SS);
l_WFIC = log(pi_c_SS*eta_a_SS^A1*eta_d_SS^A2);

l_lambda_b = -log(pi_c_SS);
//2017/08/09 Okazaki Changed
//pi_d=1;
//pi_d=pi_c_SS;
pi_d=pi_c_SS/eta_d_SS;
l_pi_d = log(pi_c_SS/eta_d_SS);
e_aa=0;
e_ad=0;
e_r=0;

//2017/09/12 Okazaki Added
e_1 = 0;
e_2 = 0;
e_3 = 0;
e_4 = 0;

//2017/09/15 Okazaki Added
EXP_R1 = log(NSR_SS);
EXP_R2 = log(NSR_SS);
EXP_R3 = log(NSR_SS);
EXP_R4 = log(NSR_SS);

//2017/11/08 Okazaki Added
e_5 = 0;
e_6 = 0;
e_7 = 0;
e_8 = 0;
e_9 = 0;
e_10 = 0;
e_11 = 0;
e_12 = 0;
EXP_R5 = log(NSR_SS);
EXP_R6 = log(NSR_SS);
EXP_R7 = log(NSR_SS);
EXP_R8 = log(NSR_SS);
EXP_R9 = log(NSR_SS);
EXP_R10 = log(NSR_SS);
EXP_R11 = log(NSR_SS);
EXP_R12 = log(NSR_SS);

e_nb=0;
e_nc=0;
eta_a=eta_a_SS; // added 2012/2/1
eta_d=eta_d_SS; // added 2012/2/1
nu_aa=0;
nu_ad=0;
nu_r=0;
nu_sb=0;
nu_sc=0;
nu_nb=0;
nu_nc=0;
nu_ea=0; // added 2012/2/1
nu_ed=0; // added 2012/2/1
NB  = NB_SS;
NEC = NEC_SS;
gec = gec_SS;
nu_gc = 0;
l_NSR = log(NSR_SS);
//2017/07/03 Okazaki Changed
//l_RSR = log(NSR_SS);
l_RSR = log(RSR_SS);

e_kc = 0;
e_pc = 0;
e_wc = 0;
nu_kc = 0;
nu_pc = 0;
nu_wc = 0;
e_d =0;
nu_d = 0;
//NLQ = NLQ_SS;
//RLQ = RLQ_SS;
//nu_l = 0;
//phi_l = phi_l_SS;
nu_pop = 0;
POP=POP_SS;

mnc_cg =(mu_ec*gcdfUEC_SS*(REC_SS)*QC_SS*KC_SS+mu_b*gcdfUB_SS*((gammacdfUEC_SS -mu_ec*gcdfUEC_SS)*(REC_SS)*QC_SS*KC_SS))/CG_SS;

l_REC=log(REC_SS);
e_ut=0;

//2017/08/09 Okazaki Changed
//UC=1;
//l_UC=0;
UC = UC_SS;
l_UC = log(UC_SS);

SP_EC = 1/(1-NECQKC_SS)*omegabarUEC_SS*REC_SS - (1-NECQKC_SS)/(1-NECQKC_SS-NBQK_SS)*omegabarUB_SS*(RB_SS);
SP_B = (1-NECQKC_SS)/(1-NECQKC_SS-NBQK_SS)*omegabarUB_SS*(RB_SS) - RSR_SS;

EX_REC = REC_SS;
SP_RE = REC_SS - RSR_SS;
SP_EX_RE = REC_SS - RSR_SS;
RE_P = log(pi_c_SS/eta_d_SS) - log(pi_c_SS);

//2017/09/07 Okazaki Changed
//GAP = 0;
GAP =1;
RGAP = 0;

pi_star = pi_c_SS;

//2017/08/08 Okazaki Added
l_C = log(eta_a_SS^A1*eta_d_SS^A2);
l_C_C = log(eta_a_SS^A1*eta_d_SS^A2);
l_Cg = log(eta_a_SS^A1*eta_d_SS^A2);
l_NB = log(eta_a_SS^A1*eta_d_SS^A2);
l_NEC = log(eta_a_SS^A1*eta_d_SS^A2);
l_IC = log(eta_a_SS^A1*eta_d_SS^(A2+1));
l_KC = log(eta_a_SS^A1*eta_d_SS^(A2+1));
l_GRKC = log(pi_c_SS*eta_d_SS^-1);
l_QC = log(1/eta_d_SS);
l_lambda_b = log((pi_c_SS*eta_a_SS^A1*eta_d_SS^A2)^-1);
l_Output = log(eta_a_SS^A1*eta_d_SS^A2);
l_Solow = log(eta_a_SS^(A1*Share_Labor)*eta_d_SS^(A2*Share_Labor-(1-Share_Labor)));

L_Output = log(eta_a_SS^A1*eta_d_SS^A2) + log(POP_SS);
L_C = log(eta_a_SS^A1*eta_d_SS^A2) + log(POP_SS);
L_IC = log(eta_a_SS^A1*eta_d_SS^(A2+1)) + log(POP_SS);
L_KC = log(eta_a_SS^A1*eta_d_SS^(A2+1)) + log(POP_SS);
L_NB = log(eta_a_SS^A1*eta_d_SS^A2) + log(POP_SS);
L_NEC = log(eta_a_SS^A1*eta_d_SS^A2) + log(POP_SS);

C_FL=C_FL_SS;
GRKC_FL=GRKC_FL_SS;
REC_FL=REC_FL_SS;
QC_FL=QC_FL_SS;
IC_FL=IC_FL_SS;
KC_FL=KC_FL_SS;
gcdfUB_FL=gcdfUB_FL_SS;
gcdfUEC_FL=gcdfUEC_FL_SS;
omegabarUB_FL=omegabarUB_FL_SS;
omegabarUEC_FL=omegabarUEC_FL_SS;
d_gcdfUB_FL=d_gcdfUB_FL_SS;
d_gcdfUEC_FL=d_gcdfUEC_FL_SS;
gammacdfUB_FL=gammacdfUB_FL_SS;
gammacdfUEC_FL=gammacdfUEC_FL_SS;
d_gammacdfUB_FL=d_gammacdfUB_FL_SS;
d_gammacdfUEC_FL=d_gammacdfUEC_FL_SS;
NECQKC_FL=NECQKC_FL_SS;
ZEC_FL = 1/(1-NECQKC_FL_SS)*omegabarUEC_FL_SS*REC_FL_SS;
ZB_FL = (1-NECQKC_FL_SS)/(1-NECQKC_FL_SS-NBQK_FL_SS)*omegabarUB_FL_SS*(RB_FL_SS);
//RB_FL = 1/(gammacdfUB_FL_SS-mu_b*gcdfUB_FL_SS)*(1-NECQKC_FL_SS-NBQK_FL_SS)/(1-NECQKC_FL_SS)*(RSR_FL_SS+RLQ_FL_SS);
RB_FL = 1/(gammacdfUB_FL_SS-mu_b*gcdfUB_FL_SS)*(1-NECQKC_FL_SS-NBQK_FL_SS)/(1-NECQKC_FL_SS)*RSR_FL_SS;
NBQK_FL=NBQK_FL_SS;
RSR_FL=RSR_FL_SS;

Cg_FL=CG_FL_SS;
MCC_FL=MCC_FL_SS;
WC_FL=WC_FL_SS;
WEC_FL=WEC_FL_SS;
WFIC_FL=WFIC_FL_SS;
HC_FL=HC_FL_SS;
C_C_FL=C_C_FL_SS;
Output_FL=Output_FL_SS;
Solow_FL=Solow_FL_SS;
lambda_b_FL=lambda_b_FL_SS;
NB_FL  = NB_FL_SS;
NEC_FL = NEC_FL_SS;
l_RSR_FL = log(RSR_FL_SS);


mnc_cg_FL =(mu_ec*gcdfUEC_FL*(REC_FL_SS)*QC_FL_SS*KC_FL_SS+mu_b*gcdfUB_FL*((gammacdfUEC_FL -mu_ec*gcdfUEC_FL)*(REC_FL_SS)*QC_FL_SS*KC_FL_SS))/CG_FL_SS;

l_REC_FL=log(REC_FL_SS);

//2017/08/09 Okazaki Changed
//UC_FL=1;
//l_UC_FL=0;
UC_FL=UC_FL_SS;
l_UC_FL=log(UC_FL_SS);

//RLQ_FL = RLQ_FL_SS;

SP_EC_FL = 1/(1-NECQKC_SS)*omegabarUEC_SS*REC_SS - (1-NECQKC_SS)/(1-NECQKC_SS-NBQK_SS)*omegabarUB_SS*(RB_SS);
SP_B_FL = (1-NECQKC_SS)/(1-NECQKC_SS-NBQK_SS)*omegabarUB_SS*(RB_SS) - RSR_SS;

EX_REC_FL = REC_FL_SS;
SP_RE_FL = REC_FL_SS - RSR_FL_SS;
SP_EX_RE_FL = REC_FL_SS - RSR_FL_SS;

//2017/08/09 Okazaki Added
l_C_FL = log(eta_a_SS^A1*eta_d_SS^A2);
l_C_C_FL = log(eta_a_SS^A1*eta_d_SS^A2);
l_Cg_FL = log(eta_a_SS^A1*eta_d_SS^A2);
l_NB_FL = log(eta_a_SS^A1*eta_d_SS^A2);
l_NEC_FL = log(eta_a_SS^A1*eta_d_SS^A2);
l_IC_FL = log(eta_a_SS^A1*eta_d_SS^(A2+1));
l_KC_FL = log(eta_a_SS^A1*eta_d_SS^(A2+1));
l_GRKC_FL = log(eta_d_SS^-1);

l_WC_FL = log(eta_a_SS^A1*eta_d_SS^A2);
l_WEC_FL = log(eta_a_SS^A1*eta_d_SS^A2);
l_WFIC_FL = log(eta_a_SS^A1*eta_d_SS^A2);

l_QC_FL = log(1/eta_d_SS);
l_lambda_b_FL = log((eta_a_SS^A1*eta_d_SS^A2)^-1);
l_Output_FL = log(eta_a_SS^A1*eta_d_SS^A2);
l_Solow_FL = log(eta_a_SS^(A1*Share_Labor)*eta_d_SS^(A2*Share_Labor-(1-Share_Labor)));

L_Output_FL = log(eta_a_SS^A1*eta_d_SS^A2) + log(POP_SS);
L_C_FL = log(eta_a_SS^A1*eta_d_SS^A2) + log(POP_SS);
L_IC_FL = log(eta_a_SS^A1*eta_d_SS^(A2+1)) + log(POP_SS);
L_KC_FL = log(eta_a_SS^A1*eta_d_SS^(A2+1)) + log(POP_SS);
L_NB_FL = log(eta_a_SS^A1*eta_d_SS^A2) + log(POP_SS);
L_NEC_FL = log(eta_a_SS^A1*eta_d_SS^A2) + log(POP_SS);

//2017/09/21 Okazaki Added
RLR = RLR_SS;
l_RLR = log(RLR_SS);
RLR_FL = RLR_FL_SS;
l_RLR_FL = log(RLR_FL_SS);

l_NLR = log(NSR_SS);
NLR = NSR_SS;

RSR_MA=RSR_SS;
RSR_FL_MA=RSR_SS;
l_RSR_MA=log(RSR_SS);
l_RSR_FL_MA=log(RSR_SS);
N1R=NSR_SS;
N5R=NSR_SS;
l_N1R=log(NSR_SS);
l_N5R=log(NSR_SS);
R1R=RSR_SS;
R5R=RSR_SS;
R1R_FL=RSR_SS;
R5R_FL=RSR_SS;
l_R1R=log(RSR_SS);
l_R1R_FL=log(RSR_SS);
l_R5R=log(RSR_SS);
l_R5R_FL=log(RSR_SS);



end;

shocks;
var nu_ad; stderr sigma_ad;
var nu_aa; stderr sigma_aa;
var nu_r ; stderr sigma_r;
//var nu_sb; stderr sigma_sb;
//var nu_sc; stderr sigma_sc;
var nu_nb; stderr sigma_nb;
var nu_nc; stderr sigma_nc;
var nu_ed; stderr sigma_ed;
var nu_ea; stderr sigma_ea;
var nu_gc; stderr sigma_gc;
var nu_kc; stderr sigma_kc;
var nu_pc; stderr sigma_pc;
var nu_wc; stderr sigma_wc;
var nu_d;  stderr sigma_d;
//var nu_l;  stderr sigma_l;

var nu_pop; stderr 0.00002;

var nu_pi;  stderr sigma_pi;

//2017/09/12 Okazaki Added
var nu_1;  stderr sigma_1;
var nu_2;  stderr sigma_2;
var nu_3;  stderr sigma_3;
var nu_4;  stderr sigma_4;

//2017/11/08 Okazaki Added
var nu_5;  stderr sigma_5;
var nu_6;  stderr sigma_6;
var nu_7;  stderr sigma_7;
var nu_8;  stderr sigma_8;
var nu_9;  stderr sigma_9;
var nu_10;  stderr sigma_10;
var nu_11;  stderr sigma_11;
var nu_12;  stderr sigma_12;

end;

resid;
steady(solve_algo = 3);
// steady(solve_algo = 0);

model_diagnostics;

check;
//check(solve_algo = 0);

// Estimation
est=1;
if est==1;

//2017/09/20 Okazaki Changed
//varobs l_pi_c l_IC l_Output l_WC l_HC l_NSR l_NB l_NEC l_pi_d l_Solow NLQ l_POP;
//varobs l_pi_c l_IC l_Output l_WC l_HC l_NSR l_NB l_NEC l_pi_d l_Solow NLQ l_POP EXP_R1 EXP_R2 EXP_R3 EXP_R4;
//varobs l_pi_c l_IC l_Output l_WC l_HC l_NSR l_NB l_NEC l_pi_d l_Solow l_POP EXP_R1 EXP_R2 EXP_R3 EXP_R4;
varobs l_pi_c l_IC l_Output l_WC l_HC l_NSR l_NB l_NEC l_pi_d l_Solow l_POP EXP_R1 EXP_R2 EXP_R3 EXP_R4 EXP_R5 EXP_R6 EXP_R7 EXP_R8 EXP_R9 EXP_R10 EXP_R11 EXP_R12;

estimated_params;

//vv,		gamma_pdf,	1,	0.1;
vv,		gamma_pdf,	0.8,	0.075;
//kappaC,	gamma_pdf,	1,	0.1;
kappaC,	gamma_pdf,	2,	0.25;
//kappa_pc,	gamma_pdf,	20,	10;
//kappa_wc,	gamma_pdf,	20,	10;
kappa_pc,	gamma_pdf,	12,	1;
kappa_wc,	gamma_pdf,	2.5,	0.5;
rpi,		normal_pdf,	2.75,	0.05;
phi_u,		gamma_pdf,	5,	1;
//rho_NSR,	normal_pdf,	0.9,	.1;
rho_NSR,	beta_pdf,	0.5,	.01; //MJEM

sigmaUEC_SS,gamma_pdf,	 0.309309721721658, .002;
//sigmaUEC_SS,gamma_pdf,	 0.7, 0.1; //2016/09/20 Okazaki Changed
sigmaUB_SS,	gamma_pdf,	0.104215265147956, .002;
mu_ec,		gamma_pdf,	 0.019618891439117, .01; //03
mu_b,		gamma_pdf,	 0.538581453181288, .01; //03
//gammaec,	beta_pdf,		 0.973943610929912, .001;
gammaec,	beta_pdf,		 0.96, .001;
//gammab,		beta_pdf,		0.922805418053119, .001;
gammab,		beta_pdf,		0.86, .001;

rho_aa,		beta_pdf,	0.5,	0.15;
rho_ad,		beta_pdf,	0.5,	0.15;

// rho_r,	beta_pdf,	0.9,	0.15;

//rho_sc,	beta_pdf,	0.85,	0.1;
//rho_sb,	beta_pdf,	0.85,	0.1;
rho_nb,		beta_pdf,	0.5,	0.15;
rho_nc,		beta_pdf,	0.5,	0.15;
rho_ea,		beta_pdf,	0.5,	0.15;
rho_gc,		beta_pdf,	0.5,	0.15;
rho_kc,		beta_pdf,	0.5,	0.15;
rho_pc,		beta_pdf,	0.5,	0.15;
rho_wc,    	beta_pdf,	0.5,	0.15;

//rho_ut,		beta_pdf,	0.5,	0.15;

rho_ed,		beta_pdf,	0.5,	0.15;
//2017/06/05 Okazaki Deleted
rho_d,		beta_pdf,	0.5,	0.15;

//2017/08/08 Okazaki Added
eta_a_SS, 	gamma_pdf,	1.0011,	0.001;
eta_d_SS,     gamma_pdf,	1.0016,	0.001;
//eta_a_SS, 	normal_pdf,	1.0011,	0.001;
//eta_d_SS,	normal_pdf,	1.0016,	0.001;


stderr nu_aa,	inv_gamma_pdf,	.05,	5;
stderr nu_ad,	inv_gamma_pdf,	.05,	5;
stderr nu_r,	inv_gamma_pdf,	.01,	5;
stderr nu_nb,	inv_gamma_pdf,	.02,	5;
stderr nu_nc,	inv_gamma_pdf,	.02,	5;
stderr nu_ea,	inv_gamma_pdf,	.01,	5;
stderr nu_gc,	inv_gamma_pdf,	.01,	5;
stderr nu_kc,	inv_gamma_pdf,	.01,	5;
stderr nu_pc,	inv_gamma_pdf,	.01,	5;
stderr nu_wc,	inv_gamma_pdf,	.01,	5;
stderr nu_ed,	inv_gamma_pdf,	.01,	5;
stderr nu_d,	inv_gamma_pdf,	.015,	5;

//rho_l,		beta_pdf,	0.5,	0.15;
//stderr nu_l,	inv_gamma_pdf,	.05,	5;

rho_pi,		beta_pdf,	0.5,	0.15;
stderr nu_pi,	inv_gamma_pdf,	.01,	5;

stderr l_NB,	inv_gamma_pdf,	0.1,	5;
stderr l_NEC,	inv_gamma_pdf,	0.1,	5;
//stderr NLQ,	inv_gamma_pdf,	0.01,	5;

//2017/08/22 Okazaki Added
//pi_c_SS,	normal_pdf,	1.000318454,	0.0001;
pi_c_SS,	normal_pdf,	1.0025,		0.001;

//2017/09/12 Okazaki Added
stderr nu_1,	inv_gamma_pdf,	0.01,	5;
stderr nu_2,	inv_gamma_pdf,	0.01,	5;
stderr nu_3,	inv_gamma_pdf,	0.01,	5;
stderr nu_4,	inv_gamma_pdf,	0.01,	5;

//2017/11/08 Okazaki Added
stderr nu_5,	inv_gamma_pdf,	0.01,	5;
stderr nu_6,	inv_gamma_pdf,	0.01,	5;
stderr nu_7,	inv_gamma_pdf,	0.01,	5;
stderr nu_8,	inv_gamma_pdf,	0.01,	5;
stderr nu_9,	inv_gamma_pdf,	0.01,	5;
stderr nu_10,	inv_gamma_pdf,	0.01,	5;
stderr nu_11,	inv_gamma_pdf,	0.01,	5;
stderr nu_12,	inv_gamma_pdf,	0.01,	5;

habit,    	beta_pdf,	0.6,	0.15;

end;


estimated_params_init;

vv, 0.9151;
kappaC, 2.0941;
kappa_pc, 6.6960;
kappa_wc, 1.0773;
rpi, 2.8715;
phi_u, 6.10475;
rho_NSR, 0.5140;
mu_ec, 0.0458;
mu_b, 0.5335;
rho_aa, 0.9194;
rho_ad, 0.6963;
rho_nb, 0.2280;
rho_nc, 0.6925;
rho_ea, 0.0997;
rho_gc, 0.9616;
rho_kc, 0.1781;
rho_pc, 0.8730;
rho_wc, 0.8617;
rho_ed, 0.5141;
//rho_l, 0.2489;
rho_pi, 0.4616;
rho_d, 0.8107;
stderr nu_aa, 0.0059;
stderr nu_ad, 0.0059;
stderr nu_r, 0.0015;
stderr nu_nb, 0.032;
stderr nu_nc, 0.051;
stderr nu_ea, 0.0023;
stderr nu_gc, 0.0466;
stderr nu_kc, 0.0302;
stderr nu_pc, 0.0277;
stderr nu_wc, 0.1010;
stderr nu_ed, 0.0024;
//stderr nu_l, 0.0059;
stderr nu_pi, 0.0038;
stderr nu_d, 0.0043;

end;


estimation(
//optim=('MaxIter',30),
//optim=('TolFun',1e-5),
//optim=('TolX',1e-5),
optim=('Display','iter')
,datafile=DATA_BJ1115, nobs=149, order=1
//,mode_check
//,mode_compute=1 //case of not using mode_file
,mode_compute=0 //case of using mode_file
,mode_file=bjem1120r1_mode
,kalman_algo=2
,mh_replic=400000 //2016/09/20 Okazaki Changed
//,mh_replic=0
//,load_mh_file
,mh_nblocks=2
,mh_jscale=0.09 //2016/09/27 Okazaki Changed
,mh_init_scale=0.01
//,conf_sig = 0.95
//,mh_conf_sig = 0.95
,smoother
,filtered_vars
//,filter_step_ahead=1
//,filter_decomposition
,bayesian_irf
,irf=41
,forecast=41
//,diffuse_filter
,graph_format = (pdf)
//,nograph
,nodisplay // 2016/09/20 Okazaki Added
)
l_Output,l_C,l_IC,l_NB,l_NEC,l_KC,l_WC,l_HC,L_Output,L_C,L_IC,L_NB,L_NEC,L_KC,L_HC,l_Output_FL,l_C_FL,l_IC_FL,l_NB_FL,l_NEC_FL,l_KC_FL,l_WC_FL,l_HC_FL,L_Output_FL,L_C_FL,L_IC_FL,L_NB_FL,L_NEC_FL,L_KC_FL,L_HC_FL,NSR,RSR,RSR_FL,ZB,ZB_FL,ZEC,ZEC_FL,REC,REC_FL,SP_B,SP_B_FL,SP_EC,SP_EC_FL,EX_REC,EX_REC_FL,SP_RE,SP_RE_FL,SP_EX_RE,SP_EX_RE_FL,l_pi_c,l_pi_d,RE_P,l_POP,GAP,RGAP,C_C,Output,NB,NEC,C_C_FL,Cg_FL,Output_FL,NB_FL,NEC_FL,KC,KC_FL,Pd,GRKC,GRKC_FL,WEC,WFIC,WC_FL,WEC_FL,WFIC_FL,MCC,MCC_FL,Solow,Solow_FL,NECQKC,NBQK,gec,HC_FL,NECQKC_FL,NBQK_FL,sigmaUB,sigmaUEC,e_aa,e_ad,e_r,e_nb,e_nc,e_pc,e_wc,e_ut,UC,UC_FL,pi_c,eta_a,eta_d,C,Cg,C_FL,IC,IC_FL,WC,QC,QC_FL,lambda_b,lambda_b_FL,HC,e_kc,e_d, l_Solow EXP_R1 EXP_R2 EXP_R3 EXP_R4 RLR RLR_FL l_RLR l_RLR_FL RSR_MA RSR_FL_MA l_RSR_MA l_RSR_FL_MA R1R R1R_FL l_R1R l_R1R_FL R5R R5R_FL l_R5R l_R5R_FL NSR l_NSR NLR l_NLR N1R l_N1R N5R l_N5R;

//shock_decomposition RSR RSR_FL l_Output_FL l_pi_c;
shock_decomposition RSR_FL RSR;

stoch_simul(
order=1
,irf=40 // Default
,graph_format = (pdf)
//,nograph
) RSR_FL e_d e_aa NSR l_pi_c RLR RLR_FL;

save BJ_all_results;

else;
end;
// Estimation ends