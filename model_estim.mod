%-------------------------------------------------------------------------------------------------%
% File Name --- model_estim.mod
%
% README (written by Ryuichiro Hashimoto, Oct. 2021)
% This code section estimates the model parameters.
% The structure of this section is:
%   1) Observed variables
%   2) Parameters to be estimated
%   3) Initial values
%   4) Estimation
%   5) Simulation (stochastic simulation using the estimated parameters)
%
% This section uses the following files/codes:
%   1) model_estim_obs.m            - loads data for observables
%   2) main_mode.mat                - posterior mode from the previous estimation
%
% Note:
%   - set mode_compute = 0 when using the mode file, otherwise whatever appropriate value based on
%       the Dynare documentation. For instance, set mode_compute = 6 when using Monte-Carlo based
%       optimization routine.
%-------------------------------------------------------------------------------------------------%

%-------------------------------------------------------------------------------------------------%
                                        % Observed variables
%-------------------------------------------------------------------------------------------------%
varobs

l_pi_c
l_IC
l_pi_d
l_Output
l_HC
l_WC
l_NSR
l_NB
l_NEC
l_Solow
fdr
edr
//fdr_p
//edr_p
;


%-------------------------------------------------------------------------------------------------%
                                        % Parameteres to be estimated
%-------------------------------------------------------------------------------------------------%
estimated_params;

vv,             gamma_pdf,      0.8,        0.075;
kappaC,         gamma_pdf,      2,          0.25;
kappa_pc,       gamma_pdf,      12,         1;
kappa_wc,       gamma_pdf,      2.5,        0.5;
rpi,            normal_pdf,     2.75,       0.05;
phi_u,		      gamma_pdf,      5,          1;
rho_NSR,        beta_pdf,       0.5,        0.01;
sigmaUEC_SS,    gamma_pdf,      0.3093,     0.002;
sigmaUB_SS,	    gamma_pdf,      0.1042,     0.002;
mu_ec,		      gamma_pdf,      0.0196,     0.01;
mu_b,           gamma_pdf,      0.5386,     0.01;
gammaec,        beta_pdf,       0.96,       0.001;
gammab,         beta_pdf,       0.86,       0.001;
rho_aa,         beta_pdf,       0.5,        0.15;
rho_ad,         beta_pdf,       0.5,        0.15;
rho_ea,	        beta_pdf,       0.5,        0.15;
rho_ed,	        beta_pdf,       0.5,        0.15;
rho_nb,	        beta_pdf,       0.5,        0.15;
rho_nc,	        beta_pdf,       0.5,        0.15;
rho_gc,	        beta_pdf,       0.5,        0.15;
rho_kc,	        beta_pdf,       0.5,        0.15;
rho_pc,	        beta_pdf,       0.5,        0.15;
rho_wc,    	    beta_pdf,       0.5,        0.15;
rho_d,		      beta_pdf,       0.5,        0.15;
rho_pi,         beta_pdf,       0.5,        0.15;
rho_ut,         beta_pdf,       0.5,        0.15;
eta_a_SS,       gamma_pdf,      1.0011,	    0.001;
eta_d_SS,       gamma_pdf,      1.0016,	    0.001;
pi_c_SS,        normal_pdf,     1.0025,     0.001;
stderr nu_aa,   inv_gamma_pdf,  0.05,       5;
stderr nu_ad,   inv_gamma_pdf,  0.05,       5;
stderr nu_ea,   inv_gamma_pdf,  0.01,       5;
stderr nu_ed,   inv_gamma_pdf,  0.01,       5;
stderr nu_r,    inv_gamma_pdf,  0.01,       5;
stderr nu_nb,   inv_gamma_pdf,  0.02,       5;
stderr nu_nc,   inv_gamma_pdf,  0.02,       5;
stderr nu_gc,   inv_gamma_pdf,  0.01,       5;
stderr nu_kc,   inv_gamma_pdf,  0.01,       5;
stderr nu_pc,   inv_gamma_pdf,  0.01,       5;
stderr nu_wc,   inv_gamma_pdf,  0.01,       5;
stderr nu_d,    inv_gamma_pdf,  0.015,      5;
stderr nu_pi,   inv_gamma_pdf,  0.01,       5;
stderr nu_ut,   inv_gamma_pdf,  0.01,       5;
stderr l_NB,    inv_gamma_pdf,  0.1,        5;
stderr l_NEC,   inv_gamma_pdf,  0.1,        5;
habit,          beta_pdf,       0.6,        0.15;

rho_fdr,        beta_pdf,       0.5,        0.15;
rho_edr,	      beta_pdf,       0.5,        0.15;
//rho_fdr_p,      beta_pdf,       0.5,        0.15;
//rho_edr_p,      beta_pdf,       0.5,        0.15;
theta_fdr,      beta_pdf,       0.5,        0.15;
//theta_fdr_p,    beta_pdf,       0.5,        0.15;
theta_edr,      beta_pdf,       0.5,        0.15;
//theta_edr_p,    beta_pdf,       0.5,        0.15;
stderr nu_fdr,  inv_gamma_pdf,  0.01,     5 ;
stderr nu_edr,  inv_gamma_pdf,  0.01,     5;
//stderr nu_fdr_p,  inv_gamma_pdf,  0.0001,     0.005;
//stderr nu_edr_p,  inv_gamma_pdf,  0.0009,     0.005;

end;

%-------------------------------------------------------------------------------------------------%
                                        % Initial values
%-------------------------------------------------------------------------------------------------%
estimated_params_init;

vv,             0.9210;
kappaC,         2.0457;
kappa_pc,       6.7397;
kappa_wc,       1.0621;
rpi,            2.8280;
phi_u,          6.1398;
rho_NSR,        0.5163;
mu_ec,          0.0314;
mu_b,           0.5350;
rho_aa,         0.9388;
rho_ad,         0.7755;
rho_nb,         0.2421;
rho_nc,         0.7772;
rho_ea,         0.1330;
rho_gc,         0.9356;
rho_kc,         0.2287;
rho_pc,         0.8319;
rho_wc,         0.8617;
rho_ed,         0.6199;
rho_pi,         0.4326;
rho_d,          0.8107;
habit,          0.136;
stderr nu_aa,   0.0059;
stderr nu_ad,   0.0059;
stderr nu_r,    0.0015;
stderr nu_nb,   0.032;
stderr nu_nc,   0.051;
stderr nu_ea,   0.0023;
stderr nu_gc,   0.0466;
stderr nu_kc,   0.0302;
stderr nu_pc,   0.0277;
stderr nu_wc,   0.1010;
stderr nu_ed,   0.0024;
stderr nu_pi,   0.0038;
stderr nu_d,    0.0043;


rho_fdr,        0.5;
//rho_fdr_p,      0.5;
theta_fdr,      0.5;
//theta_fdr_p,    0.5;
//stderr nu_fdr,  0.01;
rho_edr,        0.5;
//rho_edr_p,      0.5;
theta_edr,      0.5;
//theta_edr_p,    0.5;
//stderr nu_edr,  0.01;

end;

%-------------------------------------------------------------------------------------------------%
                                        % Estimation
%-------------------------------------------------------------------------------------------------%
estimation(
  //optim=('Display','iter'),
  datafile='../../11_Util/model_setEstObs.m',
  nobs=156,
  order=1,
  mode_compute=6,  // when computing the posterior mode
  //mode_compute=0,
  //mode_file=main_mode,
  mh_replic=200000,
  mh_nblocks=2,
  mh_jscale=0.34,
  mh_init_scale=0.01,
  smoother,
  kalman_algo=2,
  filtered_vars,
  bayesian_irf,
  irf=121,
  forecast=121,
  mode_check,
  graph_format = (pdf),
  nodisplay,
  consider_all_endogenous
);

save output/estim_AllResults.mat;
