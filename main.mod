%-------------------------------------------------------------------------------------------------%
% File Name --- main.mod
%
% README (written by Ryuichiro Hashimoto in April 2023)
% This code is for the joint research by MLIT and BOJ, which basically follows Okazaki-Sudo [2019].
% Main differences between this model and Okazaki-Sudo (and Hashimoto-Sudo) are:
%   1) This model introduces capital depreciation when the economy faces flood damage
%   2) This model also incorporates TFP level shocks as with Gallic-Vermandel [2019]
%   3) This model also explicitly introduces TFP shock by damage to public capital stock
%
% The structure of this code is:
%   1) main.mod (this file)     - controls estimation and simulation
%   2) model_desc.mod           - declares the model (Exo/Endogenous vars, SS, Equations, etc.)
%   3) model_estim.mod          - estimates model parameters and shocks given observables
%   4) model_simul.mod          - conducts stochastic simulation using estimated parameters in 3)
%
% Other useful code sections include:
%   - get_simul_output.m        - extracts simulated IRFs from 4)
%   - fdr_forecast.m            - forecasts endogenous variables up to 2050 given the bootstrapped
%                                   5000 fdr shock series
%-------------------------------------------------------------------------------------------------%


%-------------------------------------------------------------------------------------------------%
                                        % Model Setting
%-------------------------------------------------------------------------------------------------%
// shocks
//shock_list = {'nu_fdr','nu_aa','nu_ad','nu_ea','nu_ed','nu_nb','nu_nc','nu_gc','nu_kc','nu_pc','nu_wc', ...
//              'nu_d','nu_pi','nu_r','nu_ut'};
shock_list = {'nu_fdr'};
num_shocks = length(shock_list);

// global parameters
yearFrom        = 2023;
yearTo          = 2050;
nMLITscenario   = 1;
workdir         = pwd;
basedir         = fullfile(workdir, "..", "..");
runEstimation   = 0;
runStochSimul   = 1;
runMLITShock    = 0;
@#define withIns = 0
@#define bankLoss = 1
@#define noTFP   = 0

// set variables to extract irf results
variable_list = {'KKC','IC','C','Output','HC','Solow','QC','ZEC','NEC','NB','WC','LEC','KC','C_C','pi_c','NSR'};

// binary to control calibration source
loadOriginal    = 0;
loadPrevResult  = 1;

// set model type
@#if noTFP
  modelType = "noTFP";
@#elseif withIns
  modelType = "withIns";
@#elseif bankLoss
  modelType = "bankLoss";
@#else
  modelType = "BL";
@#endif

%-------------------------------------------------------------------------------------------------%
                                        % Model Declaration
%-------------------------------------------------------------------------------------------------%
@#include "model_desc.mod"


%-------------------------------------------------------------------------------------------------%
                                        % Shock Setting
%-------------------------------------------------------------------------------------------------%
shocks;
var nu_fdr;       stderr sigma_fdr;
var nu_aa;        stderr sigma_aa;
var nu_ad;        stderr sigma_ad;
var nu_ea;        stderr sigma_ea;
var nu_ed;        stderr sigma_ed;
var nu_nb;        stderr sigma_nb;
var nu_nc;        stderr sigma_nc;
var nu_gc;        stderr sigma_gc;
var nu_kc;        stderr sigma_kc;
var nu_pc;        stderr sigma_pc;
var nu_wc;        stderr sigma_wc;
var nu_d;         stderr sigma_d;
var nu_pi;        stderr sigma_pi;
var nu_r;         stderr sigma_r;
var nu_ut;        stderr sigma_ut;

end;

%-------------------------------------------------------------------------------------------------%
                                        % Steady-State Computation
%-------------------------------------------------------------------------------------------------%
resid;
steady(solve_algo = 3);
model_diagnostics;

check;

%-------------------------------------------------------------------------------------------------%
                                        % Estimation & Simulation
%-------------------------------------------------------------------------------------------------%
// Estimation (set 1 if you want to estimate parameters)
if runEstimation==1;
  @#include "model_estim.mod"
  run(fullfile(basedir, '11_Util', 'result_estIRF.m'));
  run(fullfile(basedir, '11_Util', 'write_results.m'));

end;

if runStochSimul==1;
  @#include "model_simul.mod"
  run(fullfile(basedir, '11_Util', 'result_simIRF.m'));

end;

%-------------------------------------------------------------------------------------------------%
                                        % MLIT simulations
%-------------------------------------------------------------------------------------------------%
// Simulations (shocks fed in to the model are created by MLIT)
if runMLITShock==1;

  run(fullfile(basedir, '11_Util', 'run_MLIT.m'));

end;

