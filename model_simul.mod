%-------------------------------------------------------------------------------------------------%
% File Name --- model_simul.mod
%
% README (written by Ryuichiro Hashimoto, Oct. 2021)
% This code section simulates the model using the calibrated parameters.
% The structure of this section is:
%   1) Calibration
%   2) Simulation
%-------------------------------------------------------------------------------------------------%







%-------------------------------------------------------------------------------------------------%
                                        % Simulation
%-------------------------------------------------------------------------------------------------%
stoch_simul(
  order=1,
  irf=41, // Default
  graph_format = (pdf),
  conditional_variance_decomposition = [1:41],
  nocorr,
  nodecomposition,
  nomoments,
  nofunctions,
  nograph
);

save output/simul_AllResult.mat;
