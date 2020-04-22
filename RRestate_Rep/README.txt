README FILE FOR REGIONAL REAL ESTATE ILLIQDUITIY REPLCIATION 

CODE COMES WITHOUT TECHNICAL SUPPORT OF ANY KIND USE AT YOUR OWN RISK

REQUIRES PARALLEL COMPUTING TOOLBOX IN MATLAB.

IF YOU USE THIS CODE PLEASE CITE THE PAPER 


YOU WILL NEED TO SET YOUR OWN PATHS WHERE THE CODE STATES <INSERT PATH HERE>

FOR US AND UK DATA FORECASTING RESULTS:

1. Run TVP_VAR_FORECASTING_DENSITY.m
	This will provide results provided in paper and in Supplementary Appendix pertaining to forecasting

FOR US AND UK DATA IRFs, FEVDs, and NETWORK CONNECTEDNESS

1. Run TVP_VAR_MASTER.m
2. Run TVP_results_irf_fev_network.m for regional results
3. Run TVP_results_N_irf_fev_network.m for national results

You can also get forecasting results from Bayesian VARs using a Minnesota NW prior

1. Run BVAR_FORECASTING_MASTER.m