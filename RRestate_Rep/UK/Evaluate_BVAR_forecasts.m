% evaluates forecasts from 3 models
strFolder = 'D:\Dropbox\LRA_Matlab_2020\UK';
% Simulation folders:
dirnamefcs  = strcat(strFolder, '\BVAR\Regional\');


fdirnameNoLiq = strcat(dirnamefcs , 'NO_LIQ\');
fdirnameRtoV = strcat(dirnamefcs , 'RtoV\');
fdirnameiVOL= strcat(dirnamefcs , 'invV\'); 

minfile=0; maxfile=190; % IDs of first and last forecast data files
T0=T_thres;
horizons =[1 3 6 12]; % horizons to be examined
percentiles=[10 25 50 75 90]; % Density percentiles to be stored
probs=[0 -5 -10]; % thresholds X to store prob{variable<X}

% SET RAW DATA FILES
dat_NL=RETS; dat_RtV=RTOV; dat_IV=IVOL;

T=length(dat_NL); % this is the same for all models
N1=size(dat_NL,2); % VAR model with NO LIQ
N2=size(dat_RtV,2); % VAR models with liq measures
T2=maxfile+1-max(horizons); % T2 is the # of usable forecasts. If minfile=0 we need to add +1, because
% we have [forecast0 .. forecastK] = (K+1) forecasts.

NO_LIQ = struct;
NO_LIQ.pct     = NaN(T2, N1, length(horizons), length(percentiles));
NO_LIQ.pcr     = NaN(T2, N1, length(horizons), length(percentiles)); % Prob coverage ratios
NO_LIQ.pf      = NaN(T2, N1, length(horizons));
NO_LIQ.rmse    = NaN(T2, N1, length(horizons));
NO_LIQ.ls      = NaN(T2, N1, length(horizons));
NO_LIQ.jls     = NaN(T2, length(horizons)); % joint log-score (y, pi)
NO_LIQ.probs   = NaN(T2, length(horizons), length(probs)); 
NO_LIQ.wlsL    = NaN(T2, N1, length(horizons));
NO_LIQ.wlsLR   = NaN(T2, N1, length(horizons));
%NO_LIQ.pit     = NaN(T2, N1);
NO_LIQ.pit     = NaN(T2, N1, length(horizons));
NO_LIQ.pitin   = NaN(T2, N1, length(horizons));

RtV_LIQ = struct;
RtV_LIQ.pct     = NaN(T2, N2, length(horizons), length(percentiles));
RtV_LIQ.pcr     = NaN(T2, N2, length(horizons), length(percentiles)); % Prob coverage ratios
RtV_LIQ.pf      = NaN(T2, N2, length(horizons));
RtV_LIQ.rmse    = NaN(T2, N2, length(horizons));
RtV_LIQ.ls      = NaN(T2, N2, length(horizons));
RtV_LIQ.jls     = NaN(T2, length(horizons)); % joint log-score (y, pi)
RtV_LIQ.probs   = NaN(T2, length(horizons), length(probs)); 
RtV_LIQ.wlsL    = NaN(T2, N2, length(horizons));
RtV_LIQ.wlsLR   = NaN(T2, N2, length(horizons));
%RtV_LIQ.pit     = NaN(T2, N2);
RtV_LIQ.pit     = NaN(T2, N2, length(horizons));
RtV_LIQ.pitin   = NaN(T2, N2, length(horizons));

inV_LIQ = struct;
inV_LIQ.pct     = NaN(T2, N2, length(horizons), length(percentiles));
inV_LIQ.pcr     = NaN(T2, N2, length(horizons), length(percentiles)); % Prob coverage ratios
inV_LIQ.pf      = NaN(T2, N2, length(horizons));
inV_LIQ.rmse    = NaN(T2, N2, length(horizons));
inV_LIQ.ls      = NaN(T2, N2, length(horizons));
inV_LIQ.jls     = NaN(T2, length(horizons)); % joint log-score (y, pi)
inV_LIQ.probs   = NaN(T2, length(horizons), length(probs)); 
inV_LIQ.wlsL    = NaN(T2, N2, length(horizons));
inV_LIQ.wlsLR   = NaN(T2, N2, length(horizons));
%inV_LIQ.pit     = NaN(T2, N2);
inV_LIQ.pit     = NaN(T2, N2, length(horizons));
inV_LIQ.pitin   = NaN(T2, N2, length(horizons));

% PG: temporary - to be eliminated
% checkLS = NaN(T2, N, length(horizons)); 

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for file = 0 : T2-1
   
    % Load draws
    %----------------------------------------------------------------------
	fnameNO_LIQ =strcat(fdirnameNoLiq,'forecast',num2str(file),'.mat');
    fnameRtoV =strcat(fdirnameRtoV,'forecast',num2str(file),'.mat');
    fnameIV=strcat(fdirnameiVOL,'forecast',num2str(file),'.mat');
    
    load(fnameNO_LIQ);
    fnoliq = Y_fore;

    load(fnameRtoV);
    frtov = Y_fore;
    
    load(fnameIV);
    finvv = Y_fore;
        
	% Split data into estimation sample and observations for fc/density evaluation
    %----------------------------------------------------------------------
    sample_noliq = dat_NL(1:T0+file, :);
    sample_rtov = dat_RtV(1:T0+file, :);
    sample_invv= dat_IV(1:T0+file, :);   
    
    obs_noliq = dat_NL(T0+file+1 : T0+file+max(horizons), :) ;
    obs_RtV = dat_RtV(T0+file+1 : T0+file+max(horizons), :) ; 
    obs_IV= dat_IV(T0+file+1 : T0+file+max(horizons), :) ;     
    
    % Get weights for Amisano-Giacomini weighted log-scores, based on
    % estimation sample:
    wnoliq_L  = get_LS_weights(sample_noliq, obs_noliq, 'ltail');
    wnoliq_LR = get_LS_weights(sample_noliq, obs_noliq, 'tails');
    wrtov_L  = get_LS_weights(sample_rtov, obs_RtV, 'ltail');
    wrtov_LR = get_LS_weights(sample_rtov, obs_RtV, 'tails');
    winvv_L = get_LS_weights(sample_invv, obs_IV, 'ltail');
    winvv_LR= get_LS_weights(sample_invv, obs_IV, 'tails');
    % TEMP: if the data differ, tar (tvtp) weights should be calculated
    % using sampletar (sampletvtp). If the data are the same for all
    % models, the weights are the same too and we can save function calls.
    
    
	% Get cumulative observations and forecasts 
    %----------------------------------------------------------------------    
    frtovc = cumulate_data(frtov, 2, 3, 1:N2);
    fnoliqc = cumulate_data(fnoliq, 2, 3, 1:N1);
    finvvc= cumulate_data(finvv, 2, 3, 1:N2);
    
    obs_noliqc = cumulate_data(obs_noliq, 1, 2, 1:N1);
    obs_RtVc = cumulate_data(obs_RtV, 1, 2, 1:N2);
    obs_IVc= cumulate_data(obs_IV, 1, 2, 1:N2);   
    
    % Evaluate point forecasts
    %----------------------------------------------------------------------
    pfrtov   = squeeze(mean(frtov,1));
    pfnoliq   = squeeze(mean(fnoliq,1));
    pfinvv  = squeeze(mean(finvv,1));
    
    rmsenoliq = getrmse(pfnoliq-obs_noliq,horizons);
    rmsertov = getrmse(pfrtov-obs_RtV,horizons);
    rmseinvv= getrmse(pfinvv-obs_IV,horizons);
    
    rmsenoliq = reshape(rmsenoliq, 1, N1, length(horizons));
    rmsertov = reshape(rmsertov, 1, N2, length(horizons));
    rmseinvv= reshape(rmseinvv, 1, N2, length(horizons));
    % After reshape, each of these is a 1*N*horizons array where (1,:,j)
    % are the rmse for all variables at horizon horizons(h).
	    
    % PG 7.1.14: RMSE still calculated on non-cumulated data
    
 	% Get percentiles
    %----------------------------------------------------------------------
    pctlrtov = prctile(frtovc(:,horizons,:), percentiles, 1);
    pctlnoliq = prctile(fnoliqc(:,horizons,:), percentiles, 1);
    pctlinvv = prctile(finvvc(:,horizons,:), percentiles, 1);    
    
    pctlrtov = permute(pctlrtov, [3 2 1]); 
    pctlnoliq = permute(pctlnoliq, [3 2 1]); 
    pctlinvv= permute(pctlinvv, [3 2 1]);     
    % This rearranges the dimensions as variables*horizons*percentiles
    % (purely for storage, see below)
    
    pfrtov0 = reshape(pfrtov(horizons,:)', [1 N2 length(horizons)]);
    pfnoliq0 = reshape(pfnoliq(horizons,:)', [1 N1 length(horizons)]);
    pfinvv0= reshape(pfinvv(horizons,:)',[1 N2 length(horizons)]);    
    % This rearranges the dimensions as 1*variables*horizons (purely for 
    % storage, see below)
    
    % Get Probs{output<X} COME BACK TO THIS
    %----------------------------------------------------------------------
    % These are stored directly to save one loop
    % for hh = 1:length(horizons)
        
    %   temp1 = ksdensity(frtovc(:,horizons(hh),1), probabilities, 'function', 'cdf');
    %   temp2 = ksdensity(fnoliqc(:,horizons(hh),1), probabilities, 'function', 'cdf');
    %   temp3 = ksdensity(finvvc(:,horizons(hh),1), probabilities, 'function', 'cdf');
       
    %   NO_LIQ.probs(file+1, hh, :)  = temp1;
    %   RtV_LIQ.probs(file+1, hh, :)  = temp2;
    %   inV_LIQ.probs(file+1, hh, :) = temp3;
    %end
    
    % Evaluate densities
    %----------------------------------------------------------------------
    for k = 1:N1

        scoreNOLIQ(1, k, 1:length(horizons)) = get_logscoreH( squeeze(fnoliqc(:,:,k)), obs_noliqc(:,k), horizons);
        % PG 7.1.14: all LS are now calculated using *c data and forecasts. 
               
        % PITs 
        for hh=1:length(horizons)
        pitsNOLIQ(k, hh) = ksdensity(squeeze(fnoliq(:,hh,k)), obs_noliq(hh,k), 'function', 'cdf') ;       
        % PG 7.1.14: for 1-period ahead we can indifferently use normal (*) 
        % or cumulative (*c) arrays because they coincide.  
        end 
        
       
        
  
    end
    for k=1:N2
        scoreRTOV(1, k, 1:length(horizons)) = get_logscoreH(squeeze(frtovc(:,:,k)), obs_RtVc(:,k), horizons);
        scoreIVOL(1, k, 1:length(horizons))= get_logscoreH(squeeze(finvvc(:,:,k)),obs_IVc(:,k), horizons);
        for hh=1:length(horizons)
        pitsRTOV(k, hh) = ksdensity(squeeze(frtov(:,hh,k)), obs_RtV(hh,k), 'function', 'cdf') ;
        pitsIVOL(k, hh)= ksdensity(squeeze(finvv(:,hh,k)),obs_IV(hh,k), 'function', 'cdf') ; 
        end
    end
  
% for hh=1:length(horizons)
%     kk=horizons(hh);
%     for k=1:N1
%       [scoreNOLIQ(1,k,hh),pitsNOLIQ(1,k,hh)]=lnsc(squeeze(fnoliqc(:,kk,k)),obs_noliqc(kk,k));     
%    end
%    for k=1:N2
%       [scoreRTOV(1,k,hh),pitsRTOV(1,k,hh)]=lnsc(squeeze(frtovc(:,kk,k)),obs_RtVc(kk,k));   
%       [scoreIVOL(1,k,hh),pitsIVOL(1,k,hh),grid]=lnsc(squeeze(finvvc(:,kk,k)),obs_IVc(kk,k));  
%    end  
% end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COME BACK TO THIS
    
    %jscoreVAR=get_jls(fnoliqc,obs_noliqc,horizons);
    
    % Joint log-score for (y, pi), i.e. variables #1 and #3 all horizons:
    %jscoreVAR = get_jointlogscore(squeeze(frtovc(:,:,1)), squeeze(frtovc(:,:,3)), obs_noliqc(:,1), obs_noliqc(:,3), horizons);
    %jscoreTAR = get_jointlogscore(squeeze(fnoliqc(:,:,1)), squeeze(fnoliqc(:,:,3)), obs_RtVc(:,1), obs_RtVc(:,3), horizons);
    %jscoreTVTP= get_jointlogscore(squeeze(finvvc(:,:,1)), squeeze(finvvc(:,:,3)),obs_IVc(:,1), obs_IVc(:,3), horizons);
    
    % DEBUG: break if scores go imaginary or NaN
    %if sum(~isreal([jscoreVAR jscoreTAR jscoreTVTP])) | sum(isnan([jscoreVAR jscoreTAR jscoreTVTP]))
    %    return
    %end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Calculate indicators for prob coverage ratios:
    tempdata = repmat(obs_noliqc(horizons,:)', [1 1 length(percentiles)]); 
    pcrnoliq = (tempdata < pctlnoliq);
    tempdata = repmat(obs_RtVc(horizons,:)', [1 1 length(percentiles)]);
    pcrrtov = (tempdata < pctlrtov);
    tempdata = repmat(obs_IVc(horizons,:)', [1 1 length(percentiles)]);
    pcrinvv= (tempdata < pctlinvv);
    % pcr* is a variables*horizons*percentiles matrix of logical values
    % with entry (v,h,p) = 1 if variable v at horizon h is below percentile
    % p and zero otherwise.
    % NOTE: assumes that all models are evaluated on the same data.
    
    % Weighted log-scores
    for hh=1:length(horizons)
    wlsnoliq_L(hh,:)  = squeeze(scoreNOLIQ(1,:,hh)) .* wnoliq_L(hh,:);
    wlsnoliq_LR(hh,:) = squeeze(scoreNOLIQ(1,:,hh)) .* wnoliq_LR(hh,:);
    wlsrtov_L(hh,:)  = squeeze(scoreRTOV(1,:,hh)) .* wrtov_L(hh,:);
    wlsrtov_LR(hh,:) = squeeze(scoreRTOV(1,:,hh)) .* wrtov_LR(hh,:);
    wlsinvv_L(hh,:) = squeeze(scoreIVOL(1,:,hh)).* winvv_L(hh,:);
    wlsinvv_LR(hh,:)= squeeze(scoreIVOL(1,:,hh)).* winvv_LR(hh,:);
    end
    
    
    % Save results
    %----------------------------------------------------------------------     
    NO_LIQ.pct(file+1,:,:,:)  = pctlnoliq;
    NO_LIQ.pcr(file+1,:,:,:)  = pcrnoliq;
    NO_LIQ.pf(file+1,:,:)     = pfnoliq0;
    NO_LIQ.rmse(file+1, :, :) = rmsenoliq;
    NO_LIQ.ls(file+1, :, :)   = scoreNOLIQ;
    %NO_LIQ.jls(file+1, :)     = jscoreVAR; 
    NO_LIQ.pit(file+1, :, :)     = pitsNOLIQ;
    NO_LIQ.wlsL(file+1, :, :)    = wlsnoliq_L';
    NO_LIQ.wlsLR(file+1, :, :)   = wlsnoliq_LR';
    
    RtV_LIQ.pct(file+1,:,:,:)  = pctlrtov;
	RtV_LIQ.pcr(file+1,:,:,:)  = pcrrtov;
    RtV_LIQ.pf(file+1,:,:)     = pfrtov0;
	RtV_LIQ.rmse(file+1, :, :) = rmsertov;
    RtV_LIQ.ls(file+1, :, :)   = scoreRTOV;
	% RtV_LIQ.jls(file+1, :)     = jscoreTAR; 
    RtV_LIQ.pit(file+1, :, :)     = pitsRTOV;
    RtV_LIQ.wlsL(file+1, :, :)    = wlsrtov_L';
    RtV_LIQ.wlsLR(file+1, :, :)   = wlsrtov_LR';

	inV_LIQ.pct(file+1,:,:,:)  = pctlinvv;
	inV_LIQ.pcr(file+1,:,:,:)  = pcrinvv;
    inV_LIQ.pf(file+1,:,:)     = pfinvv0;
	inV_LIQ.rmse(file+1, :, :) = rmseinvv;
    inV_LIQ.ls(file+1, :, :)   = scoreIVOL;
	% inV_LIQ.jls(file+1, :)     = jscoreTVTP; 
    inV_LIQ.pit(file+1, :, :)     = pitsIVOL;
    inV_LIQ.wlsL(file+1, :, :)    = wlsinvv_L';
    inV_LIQ.wlsLR(file+1, :, :)   = wlsinvv_LR';
    
    disp(file);
end

%% Output manipulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% De-annualize RMSE for IP growth and inflation (NO_LIQ.s 1 and 3)
%NO_LIQ.rmse(:, [1 3], :) = NO_LIQ.rmse(:, [1 3], :) / 12;
%RtV_LIQ.rmse(:, [1 3], :) = RtV_LIQ.rmse(:, [1 3], :) / 12;
%inV_LIQ.rmse(:, [1 3], :)= inV_LIQ.rmse(:, [1 3], :)/ 12;

% Get coverage ratios:
NO_LIQ.pcr = 100*squeeze(mean(NO_LIQ.pcr, 1));
RtV_LIQ.pcr = 100*squeeze(mean(RtV_LIQ.pcr, 1));
inV_LIQ.pcr= 100*squeeze(mean(inV_LIQ.pcr, 1));
% Note we get rid of the time dimesion and keep a 
%   (#horizons)*(#variables)*(#percentiles) array for each model.

% Inverse-N PITs:
NO_LIQ.pitin = icdf('Normal', NO_LIQ.pit, 0, 1);
RtV_LIQ.pitin = icdf('Normal', RtV_LIQ.pit, 0, 1);
inV_LIQ.pitin= icdf('Normal', inV_LIQ.pit, 0, 1);

% Sample averages (all in .M ):
% -------------------------------------------------------------------------
NO_LIQ.M.rmse  = squeeze(mean(NO_LIQ.rmse(:, :, 1:end)))';
NO_LIQ.M.ls    = squeeze(mean(NO_LIQ.ls(:, :, 1:end)))';
NO_LIQ.M.jls   = mean(NO_LIQ.jls, 1);
NO_LIQ.M.wlsL  = squeeze(mean(NO_LIQ.wlsL, 1));
NO_LIQ.M.wlsLR = squeeze(mean(NO_LIQ.wlsLR, 1));

RtV_LIQ.M.rmse  = squeeze(mean(RtV_LIQ.rmse(:, :, 1:end)))';
RtV_LIQ.M.ls    = squeeze(mean(RtV_LIQ.ls(:, :, 1:end)))';
RtV_LIQ.M.jls   = mean(RtV_LIQ.jls, 1);
RtV_LIQ.M.wlsL  = squeeze(mean(RtV_LIQ.wlsL, 1));
RtV_LIQ.M.wlsLR = squeeze(mean(RtV_LIQ.wlsLR, 1));

inV_LIQ.M.rmse  = squeeze(mean(inV_LIQ.rmse(:, :, 1:end)))';
inV_LIQ.M.ls    = squeeze(mean(inV_LIQ.ls(:, :, 1:end)))';
inV_LIQ.M.jls   = mean(inV_LIQ.jls, 1);
inV_LIQ.M.wlsL  = squeeze(mean(inV_LIQ.wlsL, 1));
inV_LIQ.M.wlsLR = squeeze(mean(inV_LIQ.wlsLR, 1));
%%

% Giacomini-White pvalues:
% -------------------------------------------------------------------------
GWpv      = struct;

warning off  

% UNCONDITIONAL TESTS (.UC):

% tar viz var
GWpv.UC.RtV_LIQ.rmse = get_GW_pvalues(RtV_LIQ.rmse(:,1:4,:)  - NO_LIQ.rmse(:,1:4,:) , 'uncond', horizons);
GWpv.UC.RtV_LIQ.ls   = get_GW_pvalues(RtV_LIQ.ls(:,1:4,:)    - NO_LIQ.ls(:,1:4,:), 'uncond', horizons);
GWpv.UC.RtV_LIQ.wlsL = get_GW_pvalues(RtV_LIQ.wlsL(:,1:4,:)  - NO_LIQ.wlsL(:,1:4,:) , 'uncond', horizons);
GWpv.UC.RtV_LIQ.wlsLR= get_GW_pvalues(RtV_LIQ.wlsLR(:,1:4,:) - NO_LIQ.wlsLR(:,1:4,:), 'uncond', horizons);

% tvtp viz var
GWpv.UC.inV_LIQ.rmse = get_GW_pvalues(inV_LIQ.rmse(:,1:4,:) - NO_LIQ.rmse(:,1:4,:), 'uncond', horizons);
GWpv.UC.inV_LIQ.ls   = get_GW_pvalues(inV_LIQ.ls(:,1:4,:)   - NO_LIQ.ls(:,1:4,:), 'uncond', horizons);
GWpv.UC.inV_LIQ.wlsL = get_GW_pvalues(inV_LIQ.wlsL(:,1:4,:) - NO_LIQ.wlsL(:,1:4,:) , 'uncond', horizons);
GWpv.UC.inV_LIQ.wlsLR= get_GW_pvalues(inV_LIQ.wlsLR(:,1:4,:) - NO_LIQ.wlsLR(:,1:4,:), 'uncond', horizons);

% tvtp viz tar
GWpv.UC.inV_LIQ2.rmse = get_GW_pvalues(inV_LIQ.rmse(:,1:4,:) - RtV_LIQ.rmse(:,1:4,:), 'uncond', horizons);
GWpv.UC.inV_LIQ2.ls   = get_GW_pvalues(inV_LIQ.ls(:,1:4,:)   - RtV_LIQ.ls(:,1:4,:), 'uncond', horizons);
GWpv.UC.inV_LIQ2.wlsL = get_GW_pvalues(inV_LIQ.wlsL(:,1:4,:) - RtV_LIQ.wlsL(:,1:4,:) , 'uncond', horizons);
GWpv.UC.inV_LIQ2.wlsLR= get_GW_pvalues(inV_LIQ.wlsLR(:,1:4,:)- RtV_LIQ.wlsLR(:,1:4,:), 'uncond', horizons);


% CONDITIONAL TESTS (.C):

% tar viz var
GWpv.C.RtV_LIQ.rmse = get_GW_pvalues(RtV_LIQ.rmse(:,1:4,:)   - NO_LIQ.rmse(:,1:4,:), 'cond', horizons);
GWpv.C.RtV_LIQ.ls   = get_GW_pvalues(RtV_LIQ.ls(:,1:4,:)     - NO_LIQ.ls(:,1:4,:), 'cond', horizons);
GWpv.C.RtV_LIQ.wlsL = get_GW_pvalues(RtV_LIQ.wlsL(:,1:4,:)  - NO_LIQ.wlsL(:,1:4,:), 'cond', horizons);
GWpv.C.RtV_LIQ.wlsLR= get_GW_pvalues(RtV_LIQ.wlsLR(:,1:4,:)  - NO_LIQ.wlsLR(:,1:4,:), 'cond', horizons);

% tvtp viz var
GWpv.C.inV_LIQ.rmse = get_GW_pvalues(inV_LIQ.rmse(:,1:4,:) - NO_LIQ.rmse(:,1:4,:), 'cond', horizons);
GWpv.C.inV_LIQ.ls   = get_GW_pvalues(inV_LIQ.ls(:,1:4,:)   - NO_LIQ.ls(:,1:4,:), 'cond', horizons);
GWpv.C.inV_LIQ.wlsL = get_GW_pvalues(inV_LIQ.wlsL(:,1:4,:) - NO_LIQ.wlsL(:,1:4,:), 'cond', horizons);
GWpv.C.inV_LIQ.wlsLR= get_GW_pvalues(inV_LIQ.wlsLR(:,1:4,:)- NO_LIQ.wlsLR(:,1:4,:), 'cond', horizons);

% tvtp viz tar
GWpv.C.inV_LIQ2.rmse = get_GW_pvalues(inV_LIQ.rmse(:,1:4,:) - RtV_LIQ.rmse(:,1:4,:), 'cond', horizons);
GWpv.C.inV_LIQ2.ls   = get_GW_pvalues(inV_LIQ.ls(:,1:4,:)  - RtV_LIQ.ls(:,1:4,:), 'cond', horizons);
GWpv.C.inV_LIQ2.wlsL = get_GW_pvalues(inV_LIQ.wlsL(:,1:4,:) - RtV_LIQ.wlsL(:,1:4,:), 'cond', horizons);
GWpv.C.inV_LIQ2.wlsLR= get_GW_pvalues(inV_LIQ.wlsLR(:,1:4,:)- RtV_LIQ.wlsLR(:,1:4,:), 'cond', horizons);

% warning on 

% Giacomini-White decision rules:
% -------------------------------------------------------------------------
GWdr = struct;

GWdr.rmse =  get_GW_decision(NO_LIQ.rmse(:,1:4,:) - RtV_LIQ.rmse(:,1:4,:), horizons);
GWdr.ls   = -get_GW_decision(NO_LIQ.ls(:,1:4,:) - RtV_LIQ.ls(:,1:4,:), horizons);
GWdr.jls  = -get_GW_decision(NO_LIQ.jls(:,1:4,:) - RtV_LIQ.jls(:,1:4,:), horizons);
GWdr.wlsL = -get_GW_decision(NO_LIQ.wlsL(:,1:4,:) - RtV_LIQ.wlsL(:,1:4,:), horizons);
GWdr.wlsLR= -get_GW_decision(NO_LIQ.wlsLR(:,1:4,:) - RtV_LIQ.wlsLR(:,1:4,:), horizons);

GWdr.rmse0 =  get_GW_decision(NO_LIQ.rmse(:,1:4,:) -inV_LIQ.rmse(:,1:4,:), horizons);
GWdr.ls0   = -get_GW_decision(NO_LIQ.ls(:,1:4,:) - inV_LIQ.ls(:,1:4,:), horizons);
GWdr.jls0  = -get_GW_decision(NO_LIQ.jls(:,1:4,:) - inV_LIQ.jls(:,1:4,:), horizons);
GWdr.wlsL0 = -get_GW_decision(NO_LIQ.wlsL(:,1:4,:) - inV_LIQ.wlsL(:,1:4,:), horizons);
GWdr.wlsLR0= -get_GW_decision(NO_LIQ.wlsLR(:,1:4,:) - inV_LIQ.wlsLR(:,1:4,:), horizons);


GWdr.ls1   = -get_GW_decision(inV_LIQ.ls(:,1:4,:) - RtV_LIQ.ls(:,1:4,:), horizons);
GWdr.wlsL1 = -get_GW_decision(inV_LIQ.wlsL(:,1:4,:) - RtV_LIQ.wlsL(:,1:4,:), horizons);
GWdr.wlsLR1= -get_GW_decision(inV_LIQ.wlsLR(:,1:4,:) - RtV_LIQ.wlsLR(:,1:4,:), horizons);
% NB: Signs are set up so that for both RMSE and LS:
%         "criterion>0" <=> "pick model #2" 
% So positives mean that TAR beats VAR in .crit and VAR beats TVTP in .crit0 



