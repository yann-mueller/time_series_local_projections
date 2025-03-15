function adf_obj=adf_ic(y,IC,lag_ADF,pmaxFD,deter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Function for ADF test %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% With pre-specified lag order or lag order chosen via IC %%%%%%%%
%%%% With Deterministic terms: none, constant only, consant and trend %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INPUT: 
%	y   ( T x 1 )       time series
%   IC   ( 1 x 1 )      none or information criterion (char array)
%   lag_ADF ( 1 x 1)    pre-determined lag order in case IC = 'None'
%   pmaxFD (1 x 1)      maximum number of lags in ADF test regression
%                       for chosen IC; refers to first differences
%   deter (1 x 1)       specifies deterministics (char array)


% Outputs: 
%	adf_obj        
%         .stat         ADF-t statistic
%         .pValue       p-value wrt to ADF-t statistic
%         .lag          chosen lag in test regression
%         .ICs          string: None or chosen IC
%         .deters       string: determinsitics


pmax = pmaxFD + 1; % Maximum number of lags wrt levels

T_adf = size(y,1)-pmax; % Sample size of ADF test regression

if strcmp(IC,'None')
    
    % pre-determined lag order lag_ADF will be taken
    ICs = {'None'}; % just used for output
    
else
    % dependent variable in levels (corresponds to Matlab implementation of ADF
    % test regression
    yd = y(pmax+1:end); 

    % Matrix containing all lags of FD of un
    Dy_lagm = lagmatrix(diff(y),1:pmaxFD); 
    % Note: lagmatrix produces NaN-entries but keeps number of obs. of input
    % series

    % Regressormatrix: deterministics, 1st lag of y, pmax lags of FD of y
    if strcmp(deter,'AR')
        xm = [y(pmax:end-1), Dy_lagm(pmax:end,:)];
    elseif strcmp(deter,'ARD')
        xm = [ones(T_adf,1), y(pmax:end-1), Dy_lagm(pmax:end,:)];
    else
        xm = [ones(T_adf,1), (1:1:T_adf)', y(pmax:end-1), Dy_lagm(pmax:end,:)];
    end 
       

    aic = NaN(pmax,1);
    hq = NaN(pmax,1);
    bic = NaN(pmax,1);

    % Loop: ADF test regression with zero lags up to pmax-1 lags wrt first
    % differences of y; computation of information criteria

    for lag = 0:pmax-1
    
        % ADF test regression with corresponding number of lags of FD of un
        [b,bint,res] = regress(yd,xm(:,1:lag+2));  
        aic(lag+1) = log(res'*res/(T_adf))+(lag+1)*2/(T_adf);
        hq(lag+1) = log(res'*res/(T_adf))+(lag+1)*2*log(log(T_adf))/(T_adf);
        bic(lag+1) = log(res'*res/(T_adf))+(lag+1)*log(T_adf)/(T_adf);   

    end

    % Determining "optimal lag order"
    if strcmp(IC,'AIC')
        [v,MI] = min(aic);
        % Note: MI= row index with smallest value = optimial lag order + 1
        lag_ADF = MI-1;
        ICs = {'AIC'}; % just for output
    elseif strcmp(IC,'HQ')
        [v,MI] = min(hq);
        lag_ADF = MI-1;
        ICs = {'HQ'};
    else
        [v,MI] = min(bic);
        lag_ADF = MI-1;
        ICs = {'BIC'};
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADF test with lag input based on chosen information criterion
[h,pValue,stat,cValue,reg] = adftest(y,'model',deter,'lags',lag_ADF);

% Just for output

if strcmp(deter,'AR')
    deters = {'none'};
elseif strcmp(deter,'ARD')
    deters = {'constant'};
else
    deters = {'trend'};
end

% Generating output object adf_obj and filling it with the relevant results
% Note 1: adf_adj would be called a 'structure' in Matlab
% Note 2: the attachments like '.stat' are called 'field' in Matlab

adf_obj.stat = stat;
adf_obj.pValue = pValue;
adf_obj.lag = lag_ADF;
adf_obj.ICs = ICs;
adf_obj.deters = deters;
