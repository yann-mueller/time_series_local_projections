# Time Series: Local Projections
Local Projections offer a convenient way to directly estimate cumulative impulse response functions (IRF) with a model in differences in the background. In this code example, we use an exogenous shock series of supply chain disruptions to predict the effect on German industrial production following the framework of [Jordá (2005)](#references) and [Plagborg-Møller and Wolf (2021)](#references).
This is an outtake from our project on <a href="https://github.com/yann-mueller/waterway_shocks" target="_blank" rel="noopener noreferrer">Quantifying Waterway Supply Chain Shocks</a>.

## Statistical Background
The repository contains code from our baseline model, which is defined as follows:

![grafik](https://github.com/user-attachments/assets/646170a2-30dc-408d-86f5-9c9bf507b89b)

where we estimate the cumulative response of industrial production *ip*<sub>*t+h*</sub>​ to a supply chain shock *s*<sub>*t*</sub>​ at different forecast horizons *h=0,1,…,12*. The dependent variable is the cumulative change in industrial production from time *t* to *t+h*, denoted *Δ*<sup>*h+1*</sup>*ip*<sub>*t+h*</sub> *=* *ip*<sub>*t+h*</sub>​ *-* *ip*<sub>*t-1*</sub>. The shock *s*<sub>*t*</sub>​ enters contemporaneously with a horizon-specific coefficient *⁠β⁠*<sup>*h*</sup>, capturing the impulse response at horizon *h*. The model also controls for recent dynamics in industrial production by including lagged changes *Δip*<sub>*t-1*</sub>=*ip*<sub>*t-1*</sub> - *ip*<sub>*t-2*</sub> with coefficient *γ*<sub>*h*</sub>. The error term *ν*<sub>*t*</sub>​ captures unexplained variation.

This setup follows the framework proposed by [Jordá (2005)](#references), offering flexibility in estimating impulse responses without imposing the dynamic restrictions of a full vector autoregression (VAR).


### Autocorrelation
Since the regression is run separately for each horizon *h*, the error terms *ν*<sub>*t*</sub>​ can exhibit serial correlation across time. To ensure consistent estimation of the standard errors despite this, we compute our standard errors according to **Newey-West (HAC)**, which are robust to heteroskedasticity and autocorrelation. The lag length is chosen based on the forecast horizon *h*, following standard practice in the local projections literature.

## Matlab Code Explanations

### Configuration
This section sets up the key parameters and loads the necessary data for the impulse response function (IRF) analysis.

- **Forecast Horizon:** Defines the number of months (0–12) over which IRFs will be computed.
- **Time Series Selection:** Chooses five industrial production series (*n_spec = 5*) from the production data.
- **Time Range:** Filters the data between January 1991 and December 2019.
- **Data Preparation:** Extracts and logs the relevant time series to ensure stationarity.
- **Shock Series:** Loads the estimated structural shock series and calculates a scaling factor for IRF plotting.
- **Unit Root Testing (optional):** Applies Augmented Dickey-Fuller tests to each series to check for stationarity, using AIC to select lag length and printing the results in a formatted table.

```
% Define length of the predicted horizon
H=0:12;

% Define the number of time series for which the IRF shall be computed
n_spec=5;

% Define period of available time series
ind1 = find(strcmp('01-Jan-1991',cMonths));
ind2 = find(strcmp('01-Dec-2019',cMonths));
data=nan(ind2-ind1+1,n_spec);

% Read all data for which IRF shall be computed
data(:,1)=production{1,1}.ip_bcdf(ind1:ind2);
data(:,2)=production{1,1}.ip_energy_intensive(ind1:ind2);
data(:,3)=production{1,1}.ip_b(ind1:ind2);
data(:,4)=production{1,1}.ip_c(ind1:ind2);
data(:,5)=production{1,1}.ip_d(ind1:ind2);

% Estimating log-model
data = log(data);

% Read shock series
shock = shocks{1,1}.estimates(ind1:ind2, 1);
scaler = 100*sqrt(var(shock,'omitnan'));

% Testing for unit root (Optional)
for i=1:n_spec
    % Using a pre-defined sub-routine here
    adf_obj = adf_ic(data(:,i),'AIC',0,13,'ARD');
    ColNames = {'TestStatistic', 'pValue', 'LagOrder_FD', 'IC', 'Deterministics'};
    ADF_Table_FF = table(adf_obj.stat,adf_obj.pValue,...
    adf_obj.lag,adf_obj.ICs,adf_obj.deters, 'VariableNames',ColNames)
end
```

### Estimation
This section estimates the impulse responses of each variable to the shock using local projections.
- **IRF Initialization:** Creates matrices to store IRF estimates and their standard errors for each horizon and variable.
- **Local Projection Loop:**
    - Loops through each variable and horizon.
    - Computes a differenced regression of the outcome on lagged shocks and lags of the differenced outcome.
    - Uses heteroskedasticity- or autocorrelation-consistent standard errors (HAC), with bandwidth depending on the horizon.
    - Stores the estimated coefficient and standard error of the shock for each horizon.
```
% Initialize IRF and SE matrices
irf=nan(size(H,2),n_spec);
se=nan(size(H,2),n_spec);

% Compute impulse response for every horizon
for n=1:n_spec
    for i=0:size(H,2)-1
        s_selec=lagmatrix(shock,i);
        y=data(:,n)-lagmatrix(data(:,n),i+1);  
        X=[s_selec,lagmatrix(data(:,n),i+1)-lagmatrix(data(:,n),i+2)]; 
        
        if i==0
        [~,se_est,beta]=hac(X,y,Type='HC',Intercept=true,Display="off");
        else
        [~,se_est,beta]=hac(X,y,Type='HAC',Bandwidth=i,Intercept=true,Display="off");
        end

        % Store results
        irf(i+1,n)=beta(2);
        se(i+1,n)=se_est(2);
    end
end
```

### Plotting
This section generates the IRF plots with 90% confidence intervals for the selected variables.

- **Confidence Bands:** Constructs 90% confidence intervals using ±1.64 times the standard errors.
- **Plot Range:** Sets consistent y-axis limits across panels for visual comparability.
- **IRF Plots:** For each selected series:
    - Draws shaded areas for the confidence interval.
    - Plots the IRF (scaled and sign-flipped to analyze a negative shock).
    - Adds a zero reference line.
    - Uses tiledlayout for side-by-side comparison of variables like total industrial production and energy-intensive industries.
```
cVarname={'Industrial Production', '... Energy-Intensive Sectors Only'};

f1=tiledlayout(1,2,"TileSpacing","compact");

% Define confidence intervals (90% levels)
inBetween = [-1*scaler*(irf(:,:)'+1.64*se(:,:)'), fliplr(-1*scaler*(irf(:,:)'-1.64*se(:,:)'))];
min_lim = min(inBetween(1:2,1:13),[],"all") * 1.1;
max_lim = max(inBetween(1:2,2:26),[],"all") * 1.1;

for i=1:2
    nexttile
    x1 = [H, fliplr(H)];

    fill(x1, inBetween(i,:),[0, 0.7058, 0.8470] ,'FaceAlpha',0.2);
    hold on
    plot(H,(-1*scaler*irf(:,i)), 'LineWidth', 1,'Color',[0.0117, 0.0156, 0.3686])
    plot(H, zeros(size(H,2), 1),'k-');
    title(cVarname(i),'FontWeight','normal')
    xlabel('Horizon in months')
    xlim([0,12])
    ylim([min_lim max_lim])
    set(gca, 'FontSize', 9);
    grid on
    hold off
end
```

![grafik](https://github.com/user-attachments/assets/7a453024-b828-4db1-9d33-a2ecc920efe3)


## References
Jordá, O. (2005): “<a href="https://www.aeaweb.org/articles?id=10.1257/0002828053828518" target="_blank" rel="noopener noreferrer">Estimation and Inference of Impulse Responses by Local Projections</a>”, American Economic Review, 95, 161–182.

Newey, W. K. and K. D. West (1987): “<a href="https://scholar.google.de/scholar?hl=en&as_sdt=0%2C5&q=+Simple%2C+Positive+Semi-Definite%2C+Heteroskedas-+ticity+and+Autocorrelation+Consistent+Covariance+Matrix%2C%E2%80%9D+&btnG=" target="_blank" rel="noopener noreferrer">A Simple, Positive Semi-Definite, Heteroskedasticity and Autocorrelation Consistent Covariance Matrix</a>,” Econometrica, 55, 703–708.

Plagborg-Møller, M. and C. K. Wolf (2021): “<a href="https://onlinelibrary.wiley.com/doi/full/10.3982/ECTA17813" target="_blank" rel="noopener noreferrer">Local Projections and VARs Estimate the Same Impulse Responses</a>”, Econometrica, 89, 955–980.


