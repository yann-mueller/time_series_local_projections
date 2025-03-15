# Time Series: Local Projections
Local Projections offer a convenient way to directly estimate cumulative impulse response functions (IRF) with a model in differences in the background. In this code example, we use an exogenous shock series of supply chain disruptions to predict the effect on German industrial production following the framework of Jordá (2005) and Plagborg-Møller and Wolf (2021).
This is an outtake from our project on <a href="https://github.com/yann-mueller/waterway_shocks" target="_blank" rel="noopener noreferrer">Quantifying Waterway Supply Chain Shocks</a>.

## Statistical Background
The model is defined as follows:

![grafik](https://github.com/user-attachments/assets/646170a2-30dc-408d-86f5-9c9bf507b89b)


### Autocorrelation

## Matlab Code Explanations

### Configuration
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

## References
Jordá, O. (2005): “<a href="https://www.aeaweb.org/articles?id=10.1257/0002828053828518" target="_blank" rel="noopener noreferrer">Estimation and Inference of Impulse Responses by Local Projections</a>”, American Economic Review, 95, 161–182.

Plagborg-Møller, M. and C. K. Wolf (2021): “<a href="https://onlinelibrary.wiley.com/doi/full/10.3982/ECTA17813" target="_blank" rel="noopener noreferrer">Local Projections and VARs Estimate the Same Impulse Responses</a>”, Econometrica, 89, 955–980.


