
load('station_data_yr_525');
%% calculate and plot linear trend
    prcp_regress_coeff_temp = [];
    prcp_regress_coeff = [];
    for x=1:525  % if compare full period,use commentted two lines to replace the line 'extr_regress_coeff_temp...'
        ind = 1:length(station_data_yr(:,x));
        k = ~isnan(station_data_yr(:,x)); % find index with non-NaN values
        %     y = find(k==1);
        %     pref_regress_coeff_temp = polyfit(ind(k)',36.5246.*station_data_yr(78+y,x+4),1);
        prcp_regress_coeff_temp = polyfit(ind(k)',0.1*station_data_yr(k,x),1);
        prcp_regress_coeff(x) = prcp_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
    end
    
    
    
    prcp_decade_trend = 10.*prcp_regress_coeff;  % changing yearly trends to be mm/decade
    
    
    lon_cell = stations_NE525_pro(3,:); 
    lat_cell = stations_NE525_pro(2,:);
    
    clear decadal_trend
    decadal_trend(:,1) = cell2mat(lon_cell)' ; 
    decadal_trend(:,2) = cell2mat(lat_cell)';
    decadal_trend(:,3) = prcp_decade_trend';
    
    indices = find(decadal_trend(:,3)>80 |decadal_trend(:,3)<-30); 
    decadal_trend(indices, :) = []
    
    [Z, refvec] = geoloc2grid(decadal_trend(:,2), decadal_trend(:,1), ...
        decadal_trend(:,3), 0.25);
    
    
    
    