clear all, close all

%%%%% CHANGE THIS %%%%%

river = 'Diamond'; 

%%%%% CHANGE NOTHING BELOW THIS LINE %%%%%

if(contains(river,'White'))
    load('WR_future_data_MPI.mat'); 
          
    %load('reg_snow_WRF_WR.mat'); 
    %fut_data = [WR_future_data, array2table(reg_snow_WRF_WR.snowdepth)];
    load('WR_snowdepth_avg_WRFfut_MPI.mat');
    fut_data = [WR_future_data_MPI, table(regressed_snow_WRFfut_avg)];

    acronym = 'WR';
    
    k_val = 0.74;

        
elseif(contains(river,'Shenandoah'))
     load('SR_future_data_MPI.mat'); 
          
    load('SR_snowdepth_avg_WRFfut_MPI.mat');
    fut_data = [SR_future_data_MPI, table(regressed_snow_WRFfut_avg)];    
    acronym = 'SR';
    
    k_val = 0.81;



elseif(contains(river,'Mattawamkeag'))
    load('MR_future_data_MPI.mat');         
    load('MR_snowdepth_avg_WRFfut_MPI.mat');
    fut_data = [MR_future_data_MPI, table(regressed_snow_WRFfut_avg)];
    acronym = 'MR';

     k_val = 0.98;
    
elseif(contains(river, 'Diamond'))
    load('DR_future_data_MPI.mat'); 
          
    load('DR_snowdepth_avg_WRFfut_MPI.mat');
    fut_data = [DR_future_data_MPI, table(regressed_snow_WRFfut_avg)];
    acronym = 'DR';
    
    k_val = 0.92;    


end

fut_data.Properties.VariableNames = {'date', 'av_temp', 'precip', 'snowdepth'};
%%
precip = fut_data.precip; 
av_temp = fut_data.av_temp;

snow_test = zeros(length(precip),1); 
snow_test(av_temp<0) = precip(av_temp<0); 

precip_test = zeros(length(precip),1); 
precip_test(av_temp>=0) = precip(av_temp>=0); 

%% Precipitation

antec_precip_3day = zeros(length(precip_test), 1);
antec_precip_5day = zeros(length(precip_test), 1);
antec_precip_30day = zeros(length(precip_test), 1);

for i = 1:2
    antec_precip_3day(i) = sum(precip_test(1:i));
end

for i = 3:length(precip_test)
    antec_precip_3day(i) = sum(precip_test(i-2:i));
end

for i = 1:4
    antec_precip_5day(i) = sum(precip_test(1:i));
end

for i = 5:length(precip_test)
    antec_precip_5day(i) = sum(precip_test(i-4:i));
end

for i = 1:29
    antec_precip_30day(i) = sum(precip_test(1:i));
end

for i = 30:length(precip_test)
    antec_precip_30day(i,:) = sum(precip_test(i-29:i));     
end


%% snowitation

antec_snow_3day = zeros(length(snow_test), 1);
antec_snow_5day = zeros(length(snow_test), 1);
antec_snow_30day = zeros(length(snow_test), 1);

for i = 1:2
    antec_snow_3day(i) = sum(snow_test(1:i));
end

for i = 3:length(snow_test)
    antec_snow_3day(i) = sum(snow_test(i-2:i));
end

for i = 1:4
    antec_snow_5day(i) = sum(snow_test(1:i));
end

for i = 5:length(snow_test)
    antec_snow_5day(i) = sum(snow_test(i-4:i));
end

for i = 1:29
    antec_snow_30day(i) = sum(snow_test(1:i));
end

for i = 30:length(snow_test)
    antec_snow_30day(i,:) = sum(snow_test(i-29:i));     
end


%% API and SPI 

precip = fut_data.precip; 
year_data = year(fut_data.date); 
month_data= month(fut_data.date); 
k = 1;
for i = year_data(1):year_data(end)
   
    this_data = precip(find(year_data==i));
    these_months = month_data(find(year_data==i));
    for j = 1:12
        this_month = this_data(find(these_months ==j));
        this_monthly_mean = nanmean(this_month);
        monthly_mean_precip(k) = this_monthly_mean;
        k = k+1;
    end
    
end

monthly_mean_precip = monthly_mean_precip';


precip_spi = SPI(monthly_mean_precip, 1, 12); %DOES THIS WORK? 
%%

daily_spi = zeros(length(fut_data.date), 1);
daily_spi_prvs = zeros(length(fut_data.date),1);
precip_spi_prvs = cat(2, 0, precip_spi);
precip_spi_prvs(1) = precip_spi_prvs(2); %set the 1st month to the value of itself because we don't have the month before it :/
n = 1;
m = 1;
for i = year_data(1):year_data(end)

    these_months = month_data((year_data==i));
    for j = 1:12
        this_month = these_months(these_months ==j); 
        for k = 1:length(this_month)
            daily_spi(n) = precip_spi(m);
            daily_spi_prvs(n) = precip_spi_prvs(m);
            n = n+1;
        end
        m = m +1;      
    end 
end

%% API
api = zeros(length(precip), 1);
api(1:60) = precip(1:60);
for i = 61:length(precip)
    
    decay_terms = zeros(60,1); 
    for j = 1:60 
        decay_terms(j) = k_val^j*precip(i-j); 
    end
    
    api(i) = precip(i) + sum(decay_terms);  
    
end

%% Snow depth 

snow_data = fut_data.snowdepth; 

melt = zeros(length(snow_data), 1); 
for i = 2:length(snow_data)
    snowdepth_change = snow_data(i) - snow_data(i-1);
    if(snowdepth_change<0) 
        melt(i) = snowdepth_change*-1; 
    else
        melt(i) = 0; 
    end
end

melt_5day = zeros(length(snow_data), 1);

for i =1:4
    melt_5day(i) = sum(melt(1:i));
end

for i =5:length(snow_data) 
    melt_5day(i) = sum(melt(i-4:i));
end

melt_30day = zeros(length(snow_data), 1);

for i =1:29
    melt_30day(i) = sum(melt(1:i));
end

for i =30:length(snow_data) 
    melt_30day(i) = sum(melt(i-29:i));
end





%% Combine it all! 

RF_data = [day(fut_data.date), month(fut_data.date), year(fut_data.date), ...
    fut_data.av_temp, antec_precip_3day, antec_precip_5day, antec_precip_30day, ...
    antec_snow_3day, antec_snow_5day, antec_snow_30day, ...
    api, daily_spi, daily_spi_prvs, melt_5day, melt_30day]; 
    
RF_data_table = array2table(RF_data); 
RF_data_table.Properties.VariableNames = {'Day', 'Month', 'Year', ...
    'av_temp', 'precip_3day','precip_5day','precip_30day', 'snow_3day', ...
    'snow_5day', 'snow_30day','api', 'spi', ...
    'spi_prvs','melt_5day', 'melt_30day'}; 

if(contains(river,'White'))
    RF_WRF_fut_WR = RF_data_table;
    save('RF_WRF_fut_WR_MPI.mat', 'RF_WRF_fut_WR');
elseif(contains(river,'Shenandoah'))
    RF_WRF_fut_SR = RF_data_table;
    save('RF_WRF_fut_SR_MPI.mat', 'RF_WRF_fut_SR');
elseif(contains(river,'Mattawamkeag'))
    RF_WRF_fut_MR = RF_data_table;
    save('RF_WRF_fut_MR_MPI.mat', 'RF_WRF_fut_MR');
elseif(contains(river, 'Diamond'))
    RF_WRF_fut_DR = RF_data_table;
    save('RF_WRF_fut_DR_MPI.mat', 'RF_WRF_fut_DR');
end


