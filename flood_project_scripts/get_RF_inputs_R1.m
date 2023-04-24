clear all, close all

%%%%% CHANGE THIS %%%%%

river = 'White'; 

%%%%% CHANGE NOTHING BELOW THIS LINE %%%%

if(contains(river,'White'))
    load('Livneh_data_white.mat'); 
    Livneh_data = Livneh_data_white;
    
    load('Livneh_sm_data_white.mat'); 
    soilmoist_data = Livneh_sm_data_white;
    
    load('White_discharge.mat');
    discharge_data = White_discharge; 
    
    %load('regressed_snowdepth_WR.mat'); 
    %snow_data = regressed_snowdepth_WR;
    
    load('WR_snowdepth_avg_predictions_test.mat'); 
    snow_data = regressed_snow_avg; 
    
    load('regressed_soilmoist_WR.mat'); 
    regressed_soilmoist_data = regressed_soilmoist_WR;
    
    soil_layers = readmatrix('vic.nldas.mexico.soil.txt');
    watershed = shaperead('Data/White_river/watershed_WR/globalwatershed.shp'); 
    
    min_lat = min(watershed.Y); 
    max_lat = max(watershed.Y); 
    min_lon = min(watershed.X); 
    max_lon = max(watershed.X);
    
    k_val = 0.74;
    
    two_year_count = 48; 
    
    acronym = 'WR';
    
elseif(contains(river,'Shenandoah'))
    load('Livneh_data_shenandoah.mat'); 
    Livneh_data = Livneh_data_shenandoah(4019:35429,:);
    
    load('Livneh_sm_data_shenandoah.mat'); 
    soilmoist_data = Livneh_sm_data_shenandoah;
    
    load('Shenandoah_discharge.mat');
    discharge_data = shenandoah_discharge(276:31686,:); 
    
    %load('regressed_snowdepth_SR.mat'); 
    %snow_data = regressed_snowdepth_SR;
    
    load('SR_snowdepth_avg_predictions_test.mat'); 
    snow_data = regressed_snow_avg; 
    
    load('regressed_soilmoist_SR.mat'); 
    regressed_soilmoist_data = regressed_soilmoist_SR;
    
    soil_layers = readmatrix('vic.nldas.mexico.soil.txt');
    watershed = shaperead('Data/Shenandoah_river/watershed_SR/globalwatershed.shp'); 
    min_lat = min(watershed.Y); 
    max_lat = max(watershed.Y); 
    min_lon = min(watershed.X); 
    max_lon = max(watershed.X);
    
    k_val = 0.81;
    
    two_year_count = 43; 
    
    acronym = 'SR';
    
elseif(contains(river,'Mattawamkeag'))
    load('Livneh_data_mattawamkeag.mat'); 
    Livneh_data = Livneh_data_mattawamkeag(7306:35429,:);
    
    load('Livneh_sm_data_mattawamkeag.mat'); 
    soilmoist_data = Livneh_sm_data_mattawamkeag;
    
    load('Mattawamkeag_discharge.mat');
    discharge_data = mattawamkeag_discharge(458:28216,:); 
    
    %load('regressed_snowdepth_MR.mat'); 
    %snow_data = regressed_snowdepth_MR;
    
    load('MR_snowdepth_avg_predictions_test.mat'); 
    snow_data = regressed_snow_avg; 
    
    load('regressed_soilmoist_MR.mat'); 
    regressed_soilmoist_data = regressed_soilmoist_MR;
    
    soil_layers = readmatrix('vic.nldas.mexico.soil.txt');
    watershed = shaperead('Data/Mattawamkeag_river/watershed_MR/globalwatershed.shp'); 
    min_lat = min(watershed.Y); 
    max_lat = max(watershed.Y); 
    min_lon = min(watershed.X); 
    max_lon = max(watershed.X);
    
    k_val = 0.98;
    
    two_year_count = 38; 
    
    acronym = 'MR';
    
elseif(contains(river, 'Diamond'))
    
    load('Livneh_data_diamond.mat'); 
    Livneh_data = Livneh_data_diamond(9863:35429,:);
    
    load('Livneh_sm_data_diamond.mat'); 
    soilmoist_data = Livneh_sm_data_diamond(9863:35429,:);
    
    load('diamond_discharge.mat');
    discharge_data = diamond_discharge(167:25733,:); 
    
    %load('regressed_snowdepth_DR.mat'); 
    %snow_data = regressed_snowdepth_DR;
    
    load('DR_snowdepth_avg_predictions_test.mat'); 
    snow_data = regressed_snow_avg; 
    
    load('regressed_soilmoist_DR.mat'); 
    regressed_soilmoist_data = regressed_soilmoist_DR;
    
    soil_layers = readmatrix('vic.nldas.mexico.soil.txt');
    watershed = shaperead('Data/DeadDiamond_river/watershed_DR/globalwatershed.shp'); 
    min_lat = min(watershed.Y); 
    max_lat = max(watershed.Y); 
    min_lon = min(watershed.X); 
    max_lon = max(watershed.X);
    
    k_val = 0.92;    
    two_year_count = 35; 
    
    acronym = 'DR';
    
end

first_year = year(discharge_data.Date(1));
last_year = 2011;
two_year_count = round(((last_year-first_year)+1)/2); 

%%

dates_livneh = datetime(Livneh_data(:,3), Livneh_data(:,2), ...
    Livneh_data(:,1)); 

dates_discharge = discharge_data.Date;

livneh_TT = timetable(dates_livneh, Livneh_data(:,1), ...
    Livneh_data(:,2), Livneh_data(:,3), Livneh_data(:,4), ...
    Livneh_data(:,5), Livneh_data(:,6));

discharge_TT = timetable(dates_discharge, discharge_data.Discharge); 

combined_TT = synchronize(livneh_TT, discharge_TT); 
combined_TT.Properties.VariableNames = {'day', 'month', 'year', 'precip', ...
    'max_temp','min_temp', 'discharge'};  

combined_TT(find(year(combined_TT.dates_livneh)>2011),:) = [];




%% Getting RF inputs - temperature 

av_temp = nanmean(table2array(combined_TT(:,5:6)), 2); 


cdd = zeros(length(table2array(combined_TT(:,1))),1);
for i = combined_TT.year(1):2011
   
    
    dates = table2array(combined_TT(:,1:3));
    indices = find(dates(:,3) == i-1 & (combined_TT.min_temp) <0 & ...
        (dates(:,2) == 11 | dates(:,2) == 12));
    if(i == combined_TT.year(1))
        indices = find(dates(:,3) == i & (combined_TT.min_temp) <0 & ...
        (dates(:,2) == 1));
    end
    this_cdd = sum(combined_TT.min_temp(indices));
    year_indices = find(dates(:,3) == i); 
    cdd(year_indices) = this_cdd;
end

%% Divide precip 

precip = combined_TT.precip; 

snow_test = zeros(length(precip),1); 
snow_test(av_temp<0) = precip(av_temp<0); 

precip_test = zeros(length(precip),1); 
precip_test(av_temp>=0) = precip(av_temp>=0); 


%% Precipitation

antec_precip_3day = zeros(length(combined_TT.precip), 1);
antec_precip_5day = zeros(length(combined_TT.precip), 1);
antec_precip_30day = zeros(length(combined_TT.precip), 1);

for i = 1:2
    antec_precip_3day(i) = sum(precip_test(1:i));
end

for i = 3:length(combined_TT.precip)
    antec_precip_3day(i) = sum(precip_test(i-2:i));
end

for i = 1:4
    antec_precip_5day(i) = sum(precip_test(1:i));
end

for i = 5:length(combined_TT.precip)
    antec_precip_5day(i) = sum(precip_test(i-4:i));
end

for i = 1:29
    antec_precip_30day(i) = sum(precip_test(1:i));
end

for i = 30:length(combined_TT.precip)
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

precip = combined_TT.precip; 
year_tt = combined_TT.year; 
month_tt = combined_TT.month; 
k = 1;
for i = year_tt(1):2011
   
    this_data = precip(find(year_tt==i));
    these_months = month_tt(find(year_tt==i));
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
daily_spi = zeros(length(combined_TT.year), 1);
daily_spi_prvs = zeros(length(combined_TT.year),1);
precip_spi_prvs = cat(2, 0, precip_spi);
precip_spi_prvs(1) = precip_spi_prvs(2); %set the 1st month to the value of itself because we don't have the month before it :/
n = 1;
m = 1;
for i = year_tt(1):2011

    these_months = month_tt((year_tt==i));
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

% prvsly used k = 0.92

%% Changing soil moisture to VOLUMETRIC 

watershed_indices = find(soil_layers(:,3)>min_lat & soil_layers(:,3) <max_lat & ...
    soil_layers(:,4)>min_lon & soil_layers(:,4) <-max_lon); 

layer1_depth = nanmean(soil_layers(watershed_indices, 23))*1000; 
layer2_depth = nanmean(soil_layers(watershed_indices, 24))*1000; 
layer3_depth = nanmean(soil_layers(watershed_indices, 25))*1000; 

% OG data is in mm, these layers are in m, but we multiply
% by 1000 above to get mm. So this would give mm/mm. then
% x100 to get a %

% soilmoist_data(:,2) = array2table((table2array(soilmoist_data(:,2))./layer1_depth).*100);
% soilmoist_data(:,3) = array2table((table2array(soilmoist_data(:,3))./layer2_depth).*100);
% soilmoist_data(:,4) = array2table((table2array(soilmoist_data(:,4))./layer3_depth).*100);
% 
% soilmoist_data(:,6) = array2table((table2array(soilmoist_data(:,6))./layer1_depth).*100);
% soilmoist_data(:,7) = array2table((table2array(soilmoist_data(:,7))./layer2_depth).*100);
% soilmoist_data(:,8) = array2table((table2array(soilmoist_data(:,8))./layer3_depth).*100);
% 
% soilmoist_data(:,5) = array2table(nanmean(table2array(soilmoist_data(:, 2:4)),2));
% soilmoist_data(:,9) = array2table(nanmean(table2array(soilmoist_data(:, 6:8)),2));
% 

%% Soil moisture - month 
soilmoist_30days = zeros(length(table2array(soilmoist_data(:,1))), 4);

for i = 1:29
    soilmoist_30days(i,:) = sum(table2array(soilmoist_data(1:i, 2:5)));    
end

for i = 30:length(table2array(soilmoist_data(:,1)))
    soilmoist_30days(i,:) = sum(table2array(soilmoist_data(i-29:i, 2:5)));     
end

%% Snow depth 

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

RF_data = [(Livneh_data(:,1:3)), av_temp, combined_TT.max_temp, ...
    combined_TT.min_temp, cdd, combined_TT.precip,antec_precip_3day, ...
    antec_precip_5day, antec_precip_30day,antec_snow_3day, antec_snow_5day, antec_snow_30day, api, daily_spi, daily_spi_prvs, ...
    regressed_soilmoist_data, table2array(soilmoist_data(:,2:9)), soilmoist_30days, snow_data,...
    melt, melt_5day, melt_30day, combined_TT.discharge]; 
RF_data_table = array2table(RF_data); 
RF_data_table.Properties.VariableNames = {'Day', 'Month', 'Year', ...
    'av_temp', 'max_temp', 'min_temp', 'cdd','precip', 'precip_3day', ...
    'precip_5day','precip_30day', 'snow_3day', 'snow_5day', 'snow_30day', 'api', 'spi', 'spi_prvs', 'regressed_soilmoist', ...
    'sm1', 'sm2', 'sm3', 'sm_av', 'sm1_5day', 'sm2_5day', 'sm3_5day', 'sm_av_5day', ...
    'sm1_30day', 'sm2_30day', 'sm3_30day', 'sm_av_30day', 'snow_depth', ...
    'melt', 'melt_5day', 'melt_30day','discharge'}; 

RF_data_table(isnan((RF_data_table.discharge)),:) =[];



%%

discharge = RF_data_table.discharge; 
index = 1:length(discharge); 
index = index';
dates = datetime(RF_data_table.Year, RF_data_table.Month, RF_data_table.Day);
jdays = juliandate(dates); 

discharge = [jdays discharge index]; 

[top_200, top200_indices]  = maxk(discharge(:,2), 200);

top_200_table = discharge(top200_indices,:);
%%
n=1;
for i = 1:200 
    i
    
     if(n >two_year_count) 
        break;
     end
    
     
    if(i==1)
        top_flows(n,:) = top_200_table(i,:);
        n = n+1;
        continue;
    end
    
    
    difference_count = 0;
    for j =1:n-1
        difference = abs(top_flows(j,1) - top_200_table(i,1)); 
        if(difference <10)
            difference_count = difference_count + 1;  
        end
        
    end
    
    if(difference_count == 0)
        top_flows(n,:) = top_200_table(i,:); 
        n = n+1
    end
         
end
    
dates = datetime(top_flows(:,1), 'ConvertFrom', 'juliandate');

topflow_indicator = zeros(length(RF_data_table.Year),1); 
topflow_indicator(top_flows(:,3)) = 1; 

%%

winter_indices = find(RF_data_table.Month == 1| ...
    RF_data_table.Month ==2|RF_data_table.Month ==3 | RF_data_table.Month ==4 | ...
    RF_data_table.Month ==5 | RF_data_table.Month == 11 | RF_data_table.Month ==12);
summer_indices = find(RF_data_table.Month == 6 | RF_data_table.Month ==7 | ...
    RF_data_table.Month ==8 | RF_data_table.Month == 9 | RF_data_table.Month ==10);

discharge_winter = RF_data_table.discharge(winter_indices); 
index_winter = index(winter_indices); 
dates_winter = datetime(RF_data_table.Year, RF_data_table.Month, RF_data_table.Day);
jdays_winter = juliandate(dates_winter); 
jdays_winter = jdays_winter(winter_indices); 

discharge_winter = [jdays_winter discharge_winter index_winter]; 

[top_200, top200_winter_indices]  = maxk(discharge_winter(:,2), 200);

top_200_winter_table = discharge_winter(top200_winter_indices,:);

n=1;
for i = 1:200 
    i
    
     if(n >two_year_count) 
        break;
     end
    
     
    if(i==1)
        top_flows_winter(n,:) = top_200_winter_table(i,:);
        n = n+1;
        continue;
    end
    
    
    difference_count = 0;
    for j =1:n-1
        difference = abs(top_flows_winter(j,1) - top_200_winter_table(i,1)); 
        if(difference <10)
            difference_count = difference_count + 1;  
        end
        
    end
    
    if(difference_count == 0)
        top_flows_winter(n,:) = top_200_winter_table(i,:); 
        n = n+1
    end
         
end
    
dates_winter = datetime(top_flows_winter(:,1), 'ConvertFrom', 'juliandate');

topflow_indicator_winter = zeros(length(RF_data_table.Year),1); 
topflow_indicator_winter(top_flows_winter(:,3)) = 1; 

%%

discharge_summer = RF_data_table.discharge(summer_indices); 
index_summer = index(summer_indices); 
dates_summer = datetime(RF_data_table.Year, RF_data_table.Month, RF_data_table.Day);
jdays_summer = juliandate(dates_summer); 
jdays_summer = jdays_summer(summer_indices); 

discharge_summer = [jdays_summer discharge_summer index_summer]; 

[top_200, top200_summer_indices]  = maxk(discharge_summer(:,2), 250);

top_200_summer_table = discharge_summer(top200_summer_indices,:);

n=1;
for i = 1:250 
    i
    
     if(n >two_year_count) 
        break;
     end
    
     
    if(i==1)
        top_flows_summer(n,:) = top_200_summer_table(i,:);
        n = n+1;
        continue;
    end
    
    
    difference_count = 0;
    for j =1:n-1
        difference = abs(top_flows_summer(j,1) - top_200_summer_table(i,1)); 
        if(difference <10)
            difference_count = difference_count + 1;  
        end
        
    end
    
    if(difference_count == 0)
        top_flows_summer(n,:) = top_200_summer_table(i,:); 
        n = n+1
    end
         
end
    
dates_summer = datetime(top_flows_summer(:,1), 'ConvertFrom', 'juliandate');

topflow_indicator_summer = zeros(length(RF_data_table.Year),1); 
topflow_indicator_summer(top_flows_summer(:,3)) = 1; 

%%

RF_data_table = [RF_data_table array2table(topflow_indicator) array2table(topflow_indicator_winter) ...
    array2table(topflow_indicator_summer)]; 

if(contains(river,'White'))
    RF_data_WR = RF_data_table;
    save('RF_data_WR_snowModelTest.mat', 'RF_data_WR');
    
elseif(contains(river,'Shenandoah'))
    RF_data_SR = RF_data_table;
    save('RF_data_SR_snowModelTest.mat', 'RF_data_SR');

elseif(contains(river,'Mattawamkeag'))
    RF_data_MR = RF_data_table; 
    save('RF_data_MR_snowModelTest.mat', 'RF_data_MR');

elseif(contains(river, 'Diamond'))
    RF_data_DR = RF_data_table;
    save('RF_data_DR_snowModelTest.mat', 'RF_data_DR');
end


