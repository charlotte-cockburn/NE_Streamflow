%% Comparing future WRF to Livneh obs 
clear all, close all
%Load data 

load('historical_data_d1_MPI.mat')
%load('wrf_d1_19762005.mat'); 

og_lon = lon1; 
og_lat = lat1;


%% Trim to watersheds 

S_WR = shaperead('watershed_WR/globalwatershed.shp');
[in_WR] = inpolygon(lon1, lat1, S_WR.X, S_WR.Y);
in_WR(40, 47) = true; 
in_WR(39,48) = true; 


S_SR = shaperead('watershed_SR/globalwatershed.shp');
[in_SR] = inpolygon(lon1, lat1, S_SR.X, S_SR.Y);

S_MR = shaperead('watershed_MR/globalwatershed.shp');
[in_MR] = inpolygon(lon1, lat1, S_MR.X, S_MR.Y);

DR_find = find(lat1>44.63 & lat1<45.32 & lon1>-71.54 & lon1<-70.86); 
in_DR = false(69,69);
in_DR(3493) = true;

%%
rain_WR = zeros(length(rain(1,1,:)), 1);
rain_SR = zeros(length(rain(1,1,:)), 1);
rain_MR = zeros(length(rain(1,1,:)), 1);
rain_DR = zeros(length(rain(1,1,:)), 1);

temp_WR = zeros(length(temperature(1,1,:)), 1);
temp_SR = zeros(length(temperature(1,1,:)), 1);
temp_MR = zeros(length(temperature(1,1,:)), 1);
temp_DR = zeros(length(temperature(1,1,:)), 1);


for i = 1:length(rain(1,1,:))
    
    
    rain_temp = rain(:,:,i);
    temp_temp = temperature(:,:,i);
    
    this_rain_WR = nanmean(nanmean(rain_temp(in_WR),1),2);
    rain_WR(i,:) = this_rain_WR;
    this_temp_WR = nanmean(nanmean(temp_temp(in_WR),1),2);
    temp_WR(i,:) = this_temp_WR;
    
    this_rain_SR = nanmean(nanmean(rain_temp(in_SR),1),2);
    rain_SR(i,:) = this_rain_SR;
    this_temp_SR = nanmean(nanmean(temp_temp(in_SR),1),2);
    temp_SR(i,:) = this_temp_SR;
    
    this_rain_MR = nanmean(nanmean(rain_temp(in_MR),1),2);
    rain_MR(i,:) = this_rain_MR;
    this_temp_MR = nanmean(nanmean(temp_temp(in_MR),1),2);
    temp_MR(i,:) = this_temp_MR;
    
    this_rain_DR = nanmean(nanmean(rain_temp(in_DR),1),2);
    rain_DR(i,:) = this_rain_DR;
    this_temp_DR = nanmean(nanmean(temp_temp(in_DR),1),2);
    temp_DR(i,:) = this_temp_DR;
       
end

%% 
% count = 1;
% for i = 1:94 
%     time2(count:count+364,:) = squeeze(time(:,:,i)); 
%     count = count+365; 
% end
%%
%date = datetime(time2, 'InputFormat', 'dd-MMM-yyyy'); 

WR_historical_data_MPI = table(date, temp_WR, rain_WR); 
SR_historical_data_MPI = table(date, temp_SR, rain_SR); 
MR_historical_data_MPI = table(date, temp_MR, rain_MR); 
DR_historical_data_MPI = table(date, temp_DR, rain_DR); 



%% Load observations 

load('RF_data_WR.mat');
first_date = datetime(1976,01,01);
last_date =datetime(2005,12,31);
dates_WR = datetime(RF_data_WR.Year, RF_data_WR.Month, RF_data_WR.Day);
WR_data = [RF_data_WR.av_temp, RF_data_WR.precip];
WR_data = WR_data(dates_WR>first_date & dates_WR<last_date, :); 

load('RF_data_SR.mat');
first_date = datetime(1976,01,01);
last_date =datetime(2005,12,31);
dates_SR = datetime(RF_data_SR.Year, RF_data_SR.Month, RF_data_SR.Day);
SR_data = [RF_data_SR.av_temp, RF_data_SR.precip];
SR_data = SR_data(dates_SR>first_date & dates_SR<last_date, :); 

load('RF_data_MR.mat');
first_date = datetime(1976,01,01);
last_date =datetime(2005,12,31);
dates_MR = datetime(RF_data_MR.Year, RF_data_MR.Month, RF_data_MR.Day);
MR_data = [RF_data_MR.av_temp, RF_data_MR.precip];
MR_data = MR_data(dates_MR>first_date & dates_MR<last_date, :); 


load('RF_data_DR.mat');
first_date = datetime(1976,01,01);
last_date =datetime(2005,12,31);
dates_DR = datetime(RF_data_DR.Year, RF_data_DR.Month, RF_data_DR.Day);
DR_data = [RF_data_DR.av_temp, RF_data_DR.precip];
DR_data = DR_data(dates_DR>first_date & dates_DR<last_date, :); 

%%


temp_fig = figure(2), clf 
set(temp_fig, 'Position', [300 300 1100 9000]); 

subplot(2,2,1) 
histogram(WR_data(:,1),'BinWidth',1);
hold on
histogram(temp_WR,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Temperature [^oC]'); 
title('White River', 'FontSize', 14); 

subplot(2,2,2) 
histogram(SR_data(:,1),'BinWidth',1);
hold on
histogram(temp_SR,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Temperature [^oC]'); 
title('Shenandoah River', 'FontSize', 14); 

subplot(2,2,3) 
histogram(MR_data(:,1),'BinWidth',1);
hold on
histogram(temp_MR,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Temperature [^oC]'); 
title('Mattawamkeag River', 'FontSize', 14);

subplot(2,2,4) 
histogram(DR_data(:,1),'BinWidth',1);
hold on
histogram(temp_DR,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Temperature [^oC]'); 
title('Dead Diamond River', 'FontSize', 14);

saveas(temp_fig, 'temp_comp.png'); 

%% Precipitation 

% First get wet days 
WR_data_wet = WR_data(WR_data(:,2)>1, 2); 
rain_WR_wet = rain_WR(rain_WR>1); 

SR_data_wet = SR_data(SR_data(:,2)>1, 2); 
rain_SR_wet = rain_SR(rain_SR>1); 

MR_data_wet = MR_data(MR_data(:,2)>1, 2); 
rain_MR_wet = rain_MR(rain_MR>1); 

DR_data_wet = DR_data(DR_data(:,2)>1, 2); 
rain_DR_wet = rain_DR(rain_DR>1); 

wet_fig = figure(2), clf 

set(wet_fig, 'Position', [300 300 1100 9000]); 
subplot(2,2,1) 
histogram(WR_data_wet,'BinWidth',1);
hold on
histogram(rain_WR_wet,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Precipitation [mm]'); 
title('White River', 'FontSize', 14); 

subplot(2,2,2) 
histogram(SR_data_wet,'BinWidth',1);
hold on
histogram(rain_SR_wet,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Precipitation [mm]'); 
title('Shenandoah River', 'FontSize', 14); 

subplot(2,2,3) 
histogram(MR_data_wet,'BinWidth',1);
hold on
histogram(rain_MR_wet,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Precipitation [mm]'); 
title('Mattawamkeag River', 'FontSize', 14);

subplot(2,2,4) 
histogram(DR_data_wet,'BinWidth',1);
hold on
histogram(rain_DR_wet,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Precipitation [mm]'); 
title('Dead Diamond River', 'FontSize', 14);

saveas(wet_fig, 'wet_comp.png'); 


%% Extreme precipitation 

% First get xt days 
WR_data_xt = WR_data(WR_data(:,2)>30, 2); 
rain_WR_xt = rain_WR(rain_WR>30); 

SR_data_xt = SR_data(SR_data(:,2)>30, 2); 
rain_SR_xt = rain_SR(rain_SR>30); 

MR_data_xt = MR_data(MR_data(:,2)>30, 2); 
rain_MR_xt = rain_MR(rain_MR>30); 

DR_data_xt = DR_data(DR_data(:,2)>30, 2); 
rain_DR_xt = rain_DR(rain_DR>30); 

xt_fig = figure(3), clf 
set(xt_fig, 'Position', [300 300 1100 9000]); 

subplot(2,2,1) 
histogram(WR_data_xt,'BinWidth',1);
hold on
histogram(rain_WR_xt,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Precipitation [mm]'); 
title('White River', 'FontSize', 14); 

subplot(2,2,2) 
histogram(SR_data_xt,'BinWidth',1);
hold on
histogram(rain_SR_xt,'BinWidth',1);
set(h, 'Location', 'Best');
h = legend('Observed', 'future WRF');
ylabel('Frequency') 
xlabel('Precipitation [mm]'); 
title('Shenandoah River', 'FontSize', 14); 

subplot(2,2,3) 
histogram(MR_data_xt,'BinWidth',1);
hold on
histogram(rain_MR_xt,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Precipitation [mm]'); 
title('Mattawamkeag River', 'FontSize', 14); 


subplot(2,2,4) 
histogram(DR_data_xt,'BinWidth',1);
hold on
histogram(rain_DR_xt,'BinWidth',1);
h = legend('Observed', 'future WRF');
set(h, 'Location', 'Best');
ylabel('Frequency') 
xlabel('Precipitation [mm]'); 
title('Dead Diamond River', 'FontSize', 14); 

saveas(xt_fig, 'xt_comp.png'); 

%% KS testing 

[ks_test_temp_WR, p_temp_WR] = kstest2(WR_data(:,1), temp_WR); 
[ks_test_temp_SR, p_temp_SR]= kstest2(SR_data(:,1), temp_SR); 
[ks_test_temp_MR, p_temp_MR]= kstest2(MR_data(:,1), temp_MR); 
[ks_test_temp_DR, p_temp_DR]= kstest2(DR_data(:,1), temp_DR); 

% First get wet days 
WR_data_wet = WR_data(WR_data(:,2)>1, 2); 
rain_WR_wet = rain_WR(rain_WR>1); 

SR_data_wet = SR_data(SR_data(:,2)>1, 2); 
rain_SR_wet = rain_SR(rain_SR>1); 

MR_data_wet = MR_data(MR_data(:,2)>1, 2); 
rain_MR_wet = rain_MR(rain_MR>1); 

DR_data_wet = DR_data(DR_data(:,2)>1, 2); 
rain_DR_wet = rain_DR(rain_DR>1); 


[ks_test_1_WR, p1_WR]= kstest2(WR_data_wet, rain_WR_wet); 
[ks_test_1_SR, p1_SR]= kstest2(SR_data_wet, rain_SR_wet); 
[ks_test_1_MR, p1_MR]= kstest2(MR_data_wet, rain_MR_wet); 
[ks_test_1_DR, p1_DR]= kstest2(DR_data_wet, rain_DR_wet); 



WR_data_10 = WR_data(WR_data(:,2)>10, 2); 
rain_WR_10 = rain_WR(rain_WR>10); 

SR_data_10 = SR_data(SR_data(:,2)>10, 2); 
rain_SR_10 = rain_SR(rain_SR>10); 

MR_data_10 = MR_data(MR_data(:,2)>10, 2); 
rain_MR_10 = rain_MR(rain_MR>10); 

DR_data_10 = DR_data(DR_data(:,2)>10, 2); 
rain_DR_10 = rain_DR(rain_DR>10); 

[ks_test_10_WR, p10_WR]= kstest2(WR_data_10, rain_WR_10); 
[ks_test_10_SR, p10_SR]= kstest2(SR_data_10, rain_SR_10); 
[ks_test_10_MR, p10_MR]= kstest2(MR_data_10, rain_MR_10); 
[ks_test_10_DR, p10_DR]= kstest2(DR_data_10, rain_DR_10); 



WR_data_20 = WR_data(WR_data(:,2)>20, 2); 
rain_WR_20 = rain_WR(rain_WR>20); 

SR_data_20 = SR_data(SR_data(:,2)>20, 2); 
rain_SR_20 = rain_SR(rain_SR>20); 

MR_data_20 = MR_data(MR_data(:,2)>20, 2); 
rain_MR_20 = rain_MR(rain_MR>20); 

DR_data_20 = DR_data(DR_data(:,2)>20, 2); 
rain_DR_20 = rain_DR(rain_DR>20); 

[ks_test_20_WR, p20_WR]= kstest2(WR_data_20, rain_WR_20); 
[ks_test_20_SR, p20_SR]= kstest2(SR_data_20, rain_SR_20); 
[ks_test_20_MR, p20_MR]= kstest2(MR_data_20, rain_MR_20); 
[ks_test_20_DR, p20_DR]= kstest2(DR_data_20, rain_DR_20); 



WR_data_30 = WR_data(WR_data(:,2)>30, 2); 
rain_WR_30 = rain_WR(rain_WR>30); 

SR_data_30 = SR_data(SR_data(:,2)>30, 2); 
rain_SR_30 = rain_SR(rain_SR>30); 

MR_data_30 = MR_data(MR_data(:,2)>30, 2); 
rain_MR_30 = rain_MR(rain_MR>30); 

DR_data_30 = DR_data(DR_data(:,2)>30, 2); 
rain_DR_30 = rain_DR(rain_DR>30); 

[ks_test_30_WR, p30_WR]= kstest2(WR_data_30, rain_WR_30); 
[ks_test_30_SR, p30_SR]= kstest2(SR_data_30, rain_SR_30); 
[ks_test_30_MR, p30_MR]= kstest2(MR_data_30, rain_MR_30); 
[ks_test_30_DR, p30_DR]= kstest2(DR_data_30, rain_DR_30); 


ks_WR = [ks_test_temp_WR, ks_test_1_WR, ks_test_10_WR, ...
    ks_test_20_WR, ks_test_30_WR]; 
ks_SR = [ks_test_temp_SR, ks_test_1_SR, ks_test_10_SR, ...
    ks_test_20_SR, ks_test_30_SR]; 
ks_MR = [ks_test_temp_MR, ks_test_1_MR, ks_test_10_MR, ...
    ks_test_20_MR, ks_test_30_MR]; 
ks_DR = [ks_test_temp_DR, ks_test_1_DR, ks_test_10_DR, ...
    ks_test_20_DR, ks_test_30_DR];

pvals_WR = [p_temp_WR, p1_WR, p10_WR, ...
    p20_WR, p30_WR]; 
pvals_SR = [p_temp_SR, p1_SR, p10_SR, ...
    p20_SR, p30_SR]; 
pvals_MR = [p_temp_MR, p1_MR, p10_MR, ...
    p20_MR, p30_MR]; 
pvals_DR = [p_temp_DR, p1_DR, p10_DR, ...
    p20_DR, p30_DR];

%%

rain_mean = squeeze(mean(mean(rain,1),2));
rain_mean = rain_mean(23361:end); 

load('wrf_d1_19762005.mat'); 

rain_mean_hist = squeeze(mean(mean(rain,1),2));

%%
figure(1), clf 

histogram(rain_mean_hist); 
hold on
histogram(rain_mean); 

axes('Position',[0.7 0.7 0.2 0.2]);
histogram(rain_mean_hist, 'BinLimits', [50 100], 'BinWidth', 5); 
hold on
histogram(rain_mean, 'BinLimits', [50 100], 'BinWidth', 5); 






