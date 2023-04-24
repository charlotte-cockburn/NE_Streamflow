clear all, close all

%%%%% CHANGE THIS %%%%%

river = 'Diamond';

%%%%% CHANGE NOTHING BELOW THIS LINE %%%%%

if(contains(river,'White'))
    load('RF_data_WR_R1.mat'); 
    RF_table = RF_data_WR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
    
    load('RF_WRF_hist_WR_MPI.mat'); 
    RF_hist_table = RF_WRF_hist_WR; 
    
    load('RF_WRF_fut_WR_MPI.mat'); 
    RF_fut_table = RF_WRF_fut_WR; 
  
    acronym = 'WR';
    title_name = 'c) White River'; 

  
elseif(contains(river,'Shenandoah'))
    load('RF_data_SR_R1.mat'); 
    RF_table = RF_data_SR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
    
    load('RF_WRF_hist_SR_MPI.mat'); 
    RF_hist_table = RF_WRF_hist_SR; 
    
    load('RF_WRF_fut_SR_MPI.mat'); 
    RF_fut_table = RF_WRF_fut_SR;  
    
    acronym = 'SR';
    title_name = 'd) Shenandoah River'; 
    
elseif(contains(river,'Mattawamkeag'))
    load('RF_data_MR_R1.mat'); 
    RF_table = RF_data_MR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
    
    load('RF_WRF_hist_MR_MPI.mat'); 
    RF_hist_table = RF_WRF_hist_MR; 
    
    load('RF_WRF_fut_MR_MPI.mat'); 
    RF_fut_table = RF_WRF_fut_MR;  
    
    acronym = 'MR';
    title_name = 'a) Mattawamkeag River'; 

    
elseif(contains(river,'Diamond'))
    load('RF_data_DR_R1.mat'); 
    RF_table = RF_data_DR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 

    load('RF_WRF_hist_DR_MPI.mat'); 
    RF_hist_table = RF_WRF_hist_DR; 
    
    load('RF_WRF_fut_DR_MPI.mat'); 
    RF_fut_table = RF_WRF_fut_DR;  
    
    acronym = 'DR';
    title_name = 'b) Dead Diamond River'; 
    
end

RF_table = RF_table(:,[1:4 9:17 33:38]);


winter_indices = find(RF_table.Month == 1| ...
    RF_table.Month ==2|RF_table.Month ==3 | RF_table.Month ==4 | ...
    RF_table.Month ==5 | RF_table.Month == 11  |RF_table.Month ==12);

summer_indices = find(RF_table.Month == 6 | RF_table.Month ==7 | ...
    RF_table.Month ==8 | RF_table.Month == 9 | RF_table.Month == 10);

RF_dataset = RF_table(:, 4:width(RF_table(1,:))-3); 
max_flag = logical(RF_table.topflow_indicator); 

RF_dataset_winter = RF_table(winter_indices,4:width(RF_table(1,:))-3);
[top_vals_win, top_indices_win] = maxk(RF_dataset_winter.discharge, ...
    round(length(RF_dataset_winter.discharge)*0.01)); 
%max_flag_win = logical(RF_table.topflow_indicator_winter(winter_indices)); 
max_flag_win = false(length(RF_dataset_winter.discharge),1); 
max_flag_win(top_indices_win) = true; 

RF_dataset_summer = RF_table(summer_indices,4:width(RF_table(1,:))-3);
%max_flag_sum = logical(RF_table.topflow_indicator_summer(summer_indices)); 
[top_vals_sum, top_indices_sum] = maxk(RF_dataset_summer.discharge, ...
    round(length(RF_dataset_summer.discharge)*0.01)); 
%max_flag_win = logical(RF_table.topflow_indicator_winter(winter_indices)); 
max_flag_sum = false(length(RF_dataset_summer.discharge),1); 
max_flag_sum(top_indices_sum) = true; 


dates_winter = datetime(RF_table.Year(winter_indices),...
    RF_table.Month(winter_indices), RF_table.Day(winter_indices)); 
dates_summer = datetime(RF_table.Year(summer_indices),...
    RF_table.Month(summer_indices), RF_table.Day(summer_indices));


%RF_dataset_winter.api = [];
%RF_dataset_winter.spi = [];
%RF_dataset_winter.spi_prvs = []; 
RF_dataset_winter.snow_5day = []; 
RF_dataset_winter.precip_5day = []; 


RF_dataset_summer.melt_5day = [];
RF_dataset_summer.melt_30day = [];
RF_dataset_summer.precip_5day = [];
RF_dataset_summer.snow_3day = []; 
RF_dataset_summer.snow_5day = []; 
RF_dataset_summer.snow_30day = []; 

RF_dataset_winter.melt_5day = RF_dataset_winter.melt_5day.*25.4; 
RF_dataset_winter.melt_30day = RF_dataset_winter.melt_30day*25.4;

% Change discharge from m^3/s to m3s

RF_dataset_winter.discharge = RF_dataset_winter.discharge.*0.028316847;
RF_dataset_summer.discharge = RF_dataset_summer.discharge.*0.028316847;

RF_dataset_winter_ws = RF_dataset_winter;
RF_dataset_winter_ws.melt_5day = [];
RF_dataset_winter_ws.melt_30day = [];
RF_dataset_winter_ws.snow_3day = []; 
RF_dataset_winter_ws.snow_30day = [];



%%
rng(1945,'twister')
if(contains(river, 'Mattawamkeag'))
    RF_model_summer = TreeBagger(150,RF_dataset_summer,'discharge', ...
        'OOBPredictorImportance', 'on', 'OOBPrediction', 'on', 'Method', ...
        'regression','NumPredictorstoSample',7, ...
        'MinLeafSize', 1);
    RF_model_winter = TreeBagger(150,RF_dataset_winter,'discharge', ...
        'OOBPredictorImportance', 'on', 'OOBPrediction', 'on', 'Method',...
        'regression', 'NumPredictorstoSample',7, ...
        'MinLeafSize', 1);
    RF_model_winter_NoSnow = TreeBagger(150,RF_dataset_winter_ws,'discharge', ...
        'OOBPredictorImportance', 'on', 'OOBPrediction', 'on', 'Method',...
        'regression', 'NumPredictorstoSample',7, ...
        'MinLeafSize', 1);
else 
    RF_model_summer = TreeBagger(150,RF_dataset_summer,'discharge', ...
    'OOBPredictorImportance', 'on', 'OOBPrediction', 'on', 'Method', ...
    'regression');
    RF_model_winter = TreeBagger(150,RF_dataset_winter,'discharge', ...
        'OOBPredictorImportance', 'on', 'OOBPrediction', 'on', 'Method',...
        'regression');
    RF_model_winter_NoSnow = TreeBagger(150,RF_dataset_winter_ws,'discharge', ...
        'OOBPredictorImportance', 'on', 'OOBPrediction', 'on', 'Method',...
        'regression');
end

%

%%%% Get metrics for the model %%%%

%% Summer!

y_all_sum = RF_dataset_summer.discharge; 
yhat_all_sum = oobPredict(RF_model_summer); 
[corr_all_sum, pval_all_sum] = corrcoef(y_all_sum, yhat_all_sum); 
rmse_std_all_sum = (sqrt(mean((yhat_all_sum-y_all_sum).^2)))/std(y_all_sum); 
% Winter!

y_all_win = RF_dataset_winter.discharge; 
yhat_all_win = oobPredict(RF_model_winter); 
[corr_all_win, pval_all_win] = corrcoef(y_all_win, yhat_all_win); 
rmse_std_all_win = (sqrt(mean((yhat_all_win-y_all_win).^2)))/std(y_all_win);

yhat_all_win_ws = oobPredict(RF_model_winter_NoSnow); 
[corr_all_win_ws, pval_all_win_ws] = corrcoef(y_all_win, yhat_all_win_ws); 
rmse_std_all_win_ws = (sqrt(mean((yhat_all_win_ws-y_all_win).^2)))/std(y_all_win);

[y_99_win, y_99_win_indices] = maxk(RF_dataset_winter.discharge, ...
    round(length(RF_dataset_winter.discharge)*0.01));
yhat_99_win_ws = yhat_all_win_ws(y_99_win_indices); 
[corr_99_win_ws, pval_99_win_ws] = corrcoef(y_99_win, yhat_99_win_ws); 
rmse_std_99_win_ws = (sqrt(mean((yhat_99_win_ws-y_99_win).^2)))/std(y_99_win); 


[y_99_sum, y_99_sum_indices] = maxk(RF_dataset_summer.discharge, ...
    round(length(RF_dataset_summer.discharge)*0.01));
yhat_99_sum = yhat_all_sum(y_99_sum_indices); 
[corr_99_sum, pval_99_sum] = corrcoef(y_99_sum, yhat_99_sum); 
rmse_std_99_sum = (sqrt(mean((yhat_99_sum-y_99_sum).^2)))/std(y_99_sum); 


performance = [rmse_std_all_win_ws.*100, rmse_std_all_sum.*100;...
    corr_all_win_ws(2)^2, corr_all_sum(2)^2];
performance_99 = [rmse_std_99_win_ws.*100, rmse_std_99_sum.*100; ...
    corr_99_win_ws(2)^2, corr_99_sum(2)^2]; 


%% WRF predictions - historical  

% Change this for WRF dates 
winter_indices_hist = find(RF_hist_table.Month == 1| ...
    RF_hist_table.Month ==2|RF_hist_table.Month ==3 | RF_hist_table.Month ==4 | ...
    RF_hist_table.Month ==5 | RF_hist_table.Month == 11 | RF_hist_table.Month ==12);

summer_indices_hist = find(RF_hist_table.Month == 6 | RF_hist_table.Month ==7 | ...
    RF_hist_table.Month ==8 | RF_hist_table.Month == 9 | RF_hist_table.Month ==10);

RF_dataset_hist = RF_hist_table(:, 4:end); 

WRF_dates_winter = datetime(RF_hist_table.Year(winter_indices_hist),...
    RF_hist_table.Month(winter_indices_hist), RF_hist_table.Day(winter_indices_hist));
RF_hist_winter = RF_table(winter_indices_hist,4:end);

RF_hist_summer = RF_table(summer_indices_hist,4:end);
WRF_dates_summer = datetime(RF_hist_table.Year(summer_indices_hist),...
    RF_hist_table.Month(summer_indices_hist), RF_hist_table.Day(summer_indices_hist));

%RF_hist_winter.api = [];
%RF_hist_winter.spi = [];
%RF_hist_winter.spi_prvs = []; 
RF_hist_winter.snow_5day = []; 
RF_hist_winter.precip_5day = [];

RF_hist_summer.melt_5day = [];
RF_hist_summer.melt_30day = [];
RF_hist_summer.precip_5day = [];
RF_hist_summer.snow_3day = []; 
RF_hist_summer.snow_5day = []; 
RF_hist_summer.snow_30day = []; 

RF_hist_winter.melt_5day = RF_hist_winter.melt_5day.*25.4; 
RF_hist_winter.melt_30day = RF_hist_winter.melt_30day*25.4; 

WRF_preds_win = predict(RF_model_winter,RF_hist_winter); 
WRF_preds_sum = predict(RF_model_summer,RF_hist_summer);

% WRF predictions - futorical  

% Change this for WRF dates 
winter_indices_fut = find(RF_fut_table.Month == 1| ...
    RF_fut_table.Month ==2|RF_fut_table.Month ==3 | RF_fut_table.Month ==4 | ...
    RF_fut_table.Month ==5 | RF_fut_table.Month == 11 | RF_fut_table.Month ==12);

summer_indices_fut = find(RF_fut_table.Month == 6 | RF_fut_table.Month ==7 | ...
    RF_fut_table.Month ==8 | RF_fut_table.Month == 9 | RF_fut_table.Month ==10);

RF_dataset_fut = RF_fut_table(:, 4:end); 

win_dates_temp = datetime(RF_fut_table.Year(winter_indices_fut),...
    RF_fut_table.Month(winter_indices_fut), RF_fut_table.Day(winter_indices_fut));
RF_fut_winter_temp = RF_fut_table(winter_indices_fut,4:end);

WRF_dates_fut_winter = win_dates_temp(win_dates_temp>=datetime(2070,01,01) & ...
    win_dates_temp<=datetime(2099,12,31)); 
RF_fut_winter = RF_fut_winter_temp(find(win_dates_temp>=datetime(2070,01,01) & ...
    win_dates_temp<=datetime(2099,12,31)), :); 

RF_fut_summer_temp = RF_fut_table(summer_indices_fut,4:end);
sum_dates_temp = datetime(RF_fut_table.Year(summer_indices_fut),...
    RF_fut_table.Month(summer_indices_fut), RF_fut_table.Day(summer_indices_fut));
WRF_dates_fut_summer = sum_dates_temp(sum_dates_temp>=datetime(2070,01,01) & ...
    sum_dates_temp<=datetime(2099,12,31)); 
RF_fut_summer = RF_fut_summer_temp(find(sum_dates_temp>=datetime(2070,01,01) & ...
    sum_dates_temp<=datetime(2099,12,31)),:); 

%RF_fut_winter.api = [];
%RF_fut_winter.spi = [];
%RF_fut_winter.spi_prvs = []; 
% Adding!
RF_fut_winter.snow_5day = []; 
RF_fut_winter.precip_5day = [];

RF_fut_summer.melt_5day = [];
RF_fut_summer.melt_30day = [];
RF_fut_summer.precip_5day = [];
RF_fut_summer.snow_3day = []; 
RF_fut_summer.snow_5day = []; 
RF_fut_summer.snow_30day = []; 

RF_fut_winter.melt_5day = RF_fut_winter.melt_5day.*25.4; 
RF_fut_winter.melt_30day = RF_fut_winter.melt_30day*25.4; 

WRF_preds_fut_win = predict(RF_model_winter,RF_fut_winter); 
WRF_preds_fut_sum = predict(RF_model_summer,RF_fut_summer);

%% WRF predictions - future cs data in ws model 

RF_hist_winter_ws = RF_hist_winter;
RF_fut_winter_ws = RF_fut_winter; 


RF_hist_winter_ws.melt_5day = [];
RF_hist_winter_ws.melt_30day = [];
RF_hist_winter_ws.snow_3day = []; 
RF_hist_winter_ws.snow_30day = []; 

RF_fut_winter_ws.melt_5day = [];
RF_fut_winter_ws.melt_30day = [];
RF_fut_winter_ws.snow_3day = []; 
RF_fut_winter_ws.snow_30day = []; 

WRF_preds_win_ws = predict(RF_model_winter_NoSnow,RF_hist_winter_ws);
WRF_preds_fut_win_ws = predict(RF_model_winter_NoSnow,RF_fut_winter_ws); 


%% FIND JUST DAILY 

top_historical_win = maxk(WRF_preds_win, round(length(WRF_preds_win)*0.01)); 
threshold_win = min(top_historical_win); 
top_future_win = WRF_preds_fut_win(WRF_preds_fut_win>=threshold_win);

top_historical_sum = maxk(WRF_preds_sum, round(length(WRF_preds_sum)*0.01)); 
threshold_sum = min(top_historical_sum); 
top_future_sum = WRF_preds_fut_sum(WRF_preds_fut_sum>=threshold_sum);


% no-snow model (JUST FUTURE)
top_future_win_ws = WRF_preds_fut_win_ws(WRF_preds_fut_win_ws>=threshold_win);

% no-snow model (BOTH)
top_historical_win_ws = maxk(WRF_preds_win_ws, round(length(WRF_preds_win_ws)*0.01)); 
threshold_win_ws = min(top_historical_win_ws); 
top_future_win_ws_both = WRF_preds_fut_win_ws(WRF_preds_fut_win_ws>=threshold_win_ws);


% getting percent change 

win_hist = sum(top_historical_win); 
win_fut = sum(top_future_win);
perc_change_win = (win_fut-win_hist)/win_hist;
perc_change_win = round(perc_change_win*100, 1)

sum_hist = sum(top_historical_sum); 
sum_fut = sum(top_future_sum);
perc_change_sum = (sum_fut-sum_hist)/sum_hist;
perc_change_sum = round(perc_change_sum*100, 1)

% no-snow BOTH 
win_hist_ws = sum(top_historical_win_ws); 
win_fut_ws_both = sum(top_future_win_ws_both);
perc_change_win_ws_both = (win_fut_ws_both-win_hist_ws)/win_hist_ws;
perc_change_win_ws_both = round(perc_change_win_ws_both*100, 1)

% No-snow just FUTURE
win_fut_ws = sum(top_future_win_ws);
perc_change_win_ws = (win_fut_ws-win_hist)/win_hist;
perc_change_win_ws = round(perc_change_win_ws*100, 1)

 

%% Compare all flows! 

all_win_hist = sum(WRF_preds_win); 
all_win_fut = sum(WRF_preds_fut_win); 
perc_change_all_win = (all_win_fut-all_win_hist)/all_win_hist;
perc_change_all_win = round(perc_change_all_win*100, 1)

all_sum_hist = mean(WRF_preds_sum); 
all_sum_fut = mean(WRF_preds_fut_sum); 
perc_change_all_sum = (all_sum_fut-all_sum_hist)/all_sum_hist;
perc_change_all_sum = round(perc_change_all_sum*100, 1)

all_win_hist_ws = sum(WRF_preds_win_ws); 
all_win_fut_ws = sum(WRF_preds_fut_win_ws); 
perc_change_all_win_ws_both = (all_win_fut_ws-all_win_hist_ws)/all_win_hist_ws;
perc_change_all_win_ws_both = round(perc_change_all_win_ws_both*100, 1)

perc_change_all_win_ws = (all_win_fut_ws-all_win_hist)/all_win_hist;
perc_change_all_win_ws = round(perc_change_all_win_ws*100, 1)



%% KS TEST TIME 


% [ks_win, p_win] = kstest2(WRF_preds_win, WRF_preds_fut_win); 
% %[ks_win_top, p_win_top] = kstest2(top_historical_win, top_future_win)
% 
% [ks_sum, p_sum] = kstest2(WRF_preds_sum, WRF_preds_fut_sum); 
% %[ks_sum_top, p_sum_top] = kstest2(top_historical_sum, top_future_sum); 
% 
% [ks_win, p_win] = kstest2(WRF_preds_win, WRF_preds_fut_win); 
% [ks_win_top, p_win_top] = kstest2(top_historical_win, top_future_win)
% 
% ks_tests = [ks_win, p_win, ks_sum, p_sum; ks_win_top, p_win_top, ...
%     ks_sum_top, p_sum_top]; 

%% LOG FLOW DIAGRAMS 


flows_fig_cs = figure(2), clf;
set(flows_fig_cs, 'Position', [200 200 500 400]);
histogram(log10(WRF_preds_win), 'BinWidth', 0.15); 
hold on;
histogram(log10(WRF_preds_fut_win), 'BinWidth', 0.15); 
leg = legend('Historical', 'Future', 'FontSize', 16) ;
set(gca,'FontSize',16)
xlabel('log(Discharge)', 'FontSize', 18); 
ylabel('Days', 'FontSize', 18);

ks_all = kstest2(WRF_preds_win, WRF_preds_fut_win); 
if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end

if(perc_change_all_win>=0)
    title(strcat(title_name, ' (+',num2str(perc_change_all_win), '%)', asterisk), 'FontSize', 22);
else
    title(strcat(title_name, ' (',num2str(perc_change_all_win), '%)', asterisk), 'FontSize', 22);
end

if(contains(river, 'White'))
    xlim([0.75 3.4]); 
    axes('Position',[0.58 0.48 0.3 0.25]);
elseif(contains(river, 'Shenandoah'))
    xlim([-1 3.2])
    axes('Position',[0.58 0.48 0.3 0.25]);
elseif(contains(river, 'Mattawamkeag'))
    xlim([0.25 4])
    axes('Position',[0.6 0.48 0.29 0.25]);
elseif(contains(river, 'Diamond'))
    xlim([0 2.6])
    axes('Position',[0.6 0.48 0.29 0.25]);
end

histogram(log10(top_historical_win), 'BinWidth', 0.05); 
hold on;
histogram(log10(top_future_win), 'BinWidth', 0.05); 


if(isempty(top_future_win))
    ks_top = 0;
    xlim([min(log10(top_historical_win)) max(log10(top_historical_win))]);
else
    ks_top = kstest2(top_historical_win, top_future_win); 
end

if(ks_top ==1) 
    asterisk = '';
elseif(ks_top ==0);
    asterisk = '';
end

if(perc_change_win>=0)
    title(strcat('99th % flows (+',num2str(perc_change_win), '%)', asterisk), ...
        'FontSize', 14);
else
    title(strcat('99th % flows (',num2str(perc_change_win), '%)', asterisk), ...
        'FontSize', 14);
end

saveas(flows_fig_cs, string(strcat(acronym, '_log_flowchange_cs_MPI.png')));

%% SUMMER LOG FLOWS
flows_fig_ws = figure(3), clf ;
set(flows_fig_ws, 'Position', [200 200 500 400]);
histogram(log10(WRF_preds_sum), 'BinWidth', 0.2); 
hold on;
histogram(log10(WRF_preds_fut_sum), 'BinWidth', 0.2); 

leg2 = legend('Historical', 'Future', 'FontSize', 16) ;
set(gca,'FontSize',16)

ks_all = kstest2(WRF_preds_sum, WRF_preds_fut_sum); 
if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end

if(perc_change_all_sum>=0)
    title(strcat(title_name, ' (+',num2str(perc_change_all_sum), '%)', asterisk), 'FontSize', 22);
else
    title(strcat(title_name, ' (',num2str(perc_change_all_sum), '%)', asterisk), 'FontSize', 22);
end
xlabel('log(Discharge)', 'FontSize', 18); 
ylabel('Days', 'FontSize', 18); 

if(contains(river, 'White'))
    %xlim([4.5 10]); 
    axes('Position',[0.5 0.42 0.4 0.3]);
elseif(contains(river, 'Shenandoah'))
    xlim([-1.5 4]); 
    axes('Position',[0.57 0.42 0.3 0.3]);
elseif(contains(river, 'Mattawamkeag'))
    xlim([0.25 3.75]); 
    axes('Position',[0.575 0.42 0.3 0.3]);
elseif(contains(river, 'Diamond'))
    %xlim([2.5 9]); 
    axes('Position',[0.575 0.42 0.3 0.3]);
end

histogram(log10(top_historical_sum), 'BinWidth', 0.1); 
hold on;
histogram(log10(top_future_sum), 'BinWidth', 0.1); 

ks_top = kstest2(top_historical_sum, top_future_sum); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(perc_change_sum>=0)
    title(strcat('99th % flows (+',num2str(perc_change_sum), '%)', asterisk), ...
        'FontSize', 16);
else
    title(strcat('99th % flows (',num2str(perc_change_sum), '%)', asterisk),...
        'FontSize', 16);
end

saveas(flows_fig_ws, string(strcat(acronym, '_log_flowchange_ws_MPI.png')));

%% LOG FLOW DIAGRAMS - cold season w warm season model


flows_fig_cs_ws = figure(2), clf;
set(flows_fig_cs_ws, 'Position', [200 200 500 400]);
histogram(log10(WRF_preds_win), 'BinWidth', 0.15); 
hold on;
histogram(log10(WRF_preds_fut_win_ws), 'BinWidth', 0.15); 
leg = legend('Historical', 'Future', 'FontSize', 16) ;
set(gca,'FontSize',16)
xlabel('log(Discharge)', 'FontSize', 18); 
ylabel('Days', 'FontSize', 18);

ks_all = kstest2(WRF_preds_win, WRF_preds_fut_win_ws); 
if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end


if(perc_change_all_win_ws>=0)
    title(strcat(title_name, ' (+',num2str(perc_change_all_win_ws), '%)', asterisk), 'FontSize', 22);
else
    title(strcat(title_name, ' (',num2str(perc_change_all_win_ws), '%)', asterisk), 'FontSize', 22);
end

if(contains(river, 'White'))
    xlim([0.75 3.4]); 
    axes('Position',[0.58 0.48 0.3 0.25]);
elseif(contains(river, 'Shenandoah'))
    xlim([-1 3.2])
    axes('Position',[0.58 0.48 0.3 0.25]);
elseif(contains(river, 'Mattawamkeag'))
    xlim([0.25 4.3])
    axes('Position',[0.6 0.48 0.29 0.25]);
elseif(contains(river, 'Diamond'))
    xlim([0 3])
    axes('Position',[0.6 0.48 0.29 0.25]);
end

histogram(log10(top_historical_win), 'BinWidth', 0.05); 
hold on;
histogram(log10(top_future_win_ws), 'BinWidth', 0.05);

ks_top = kstest2(top_historical_win, top_future_win_ws); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(perc_change_win_ws>=0)
    title(strcat('99th % flows (+',num2str(perc_change_win_ws), '%)', asterisk), ...
        'FontSize', 14);
else
    title(strcat('99th % flows (',num2str(perc_change_win_ws), '%)', asterisk), ...
        'FontSize', 14);
end

saveas(flows_fig_cs_ws, string(strcat(acronym, '_log_flowchange_cs_ws_MPI.png')));

%% LOG FLOW DIAGRAMS - cold season w warm season model in hist AND fut


flows_fig_cs_ws2 = figure(3), clf;
set(flows_fig_cs_ws2, 'Position', [200 200 500 400]);
histogram(log10(WRF_preds_win_ws), 'BinWidth', 0.15); 
hold on;
histogram(log10(WRF_preds_fut_win_ws), 'BinWidth', 0.15); 
leg = legend('Historical', 'Future', 'FontSize', 16) ;
set(gca,'FontSize',16)
xlabel('log(Discharge)', 'FontSize', 18); 
ylabel('Days', 'FontSize', 18);

ks_all = kstest2(WRF_preds_win_ws, WRF_preds_fut_win_ws); 
if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end


if(perc_change_all_win_ws_both>=0)
    title(strcat(title_name, ' (+',num2str(perc_change_all_win_ws_both), '%)', asterisk), 'FontSize', 22);
else
    title(strcat(title_name, ' (',num2str(perc_change_all_win_ws_both), '%)', asterisk), 'FontSize', 22);
end

if(contains(river, 'White'))
    xlim([0.75 3.4]); 
    axes('Position',[0.58 0.48 0.3 0.25]);
elseif(contains(river, 'Shenandoah'))
    xlim([-1 3.2])
    axes('Position',[0.58 0.48 0.3 0.25]);
elseif(contains(river, 'Mattawamkeag'))
    xlim([0.25 4.3])
    axes('Position',[0.6 0.48 0.29 0.25]);
elseif(contains(river, 'Diamond'))
    xlim([0 3])
    axes('Position',[0.6 0.48 0.29 0.25]);
end

histogram(log10(top_historical_win_ws), 'BinWidth', 0.05); 
hold on;
histogram(log10(top_future_win_ws_both), 'BinWidth', 0.05);

ks_top = kstest2(top_historical_win_ws, top_future_win_ws_both); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(perc_change_win_ws_both>=0)
    title(strcat('99th % flows (+',num2str(perc_change_win_ws_both), '%)', asterisk), ...
        'FontSize', 14);
else
    title(strcat('99th % flows (',num2str(perc_change_win_ws_both), '%)', asterisk), ...
        'FontSize', 14);
end

saveas(flows_fig_cs_ws2, string(strcat(acronym, '_log_flowchange_noSnowEither_MPI.png')));


%% SEASONAL CYCLES 

winter_indices = find(RF_table.Month == 1| ...
    RF_table.Month ==2|RF_table.Month ==3 | RF_table.Month ==4 | ...
    RF_table.Month ==5 | RF_table.Month == 11 | RF_table.Month ==12);

summer_indices = find(RF_table.Month == 6 | RF_table.Month ==7 | ...
    RF_table.Month ==8 | RF_table.Month == 9 | RF_table.Month ==10);


obs_dis_win = RF_dataset_winter.discharge; 
WRF_dis_hist_win = WRF_preds_win; 
WRF_dis_fut_win = WRF_preds_fut_win; 
rf_dis_hist_win = oobPredict(RF_model_winter);

obs_dis_sum = RF_dataset_summer.discharge; 
WRF_dis_hist_sum = WRF_preds_sum; 
WRF_dis_fut_sum = WRF_preds_fut_sum;
rf_dis_hist_sum = oobPredict(RF_model_summer);


WRF_dates_fut_winter;
WRF_dates_fut_summer;

dates_win = RF_table(winter_indices,1:3); 
dates_win = datetime(dates_win.Year, dates_win.Month, dates_win.Day);
dates_indices_win = find(dates_win>=datetime(1976,01,01) & ...
    dates_win<= datetime(2005,12,31)); 
dates_win = dates_win(dates_indices_win); 
obs_dis_win = obs_dis_win(dates_indices_win); 

dates_sum = RF_table(summer_indices,1:3); 
dates_sum = datetime(dates_sum.Year, dates_sum.Month, dates_sum.Day);
dates_indices_sum = find(dates_sum>=datetime(1976,01,01) & ...
    dates_sum<= datetime(2005,12,31)); 
dates_sum = dates_sum(dates_indices_sum); 
obs_dis_sum = obs_dis_sum(dates_indices_sum); 

monthly_mean_obs_win = zeros(12,1); 
monthly_mean_wrf_win_hist = zeros(12,1); 
monthly_mean_wrf_win_fut = zeros(12,1); 
monthly_mean_rf_win_hist = zeros(12,1);

for i = [1:5 11:12] 
    
    monthly_mean_obs_win(i) = nanmean(obs_dis_win(month(dates_win)==i));
    monthly_mean_wrf_win_hist(i) = nanmean(WRF_dis_hist_win(month(WRF_dates_fut_winter)==i)); 
    monthly_mean_wrf_win_fut(i) = nanmean(WRF_dis_fut_win(month(WRF_dates_fut_winter)==i)); 
    
    monthly_mean_rf_win_hist(i) = nanmean(rf_dis_hist_win(month(dates_win)==i));
    
end

monthly_mean_obs_sum = zeros(12,1); 
monthly_mean_wrf_sum_hist = zeros(12,1); 
monthly_mean_wrf_sum_fut = zeros(12,1); 
monthly_mean_rf_sum_hist = zeros(12,1);

for i = 6:10
    monthly_mean_obs_sum(i) = nanmean(obs_dis_sum(month(dates_sum)==i)); 
    monthly_mean_wrf_sum_hist(i) = nanmean(WRF_dis_hist_sum(month(WRF_dates_fut_summer)==i)); 
    monthly_mean_wrf_sum_fut(i) = nanmean(WRF_dis_fut_sum(month(WRF_dates_fut_summer)==i));  
    monthly_mean_rf_sum_hist(i) = nanmean(rf_dis_hist_sum(month(dates_sum)==i));

end

monthly_mean_wrf_hist = monthly_mean_wrf_win_hist + monthly_mean_wrf_sum_hist;
monthly_mean_wrf_fut = monthly_mean_wrf_win_fut + monthly_mean_wrf_sum_fut; 
monthly_mean_obs = monthly_mean_obs_win + monthly_mean_obs_sum; 
monthly_mean_rf = monthly_mean_rf_win_hist + monthly_mean_rf_sum_hist; 

%%
seasonal_fig= figure(12), clf; 

plot([1:12], monthly_mean_wrf_hist, 'LineWidth', 2, 'Color', 'blue'); 
hold on;
plot([1:12], monthly_mean_wrf_fut, 'LineWidth', 2, 'Color', 'red');
plot([1:12], monthly_mean_obs,'--', 'LineWidth', 2, 'Color', [0.5 0.5 0.5 0.6]);
grid on

xticks([1:12]); 
xticklabels({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', ...
    'Oct', 'Nov', 'Dec', 'Jan', 'Feb'}); 
ax = gca;
ax.FontSize = 18;
xlim([1 12]); 
legend('Historical simulated', 'Future Simulated','Observed','FontSize', 18); 
legend boxoff;
ylabel('Discharge [m^3/s]', 'FontSize', 18);
title(title_name, 'FontSize', 22); 
if(contains(river, 'Shenandoah'))
    ylim([0 20]);
end

saveas(seasonal_fig, strcat(acronym, '_fut_hist_seasonal_cycles_MPI.png')); 




