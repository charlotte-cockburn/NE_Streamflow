%% FINAL RANDOM FOREST!!!
clear all, close all

%%%%% CHANGE THIS %%%%%

river = 'Diamond';

%%%%% CHANGE NOTHING BELOW THIS LINE %%%%%

if(contains(river,'White'))
    load('RF_data_WR_R1.mat'); 
    RF_table = RF_data_WR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
  
    acronym = 'WR';
    title_name = 'c) White River'; 
    
    pdp_win_lim =  [20 140];
    pdp_sum_lim =  [0 150];
   

elseif(contains(river,'Shenandoah'))
    load('RF_data_SR_R1.mat'); 
    RF_table = RF_data_SR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
    
    acronym = 'SR';
    title_name = 'd) Shenandoah River'; 
    
    pdp_win_lim =  [0 120];
    pdp_sum_lim =  [0 150];

    
elseif(contains(river,'Mattawamkeag'))
    load('RF_data_MR_R1.mat'); 
    RF_table = RF_data_MR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 

    acronym = 'MR';
    title_name = 'a) Mattawamkeag River'; 
    
    pdp_win_lim =  [0 200];
    pdp_sum_lim =  [0 350];
    
elseif(contains(river,'Diamond'))
    load('RF_data_DR_R1.mat'); 
    RF_table = RF_data_DR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 

    acronym = 'DR';
    title_name = 'b) Dead Diamond River'; 
    
    pdp_win_lim =  [0 50];
    pdp_sum_lim =  [0 35];



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
%RF_dataset_winter.discharge = log(RF_dataset_winter.discharge);
%RF_dataset_summer.discharge = log(RF_dataset_summer.discharge);


% %Bayseian Optimization
% maxMinLS = 20;
% minLS = optimizableVariable('minLS',[1,maxMinLS],'Type','integer');
% numPTS = optimizableVariable('numPTS',[1,8],'Type','integer');
% hyperparametersRF = [minLS; numPTS];
% params = [minLS; numPTS];
% y = 'discharge';
% 
% %Return results 
% results = bayesopt(@(params)oobErrRF(params,RF_dataset, y),hyperparametersRF,...
%     'AcquisitionFunctionName','expected-improvement-plus','Verbose',0);
% bestOOBErr = results.MinObjective;
% bestHyperparameters = results.XAtMinObjective;


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
        'regression');

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

% Summer!
%
y_all_sum = RF_dataset_summer.discharge; 
yhat_all_sum = oobPredict(RF_model_summer); 
[corr_all_sum, pval_all_sum] = corrcoef(y_all_sum, yhat_all_sum); 
rmse_std_all_sum = (sqrt(mean((yhat_all_sum-y_all_sum).^2)))/std(y_all_sum); 

[y_99_sum, y_99_sum_indices] = maxk(RF_dataset_summer.discharge, ...
    round(length(RF_dataset_summer.discharge)*0.01));
yhat_99_sum = yhat_all_sum(y_99_sum_indices); 
[corr_99_sum, pval_99_sum] = corrcoef(y_99_sum, yhat_99_sum); 
rmse_std_99_sum = (sqrt(mean((yhat_99_sum-y_99_sum).^2)))/std(y_99_sum); 

% Winter!

y_all_win = RF_dataset_winter.discharge; 
yhat_all_win = oobPredict(RF_model_winter); 
[corr_all_win, pval_all_win] = corrcoef(y_all_win, yhat_all_win); 
rmse_std_all_win = (sqrt(mean((yhat_all_win-y_all_win).^2)))/std(y_all_win); 

[y_99_win, y_99_win_indices] = maxk(RF_dataset_winter.discharge, ...
    round(length(RF_dataset_winter.discharge)*0.01));
yhat_99_win = yhat_all_win(y_99_win_indices); 
[corr_99_win, pval_99_win] = corrcoef(y_99_win, yhat_99_win); 
rmse_std_99_win = (sqrt(mean((yhat_99_win-y_99_win).^2)))/std(y_99_win); 


% Winter NO SNOW 

y_all_win_ws = RF_dataset_winter_ws.discharge; 
yhat_all_win_ws = oobPredict(RF_model_winter_NoSnow); 
[corr_all_win_ws, pval_all_win_ws] = corrcoef(y_all_win_ws, yhat_all_win_ws); 
rmse_std_all_win_ws = (sqrt(mean((yhat_all_win_ws-y_all_win_ws).^2)))/std(y_all_win_ws); 

[y_99_win_ws, y_99_win_indices] = maxk(RF_dataset_winter_ws.discharge, ...
    round(length(RF_dataset_winter_ws.discharge)*0.01));
yhat_99_win_ws = yhat_all_win_ws(y_99_win_indices); 
[corr_99_win_ws, pval_99_win_ws] = corrcoef(y_99_win_ws, yhat_99_win_ws); 
rmse_std_99_win_ws = (sqrt(mean((yhat_99_win_ws-y_99_win_ws).^2)))/std(y_99_win_ws); 



summer_correlations = [corr_all_sum(2), corr_99_sum(2)];
summer_correlations = summer_correlations.^2; 
winter_correlations = [corr_all_win(2), corr_99_win(2)];
winter_correlations = winter_correlations.^2; 
winter_correlations_ws = [corr_all_win_ws(2), corr_99_win_ws(2)];
winter_correlations_ws = winter_correlations_ws.^2;

rmse_std = [rmse_std_all_win, rmse_std_99_win; ...
    rmse_std_all_win_ws, rmse_std_99_win_ws; ...
    rmse_std_all_sum, rmse_std_99_sum]; 

correlations_all = [corr_all_win(2), corr_all_win_ws(2), corr_all_sum(2)];
correlations_all = correlations_all.^2;

rmse_std_all = [rmse_std_all_win, rmse_std_all_win_ws, rmse_std_all_sum]; 

correlations_99 = [corr_99_win(2), corr_99_win_ws(2), corr_99_sum(2)]; 
correlations_99 = correlations_99.^2; 

rmse_std_99 = [rmse_std_99_win,rmse_std_99_win_ws, rmse_std_99_sum];
rmse_std_99 = rmse_std_99*100; 

% NSE 
 
y_win = [datenum(dates_winter), y_all_win];
yhat_win = [datenum(dates_winter), yhat_all_win]; 
y_sum = [datenum(dates_summer), y_all_sum];
yhat_sum = [datenum(dates_summer), yhat_all_sum]; 

NSE_win = nashsutcliffe(y_win, yhat_win); 
NSE_sum = nashsutcliffe(y_sum, yhat_sum); 

%NSE 99th 

y_win_99 = [datenum(dates_winter(y_99_win_indices)), y_99_win]; 
yhat_win_99 = [datenum(dates_winter(y_99_win_indices)), yhat_99_win]; 
y_sum_99 = [datenum(dates_summer(y_99_sum_indices)), y_99_sum]; 
yhat_sum_99 = [datenum(dates_summer(y_99_sum_indices)), yhat_99_sum]; 

NSE_win_99 = nashsutcliffe(y_win_99, yhat_win_99); 
NSE_sum_99 = nashsutcliffe(y_sum_99, yhat_sum_99); 


performance = [rmse_std_all.*100;correlations_all;...
    rmse_std_99;correlations_99];

 %% PDP winter- stacking some lines 

pdp_fig_win = figure(10), clf;
set(pdp_fig_win, 'Position', [150, 150, 1200, 400]);
sgtitle(title_name, 'FontSize', 26); 

subplot(1,4,1);
plotPartialDependence(RF_model_winter, 1); 
grid on
%ax = axes("Color","none")
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlabel('average temp [^oC]','FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 
ylim(pdp_win_lim);
set(gca,'FontSize',18)
title(''); 

subplot(1,4,2);
plotPartialDependence(RF_model_winter,2); 
hold on;
plotPartialDependence(RF_model_winter,3); 
plotPartialDependence(RF_model_winter,4); 
plotPartialDependence(RF_model_winter,5); 
plotPartialDependence(RF_model_winter,6); 
grid on
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
ylim(pdp_win_lim);
title('');
leg = legend('3-day rain', '30-day rain', '3-day snow',...
    '30-day snow', 'API', 'FontSize', 14); 
legend boxoff
set(leg, 'Location', 'best'); 
set(gca,'FontSize',18)
xlabel('precipitation [mm]','FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 


subplot(1,4,3); 
plotPartialDependence(RF_model_winter,7); 
hold on;
plotPartialDependence(RF_model_winter,8); 

grid on
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
ylim(pdp_win_lim);
title('');
leg = legend('SPI', 'SPI (month-1)','FontSize', 20); 
legend boxoff
set(leg, 'Location', 'best'); 
set(gca,'FontSize',18)
xlabel('SPI', 'FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 

subplot(1,4,4); 
plotPartialDependence(RF_model_winter,9); 
hold on;
plotPartialDependence(RF_model_winter,10); 
grid on
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
ylim(pdp_win_lim);
title('');
leg = legend('5-day melt', '30-day melt','FontSize', 20); 
legend boxoff
set(leg, 'Location', 'best'); 
set(gca,'FontSize',18)
xlabel('snow melt [mm]', 'FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 


saveas(pdp_fig_win, strcat(acronym, '_pdp_win_R1.png')); 

%% PDP winter NO SNOW- stacking some lines 

pdp_fig_win = figure(12), clf;
set(pdp_fig_win, 'Position', [150, 150, 1200, 400]);
sgtitle(title_name, 'FontSize', 26); 

subplot(1,3,1);
plotPartialDependence(RF_model_winter_NoSnow, 1); 
grid on
%ax = axes("Color","none")
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlabel('average temp [^oC]','FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 
ylim(pdp_win_lim);
set(gca,'FontSize',18)
title(''); 

subplot(1,3,2);
plotPartialDependence(RF_model_winter_NoSnow,2); 
hold on;
plotPartialDependence(RF_model_winter_NoSnow,3); 
api_plot = plotPartialDependence(RF_model_winter_NoSnow,4); 
ax = gca
grid on
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
ylim(pdp_win_lim);
title('');
leg = legend('3-day rain', '30-day rain','API', 'FontSize', 14); 
legend boxoff
set(leg, 'Location', 'best'); 
set(gca,'FontSize',18)
xlabel('precipitation [mm]','FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 


subplot(1,3,3); 
plotPartialDependence(RF_model_winter_NoSnow,5); 
hold on;
plotPartialDependence(RF_model_winter_NoSnow,6); 

grid on
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
ylim(pdp_win_lim);
title('');
leg = legend('SPI', 'SPI prvs','FontSize', 20); 
legend boxoff
set(leg, 'Location', 'best'); 
set(gca,'FontSize',18)
xlabel('SPI', 'FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 

saveas(pdp_fig_win, strcat(acronym, '_pdp_win_ws_R1.png')); 

%% PDP summer- stacking some lines 

pdp_fig_sum = figure(12), clf;
set(pdp_fig_sum, 'Position', [150, 150, 1200, 400]);
sgtitle(title_name, 'FontSize', 26); 

subplot(1,3,1);
plotPartialDependence(RF_model_summer, 1); 
grid on;
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
xlabel('average temp [^oC]','FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 
ylim(pdp_sum_lim);
set(gca,'FontSize',18)
title(''); 

subplot(1,3,2);
plotPartialDependence(RF_model_summer,2); 
hold on;
plotPartialDependence(RF_model_summer,3); 
plotPartialDependence(RF_model_summer,4);
ax= gca;
grid on
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
ylim(pdp_sum_lim);
title('');
leg = legend('3-day pre','30-day pre', 'API', 'FontSize', 20); 
legend boxoff
set(leg, 'Location', 'best'); 
set(gca,'FontSize',18)
xlabel('precipitation [mm]','FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 

subplot(1,3,3); 
plotPartialDependence(RF_model_summer,5); 
hold on;
plotPartialDependence(RF_model_summer,6); 
grid on
set(findall(gca, 'Type', 'Line'),'LineWidth',2);
ylim(pdp_sum_lim);
title('');
leg = legend('SPI', 'SPI (month-1)','FontSize', 18); 
legend boxoff
set(leg, 'Location', 'best'); 
set(gca,'FontSize',18)
xlabel('Relative soil moisture', 'FontSize', 20); 
ylabel('discharge [m^3/s]', 'FontSize', 20); 

saveas(pdp_fig_sum, strcat(acronym, '_pdp_sum_R1.png')); 
%%
% %% VARIABLE IMPORTANCE 
% 
% imp = RF_model_summer.OOBPermutedPredictorDeltaError;
% 
% predictor_names = {'av_temp', 'precip_3day', 'precip_30day', 'api', 'spi', ...
%     'spi_prvs'}; 
% 
% figure(10), clf
% bar(imp); 
% 
% title('Variable importance in flow prediction - summer', 'FontSize', 20);
% ylabel('Predictor importance estimates', 'FontSize', 20);
% xlim([0.5 6.5])
% xticks([1:6])
% h = gca;
% set(gca,'xticklabel',predictor_names, 'FontSize', 14);
% h.XTickLabelRotation = 45;
% h.TickLabelInterpreter = 'none';
% h.TitleFontSizeMultiplier = 1.7;
% h.LabelFontSizeMultiplier = 1.4;
% 
% %% WINTER 
% 
% imp = RF_model_winter.OOBPermutedPredictorDeltaError;
% 
% predictor_names = {'av_temp', 'precip_3day', 'precip_30day', 'snow_3day', 'snow_30day', ...
%     'melt_5day', 'melt_30day'}; 
% 
% figure(11), clf
% bar(imp); 
% 
% title('Variable importance in flow prediction - winter', 'FontSize', 20);
% ylabel('Predictor importance estimates', 'FontSize', 20);
% xlim([0.5 7.5])
% xticks([1:7])
% h = gca;
% set(gca,'xticklabel',predictor_names, 'FontSize', 14);
% h.XTickLabelRotation = 45;
% h.TickLabelInterpreter = 'none';
% h.TitleFontSizeMultiplier = 1.7;
% h.LabelFontSizeMultiplier = 1.4;
% 
% 
% %% Variable importances - summer
% 
% imp = RF_model_summer.OOBPermutedPredictorDeltaError;
% 
% predictor_names = {'precip_3day', 'spi', 'av_temp', 'api', ...
%     'spi_prvs', 'precip_30day'}; 
% 
% 
% imp_p30 = imp(3); 
% imp_spiprvs = imp(6); 
% imp_api = imp(4); 
% imp_temp = imp(1); 
% imp_spi = imp(5);
% imp_p3 = imp(2);
% 
% imp_ordered = [imp_p3, imp_spi, imp_temp, imp_api, ...
%     imp_spiprvs, imp_p30]; 
% 
% importance_fig_sum = figure(1),clf;
% set(importance_fig_sum, 'Position', [440, 378, 1000, 400]);
% 
% b = barh(imp_ordered);
% %b.FaceColor = 'flat';
% %b.CData(:,:) = [.09 0.56 .87;.09 0.56 .87;.09 0.56 .87;.09 0.56 .87];
% hold on
% 
% %title('Random Forest variable importance in flow prediction - summer', 'FontSize', 20);
% xlabel('Predictor importance estimates', 'FontSize', 20);
% ylim([0.5 6.5])
% yticks([1:6])
% h = gca;
% set(gca,'yticklabel',predictor_names, 'FontSize', 14);
% %h.YTickLabelRotation = 45;
% h.TickLabelInterpreter = 'none';
% h.TitleFontSizeMultiplier = 1.7;
% h.LabelFontSizeMultiplier = 1.4;
% %saveas(importance_fig_sum, strcat('importance_summer_V2_', acronym, '.png'));
% 
% %% Variable importances - winter
% 
% imp = RF_model_winter.OOBPermutedPredictorDeltaError;
% 
% predictor_names = {'snow_3day', 'melt_5day', 'snow_30day', ...
%     'precip_3day', 'av_temp', 'precip_30day', 'melt_30day'};
%     
% 
% 
% imp_m30 = imp(7); 
% imp_p30 = imp(3); 
% imp_temp = imp(1); 
% imp_p3 = imp(2); 
% imp_s30 = imp(5);
% imp_m5 = imp(6);
% imp_s3 = imp(4);
% 
% 
% imp_ordered = [imp_s3, imp_m5, imp_s30, imp_p3, ...
%     imp_temp, imp_p30, imp_m30]; 
% 
% importance_fig_sum = figure(1),clf;
% set(importance_fig_sum, 'Position', [440, 378, 1000, 400]);
% 
% b = barh(imp_ordered);
% %b.FaceColor = 'flat';
% %b.CData(:,:) = [.09 0.56 .87;.09 0.56 .87;.09 0.56 .87;.09 0.56 .87];
% hold on
% 
% %title('Random Forest variable importance in flow prediction - cold season', 'FontSize', 20);
% xlabel('Predictor importance estimates', 'FontSize', 20);
% ylim([0.5 7.5])
% yticks([1:7])
% h = gca;
% set(gca,'yticklabel',predictor_names, 'FontSize', 14);
% %h.YTickLabelRotation = 45;
% h.TickLabelInterpreter = 'none';
% h.TitleFontSizeMultiplier = 1.7;
% h.LabelFontSizeMultiplier = 1.4;
% %saveas(importance_fig_sum, strcat('importance_summer_V2_', acronym, '.png'));
% 
% 
% 
% 
%  