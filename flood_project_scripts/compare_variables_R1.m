clear all, close all

%%%%% CHANGE THIS %%%%%

river = 'Diamond'; 

%%%%% CHANGE NOTHING BELOW THIS LINE %%%%%

if(contains(river,'White'))
    load('RF_WRF_hist_WR_R1.mat'); 
    RF_historical_table = RF_WRF_hist_WR; 
    
    load('RF_WRF_fut_WR_R1.mat'); 
    RF_future_table = RF_WRF_fut_WR; 
   
    title_name = 'c) White River';
    acronym = 'WR';
    
elseif(contains(river,'Shenandoah'))
    load('RF_WRF_hist_SR_R1.mat'); 
    RF_historical_table = RF_WRF_hist_SR; 
    
    load('RF_WRF_fut_SR_R1.mat'); 
    RF_future_table = RF_WRF_fut_SR; 
    
    title_name = 'd) Shenandoah River';
    acronym = 'SR';
    
elseif(contains(river,'Mattawamkeag'))
 
    load('RF_WRF_hist_MR_R1.mat'); 
    RF_historical_table = RF_WRF_hist_MR; 
    
    load('RF_WRF_fut_MR_R1.mat'); 
    RF_future_table = RF_WRF_fut_MR; 
    
    title_name = 'a) Mattawamkeag River';
    acronym = 'MR';

elseif(contains(river,'Diamond'))
    load('RF_WRF_hist_DR_R1.mat'); 
    RF_historical_table = RF_WRF_hist_DR; 
    
    load('RF_WRF_fut_DR_R1.mat'); 
    RF_future_table = RF_WRF_fut_DR; 
   
    title_name = 'b) Dead Diamond River';
    acronym = 'DR';

end

RF_future_table = RF_future_table(23361:end,:); 

%% All RF variables 

%hist_data_win = RF_historical_table(:, [4 5 6 7 11 12]); 
%hist_data_sum = RF_historical_table(:, [4 5 7:10]);
hist_data_win = RF_historical_table(:, [4 5 7 8 10 11 12 13 14 15]); 
hist_data_sum = RF_historical_table(:, [4 5 7 11 12 13]);


winter_indices_hist = find(RF_historical_table.Month == 1| ...
    RF_historical_table.Month ==2|RF_historical_table.Month ==3 |...
    RF_historical_table.Month ==4 | RF_historical_table.Month ==5 |...
    RF_historical_table.Month == 11 | RF_historical_table.Month ==12);

summer_indices_hist = find(RF_historical_table.Month == 6 | RF_historical_table.Month ==7 | ...
    RF_historical_table.Month ==8 | RF_historical_table.Month == 9 | ...
    RF_historical_table.Month ==10);

hist_data_win = hist_data_win(winter_indices_hist,:); 
hist_data_sum = hist_data_sum(summer_indices_hist,:); 


fut_data_win = RF_future_table(:, [4 5 7 8 10 11 12 13 14 15]); 
fut_data_sum = RF_future_table(:, [4 5 7 11 12 13]); 

winter_indices_fut = find(RF_future_table.Month == 1| ...
    RF_future_table.Month ==2|RF_future_table.Month ==3 |...
    RF_future_table.Month ==4 | RF_future_table.Month ==5 |...
    RF_future_table.Month == 11 | RF_future_table.Month ==12);

summer_indices_fut = find(RF_future_table.Month == 6 | RF_future_table.Month ==7 | ...
    RF_future_table.Month ==8 | RF_future_table.Month == 9 | ...
    RF_future_table.Month ==10);

fut_data_win = fut_data_win(winter_indices_fut,:); 
fut_data_sum = fut_data_sum(summer_indices_fut,:); 

xlabel_vec_win = {'Avg temp [^oC]','3-day pre [mm]', ...
'30-day pre [mm]', '3-day snow [mm]', '30-day snow [mm]','API [mm]', 'SPI', ...
    'SPI month-1','5-day melt [mm]', '30-day melt [mm]'}; 


% xlabel_vec_win = {'Avg temp','3-day rain', '30-day rain', ...
%     '3-day snow', '30-day snow','API', 'SPI', 'SPI month-1', ...
%     '5-day melt', '30-day melt'}; 
xlabel_vec_sum = {'Avg temp [^oC]', '3-day rain [mm]', '30-day rain [mm]', ...
    'API [mm]', 'SPI', 'SPI (month-1)'}; 
% xlabel_vec_sum = {'Avg temp', '3-day pre', '30-day pre', ...
%     'API', 'SPI', 'SPI month-1'}; 


fut_data_win.melt_5day = fut_data_win.melt_5day.*25.4; 
fut_data_win.melt_30day = fut_data_win.melt_30day.*25.4; 

hist_data_win.melt_5day = hist_data_win.melt_5day.*25.4;
hist_data_win.melt_30day = hist_data_win.melt_30day.*25.4;


%% WINTER COMPARISON PLOTS 

winter_hist_fig = figure(1); clf; 
set(winter_hist_fig, 'Position', [100 100 433 724]); 
sgtitle(strcat(title_name), 'FontSize', 22); 


%TEMPERATURE 
subplot(5,2,1); 
bin_width = 5; 
this_data_hist = table2array((hist_data_win(:,1)));
this_data_fut = table2array((fut_data_win(:,1)));
perc_change = (nanmean(this_data_fut)-nanmean(this_data_hist));%/...
    %abs(nanmean(this_data_hist));
perc_change = round(perc_change,2);%*100;
p = histogram((this_data_hist), 'BinWidth', bin_width); 
hold on;    
p2 = histogram((this_data_fut), 'BinWidth', bin_width); 
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance(1,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end

set(gca,'FontSize',13);

xlabel(xlabel_vec_win(1)); 
if(perc_change>0)
    title(strcat('+',num2str(perc_change), ' ^oC', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), ' ^oC', asterisk),...
        'FontSize', 14); 
end

%3DAY PRECIP 
test = subplot(5,2,2); 
this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,2))>=3),2));
this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,2))>=3),2));
hist_omitted = length(find(table2array(hist_data_win(:,2))<3)); 
fut_omitted = length(find(table2array(fut_data_win(:,2))<3)); 
omitted_win(1,:) = [hist_omitted, fut_omitted]; 
perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;
histogram(this_data_hist, 'BinWidth', bin_width);
hold on;
histogram(this_data_fut, 'BinWidth', bin_width);
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance(2,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13);
xlabel(xlabel_vec_win(2)); 
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),'FontSize', 16); 
else
    title(strcat(num2str(perc_change), '%', asterisk),'FontSize', 16); 

end
axes('Position',[0.76 0.823 0.08 0.035]);

max_hist = max(this_data_hist); 
max_fut = max(this_data_fut); 
max_val = max(max_hist, max_fut); 

top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',10);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
histogram(top_data_fut, 'BinWidth', 10);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 

ks_top = kstest2(top_data_hist, top_data_fut); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
end

%30 DAY PRECIP
test = subplot(5,2,3); 
bin_width = 20;
this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,3))>=30),3));
this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,3))>=30),3));
hist_omitted = length(find(table2array(hist_data_win(:,3))<3)); 
fut_omitted = length(find(table2array(fut_data_win(:,3))<3)); 
omitted_win(2,:) = [hist_omitted, fut_omitted]; 
perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;
histogram(this_data_hist, 'BinWidth', bin_width);
hold on
histogram(this_data_fut, 'BinWidth', bin_width);
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance(3,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end

set(gca,'FontSize',13)

xlabel(xlabel_vec_win(3));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
end

axes('Position',[0.35 0.66 0.065 0.036]);
max_hist = max(this_data_hist); 
max_fut = max(this_data_fut); 
max_val = max(max_hist, max_fut);
top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',10);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
histogram(top_data_fut, 'BinWidth', 10);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
ks_top = kstest2(top_data_hist, top_data_fut); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
end

% 3 DAY SNOW
test = subplot(5,2,4);
bin_width = 5; 
this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,4))>=3),4));
this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,4))>=3),4));
hist_omitted = length(find(table2array(hist_data_win(:,4))<3)); 
fut_omitted = length(find(table2array(fut_data_win(:,4))<3));
omitted_win(3,:) = [hist_omitted, fut_omitted]; 
perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;
histogram(this_data_hist, 'BinWidth', bin_width);
hold on
histogram(this_data_fut, 'BinWidth', bin_width);
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance(4,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)

xlabel(xlabel_vec_win(4));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
end

axes('Position',[0.785 0.66 0.1 0.036]);
top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',5);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
histogram(top_data_fut, 'BinWidth', 5);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
ks_top = kstest2(top_data_hist, top_data_fut); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
end

% 30 DAY SNOW
test = subplot(5,2,5);
bin_width = 10;
this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,5))>=30),5));
this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,5))>=30),5));
hist_omitted = length(find(table2array(hist_data_win(:,5))<3)); 
fut_omitted = length(find(table2array(fut_data_win(:,5))<3));
omitted_win(4,:) = [hist_omitted, fut_omitted]; 
perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;
histogram((this_data_hist), 'BinWidth', bin_width); 
hold on    
histogram((this_data_fut), 'BinWidth', bin_width); 
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance(5,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)

xlabel(xlabel_vec_win(5));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
end

axes('Position',[0.35 0.5 0.08 0.032]);
top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',10);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
histogram(top_data_fut, 'BinWidth', 10);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
ks_top = kstest2(top_data_hist, top_data_fut); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
end

% 5 DAY MELT 
test = subplot(5,2,6); 
if(contains(river, 'White'))
    bin_width = 2; 
else
    bin_width =0.2;
end

if(contains(river, 'Shenandoah'))
    threshold = 0.5;
else
    threshold = 3;
end
this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,9))>=threshold),9));
this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,9))>=threshold),9));
hist_omitted = length(find(table2array(hist_data_win(:,9))<3)); 
fut_omitted = length(find(table2array(fut_data_win(:,9))<3)); 
omitted_win(5,:) = [hist_omitted, fut_omitted]; 
histogram((this_data_hist), 'BinWidth', bin_width); 
hold on    
histogram((this_data_fut), 'BinWidth', bin_width);
set(gca,'FontSize',13)
perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;

[ks_all, p_all] = kstest2(this_data_hist, this_data_fut);
significance(6,:) = [perc_change, ks_all, p_all]; 
if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end

xlabel(xlabel_vec_win(9));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
end


axes('Position',[0.78 0.495 0.08 0.035]);
max_hist = max(this_data_hist); 
max_fut = max(this_data_fut); 
max_val = max(max_hist, max_fut);
if(contains(river, 'White')); 
    inset_size = 200;
    bin_width = 3;
else
    inset_size = 50;
    bin_width = 2;
end
top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
if(isempty(top_data_fut))
    ks_top = 0; 
else
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
end

if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
end


% 30 DAY MELT 
test = subplot(5,2,7); 
if(contains(river, 'White'))
    bin_width = 5; 
else
    bin_width =0.5;
end

if(contains(river, 'Shenandoah'))
    threshold = 2;
else
    threshold = 5;
end
this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,10))>=threshold),10));
this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,10))>=threshold),10));
hist_omitted = length(find(table2array(hist_data_win(:,10))<3)); 
fut_omitted = length(find(table2array(fut_data_win(:,10))<3)); 
omitted_win(6,:) = [hist_omitted, fut_omitted]; 
perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;
p = histogram((this_data_hist), 'BinWidth', bin_width); 
hold on    
p2 = histogram((this_data_fut), 'BinWidth', bin_width); 
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut);
significance(7,:) = [perc_change, ks_all, p_all]; 
if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)

xlabel(xlabel_vec_win(10));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
end

axes('Position',[0.385 0.325 0.06 0.03]);
if(contains(river, 'White')); 
    inset_size = 1000;
    bin_width = 5;
elseif(contains(river, 'Diamond'))
    inset_size = 250;
    bin_width = 5;
else
    inset_size = 50;
    bin_width = 1;
end
top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
if(isempty(top_data_fut))
    ks_top = 0; 
else
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
end
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
end

% API
test = subplot(5,2,8);
if(contains(river, 'Mattawamkeag'))
    bin_width = 20; 
else
    bin_width =5;
end

this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,6))),6));
this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,6))),6));


perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;
p = histogram((this_data_hist), 'BinWidth', bin_width); 
hold on    
p2 = histogram((this_data_fut), 'BinWidth', bin_width); 
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance(8,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)

xlabel(xlabel_vec_win(6));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
end

legend([p, p2],'Historical', 'Future','Orientation',...
    'horizontal', 'FontSize', 16, 'Position',[0.42 0.04 0.15 0.03]);

axes('Position',[0.82 0.325 0.06 0.03]);
if(contains(river, 'White')); 
    bin_width = 10;
elseif(contains(river, 'Diamond'))
    bin_width = 10;
else
    bin_width = 10;
end
top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
ks_top = kstest2(top_data_hist, top_data_fut); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end
if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 10); 
end

% SPI
test = subplot(5,2,9); 
bin_width = 0.25; 
this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,7))),7));
this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,7))),7));


perc_change = (nanmean(this_data_fut)-nanmean(this_data_hist));%/...
    %abs(nanmean(this_data_hist));
perc_change = round(perc_change,2);%*100;
p = histogram((this_data_hist), 'BinWidth', bin_width); 
hold on    
p2 = histogram((this_data_fut), 'BinWidth', bin_width); 
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance(9,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)

xlabel(xlabel_vec_win(7));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), asterisk),...
        'FontSize', 14); 
end


% SPI (month-1)
test = subplot(5,2,10); 
bin_width = 0.25; 
this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,8))),8));
this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,8))),8));


perc_change = (nanmean(this_data_fut)-nanmean(this_data_hist));%/...
    %abs(nanmean(this_data_hist));
perc_change = round(perc_change,2);%*100;
p = histogram((this_data_hist), 'BinWidth', bin_width); 
hold on    
p2 = histogram((this_data_fut), 'BinWidth', bin_width); 
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance(10,:) = [perc_change, ks_all, p_all]; 
if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)

xlabel(xlabel_vec_win(8));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), asterisk),...
        'FontSize', 14); 
end


legend([p, p2],'Historical', 'Future','Orientation',...
    'horizontal', 'FontSize', 16, 'Position',[0.42 0.04 0.15 0.03]);



saveas(winter_hist_fig, strcat(acronym, '_RFvars_percchange_win_R1.png')); 

%% 
summer_hist_fig = figure(4);
clf; 
set(summer_hist_fig, 'Position', [100 100 433 724]); 
sgtitle(strcat(title_name), 'FontSize', 22);

%TEMPERATURE 
subplot(3,2,1); 
bin_width =5;
this_data_hist = table2array(hist_data_sum(:,1));
this_data_fut = table2array(fut_data_sum(:,1));
perc_change = (nanmean(this_data_fut)-nanmean(this_data_hist));%/...
    %abs(nanmean(this_data_hist));
perc_change = round(perc_change,2);%*100;
p = histogram(this_data_hist, 'BinWidth', bin_width); 
hold on;
p2 = histogram(this_data_fut, 'BinWidth', bin_width); 
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance_sum(1,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)

xlabel(xlabel_vec_sum(1));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), ' ^oC', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), ' ^oC', asterisk),...
        'FontSize', 14); 
end

% 3DAY PRECIP
subplot(3,2,2); 
bin_width = 5;
this_data_hist = table2array(hist_data_sum(find(table2array(hist_data_sum(:,2))>=3),2));
this_data_fut = table2array(fut_data_sum(find(table2array(fut_data_sum(:,2))>=3),2));
hist_omitted = length(find(table2array(hist_data_sum(:,2))<3)); 
fut_omitted = length(find(table2array(fut_data_sum(:,2))<3));
omitted_sum(1,:) = [hist_omitted, fut_omitted]; 
perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;
histogram(this_data_hist, 'BinWidth', bin_width);
hold on
histogram(this_data_fut, 'BinWidth', bin_width);
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance_sum(2,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)
xlabel(xlabel_vec_sum(2));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
end

axes('Position',[0.75 0.73 0.15 0.1]);
max_hist = max(this_data_hist); 
max_fut = max(this_data_fut); 
max_val = max(max_hist, max_fut); 
if(contains(river, 'White')|contains(river, 'Diamond')); 
    inset_size = 125;
    bin_width = 10;
else
    inset_size = 50;
    bin_width = 10;
end
top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
ks_top = kstest2(top_data_hist, top_data_fut); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 12); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 12); 
end


% 30DAY PRECIP
subplot(3,2,3); 
bin_width = 20;
this_data_hist = table2array(hist_data_sum(find(table2array(hist_data_sum(:,3))>=30),3));
this_data_fut = table2array(fut_data_sum(find(table2array(fut_data_sum(:,3))>=30),3));
hist_omitted = length(find(table2array(hist_data_sum(:,3))<3)); 
fut_omitted = length(find(table2array(fut_data_sum(:,3))<3));
omitted_sum(2,:) = [hist_omitted, fut_omitted]; 
perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;
histogram(this_data_hist, 'BinWidth', bin_width);
hold on
histogram(this_data_fut, 'BinWidth', bin_width);
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance_sum(3,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end

set(gca,'FontSize',13)
xlabel(xlabel_vec_sum(3));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
end

axes('Position',[0.33 0.46 0.12 0.1]);
top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
ks_top = kstest2(top_data_hist, top_data_fut); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 12); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 12); 
end


% API 
subplot(3,2,4); 
bin_width = 8; 
this_data_hist = table2array(hist_data_sum(:,4));
this_data_fut = table2array(fut_data_sum(:,4));
histogram(this_data_hist, 'BinWidth', bin_width);
hold on
histogram(this_data_fut, 'BinWidth', bin_width);
set(gca,'FontSize',13)
perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
    abs(nansum(this_data_hist));
perc_change = round(perc_change,2)*100;
if(contains(river, 'Mattawamkeag'))
    xlim([0 300])
end
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance_sum(4,:) = [perc_change, ks_all, p_all]; 
if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)
xlabel(xlabel_vec_sum(4));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), '%', asterisk),...
        'FontSize', 14); 
end
max_hist = max(this_data_hist); 
max_fut = max(this_data_fut); 
max_val = max(max_hist, max_fut); 
if(contains(river, 'White') | contains(river,'Diamond')); 
    inset_size = 125;
    bin_width = 10;
    axes('Position',[0.74 0.46 0.15 0.1]);
elseif(contains(river, 'Mattawamkeag')) 
    axes('Position',[0.785 0.47 0.1 0.08]);
else
    axes('Position',[0.74 0.46 0.15 0.1]);
    inset_size = 50;
    bin_width = 10;
end
top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
    abs(nansum(top_data_hist))*100;
top_perc_change = round(top_perc_change,1)
histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
hold on
histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
ks_top = kstest2(top_data_hist, top_data_fut); 
if(ks_top ==1) 
    asterisk = '*';
elseif(ks_top ==0);
    asterisk = '';
end

if(top_perc_change>0)
    title(strcat('+',num2str(top_perc_change), '%', asterisk),...
        'FontSize', 12); 
else
    title(strcat(num2str(top_perc_change), '%', asterisk),...
        'FontSize', 12); 
end

% SPI 
subplot(3,2,5); 
bin_width = 0.5;
this_data_hist = table2array(hist_data_sum(:,5));
this_data_fut = table2array(fut_data_sum(:,5));
perc_change = (nanmean(this_data_fut)-nanmean(this_data_hist));%/...
    %abs(nanmean(this_data_hist));
perc_change = round(perc_change,3);%*100;
p = histogram(this_data_hist, 'BinWidth', bin_width); 
hold on;
p2 = histogram(this_data_fut, 'BinWidth', bin_width); 
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance_sum(5,:) = [perc_change, ks_all, p_all]; 
if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)
xlabel(xlabel_vec_sum(5));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), asterisk),...
        'FontSize', 14); 
end


% SPI PRVS
subplot(3,2,6); 
bin_width = 0.5;
this_data_hist = table2array(hist_data_sum(:,6));
this_data_fut = table2array(fut_data_sum(:,6));
perc_change = (nanmean(this_data_fut)-nanmean(this_data_hist));%/...
    %abs(nanmean(this_data_hist));
perc_change = round(perc_change,2);%*100;
p = histogram(this_data_hist, 'BinWidth', bin_width); 
hold on;
p2 = histogram(this_data_fut, 'BinWidth', bin_width); 
[ks_all, p_all] = kstest2(this_data_hist, this_data_fut); 
significance_sum(6,:) = [perc_change, ks_all, p_all]; 

if(ks_all ==1) 
    asterisk = '*';
elseif(ks_all ==0);
    asterisk = '';
end
set(gca,'FontSize',13)
xlabel(xlabel_vec_sum(6));
if(perc_change>0)
    title(strcat('+',num2str(perc_change), asterisk),...
        'FontSize', 14); 
else
    title(strcat(num2str(perc_change), asterisk),...
        'FontSize', 14); 
end
legend([p, p2],'Historical', 'Future','Orientation','horizontal', ...
    'FontSize', 16, 'Position',[0.42 0.02 0.15 0.03]);

saveas(summer_hist_fig, strcat(acronym, '_RFvars_percchange_sum_R1.png')); 

%% KS TESTS 

%test_fig = figure(3), clf
%set(test_fig, 'Position', [100 100 433 724]); 

for i = 1:7 
    
    [ks_test_winter(i, 1), pval_win(i,1)] = kstest2(table2array(hist_data_win(:,i)), ...
        table2array(fut_data_win(:,i))); 
    [t_test_win_hist(i,1), p_win_hist(i,1)] = kstest(table2array(hist_data_win(:,i)));
    [t_test_win_fut(i,1), p_win_fut(i,1)] = kstest(table2array(fut_data_win(:,i)));

    if i ==7 
        break;
    end
    
    [ks_test_summer(i, 1), pval_sum(i,1)] = kstest2(table2array(hist_data_sum(:,i)), ...
        table2array(fut_data_sum(:,i))); 
    [t_test_sum_hist(i,1), p_sum_hist(i,1)] = kstest(table2array(hist_data_sum(:,i)));
    [t_test_sum_fut(i,1), p_sum_fut(i,1)] = kstest(table2array(fut_data_sum(:,i)));

%     subplot(3,2,i); 
%     this_data_hist = table2array(hist_data_sum(:,i));
%     this_data_fut = table2array(fut_data_sum(:,i));
%     histogram(this_data_hist); 
%     hold on 
%     histogram(this_data_fut); 
%     

end
%%
ks_test_summer(7,1) = 0; 
pval_sum(7,1) = 0;


ks_tests = [ks_test_winter(:,1), pval_win(:,1), ks_test_summer, pval_sum ]; 

t_test_win = [t_test_win_hist, p_win_hist, t_test_win_fut, p_win_fut];
t_test_sum = [t_test_sum_hist, p_sum_hist, t_test_sum_fut, p_sum_fut];

%% PICKING TOP 2 PER WATERSHED - cold season

% WR and SR- 3 day rain and API
if(contains(river, 'White')|contains(river, 'Shenandoah'))
    cs_top2_fig = figure(14), clf; 
    set(cs_top2_fig, 'Position', [150 150 800 350]); 
    sgtitle(title_name, 'FontSize', 22);
    
    %3 DAY PRECIP
    subplot(1,2,1) 
    bin_width = 5;
    this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,2))>=3),2));
    this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,2))>=3),2));
    perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
        abs(nansum(this_data_hist));
    perc_change = round(perc_change,2)*100;
    histogram(this_data_hist, 'BinWidth', bin_width);
    hold on
    histogram(this_data_fut, 'BinWidth', bin_width);
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0);
        asterisk = '';
    end
    set(gca,'FontSize',18)
    
    xlabel(xlabel_vec_win(2));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    end
    
    ylabel('Days');
    legend('Historical', 'Future', 'FontSize', 18);
    axes('Position',[0.285 0.22 0.17 0.34]);
    top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
    top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
    top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
        abs(nansum(top_data_hist))*100;
    top_perc_change = round(top_perc_change,1)
    histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    hold on
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
    if(ks_top ==1) 
        asterisk = '*';
    elseif(ks_top ==0);
        asterisk = '';
    end

    if(top_perc_change>0)
        title(strcat('+',num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    else
        title(strcat(num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    end

     % API 
    subplot(1,2,2); 
    bin_width = 5; 
    this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,6))),6));
    this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,6))),6));
    perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
        abs(nansum(this_data_hist));
    perc_change = round(perc_change,2)*100;
    p = histogram((this_data_hist), 'BinWidth', bin_width); 
    hold on    
    p2 = histogram((this_data_fut), 'BinWidth', bin_width);
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0);
        asterisk = '';
    end
    set(gca,'FontSize',18);
    xlabel(xlabel_vec_win(6));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    end
    ylabel('Days'); 
    legend('Historical', 'Future', 'FontSize', 18);
    
    axes('Position',[0.72 0.29 0.175 0.25]);
    max_hist = max(this_data_hist); 
    max_fut = max(this_data_fut); 
    bin_width = 5;
    top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
    top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
    top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
        abs(nansum(top_data_hist))*100;
    top_perc_change = round(top_perc_change,1)
    histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    hold on
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
    if(ks_top ==1) 
        asterisk = '*';
    elseif(ks_top ==0);
        asterisk = '';
    end

    if(top_perc_change>0)
        title(strcat('+',num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    else
        title(strcat(num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    end

elseif(contains(river, 'Mattawamkeag'))
    cs_top2_fig = figure(14), clf 
    set(cs_top2_fig, 'Position', [150 150 800 350]); 
    sgtitle(title_name, 'FontSize', 22);
    
    %30 DAY PRECIP
    subplot(1,2,1); 
    bin_width = 20;
    this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,3))>=30),3));
    this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,3))>=30),3));
    hist_omitted = length(find(table2array(hist_data_win(:,3))<3)); 
    fut_omitted = length(find(table2array(fut_data_win(:,3))<3)); 
    omitted_win(2,:) = [hist_omitted, fut_omitted]; 
    perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
        abs(nansum(this_data_hist));
    perc_change = round(perc_change,2)*100;
    histogram(this_data_hist, 'BinWidth', bin_width);
    hold on
    histogram(this_data_fut, 'BinWidth', bin_width);
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0);
        asterisk = '';
    end
    set(gca,'FontSize',18)
    xlabel(xlabel_vec_win(3));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    end
    ylabel('Days');
    legend('Historical', 'Future', 'FontSize', 18);
    
    axes('Position',[0.315 0.28 0.14 0.28]);
    top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
    top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
    top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
        abs(nansum(top_data_hist))*100;
    top_perc_change = round(top_perc_change,1)
    histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    hold on
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
    if(ks_top ==1) 
        asterisk = '*';
    elseif(ks_top ==0);
        asterisk = '';
    end

    if(top_perc_change>0)
        title(strcat('+',num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    else
        title(strcat(num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    end


    % API 
    subplot(1,2,2); 
    bin_width = 20; 
    this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,6))),6));
    this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,6))),6));
    perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
        abs(nansum(this_data_hist));
    perc_change = round(perc_change,2)*100;
    p = histogram((this_data_hist), 'BinWidth', bin_width); 
    hold on    
    p2 = histogram((this_data_fut), 'BinWidth', bin_width);
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0);
        asterisk = '';
    end
    set(gca,'FontSize',18);
    xlabel(xlabel_vec_win(6));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    end
    ylabel('Days'); 
    legend('Historical', 'Future', 'FontSize', 18);
    
    axes('Position',[0.8 0.38 0.1 0.2]);
    bin_width = 10;
    top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
    top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
    top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
        abs(nansum(top_data_hist))*100;
    top_perc_change = round(top_perc_change,1)
    histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    hold on
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
    if(ks_top ==1) 
        asterisk = '*';
    elseif(ks_top ==0);
        asterisk = '';
    end

    if(top_perc_change>0)
        title(strcat('+',num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    else
        title(strcat(num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    end

    
elseif(contains(river, 'Diamond'))
    cs_top2_fig = figure(14), clf 
    set(cs_top2_fig, 'Position', [150 150 800 350]); 
    sgtitle(title_name, 'FontSize', 22);
    %3 DAY PRECIP
    subplot(1,2,1);
    bin_width = 5;
    this_data_hist = table2array(hist_data_win(find(table2array(hist_data_win(:,2))>=3),2));
    this_data_fut = table2array(fut_data_win(find(table2array(fut_data_win(:,2))>=3),2));
    perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
        abs(nansum(this_data_hist));
    perc_change = round(perc_change,2)*100;
    histogram(this_data_hist, 'BinWidth', bin_width);
    hold on
    histogram(this_data_fut, 'BinWidth', bin_width);
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0);
        asterisk = '';
    end
    set(gca,'FontSize',18)
    xlabel(xlabel_vec_win(2));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    end
    ylabel('Days');
    legend('Historical', 'Future', 'FontSize', 18);
    
    axes('Position',[0.305 0.25 0.15 0.3]);
    top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
    top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
    top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
        abs(nansum(top_data_hist))*100;
    top_perc_change = round(top_perc_change,1)
    histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    hold on
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
    if(ks_top ==1) 
        asterisk = '*';
    elseif(ks_top ==0);
        asterisk = '';
    end

    if(top_perc_change>0)
        title(strcat('+',num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    else
        title(strcat(num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    end
   


    %TEMPERATURE 
    subplot(1,2,2); 
    bin_width = 5; 
    this_data_hist = table2array((hist_data_win(:,1)));
    this_data_fut = table2array((fut_data_win(:,1)));
    perc_change = (nanmean(this_data_fut)-nanmean(this_data_hist));%/...
        %abs(nanmean(this_data_hist));
    perc_change = round(perc_change,2);%*100;
    p = histogram((this_data_hist), 'BinWidth', bin_width); 
    hold on;    
    p2 = histogram((this_data_fut), 'BinWidth', bin_width); 
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0)
        asterisk = '';
    end
    set(gca,'FontSize',18)
    legend('Historical', 'Future', 'FontSize', 18);
    xlabel(xlabel_vec_win(1));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change), '%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change), '%',asterisk),...
            'FontSize', 20); 
    end 
    ylabel('Days');
end

saveas(cs_top2_fig, strcat(acronym, '_RFvars_top2_win_R1.png')); 

%% WARM SEASON TOP 2 

if(contains(river, 'Mattawamkeag'))
    
    ws_top2_fig = figure(14), clf; 
    set(ws_top2_fig, 'Position', [150 150 800 350]); 
    sgtitle(title_name, 'FontSize', 22);
    
    % API 
    subplot(1,2,1); 
    bin_width = 8; 
    this_data_hist = table2array(hist_data_sum(:,4));
    this_data_fut = table2array(fut_data_sum(:,4));
    histogram(this_data_hist, 'BinWidth', bin_width);
    hold on
    histogram(this_data_fut, 'BinWidth', bin_width);
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0);
        asterisk = '';
    end
    set(gca,'FontSize',18)
    perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
        abs(nansum(this_data_hist));
    perc_change = round(perc_change,2)*100;

    xlabel(xlabel_vec_sum(4));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    end
    
    legend('Historical', 'Future', 'FontSize', 18); 
    ylabel('Days');
    max_hist = max(this_data_hist); 
    max_fut = max(this_data_fut); 
    max_val = max(max_hist, max_fut); 
    if(contains(river, 'White') | contains(river,'Diamond')); 
        inset_size = 125;
        bin_width = 10;
    else
        inset_size = 50;
        bin_width = 10;
    end
    
    axes('Position',[0.345 0.35 0.11 0.2]);
    top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
    top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
    top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
        abs(nansum(top_data_hist))*100;
    top_perc_change = round(top_perc_change,1)
    histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    hold on
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
    if(ks_top ==1) 
        asterisk = '*';
    elseif(ks_top ==0);
        asterisk = '';
    end

    if(top_perc_change>0)
        title(strcat('+',num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    else
        title(strcat(num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    end



    % 30DAY PRECIP
    subplot(1,2,2); 
    bin_width = 20;
    this_data_hist = table2array(hist_data_sum(find(table2array(hist_data_sum(:,3))>=30),3));
    this_data_fut = table2array(fut_data_sum(find(table2array(fut_data_sum(:,3))>=30),3));
    perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
        abs(nansum(this_data_hist));
    perc_change = round(perc_change,2)*100;
    histogram(this_data_hist, 'BinWidth', bin_width);
    hold on
    histogram(this_data_fut, 'BinWidth', bin_width);
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0);
        asterisk = '';
    end
    set(gca,'FontSize',18)
    xlabel(xlabel_vec_sum(3));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    end
    
    legend('Historical', 'Future', 'FontSize', 18);
    ylabel('Days');
    axes('Position',[0.775 0.29 0.12 0.28]);
    max_hist = max(this_data_hist); 
    max_fut = max(this_data_fut); 
    max_val = max(max_hist, max_fut); 
    top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
    top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
    top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
        abs(nansum(top_data_hist))*100;
    top_perc_change = round(top_perc_change,1)
    histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    hold on
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
    if(ks_top ==1) 
        asterisk = '*';
    elseif(ks_top ==0);
        asterisk = '';
    end

    if(top_perc_change>0)
        title(strcat('+',num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    else
        title(strcat(num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    end


    
else
        
    ws_top2_fig = figure(14), clf; 
    set(ws_top2_fig, 'Position', [150 150 800 350]); 
    sgtitle(title_name, 'FontSize', 22);
    
    % API 
    subplot(1,2,1); 
    bin_width = 8; 
    this_data_hist = table2array(hist_data_sum(:,4));
    this_data_fut = table2array(fut_data_sum(:,4));
    histogram(this_data_hist, 'BinWidth', bin_width);
    hold on
    histogram(this_data_fut, 'BinWidth', bin_width);
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0);
        asterisk = '';
    end
    set(gca,'FontSize',18)
    perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
        abs(nansum(this_data_hist));
    perc_change = round(perc_change,2)*100;
    xlabel(xlabel_vec_sum(4));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    end
    
    legend('Historical', 'Future', 'FontSize', 18); 
    ylabel('Days');
    max_hist = max(this_data_hist); 
    max_fut = max(this_data_fut); 
    max_val = max(max_hist, max_fut); 
    if(contains(river, 'White') | contains(river,'Diamond')); 
        inset_size = 125;
        bin_width = 10;
    else
        inset_size = 50;
        bin_width = 10;
    end
    if(contains(river, 'White')|contains(river, 'Diamond'))
        axes('Position',[0.285 0.22 0.17 0.34]);
    else
        axes('Position',[0.32 0.22 0.135 0.27]);
    end
    top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
    top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
    top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
        abs(nansum(top_data_hist))*100;
    top_perc_change = round(top_perc_change,1)
    histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    hold on
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
    if(ks_top ==1) 
        asterisk = '*';
    elseif(ks_top ==0);
        asterisk = '';
    end

    if(top_perc_change>0)
        title(strcat('+',num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    else
        title(strcat(num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    end


    %3 DAY PRECIP
    subplot(1,2,2) 
    bin_width = 5;
    this_data_hist = table2array(hist_data_sum(find(table2array(hist_data_sum(:,2))>3),2));
    this_data_fut = table2array(fut_data_sum(find(table2array(fut_data_sum(:,2))>3),2));
    perc_change = (nansum(this_data_fut)-nansum(this_data_hist))/...
        abs(nansum(this_data_hist));
    perc_change = round(perc_change,2)*100;
    histogram(this_data_hist, 'BinWidth', bin_width);
    hold on
    histogram(this_data_fut, 'BinWidth', bin_width);
    ks_all = kstest2(this_data_hist, this_data_fut); 
    if(ks_all ==1) 
    asterisk = '*';
    elseif(ks_all ==0);
        asterisk = '';
    end
    set(gca,'FontSize',18)
    xlabel(xlabel_vec_sum(2));
    if(perc_change>0)
        title(strcat('+',num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    else
        title(strcat(num2str(perc_change),'%', asterisk),...
            'FontSize', 20); 
    end
    
    ylabel('Days');
    legend('Historical', 'Future', 'FontSize', 18);
    axes('Position',[0.725 0.22 0.17 0.34]);
    max_hist = max(this_data_hist); 
    max_fut = max(this_data_fut); 
    max_val = max(max_hist, max_fut); 
    top_data_hist = maxk(this_data_hist, round(length(this_data_hist)*0.01)); 
    top_data_fut = this_data_fut(this_data_fut>min(top_data_hist)); 
    top_perc_change = (nansum(top_data_fut)-nansum(top_data_hist))/...
        abs(nansum(top_data_hist))*100;
    top_perc_change = round(top_perc_change,1)
    histogram(top_data_hist, 'BinWidth',bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    hold on
    histogram(top_data_fut, 'BinWidth', bin_width);%, 'BinLimits', [max_val-100, max_val], 'BinWidth', 10); 
    ks_top = kstest2(top_data_hist, top_data_fut); 
    if(ks_top ==1) 
        asterisk = '*';
    elseif(ks_top ==0);
        asterisk = '';
    end

    if(top_perc_change>0)
        title(strcat('+',num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    else
        title(strcat(num2str(top_perc_change), '%', asterisk),...
            'FontSize', 16); 
    end


end

saveas(ws_top2_fig, strcat(acronym, '_RFvars_top2_sum_R1.png')); 
