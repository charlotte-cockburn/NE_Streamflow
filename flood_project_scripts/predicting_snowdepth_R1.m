%%  Section 1: Settings based on which river we pick 
clear all, close all
river = 'White'; %SET THIS VALUE!!! 

if(contains(river,'White'))

    min_year = 1915;
    load('White_River/Livneh_data_white.mat');
    Livneh_data = Livneh_data_white;
    load('WR_historical_data.mat'); 
    historical_data = WR_historical_data; 
    load('WR_future_data.mat'); 
    future_data = WR_future_data; 
    load('WR_all_snow.mat'); 
    load('WR_avg_snow.mat');
    
    acronym = 'WR';
    
elseif(contains(river,'Shenandoah')) 
    min_year = 1926;
    load('Shenandoah_River/Livneh_data_shenandoah.mat');
    Livneh_data = Livneh_data_shenandoah(4019:35429,:); 
    load('SR_historical_data.mat'); 
    historical_data = SR_historical_data; 
    load('SR_future_data.mat'); 
    future_data = SR_future_data; 
    acronym = 'SR';
    load('SR_all_snow.mat'); 
    load('SR_avg_snow.mat');
elseif(contains(river,'Mattawamkeag'))
    min_year = 1935;
    load('Mattawamkeag_River/Livneh_data_mattawamkeag.mat');
    Livneh_data = Livneh_data_mattawamkeag(7306:35429,:);
    load('MR_historical_data.mat'); 
    historical_data = MR_historical_data; 
    load('MR_future_data.mat'); 
    future_data = MR_future_data;     
    acronym = 'MR'; 
    load('MR_all_snow.mat'); 
    load('MR_avg_snow.mat');
elseif(contains(river, 'Diamond'))
    min_year = 1942;
    load('DeadDiamond_River/Livneh_data_diamond.mat');
    Livneh_data = Livneh_data_diamond(9863:35429,:);
    load('DR_historical_data.mat'); 
    historical_data = DR_historical_data; 
    load('DR_future_data.mat'); 
    future_data = DR_future_data; 
    acronym = 'DR';
    load('DR_all_snow.mat'); 
    load('DR_avg_snow.mat');
elseif(contains(river, 'Ware'))
    directory = 'Ware_River/raw_snow_data'; 
    range = 4:13; 
    badfile = 0; 
    min_year = 2005; 
    
    load('Ware_River/Livneh_data_ware.mat');
    Livneh_data = Livneh_data_ware(32874:35429,:);    
    
    acronym = 'WR2';

    
else
    print('Error! Must pick a river') 

end

historical_data.Properties.VariableNames = {'Date', 'av_temp', 'precip'};
future_data.Properties.VariableNames = {'Date', 'av_temp', 'precip'};

av_temp = nanmean(Livneh_data(:, 5:6), 2);
Livneh_data = [Livneh_data(:,1:4), av_temp, Livneh_data(:,5:6)];
Livneh_data = array2table(Livneh_data); 
Livneh_data.Properties.VariableNames = {'Day', 'Month', 'Year', 'precip', ...
  'av_temp', 'max_temp', 'min_temp'};% %'wind'};   



%% Trimming to the correct size 

% if(contains(river,'White')) 
%     out_of_range = find(year(all_snowdepths.new_date)>2011 | ...
%         year(all_snowdepths.new_date)<min_year | year(all_snowdepths.new_date) ...
%         == 1927 | year(all_snowdepths.new_date) ==1928); 
%     
% else
     out_of_range = find(year(all_snowdepths.new_date)>2011 | ...
        year(all_snowdepths.new_date)<min_year); 
%end

average_snowdepth_trimmed = average_snowdepth; 
average_snowdepth_trimmed(out_of_range) = []; 
dates_obs = all_snowdepths.new_date; 
dates_obs(out_of_range) = [];


%% Section 4: Getting accumulated winter precip 
livneh_precip = Livneh_data.precip;
livneh_temp = Livneh_data.av_temp;

precip_accumulation = zeros(length(livneh_precip),1);
snow_accumulation = zeros(length(livneh_precip),1);
for i = min_year+1:2011 
    
%     if(contains(river,'White'))
%         if(i == 1927 | i == 1928)
%             continue
%         end   
%     end
    
    dec_index = find(Livneh_data.Day ==1 & Livneh_data.Month==12 ...
        & Livneh_data.Year == i-1);
    may_index = find(Livneh_data.Day ==1 & Livneh_data.Month==5 ...
        & Livneh_data.Year == i);
    
%     if(contains(river,'White')) 
%         if(i == 1929) 
%         dec_index = find(Livneh_data.Day ==1 & Livneh_data.Month==1 ...
%             & Livneh_data.Year == i);  
%         end
%     end
    
      
    for j = dec_index:may_index 
        precip_accumulation(j) = sum(livneh_precip(dec_index:j)); 
        if(livneh_temp(j)<0)  
          snow_accumulation(j) = snow_accumulation(j-1)+ livneh_precip(j); 
        else
          snow_accumulation(j) = snow_accumulation(j-1);
        end  
    end
end

%% Section 5: calculate growing degree days 

gdd = zeros(length(livneh_temp), 1); 
for i = min_year:2011  
    
%     if(contains(river,'White'))
%         if(i ==1927 | i ==1928)
%             continue;
%         elseif(i == 1929 | i == min_year) 
%              first_index = find((Livneh_data.Year) == i ...
%                 & (Livneh_data.Month) ==1 & (Livneh_data.Day) == 1);
%              last_index = find(Livneh_data.Year == i & Livneh_data.Month ==5 ...
%                 & Livneh_data.Day == 31);
%         else
%             first_index = find(Livneh_data.Year == i-1 & Livneh_data.Month ==10 ...
%                 & Livneh_data.Day == 1);
%             last_index = find(Livneh_data.Year == i & Livneh_data.Month ==5 ...
%                 & Livneh_data.Day == 31);      
%         end     
%     else 
        if(i == min_year) 
            first_index = find(Livneh_data.Year == i ...
                & Livneh_data.Month ==1 & Livneh_data.Day == 1);
             last_index = find(Livneh_data.Year == i & Livneh_data.Month ==5 ...
                & Livneh_data.Day == 31);   
        else
            first_index = find(Livneh_data.Year == i-1 & Livneh_data.Month ==10 ...
                & Livneh_data.Day == 1);
            last_index = find(Livneh_data.Year == i & Livneh_data.Month ==5 ...
                & Livneh_data.Day == 31);     
        end
     %end
    
    gdd(first_index) = 0;
    for j = first_index+1:last_index 
        if(livneh_temp(j)>0)
            this_gdd = livneh_temp(j);
            gdd(j) = gdd(j-1) + this_gdd;

        else 
            this_gdd = 0;
            gdd(j) = gdd(j-1) + this_gdd;

        end
    end
end

%% Section 6: getting same inputs, but for WRF

precip_hist = historical_data.precip;
temp_hist = historical_data.av_temp;
dates_hist = historical_data.Date;

precip_accumulation_hist = zeros(length(precip_hist),1);
snow_accumulation_hist = zeros(length(precip_hist),1);
for i = 1976:2005
    
    if(i ==1976)
        dec_index = find(day(dates_hist) ==1 & month(dates_hist) ==1 ...
        & year(dates_hist) == i);
    else
        dec_index = find(day(dates_hist) ==1 & month(dates_hist) ==12 ...
        & year(dates_hist) == i-1);
    end
    

    may_index = find(day(dates_hist) ==1 & month(dates_hist)==5 ...
        & year(dates_hist) == i);
    
    for j = dec_index:may_index 
        precip_accumulation_hist(j) = sum(precip_hist(dec_index:j)); 
        if(temp_hist(j)<0)  
            if(i ==1976 & j == dec_index)
                snow_accumulation_hist(j) = precip_hist(j);    
            else
                snow_accumulation_hist(j) = snow_accumulation_hist(j-1)+ precip_hist(j); 
            end
            
        else
            if(i ==1976 & j == dec_index)
                snow_accumulation_hist(j) = 0;    
            else
                snow_accumulation_hist(j) = snow_accumulation_hist(j-1);
            end
        end  
    end
end

%growing degree days 

gdd_hist = zeros(length(temp_hist), 1); 
for i = 1976:2005 

        if(i == 1976) 
            first_index = find(year(dates_hist) == i ...
                & month(dates_hist) ==1 & day(dates_hist) == 1);
             last_index = find(year(dates_hist) == i & month(dates_hist) ==5 ...
                & day(dates_hist) == 31);   
        else
            first_index = find(year(dates_hist) == i-1 & month(dates_hist) ==10 ...
                & day(dates_hist) == 1);
            last_index = find(year(dates_hist) == i & month(dates_hist) ==5 ...
                & day(dates_hist) == 31);     
        end
%     end
    
    gdd_hist(first_index) = 0;
    for j = first_index+1:last_index 
        if(temp_hist(j)>0)
            this_gdd_hist = temp_hist(j);
            gdd_hist(j) = gdd_hist(j-1) + this_gdd_hist;

        else 
            this_gdd_hist = 0;
            gdd_hist(j) = gdd_hist(j-1) + this_gdd_hist;

        end
    end
end

%% Now for the future

precip_fut = future_data.precip;
temp_fut = future_data.av_temp;
dates_fut = future_data.Date;

precip_accumulation_fut = zeros(length(precip_fut),1);
snow_accumulation_fut = zeros(length(precip_fut),1);
for i = 2006:2099
    
    if(i ==2006)
        dec_index = find(day(dates_fut) ==1 & month(dates_fut) ==1 ...
        & year(dates_fut) == i);
    else
        dec_index = find(day(dates_fut) ==1 & month(dates_fut) ==12 ...
        & year(dates_fut) == i-1);
    end
    

    may_index = find(day(dates_fut) ==1 & month(dates_fut)==5 ...
        & year(dates_fut) == i);
    
    for j = dec_index:may_index 
        precip_accumulation_fut(j) = sum(precip_fut(dec_index:j)); 
        if(temp_fut(j)<0)  
            if(i ==2006 & j == dec_index)
                snow_accumulation_fut(j) = precip_fut(j);    
            else
                snow_accumulation_fut(j) = snow_accumulation_fut(j-1)+ precip_fut(j); 
            end
            
        else
            if(i ==2006 & j == dec_index)
                snow_accumulation_fut(j) = 0;    
            else
                snow_accumulation_fut(j) = snow_accumulation_fut(j-1);
            end
        end  
    end
end

%growing degree days 

gdd_fut = zeros(length(temp_fut), 1); 
for i = 2006:2099

        if(i == 1976) 
            first_index = find(year(dates_fut) == i ...
                & month(dates_fut) ==1 & day(dates_fut) == 1);
             last_index = find(year(dates_fut) == i & month(dates_fut) ==5 ...
                & day(dates_fut) == 31);   
        else
            first_index = find(year(dates_fut) == i-1 & month(dates_fut) ==10 ...
                & day(dates_fut) == 1);
            last_index = find(year(dates_fut) == i & month(dates_fut) ==5 ...
                & day(dates_fut) == 31);     
        end
%     end
    
    gdd_fut(first_index) = 0;
    for j = first_index+1:last_index 
        if(temp_fut(j)>0)
            this_gdd_fut = temp_fut(j);
            gdd_fut(j) = gdd_fut(j-1) + this_gdd_fut;

        else 
            this_gdd_fut = 0;
            gdd_fut(j) = gdd_fut(j-1) + this_gdd_fut;

        end
    end
end

%% Section 6: regression! - let's do a loop!

avg_snowdepth_prvs = [0 average_snowdepth_trimmed']';
avg_snowdepth_prvs = avg_snowdepth_prvs(1:end-1); 
%%
for i = 1%:50%:500
    
i
    winter_indices2 = find((month(dates_obs) == 1 |...
             month(dates_obs) == 2 |...
             month(dates_obs) == 3 | ...
             month(dates_obs) == 4 | ...
             month(dates_obs) == 5 | ...
             month(dates_obs) == 11 |...
             month(dates_obs) == 12));

    % Get only the winter values, only the places where average snowdepth is not NaN     
    X = table(livneh_temp, livneh_precip,snow_accumulation,gdd, ...
        avg_snowdepth_prvs, average_snowdepth_trimmed);
    X_winter = X(winter_indices2, :);
    
    index1 = ~isnan(average_snowdepth_trimmed(winter_indices2));
    index2 = ~isnan(avg_snowdepth_prvs(winter_indices2));
    indices = false(size(winter_indices2,1),1); 
    indices(index1==1 & index2 ==1) = true; 
    X_nonans = X_winter(indices, :);

    %Get a training dataset, take 10,000 values to train on, roughly 1/2 of the
    %dataset 
    [training_set, training_indices] = datasample(X_nonans,800);


    tree = fitrtree(training_set,'average_snowdepth_trimmed','MinLeafSize', 25); 

    % Get testing set (aka those that were not used to train) to see how well
    % it does 
    test_indices = true(size(X_nonans(:,1)));
    test_indices(training_indices) = false; 
    testing_set = X_nonans(test_indices, :);

    %Make predictions for the testing set, see how well it correlates  
    predictions = predict(tree, testing_set(:, 1:5)); 
    [correlation_obs_predictions, corr_obs_pre_pval] = corrcoef(predictions, table2array(X_nonans(test_indices, 5)));
    correlations(i) = correlation_obs_predictions(2);
    %
    % Use the tree to make predictions for the whole winter dataset
    all_predictions = predict(tree, X_winter(:, 1:5));


    % Timetable of just the winter values 
    TT_obs = table2timetable(table(all_predictions),'RowTimes',dates_obs(winter_indices2));

    %Synchronize it to a full time length, non-winter values will be NaN
    regressed_snowdepthTT = synchronize(timetable(dates_obs), TT_obs);
    %Set regressed snowdepth to be 0 in non-winter (when it is NaN, no
    %predictions were made for these times) 
    regressed_snowdepth(:,i) = regressed_snowdepthTT.all_predictions;
    regressed_snowdepth(isnan(regressed_snowdepth(:,i)),i) = 0; 

    observations = table2array(X_nonans(test_indices, 5));

    rmse = sqrt(mean((predictions-observations).^2));

    % WRF HISTORICAL
    WRF_predictors_hist = table(temp_hist, precip_hist,snow_accumulation_hist,gdd_hist);
    WRF_predictors_hist.Properties.VariableNames = ...
        X_nonans.Properties.VariableNames(1:4); 

    WRF_winter_indices_hist = find(month(dates_hist) ==1 |month(dates_hist) ==2 | ...
        month(dates_hist) ==3 | month(dates_hist) ==4 | month(dates_hist) ==5 ...
        | month(dates_hist) ==11 | month(dates_hist) ==12); 
    WRF_preds_winter_hist = WRF_predictors_hist(WRF_winter_indices_hist,:); 
    WRF_dates_winter_hist = dates_hist(WRF_winter_indices_hist); 


    WRF_snow_predictions_hist = zeros(size(WRF_preds_winter_hist,1),1); 
    
    for j = 1:size(WRF_preds_winter_hist,1)
        
        if(j ==1)
            precip = WRF_preds_winter_hist(j,2);
            precip.Properties.VariableNames = {'avg_snowdepth_prvs'}; 
            wrf_table_hist = [WRF_preds_winter_hist(j,:), precip];
            WRF_snow_predictions_hist(j) = predict(tree, wrf_table_hist);
        else
            WRF_table_hist = [WRF_preds_winter_hist(j,:), ...
                table(WRF_snow_predictions_hist(j-1))];
            WRF_table_hist.Properties.VariableNames(5) = {'avg_snowdepth_prvs'};
            WRF_snow_predictions_hist(j) = predict(tree, WRF_table_hist);
            

        end
         
    end
    
    TT_hist = table2timetable(table(WRF_snow_predictions_hist),'RowTimes',WRF_dates_winter_hist);

    %Synchronize it to a full time length, non-winter values will be NaN
    regressed_snowdepth_WRF_TT_hist = synchronize(timetable(dates_hist), TT_hist);
    regressed_snowdepth_WRF_hist(:,i) = (regressed_snowdepth_WRF_TT_hist.WRF_snow_predictions_hist);
    %regressed_snowdepth_WRF_hist.Properties.VariableNames = {'Dates', 'snowdepth'};

    % WRF FUTURE
    WRF_predictors_fut = table(temp_fut, precip_fut,snow_accumulation_fut,gdd_fut);
    WRF_predictors_fut.Properties.VariableNames = ...
        X_nonans.Properties.VariableNames(1:4); 

    WRF_winter_indices_fut = find(month(dates_fut) ==1 |month(dates_fut) ==2 | ...
        month(dates_fut) ==3 | month(dates_fut) ==4 | month(dates_fut) ==5 ...
        | month(dates_fut) ==11 | month(dates_fut) ==12); 
    WRF_preds_winter_fut = WRF_predictors_fut(WRF_winter_indices_fut,:); 
    WRF_dates_winter_fut = dates_fut(WRF_winter_indices_fut); 
    

    WRF_snow_predictions_fut = zeros(size(WRF_preds_winter_fut,1),1); 

    for j = 1:size(WRF_preds_winter_fut,1)
        size_vec(j,:) = size(WRF_snow_predictions_fut); 
        if(j ==1)
            precip = WRF_preds_winter_fut(j,2);
            precip.Properties.VariableNames = {'avg_snowdepth_prvs'}; 
            wrf_table_fut = [WRF_preds_winter_fut(j,:), precip];
            WRF_snow_predictions_fut(j) = predict(tree, wrf_table_fut);
        else
            WRF_table_fut = [WRF_preds_winter_fut(j,:), ...
                table(WRF_snow_predictions_fut(j-1))];
            WRF_table_fut.Properties.VariableNames(5) = {'avg_snowdepth_prvs'};
            WRF_snow_predictions_fut(j) = predict(tree, WRF_table_fut);
            
        end     
    end
    
    TT_fut = table2timetable(table(WRF_snow_predictions_fut),'RowTimes',WRF_dates_winter_fut);


    %Synchronize it to a full time length, non-winter values will be NaN
    regressed_snowdepth_WRF_TT_fut = synchronize(timetable(dates_fut), TT_fut);
    regressed_snowdepth_WRF_fut(:,i) = (regressed_snowdepth_WRF_TT_fut.WRF_snow_predictions_fut);
    %regressed_snowdepth_WRF_fut.Properties.VariableNames = {'Dates', 'snowdepth'};


end

%% Getting avg 


regressed_snow_avg = nanmean(regressed_snowdepth, 2); 
regressed_snow_std = std(regressed_snowdepth(:,2:end),0, 2); 
regressed_snow_se = regressed_snow_std./sqrt(500); 
mean_std = nanmean(regressed_snow_std);
mean_se = nanmean(regressed_snow_se); 
save(strcat(acronym, '_snowdepth_avg_predictions_test.mat'), 'regressed_snow_avg');

regressed_snow_WRFhist_avg = nanmean(regressed_snowdepth_WRF_hist,2); 
regressed_snow_WRFhist_std = std(regressed_snowdepth_WRF_hist,0,2);
regressed_snow_se = regressed_snow_WRFhist_std./sqrt(500); 
mean_std_hist = nanmean(regressed_snow_WRFhist_std);
mean_se_hist = nanmean(regressed_snow_se); 


regressed_snow_WRFfut_avg = nanmean(regressed_snowdepth_WRF_fut,2); 
regressed_snow_WRFfut_std = std(regressed_snowdepth_WRF_fut,0,2);
regressed_snow_se = regressed_snow_WRFfut_std./sqrt(500); 
mean_std_fut = nanmean(regressed_snow_WRFfut_std);
mean_se_fut = nanmean(regressed_snow_se); 

save(strcat(acronym, '_snowdepth_avg_WRFhist_test.mat'), 'regressed_snow_WRFhist_avg'); 
save(strcat(acronym, '_snowdepth_avg_WRFfut_test.mat'), 'regressed_snow_WRFfut_avg'); 


avg_corr = nanmean(correlations)

mean_se_table = [nanmean(regressed_snow_avg), mean_se; ...
    nanmean(regressed_snow_WRFhist_avg), mean_se_hist; ...
    nanmean(regressed_snow_WRFfut_avg), mean_se_fut];

% %%
% for i = 1916:2011
%     i
% %     if(i == 1927 | i == 1928)
% %         continue;
% %     end
%     this_year_indices = find(year(livneh_dates) == i); 
%     these_snows = average_snowdepth_trimmed(this_year_indices); 
%     max_snows(i-1915) = max(these_snows);     
% end 
% 
% mean_max_snowdepth= nanmean(max_snows);
% 
% 

% 
% %%
% if(contains(river,'White')) 
%     regressed_snowdepth_WR = regressed_snowdepth; 
% elseif(contains(river,'Shenandoah')) 
%     regressed_snowdepth_SR = regressed_snowdepth;
% elseif(contains(river,'Mattawamkeag'))
%     regressed_snowdepth_MR = regressed_snowdepth;
% elseif(contains(river, 'Diamond'))
%     regressed_snowdepth_DR = regressed_snowdepth;
% elseif(contains(river, 'Ware'))
%     regressed_snowdepth_WR2 = regressed_snowdepth;
% end
% 
% 
% 

%% Getting Livneh for comparison
% 
% 
% swe_files = dir('swe_Livneh'); 
% swe_filenames = {swe_files.name};
% 
% swe_info = ncinfo('swe.2008.nc');
% 
% lat_vec = ncread(string(swe_filenames(3)),'lat');
% lon_vec = ncread(string(swe_filenames(3)), 'lon');
% lon_vec = (360-lon_vec)*-1;
% 
% lat = zeros(922, 444);
% for m= 1:922
%     lat(m, :) = lat_vec;
% end
% 
% lon = zeros(922, 444);
% for m= 1:444
%     lon(:,m) = lon_vec;
% end
% 
% if(contains(river, 'White'))
%     S = shaperead('White_River/layers/globalwatershed.shp');
% elseif(contains(river, 'Shenandoah'))
%     S = shaperead('Shenandoah_River/watershed_shape_shenandoah/globalwatershed.shp');
% elseif(contains(river, 'Mattawamkeag'))
%     S = shaperead('Mattawamkeag_River/watershed_shape_mattawamkeag/globalwatershed.shp');
% end
% 
% in = inpolygon(lon, lat, S.X, S.Y);
% 
% clear Livneh_swe_data
% Livneh_swe_data = table;
% for i = 1:97
%     i
%     this_filename = string(swe_filenames(i+2));
%     time = ncread(this_filename, 'time');
%     swe = ncread(this_filename, 'swe'); 
%     
%     clear storage
%     for j = 1:length(time)
%         
%        swe_temp = swe(:,:,j);
%        this_swe = nanmean(nanmean(swe_temp(in),1),2)*0.03937; %Take average SWE over watershed, convert to inch
%        temp_date = convertCharsToStrings(datestr(time(j)+datenum('1915-01-01','yyyy-mm-dd'),1));
%        this_date = datetime(temp_date, 'InputFormat', 'dd-MMM-yyyy');
%        
%        storage(j, :) = table(this_swe, this_date); 
% 
%     end
%     
%     Livneh_swe_data = cat(1, Livneh_swe_data, storage);
% end
% 
% %Get datetime values from temp and precip
% Livneh_datetime = datetime(Livneh_data(:,3), Livneh_data(:,2), Livneh_data(:,1));
% 
% %% Trimming to correct size 
% 
% if(contains(river,'White')) 
%     out_of_range = find(year(Livneh_datetime)>2011 | ...
%         year(Livneh_datetime)<min_year | year(Livneh_datetime) ...
%         == 1927 | year(Livneh_datetime) ==1928); 
%     
% else
%      out_of_range = find(year(Livneh_datetime)>2011 | ...
%         year(Livneh_datetime)<min_year); 
% end
% 
% livneh_swe_trimmed = Livneh_swe_data.this_swe; 
% livneh_swe_trimmed(out_of_range) = []; 
% livneh_dates = Livneh_datetime;
% livneh_dates(out_of_range) = [];
% 
% 
% %% Getting correlations 
% 
% snowdepth_nan_indices = ~isnan(average_snowdepth);
% 
% winter_indices = find(month(livneh_dates) == 1 |...
%         month(livneh_dates) == 2 |...
%         month(livneh_dates) == 3 | ...
%         month(livneh_dates) == 4 | ...
%         month(livneh_dates) == 5 | ...
%         month(livneh_dates) == 10 | ...
%         month(livneh_dates) == 11 | ...
%         month(livneh_dates) == 12);
%     
% winter_snowdepths_obs = average_snowdepth_trimmed(winter_indices); 
% winter_swe_livneh = (livneh_swe_trimmed(winter_indices));
% 
% nan_indices = ~isnan(winter_snowdepths_obs); 
% 
% [correlation_livneh_winter livneh_pval_winter] = corrcoef(winter_snowdepths_obs(nan_indices), winter_swe_livneh(nan_indices));
% 
% nan_indices = ~isnan(average_snowdepth_trimmed);
% [correlation_livneh, livneh_pval] = corrcoef(average_snowdepth_trimmed(nan_indices), livneh_swe_trimmed(nan_indices));
% 
% 
% %% Plotting the three against each other 
% 
% 
% 
% comparison_fig = figure(1), clf 
% plot(livneh_dates, average_snowdepth_trimmed, 'LineWidth',2)
% hold on 
% plot(livneh_dates, livneh_swe_trimmed, 'LineWidth',2)
% %plot(livneh_dates, regressed_snowdepth, 'LineWidth',1, 'Color','g')
% title('Observed snow depth vs. Livneh SWE', 'FontSize', 22) 
% ylabel('Snow depth and SWE [inch]', 'FontSize', 14)
% xlim([datetime(1981, 9, 1) datetime(1984, 9,1)]);
% legend('Snow depth observations', 'Livneh data', 'Predicted depths', 'FontSize', 15)
% %saveas(comparison_fig, 'Figures/comparison_fig.png')


%%
figure(1), clf 
time_vec = datetime(1915,01,01):datetime(2011,12,31); 
plot(time_vec, regressed_snow_avg*25.4); 
title('White River snow depth over time', 'FontSize', 20); 
ylabel('Depth [mm]'); 

%%
figure(2), clf 
time_vec = datetime(1976,01,01):datetime(2005,12,31); 
time_vec(month(time_vec)==2 & day(time_vec)==29) = []; 
plot(time_vec, regressed_snow_WRFhist_avg*25.4);
title('White River snow depth over time- WRF', 'FontSize', 20); 
ylabel('Depth [mm]'); 


%%
figure(1), clf 

histogram(regressed_snow_WRFhist_avg); 

