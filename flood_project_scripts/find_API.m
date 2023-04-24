%% Getting API coefficient 

%Random forest! 
clear all, close all

%%%%% CHANGE THIS %%%%%

river = 'Diamond'; 

%%%%% CHANGE NOTHING BELOW THIS LINE %%%%%

if(contains(river,'White'))
    load('RF_data_WR.mat'); 
    RF_table = RF_data_WR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
    
    load('Livneh_data_white.mat'); 
    Livneh_data = Livneh_data_white;
    load('White_discharge.mat');
    discharge_data = White_discharge; 
    
    
    acronym = 'WR';
    title_name = 'White River';


elseif(contains(river,'Shenandoah'))
    load('RF_data_SR.mat'); 
    RF_table = RF_data_SR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
    
    load('Livneh_data_shenandoah.mat'); 
    Livneh_data = Livneh_data_shenandoah;
    load('Shenandoah_discharge.mat');
    discharge_data = shenandoah_discharge; 

    acronym = 'SR';
    title_name = 'Shenandoah River';
    
elseif(contains(river,'Mattawamkeag'))
    load('RF_data_MR.mat'); 
    RF_table = RF_data_MR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
    
    load('Livneh_data_mattawamkeag.mat'); 
    Livneh_data = Livneh_data_mattawamkeag;
    load('Mattawamkeag_discharge.mat');
    discharge_data = mattawamkeag_discharge; 

    acronym = 'MR';
    title_name = 'Mattawamkeag River';
    
elseif(contains(river,'Diamond'))
    load('RF_data_DR.mat'); 
    RF_table = RF_data_DR;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
   
    load('Livneh_data_diamond.mat'); 
    Livneh_data = Livneh_data_diamond;
    load('Diamond_discharge.mat');
    discharge_data = diamond_discharge; 
    
    acronym = 'DR';
    title_name = 'Dead Diamond River';
    
 elseif(contains(river,'Ware'))
    load('RF_data_WR2.mat'); 
    RF_table = RF_data_table;
    RF_table(find(isnan(RF_table.discharge)), :) = []; 
   
    load('Livneh_data_ware.mat'); 
    Livneh_data = Livneh_data_ware;
    load('ware_discharge.mat');
    discharge_data = ware_discharge; 
    
    acronym = 'WR2';
    title_name = 'Ware River';
end

%RF_table = RF_table(:,[1:4 9:14 30:35]);

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
%%
if(contains(river, 'White'))
    combined_TT(year(combined_TT.dates_livneh) == 1927 |...
        year(combined_TT.dates_livneh) == 1928 | ...
        year(combined_TT.dates_livneh) == 1915, :) = [];
end



%%
winter_indices = find(RF_table.Month == 1| ...
    RF_table.Month ==2|RF_table.Month ==3 | RF_table.Month ==4 | ...
    RF_table.Month ==5 | RF_table.Month == 11 | RF_table.Month ==12);

summer_indices = find(RF_table.Month == 6 | RF_table.Month ==7 | ...
    RF_table.Month ==8 | RF_table.Month == 9 | RF_table.Month ==10);

RF_dataset = RF_table(:, 4:width(RF_table(1,:))-3); 
%max_flag = logical(RF_table.topflow_indicator); 

RF_dataset_winter = RF_table(winter_indices,4:width(RF_table(1,:)));
%max_flag_win = logical(RF_table.topflow_indicator_winter(winter_indices)); 

RF_dataset_summer = RF_table(summer_indices,4:width(RF_table(1,:)));
%max_flag_sum = logical(RF_table.topflow_indicator_summer(summer_indices)); 

RF_dataset_winter.api = [];
RF_dataset_winter.spi = [];
RF_dataset_winter.spi_prvs = []; 

RF_dataset_summer.melt_5day = [];
RF_dataset_summer.melt_30day = [];
RF_dataset_summer(:,3) = [];

%%

precip = combined_TT.precip; 
precip_sum = precip(summer_indices); 



discharge_sum = RF_dataset_summer.discharge;
count = 1;
for i = 2:length(precip_sum) 
    if(precip_sum(i) == 0& precip_sum(i-1) ==0) 
        discharge_d(count) = discharge_sum(i); 
        discharge_d_1(count) = discharge_sum(i-1);
        count = count+1; 
    end
    
end



indices = find(discharge_d_1>=1800);
discharge_d_1(indices) = [];
discharge_d(indices) = [];


P = polyfit(discharge_d_1, discharge_d, 1);
slope = P(1)
intercept = P(2)
f = polyval(P,discharge_d_1);  % P(1) is the slope and P(2) is the intercept
figure(1), clf 
plot(discharge_d_1, discharge_d, '.') 
hold on;
plot(discharge_d_1,f,'r-.')
ylim([0 5000]);
xlim([0 5000]);








