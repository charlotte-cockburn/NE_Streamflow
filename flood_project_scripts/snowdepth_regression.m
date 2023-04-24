%% GETTING REGRESSED SNOWDEPTH INPUT %%

%%  Section 1: Settings based on which river we pick 
clear all, close all
river = 'White'; %SET THIS VALUE!!! 

if(contains(river,'White'))
    directory = 'White_River/raw_snow_data'; 
    range = 4:24;
    badfile = 19; 
    min_year = 1915;  
    load('White_River/Livneh_data_white.mat');
    Livneh_data = Livneh_data_white;
    acronym = 'WR';
    
elseif(contains(river,'Shenandoah')) 
    directory = 'Shenandoah_River/raw_snow_data'; 
    range = 4:17;
    badfile = [8 16]; 
    min_year = 1926;
    load('Shenandoah_River/Livneh_data_shenandoah.mat');
    Livneh_data = Livneh_data_shenandoah(4019:35429,:); 
    acronym = 'SR';
elseif(contains(river,'Mattawamkeag'))
    directory = 'Mattawamkeag_River/raw_snow_data'; 
    range = 4:16;
    badfile = [6]; 
    min_year = 1935;
    load('Mattawamkeag_River/Livneh_data_mattawamkeag.mat');
    Livneh_data = Livneh_data_mattawamkeag(7306:35429,:);
    acronym = 'MR';
    
elseif(contains(river, 'Diamond'))
    directory = 'DeadDiamond_River/raw_snow_data'; 
    range = 4:10;
    badfile = 0; 
    min_year = 1942;
  
    load('DeadDiamond_River/Livneh_data_diamond.mat');
    Livneh_data = Livneh_data_diamond(9863:35429,:);
    acronym = 'DR';
   
elseif(contains(river, 'Ware'))
    directory = 'Ware_River/raw_snow_data'; 
    range = 4:13; 
    badfile = 0; 
    min_year = 2005; 
    
    load('Ware_River/Livneh_data_ware.mat');
    Livneh_data = Livneh_data_ware(32874:35429,:);    
        
else
    print('Error! Must pick a river') 

end

av_temp = nanmean(Livneh_data(:, 5:6), 2);
Livneh_data = [Livneh_data(:,1:4), av_temp, Livneh_data(:,5:6)];
Livneh_data = array2table(Livneh_data); 
Livneh_data.Properties.VariableNames = {'Day', 'Month', 'Year', 'precip', ...
  'av_temp', 'max_temp', 'min_temp'};% %'wind'};   

%% Section 2: getting the snow depth data in a form that we can use  

%Load snow depth data, get it in a nice vector format 
file_directory = dir(directory);
filenames = {file_directory.name};
for i = range
    i
    %Load the raw data
    this_filename = filenames(i);
    this_data = readtable(string(this_filename));
    raw_date = char(this_data{:,1});
    raw_temp = string(this_data{:,2}); 
    raw_precip = string(this_data{:,3});
    raw_snowfall = string(this_data{:,4});
    raw_snowdepth = char(this_data{:,5});
    %Set up empty vectors 
    new_temp = zeros(length(raw_temp)-1, 1);
    new_precip = zeros(length(raw_precip)-1, 1);
    new_snowfall = zeros(length(raw_snowfall)-1, 1);
    new_snowdepth = zeros(length(raw_snowdepth)-1, 1); 
    clear new_date; 
    clear raw_new_snowdepth;
    for j = 1:(length(raw_date)-1)
        %Convert the date to datetime format 
        if(ismember(i, badfile))
            %One station had a different format, so here we are 
            new_date(j) = datetime(raw_date(j, :), 'InputFormat', 'yyyy-MM-dd');
        else
            new_date(j) = datetime(raw_date(j, :), 'InputFormat', 'dd-MMM-yyyy');
        end       
        %Data have 'M' rather than NaN value, so convert to NaN
        %Temperature 
        if(raw_temp(j) == 'M')
            new_temp(j) = NaN;
        else 
            new_temp(j) = str2double(raw_temp(j));
        end    
        %Precip
        if(raw_precip(j) == 'M')
            new_precip(j) = NaN;
        else 
            new_precip(j) = str2double(raw_precip(j)); 
        end 
        %Snowfall 
        if(raw_snowfall(j) == 'M') 
            new_snowfall(j) = NaN; 
        else 
            new_snowfall(j) = str2double(raw_snowfall(j)); 
        end   
        
        %Snow depth -- trickier due to those annoying backslashes in WR 
        if(contains(river,'White')) 
            index = find(raw_snowdepth(j,:) == '\');
            raw_new_snowdepth(j,:) = string(raw_snowdepth(j,1:index-1));
        else
            raw_new_snowdepth(j,:) = raw_snowdepth(j,:);
        end
        
        
        if(raw_new_snowdepth(j) == 'M' | raw_new_snowdepth == 'T')
            new_snowdepth(j) = NaN; 
        else 
            new_snowdepth(j) = str2double(raw_new_snowdepth(j)); 
        end
        
    end
    new_date = new_date';
    max_date(i) = max(new_date);
    recombined = table(new_date, new_temp, new_precip, new_snowfall, new_snowdepth);
    this_timetable = table2timetable(recombined);
    %for the first entry, set up new table
    if(i==range(1)) 
        all_temps = this_timetable(:,1);
        all_precips = this_timetable(:,2);
        all_snowfalls = this_timetable(:,3);
        all_snowdepths = this_timetable(:,4);
    else 
        %For the rest, use the synchronize function to line them up
        %according to datetime
        all_temps = synchronize(all_temps, this_timetable(:,1));  
        all_precips = synchronize(all_precips, this_timetable(:,2)); 
        all_snowfalls = synchronize(all_snowfalls, this_timetable(:,3));
        all_snowdepths = synchronize(all_snowdepths, this_timetable(:,4));
    end       
end
max_date = max_date';
%%
all_temps = retime(all_temps, 'daily', 'fillwithmissing'); 
all_precips = retime(all_precips, 'daily', 'fillwithmissing'); 
all_snowfalls = retime(all_snowfalls, 'daily', 'fillwithmissing'); 
all_snowdepths = synchronize(all_snowdepths, 'daily', 'fillwithmissing'); 

if(contains(river, 'Ware'))
    all_precips(2556, :) = all_precips(2555,:); 
    all_temps(2556, :) = all_temps(2555,:); 
    all_snowfalls(2556, :) = all_snowfalls(2555,:); 
    all_snowdepths(2556, :) = all_snowdepths(2555,:); 
    
    all_precips.new_date(2556) = datetime(2011,12,31); 
    all_temps.new_date(2556) = datetime(2011,12,31); 
    all_snowdepths.new_date(2556) = datetime(2011,12,31); 
    all_snowdepths.new_date(2556) = datetime(2011,12,31);
end


if(contains(river,'White'))
    var_names = {'Barre', 'Bethel', 'Bethel4N', 'WR2ookfield', 'Chelsea', 'Chelsea2', 'Corinth',...
    'Groton', 'Hanover', 'Lebanon', 'NewHaven', 'Northfield', 'Pomfret', 'Rochester', ...
    'SouthLincoln', 'Topsham', 'UnionVillage', 'Waitsfield', 'WestHartford',...
    'Woodstock', 'WRJunction'};
elseif(contains(river,'Mattawamkeag'))
    var_names = {'Danforth', 'Enfield', 'Haynesville', 'Houlton', 'Houlton2', ...
    'Howe', 'Knowles', 'Lincoln', 'Millinocket_airport', 'Millinocket', ...
    'Millinocket2', 'Prentiss', 'Springfield'};
    
elseif(contains(river,'Shenandoah'))
    var_names = {'Bartow', 'WR2andyWine', 'WR2ushyRun', 'Cootes', 'Dale', 'Dale2', ...
    'FortSeybert', 'FortSyevert2', 'Franklin', 'LostRiver', 'SugarGrove', ...
    'Timberville', 'UpperTract', 'UpperTract2'}; 
elseif(contains(river, 'Diamond'))
    var_names = {'Canaan', 'ConnecticutLake', 'DiamondPond', 'Dixville', 'Errol', ...
        'Lincoln', 'Pittsburg'}; 
    
elseif(contains(river, 'Ware'))
    var_names = {'Ashburnham', 'AshburnhamNorth', 'Barre', 'BirchHill', ...
        'EastMilford', 'Greenville', 'Hardwicke', 'Milford', 'TullyLake', ...
        'Worcester'};
end

all_temps.Properties.VariableNames = var_names;
all_precips.Properties.VariableNames = var_names;
all_snowfalls.Properties.VariableNames = var_names;
all_snowdepths.Properties.VariableNames = var_names;


%% Section 3: Get average snowdepth -- must have at least 2 stations 
raw_size = height(all_snowdepths(:,2)); 
average_snowdepth = zeros(raw_size, 1); 
num_stations = zeros(raw_size, 1);
for i = 1:raw_size
    this_row = table2array(all_snowdepths(i, :)); 
      if(sum(~isnan(this_row))<1)
          average_snowdepth(i) = NaN;
          num_stations(i) =0;
      else
        average_snowdepth(i) = nanmean(this_row);
        num_stations(i) = sum(~isnan(this_row));
      end    
end
sum(isnan(average_snowdepth)) + sum(average_snowdepth ==0)

%%

save(strcat(acronym, '_avg_snow.mat'), 'average_snowdepth'); 
save(strcat(acronym, '_all_snow.mat'), 'all_snowdepths');