%% Comparing historical WRF to Livneh obs 
clear all, close all
%Load data 

load('historical_data_d1_MPI.mat'); 
rain_h = rain; 
temperature_h = temperature; 
time_h = date;
load('future_data_d1_MPI.mat');
load('lat_d1_wrf.mat'); 
load('lon_d1_wrf.mat'); 
load('landmask_d1_wrf.mat'); 

rain_f = rain(:,:,end-10949:end); 

white_watershed = shaperead('watershed_WR/globalwatershed.shp'); 
shen_watershed = shaperead('watershed_SR/globalwatershed.shp');
matt_watershed = shaperead('watershed_MR/globalwatershed.shp');
diam_watershed = shaperead('watershed_DR/globalwatershed.shp');

%% Get 3 day precip 

rain_3day_h = zeros(size(rain_h)); 
rain_3day_f = zeros(size(rain_f)); 

rain_3day_h(:,:,1:2) = rain_h(:,:,1:2); 
rain_3day_f(:,:,1:2) = rain_f(:,:,1:2); 



for i= 3:10950 
    rain_3day_h(:,:,i) = sum(rain_h(:,:,i-2:i),3); 
    rain_3day_f(:,:,i) = sum(rain_f(:,:,i-2:i),3); 
end

%%
date1 = datetime(1976, 01,01); 
date2 = datetime(2005, 12, 31); 

dates = date1:date2;
dates(month(dates)==2 & day(dates)==29) = []; 

winter_indices = find(month(dates)==1 |month(dates)==2 | month(dates)==3 | ...
    month(dates)==4 | month(dates)==5 | month(dates)==11 |month(dates)==12); 

summer_indices = find(month(dates)==6 | month(dates)==7 | month(dates)==8 | ...
    month(dates)==9 | month(dates)==10); 

%rain_50_h_temp = zeros(size(rain_3day_h)); 
%rain_50_h_temp(rain_3day_h>50) = 1; 

%rain_50_f_temp = zeros(size(rain_3day_f)); 
%rain_50_f_temp(rain_3day_f>50) = 1; 

%rain_50_h_win = sum(rain_50_h_temp(:,:,winter_indices),3);
%rain_50_h_sum = sum(rain_50_h_temp(:,:,summer_indices),3);

%rain_50_f_win = sum(rain_50_f_temp(:,:,winter_indices),3); 
%rain_50_f_sum = sum(rain_50_f_temp(:,:,summer_indices),3); 


%diff_precip_win = rain_50_f_win-rain_50_h_win; 
%diff_precip_sum = rain_50_f_sum-rain_50_h_sum; 

rain_total_h_win = sum(rain_h(:,:,winter_indices),3)./30; 
rain_total_h_sum = sum(rain_h(:,:,summer_indices),3)./30; 

rain_total_f_win = sum(rain_f(:,:,winter_indices),3)./30; 
rain_total_f_sum = sum(rain_f(:,:,summer_indices),3)./30; 

diff_precip_win = rain_total_f_win-rain_total_h_win; 
diff_precip_sum = rain_total_f_sum-rain_total_h_sum; 

%% Plotting this change 

LANDMASK(40:69,60:69)=0; %North of Maine
LANDMASK(40:43,53:69)=0; %North of NH
LANDMASK(43,53:69)=0;
LANDMASK(44,54:69)=0;%Northeast of NH, west of ME
LANDMASK(45,56:69)=0; %west of Maine
LANDMASK(46,58:69)=0; %west of Maine
LANDMASK(34:36,52:69)=0; %North of NY

LANDMASK(49,60:69)=0;
LANDMASK(50,60:69)=0;
LANDMASK(51,56:69)=0; %East of maine
LANDMASK(53:69,40:69)=0; %east of Maine
LANDMASK(52:53,54:69)=0; %east of Maine
LANDMASK(1:40,53:69)=0; %North of NY

LANDMASK(25:27,48:69)=0;
LANDMASK(27:30,49:69)=0; %west of NY---
LANDMASK(30:32,51:69)=0;
LANDMASK(32:33,52:69)=0;
LANDMASK(33:34,54:69)=0;
LANDMASK(34:35,55:69)=0;
LANDMASK(36:37,56:69)=0;

LANDMASK(15:25,50:59)=0;
LANDMASK(15:20,43:69)=0; %west of NY
LANDMASK(15:24,44:69)=0; %west of NY (Lake Ontario area)

LANDMASK(33,51:69)=0;%straight north of NY
LANDMASK(34:35,53:69)=0; %straight north of NY
LANDMASK(28,49:69)=0;
LANDMASK(29,47:69)=0;
LANDMASK(30,47:69)=0; %Lake Ontario area -- northeawest of NY
LANDMASK(31,49:69)=0;
LANDMASK(32,49:69)=0; %Lake Ontario area
LANDMASK(24,45:69)=0;
LANDMASK(26,47:69)=0; %
LANDMASK(25,46:69)=0;

LANDMASK(20:25,48:53)=0; %west of NY
LANDMASK(25:30,52:69)=0;%west of NY
LANDMASK(31,53:69)=0;
LANDMASK(30:35,52:69)=0; %west of NY
LANDMASK(1:20,36:69)=0;%west of WV
LANDMASK(20:21,37:69)=0; %west of WV
LANDMASK(5:15,25:69)=0; %west of WV


LANDMASK(16,27:30)=0;
LANDMASK(16,33:69)=0; %southwest of WV
LANDMASK(17,34:69)=0;
LANDMASK(17,27:29)=0;
LANDMASK(18,35:69)=0;
LANDMASK(19,35:69)=0;

%Trying to get rid of VA
LANDMASK(1:27,15:27)=0; %eliminate VA---trying to include ML-----
LANDMASK(33:41,15:29)=0;
LANDMASK(21,14:28)=0; %southeast of WV


LANDMASK(22,14:29)=0;
LANDMASK(23,14:30)=0; %--
LANDMASK(24,14:31)=0;
LANDMASK(25,14:32)=0;
LANDMASK(26,14:32)=0;
LANDMASK(27,14:33)=0;
LANDMASK(28,14:32)=0;
LANDMASK(29,10:33)=0; 
LANDMASK(30,10:30)=0;
LANDMASK(31,14:30)=0;
LANDMASK(32,14:28)=0;

LANDMASK(25:26,27:31)=0; %eliminate VA
LANDMASK(27,27:32)=0;
LANDMASK(28,27:33)=0;
LANDMASK(29,27:31)=0;
LANDMASK(19:20,15:28)=0; %southwest/west of WV
LANDMASK(18,15:28)=0; %^^^^
LANDMASK(22,27:28)=0; %going up Northeast of WV
LANDMASK(23,28:29)=0; %going up northeast of WV
LANDMASK(24,29:30)=0; %going up Northeast of WV
%%
%Bringing in all states in D1
pcs = {'Vermont','New York','New Hampshire','Rhode Island','Connecticut','Massachusetts','Maine','Pennsylvania','New Jersey','Delaware','West Virginia','Maryland'};
VT = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_VT = [VT.Lat];
lon_VT = [VT.Lon];
NY = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_NY = [NY.Lat];
lon_NY = [NY.Lon];

NH = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_NH = [NH.Lat];
lon_NH = [NH.Lon];

RI = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_RI = [RI.Lat];
lon_RI = [RI.Lon];

CT = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_CT = [CT.Lat];
lon_CT = [CT.Lon];

MA = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_MA = [MA.Lat];
lon_MA = [MA.Lon];

ME = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_ME = [ME.Lat];
lon_ME = [ME.Lon];

PA = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_PA = [PA.Lat];
lon_PA = [PA.Lon];

NJ = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_NJ = [NJ.Lat];
lon_NJ = [NJ.Lon];

DE = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_DE = [DE.Lat];
lon_DE = [DE.Lon];

WV = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_WV = [WV.Lat];
lon_WV = [WV.Lon];

MD = shaperead('usastatelo.shp',...
    'UseGeoCoords', true,...
    'Selector',{@(name)any(strcmpi(name,pcs),2), 'Name'});
lat_MD = [MD.Lat];
lon_MD = [MD.Lon];
%%
hFig_win = figure(1), clf; 
    set(hFig_win, 'Color','white','Units','centimeters','Position', [3 3 24 24]);
    ha = tight_subplot(1,1,[0.09 0],[0.09 0.01],[.09 .01]); %[0.05 0],[0.07 0.01],[.04 .01]
    
   latitude = double(LAT);
   longitude = double(LON);
   
    for i = 1:69
         for j = 1:69
             for k = 1:10950
       
           
            if latitude(i,j) >= 36 & latitude(i,j)<= 47.5
       
                lat_r(i,j) = latitude(i,j);
       
            else 
                lat_r(i,j) = NaN;
       
             end
   
     
            if longitude(i,j) >= -84 & longitude(i,j) <= -66
       
                lon_r(i,j) = longitude(i,j);
       
            else 
                lon_r(i,j) = NaN;
                
            end
       end
         end
    end
 
diff_precip_win(isnan(lat_r) & isnan(lon_r)) = NaN; 
%Setting lat/lon limits
    latlim = [min(min(lat_r)) max(max(lat_r))];
    lonlim = [min(min(lon_r)) max(max(lon_r))];
   
   axes(ha(1));
    axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'Frame','off','Grid','off',...
        'MeridianLabel','off','ParallelLabel','off','FontSize',8,'MLabelLocation',1,'PLabelLocation',1);
    %longitude labels
%     text(0.16,0.037,'81° W','FontSize',8,'FontWeight','normal','Units','normalized'); 
%     text(0.265,0.037,'78° W','FontSize',8,'FontWeight','normal','Units','normalized');
%     text(0.37,0.037,'74° W','FontSize',8,'FontWeight','normal','Units','normalized');
%     text(0.49,0.037,'70° W','FontSize',8,'FontWeight','normal','Units','normalized');

    text(0.22,0.005,'81° W','FontSize',16,'FontWeight','normal','Units','normalized'); 
    text(0.36,0.005,'78° W','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.50,0.005,'74° W','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.69,0.005,'70° W','FontSize',16,'FontWeight','normal','Units','normalized');


    %latitude labels
    text(0.08,0.86,'47° N','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.08,0.63,'44° N','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.08,0.43,'41° N','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.07,0.20,'38° N','FontSize',16,'FontWeight','normal','Units','normalized');
    hold on
    geo_map = geoshow(lat_r,lon_r,diff_precip_win, 'DisplayType','texturemap'); % show annual total
    
    geo_map.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
    geo_map.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
    alpha(geo_map,double(~isnan(diff_precip_win)));
   
    %Latitude lines
    hold on
    plot3m([47 47],[min(min(lon_r))+0.05 max(max(lon_r))-0.1],[1000 1000],':k','linewidth',0.6);% plot lat lines
    plot3m([44 44],[min(min(lon_r))+0.05 max(max(lon_r))-0.1],[1000 1000],':k','linewidth',0.6);
    plot3m([41 41],[min(min(lon_r))+0.05 max(max(lon_r))-0.07],[1000 1000],':k','linewidth',0.6);
    plot3m([38 38],[min(min(lon_r))+0.05 max(max(lon_r))-0.07],[1000 1000],':k','linewidth',0.6);
    %longitude lines
    plot3m([min(min(lat_r))+0.05 max(max(lat_r))-0.05],[-81 -81],[1000 1000],':k','linewidth',0.6);% plot lon line
    plot3m([min(min(lat_r))+0.05 max(max(lat_r))-0.05],[-78 -78],[1000 1000],':k','linewidth',0.6);
    plot3m([min(min(lat_r))+0.05 max(max(lat_r))-0.05],[-74 -74],[1000 1000],':k','linewidth',0.6);
    plot3m([min(min(lat_r))+0.05 max(max(lat_r))-0.05],[-70 -70],[1000 1000],':k','linewidth',0.6);
   
    Z = ones(1,size(lat_NY,2))*1000;% to plot 3-D state line
    lon_NY(1,341)=NaN; lon_NY(1,341)=NaN; lon_NY(1,342)=lon_NY(1,342)-0.02;% process out-of-domain line
    plot3m(lat_NY,lon_NY,Z,'k','linewidth',1.2);
    
    p1 = geoshow(white_watershed.Y, white_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
    geoshow(shen_watershed.Y, shen_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
    geoshow(matt_watershed.Y, matt_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
    geoshow(diam_watershed.Y, diam_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
    
    upperlimit = max(max(diff_precip_win));
    lowerlimit = min(min(diff_precip_win));
    
   
    %Setting color bar specifics--do not need to change this
    %range = floor(-256*lowerlimit/(upperlimit-lowerlimit)); 
    range = floor(-256*-200/(200-(-200))); 

    blueColorMap = [ones(1,range), linspace(1, 0, 256-range)]; 
    redColorMap = [linspace(0, 1, range), linspace(1, 0, 256-range)];
    whiteColorMap = [linspace(0, 1, range), ones(1,256-range)]; 
    colorMap = [blueColorMap; redColorMap; whiteColorMap]'; 
    colormap(gca,colorMap);
    
    hcb=colorbar('southoutside'); %this sets the horizontal color bar below (i.e. south outside) of the figure
    hcb.FontSize = 18;
    hcb.Position = [hcb.Position(1)+0.20 hcb.Position(2)-0.06 0.5 0.015];
    hcb.Label.String = 'Change in average annual precipitation [mm/year]'; 
    set(ha(1),'clim',[-200 200]);
       
    htt=title(['Difference in average annual precipitation - cold season'],'FontWeight','normal');
    %htt.Position = [htt.Position(1) htt.Position(2)-.024 htt.Position(3)];
    set(gca,'FontSize',18);
   
    axis off
    
    saveas(hFig_win, 'EP_change_D1_win_MPI.png'); 
    
%%

hFig_sum = figure(2), clf; 
    set(hFig_sum, 'Color','white','Units','centimeters','Position', [3 3 24 24]);
    ha = tight_subplot(1,1,[0.09 0],[0.09 0.01],[.09 .01]); %[0.05 0],[0.07 0.01],[.04 .01]
    
   latitude = double(LAT);
   longitude = double(LON);
   
    for i = 1:69
         for j = 1:69
             for k = 1:10950
       
           
            if latitude(i,j) >= 36 & latitude(i,j)<= 47.5
       
                lat_r(i,j) = latitude(i,j);
       
            else 
                lat_r(i,j) = NaN;
       
             end
   
     
            if longitude(i,j) >= -84 & longitude(i,j) <= -66
       
                lon_r(i,j) = longitude(i,j);
       
            else 
                lon_r(i,j) = NaN;
                
            end
       end
         end
    end
 
diff_precip_win(isnan(lat_r) & isnan(lon_r)) = NaN; 
%Setting lat/lon limits
    latlim = [min(min(lat_r)) max(max(lat_r))];
    lonlim = [min(min(lon_r)) max(max(lon_r))];
   
   axes(ha(1));
    axesm('lambert','MapLatLimit',latlim,'MapLonLimit',lonlim,'Frame','off','Grid','off',...
        'MeridianLabel','off','ParallelLabel','off','FontSize',8,'MLabelLocation',1,'PLabelLocation',1);
    %longitude labels 
    text(0.22,0.005,'81° W','FontSize',16,'FontWeight','normal','Units','normalized'); 
    text(0.36,0.005,'78° W','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.50,0.005,'74° W','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.69,0.005,'70° W','FontSize',16,'FontWeight','normal','Units','normalized');


    %latitude labels
    text(0.08,0.86,'47° N','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.08,0.63,'44° N','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.08,0.43,'41° N','FontSize',16,'FontWeight','normal','Units','normalized');
    text(0.07,0.20,'38° N','FontSize',16,'FontWeight','normal','Units','normalized');
    hold on
    geo_map = geoshow(lat_r,lon_r,diff_precip_sum, 'DisplayType','texturemap'); % show annual total
    
    geo_map.AlphaDataMapping = 'none'; % interpet alpha values as transparency values
    geo_map.FaceAlpha = 'texturemap'; % Indicate that the transparency can be different each pixel
    alpha(geo_map,double(~isnan(diff_precip_sum)));
   
    %Latitude lines
    hold on
    plot3m([47 47],[min(min(lon_r))+0.05 max(max(lon_r))-0.1],[1000 1000],':k','linewidth',0.6);% plot lat lines
    plot3m([44 44],[min(min(lon_r))+0.05 max(max(lon_r))-0.1],[1000 1000],':k','linewidth',0.6);
    plot3m([41 41],[min(min(lon_r))+0.05 max(max(lon_r))-0.07],[1000 1000],':k','linewidth',0.6);
    plot3m([38 38],[min(min(lon_r))+0.05 max(max(lon_r))-0.07],[1000 1000],':k','linewidth',0.6);
    %longitude lines
    plot3m([min(min(lat_r))+0.05 max(max(lat_r))-0.05],[-81 -81],[1000 1000],':k','linewidth',0.6);% plot lon line
    plot3m([min(min(lat_r))+0.05 max(max(lat_r))-0.05],[-78 -78],[1000 1000],':k','linewidth',0.6);
    plot3m([min(min(lat_r))+0.05 max(max(lat_r))-0.05],[-74 -74],[1000 1000],':k','linewidth',0.6);
    plot3m([min(min(lat_r))+0.05 max(max(lat_r))-0.05],[-70 -70],[1000 1000],':k','linewidth',0.6);
   
    Z = ones(1,size(lat_NY,2))*1000;% to plot 3-D state line
    lon_NY(1,341)=NaN; lon_NY(1,341)=NaN; lon_NY(1,342)=lon_NY(1,342)-0.02;% process out-of-domain line
    plot3m(lat_NY,lon_NY,Z,'k','linewidth',1.2);
    
    p1 = geoshow(white_watershed.Y, white_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
    geoshow(shen_watershed.Y, shen_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
    geoshow(matt_watershed.Y, matt_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
    geoshow(diam_watershed.Y, diam_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
    
    upperlimit = max(max(diff_precip_sum));
    lowerlimit = min(min(diff_precip_sum));
    
   
    %Setting color bar specifics--do not need to change this
    %range = floor(-256*lowerlimit/(upperlimit-lowerlimit)); 
    range = floor(-256*-200/(200-(-200))); 

    blueColorMap = [ones(1,range), linspace(1, 0, 256-range)]; 
    redColorMap = [linspace(0, 1, range), linspace(1, 0, 256-range)];
    whiteColorMap = [linspace(0, 1, range), ones(1,256-range)]; 
    colorMap = [blueColorMap; redColorMap; whiteColorMap]'; 
    colormap(gca,colorMap);
    
    hcb=colorbar('southoutside'); %this sets the horizontal color bar below (i.e. south outside) of the figure
    hcb.FontSize = 18;
    hcb.Position = [hcb.Position(1)+0.20 hcb.Position(2)-0.06 0.5 0.015];
    set(ha(1),'clim',[lowerlimit upperlimit]);
    hcb.Label.String = 'Change in average annual precipitation [mm/year]'; 
    set(ha(1),'clim',[-200 200]);

       
    htt=title(['Difference in average annual precipitation - warm season'],'FontWeight','normal');
    %htt.Position = [htt.Position(1) htt.Position(2)-.024 htt.Position(3)];
    set(gca,'FontSize',18);
   
    axis off
    
    saveas(hFig_sum, 'EP_change_D1_sum_MPI.png'); 
