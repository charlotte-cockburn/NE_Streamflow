clear all, close all
load('decadal_trends_EP.mat');
load('Livneh_Prec_NE.mat');
load('prcp_annual_NE525.mat');
load('prec_year_extr_NE.mat');
load('station_data_525.mat');
load('station_data_yr_525.mat');
load('station_year_extr_NE525.mat'); 
load('stations_NE525_pro.mat'); 
load('stations_NE525.mat');



%%
load('prec_year_extr_NE'); 

white_watershed = shaperead('watershed_WR/globalwatershed.shp'); 
shen_watershed = shaperead('watershed_SR/globalwatershed.shp');
matt_watershed = shaperead('watershed_MR/globalwatershed.shp');
diam_watershed = shaperead('watershed_DR/globalwatershed.shp');



%% Calculating Livneh trend 

    prec_extr_regress_coeff_temp = [];
    prec_extr_regress_coeff = [];
    for x=1:size(prec_year_extr_NE,1)
        x
        for y=1:size(prec_year_extr_NE,2)
            if isnan(prec_year_extr_NE(x,y,size(prec_year_extr_NE,3)))
                prec_extr_regress_coeff(x,y) = NaN;
            else
                datain = [];
                datain = [linspace(1915,2011,33)',squeeze(prec_year_extr_NE(x,y,65:97))];
                [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
                prec_extr_regress_coeff(x,y) = sen;
%                 prec_regress_coeff_temp = polyfit([1:1:size(prec_year_extr_NE(:,:,65:97),3)],squeeze(prec_year_extr_NE(x,y,65:97))',1);
%                 prec_regress_coeff(x,y) = prec_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
            end
        end
    end
    prec_extr_decade_trend(:,:) = 10.*prec_extr_regress_coeff;
    
    
    %%
    
    pre_96 = nanmean(prec_year_extr_NE(:,:,1:81),3); 
    post_96 = nanmean(prec_year_extr_NE(:,:,82:end),3);
    period_diff = post_96-pre_96; 
    

%% LINEAR TREND
figure(1), clf
set(gcf, 'position', [10 10 1000 1000]);
latlim = [37, 47.5];
lonlim = [-82, -67];
ax = usamap(latlim, lonlim);

geoshow(prec_extr_decade_trend(:,:),prec_refvec,'DisplayType','surface',...
        'ZData',zeros(size(prec_extr_decade_trend(:,:))),'CData',prec_extr_decade_trend(:,:));
shading flat;
hold on
geoshow(lat_NE_map,lon_NE_map,'Color','k');
setm(ax,'FontSize',20)

states = shaperead('usastatehi',...
        'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
geoshow(ax, states, 'FaceColor', [1 1 1], 'facealpha', 0)


p1 = geoshow(white_watershed.Y, white_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
textm(44.2, -73.2, 'C', 'FontSize', 20);
geoshow(shen_watershed.Y, shen_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
textm(38.6, -78.7, 'D', 'FontSize', 20)
geoshow(matt_watershed.Y, matt_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
textm(46, -69, 'A', 'FontSize', 20);
geoshow(diam_watershed.Y, diam_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
textm(45, -71, 'B', 'FontSize', 20)

hcb=colorbar('eastoutside');
hcb.Label.String = strcat('Extreme precipitation trend [mm/decade]');
hcb.FontSize = 22;
% cmapb = redblue(20);
% colormap(gca,flipud(cmapb)); % change default color to jet:blue-->red
% caxis([-50 50]);

% leg = legend([p1], 'watershed outline', 'Position', [0.75 0.22 0.01 0.02], 'FontSize', 20);
% set(leg, 'Location', 'Southeast')
%legend({'A','B'},'Position',[0.8 0.1 0.1 0.2])







%% 1996 CHANGEPOINT
map_fig = figure(2), clf
set(gcf, 'position', [10 10 1000 1000]);
latlim = [37, 47.5];
lonlim = [-82, -67];
ax = usamap(latlim, lonlim);

geoshow(period_diff(:,:),prec_refvec,'DisplayType','surface',...
        'ZData',zeros(size(period_diff(:,:))),'CData',period_diff(:,:));
shading flat;
hold on
%geoshow(lat_NE_map,lon_NE_map, 'Color', 'k');
setm(ax,'FontSize',30)

states = shaperead('usastatehi',...
        'UseGeoCoords', true, 'BoundingBox', [lonlim', latlim']);
geoshow(ax, states, 'FaceColor', [1 1 1], 'facealpha', 0)


p1 = geoshow(white_watershed.Y, white_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
textm(44.4, -73.2, 'C', 'FontSize', 30, 'FontWeight', 'Bold');
geoshow(shen_watershed.Y, shen_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
textm(38.6, -78.7, 'D', 'FontSize', 30, 'FontWeight', 'Bold')
geoshow(matt_watershed.Y, matt_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
textm(46, -69.1, 'A', 'FontSize', 30, 'FontWeight', 'Bold');
geoshow(diam_watershed.Y, diam_watershed.X, 'DisplayType', 'Line', 'LineWidth', 2, 'Color', 'Black')
textm(45, -71, 'B', 'FontSize', 30, 'FontWeight', 'Bold')


upperlimit = max(max(period_diff)); 
lowerlimit = min(min(period_diff));

range = floor(-256*lowerlimit/(upperlimit-lowerlimit)); 
blueColorMap = [ones(1,range), linspace(1, 0, 256-range)]; 
redColorMap = [linspace(0, 1, range), linspace(1, 0, 256-range)];
whiteColorMap = [linspace(0, 1, range), ones(1,256-range)]; 
colorMap = [blueColorMap; redColorMap; whiteColorMap]'; 
colormap(gca,colorMap);


hcb=colorbar('eastoutside');
hcb.Label.String = strcat('Change in mean extreme precipitation after 1996 [mm]');
hcb.FontSize = 24;
% cmapb = bluewhitered(255);
% colormap(gca,flipud(cmapb)); % change default color to jet:blue-->red
%caxis([-50 200]);

% leg = legend([p1], 'watershed outline', 'Position', [0.75 0.22 0.01 0.02], 'FontSize', 20);
% set(leg, 'Location', 'Southeast')
%legend({'A','B'},'Position',[0.8 0.1 0.1 0.2])

saveas(map_fig, 'EP_change_map_watersheds.png'); 




