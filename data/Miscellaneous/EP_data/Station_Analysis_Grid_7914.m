clear all
close all

path(path,'/ahome/hhuang/Matlab/Scripts');
path(path,'/ahome/hhuang/Matlab/Scripts/altmany-export_fig-e1b8666');

plot_annual_mean_G = 0; % control running, 0 means off,1 is on
plot_annual_trend_G = 0;  
plot_extr_mean_G = 0; 
plot_extr_trend_G = 0;
plot_season_mean_G = 1;
plot_season_trend_G = 0;
plot_season_extr_mean_G = 0;
plot_season_extr_trend_G = 0;

sav = 0;
load('/ahome/hhuang/Matlab/NE_PRECIP/Proc_Data/NARR_Prec_NE.mat'); 
load '/adata/Data/Observations/GHCN-D/station_data_525_adjust.mat';  % 6*525

yr_start_G = 1979;
yr_end_G = 2014;
nyr_G = yr_end_G - yr_start_G +1;
nmon_G = (yr_end_G-yr_start_G+1)*12;
datenum_start_G = 722816; %19790101
day_mon_G = datenum(yr_start_G*ones(nmon_G+1,1),[1:nmon_G+1]',1);
tday_G = [day_mon_G(1):day_mon_G(1)+length(station_data_525_adjust)-1]';

load '/adata/Data/Observations/GHCN-D/stations_NE525_pro.mat';  % 6*525
geo = []; % 2*525
geo(1,:) = cell2mat(stations_NE525_pro(2,:)); % latitude,convert cell array element to numeric number
geo(2,:) = cell2mat(stations_NE525_pro(3,:)); % longitude

hFig = figure(1); % get the name of entire figure
set(hFig, 'Position', [200 200 720 350]); % [left,bottom,width,height]:define the size of the entire figure. Very critical,since deciding the proportion of X/Y and then gap between subplots
ax0 = gca; % get the handle of entire figure
set(ax0,'visible','off'); % set the whole background axes to be invisible
text(-0.15,0.45,'Latitude(°)','FontSize',8,'FontWeight','normal','Rotation',90); % add Y label to entire figure,[left,bottom] distance to the entire axis point 0
text(0.45,0.05,'Longitude(°)','FontSize',8,'FontWeight','normal'); % add X label to entire figure,[left,bottom] distance to the entire axis point 0

%ha = tight_subplot(3,2,[.06 .04],[.12 .04],[.06 .01]); % outside funtion:tight_subplot. limit the gaps and margins
ha = tight_subplot(1,2,[.06 .06],[.12 .04],[.06 .01]); % outside funtion:tight_subplot. limit the gaps and margins


           % 1.GHCN-D:plot 1979-2014 annual mean PRCP (mm/day)
if plot_annual_mean_G==1    
    load '/adata/Data/Observations/GHCN-D/station_data_yr_525.mat'; % 36*525
    axes(ha(1)); % use subplot 1 by naming axes
    scatter(geo(2,:),geo(1,:),7,0.1*nanmean(station_data_yr,1),'filled','MarkerEdgeColor','k','LineWidth',0.1);
    geoshow(lat_NE_map,lon_NE_map,'Color','k');   
    ax = gca; % get the handle for the current axes
    title(ax,'a) GHCN-D 1979-2014 mean','FontSize',8,'FontWeight','normal');
    set(ax,'FontSize',8);
    xlim([-84 -66]);
    ylim([36 48]);
    set(gca,'Xtick',linspace(-84,-66,7));
    set(gca,'Ytick',linspace(36,48,5));
    colormap(ax,flipud(jet)); % change default color to jet:blue-->red
    
    ax.CLim=[750 1450];
    A = colorbar('southoutside'); % put colorbar at the bottom
    set(A, 'Position', [.085 .1 .38 .03]); %[left,bottom,width,height]
    A.Label.String = 'mm/yr'; % name colorbar title
    A.Label.FontSize = 8;

    if(sav==1)
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
        saveas(gcf,'Figures/annual_mean_G.jpg'); % need to manually save .jpg by replacing with same name
% saveas(gcf,'Figures/annual_mean_G.epsc'); % need to manually save .jpg by replacing with same name
        %set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
%         print('-dpdf', 'Figures/annual_mean_G.pdf');
    end
    hold on
    %figure
end        

           % 1'.GHCN-D:plot trend in 1979-2014 annual mean PRCP (mm/decade)
if plot_annual_trend_G==1
    load '/adata/Data/Observations/GHCN-D/prcp_annual_NE525.mat';  % 11*16*36
    % calculate regional trend
    %   M-K test and Thei-sen slope
    datain = [linspace(1979,2014,36)',0.1*squeeze(nanmean(nanmean(prcp_annual_NE,2)))];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
    slope_MK = 10*sen
    p_value_MK = sig
    stop
    H_MK = h
    %   linear trend estimate
    disp('    Hereafter is linear trend estimate');
    mdl = fitlm(linspace(1979,2014,36)',0.1*squeeze(nanmean(nanmean(prcp_annual_NE,2)))); % provide linear trend and p-value
%     figure
%     plotResiduals(mdl)
    [H1,p] = lillietest(mdl.Residuals.Raw) % H1=0 means a normal distribution at default alpha level 0.05
    h1 = lbqtest(mdl.Residuals.Raw) % h1=0 means residuals are independent/non-autocorrelated
    slope = 10*table2array(mdl.Coefficients(2,1)) % convert to mm/decade
    p_value = table2array(mdl.Coefficients(2,4)) % for linear regression model

    B2 = linspace(1979,2014,36)';
    C2 = 0.1*squeeze(nanmean(nanmean(prcp_annual_NE(:,:,:),2)));
    my_poly2=polyfit(B2,C2,1)
    X2 = 1979 % X data range
    Y2 = polyval(my_poly2,X2)
    PE=10*my_poly2(1,1)/Y2
stop    


    % calculate and plot linear trend
    prcp_regress_coeff_temp = [];
    prcp_regress_coeff = [];
    for x=1:525  % if compare full period,use commentted two lines to replace the line 'extr_regress_coeff_temp...'
        ind = 1:length(station_data_yr(:,x));
        k = ~isnan(station_data_yr(:,x)); % find index with non-NaN values
        %     y = find(k==1);
        %     pref_regress_coeff_temp = polyfit(ind(k)',36.5246.*station_data_yr(78+y,x+4),1);
        prcp_regress_coeff_temp = polyfit(ind(k)',0.1*station_data_yr(k,x),1);
        prcp_regress_coeff(x) = prcp_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
    end
    
    % % analyze 1979-2014 trend only
    % for x=1:525
    %     ind = 1:length(station_data_yr(15:111,x+4)); % analyze 1979-2014 trend only
    %     k = ~isnan(station_data_yr(15:111,x+4)); % find index with non-NaN values
    %     %pref_regress_coeff_temp = polyfit(ind(k)',36.5246.*station_data_yr(k,x+4),1);
    %     y = find(k==1);
    %     pref_regress_coeff_temp = polyfit(ind(k)',36.5246.*station_data_yr(14+y,x+4),1);
    %     pref_regress_coeff(x) = pref_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
    % end
    
    prcp_decade_trend = 10.*prcp_regress_coeff;  % changing yearly trends to be mm/decade
    
    axes(ha(2)); % use subplot 2 by naming axes
    scatter(geo(2,:),geo(1,:),7,prcp_decade_trend,'filled','MarkerEdgeColor','k','LineWidth',0.1);
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    ax = gca; % get the handle for the current axes
    title(ax,'b) GHCN-D 1979-2014 trend','FontSize',8,'FontWeight','normal');
    set(ax,'FontSize',8);
    xlim([-84 -66]);
    ylim([36 48]);
    set(gca,'Xtick',linspace(-84,-66,7));
    set(gca,'Ytick',linspace(36,48,5));
    blueColorMap = [ones(1,128), linspace(1, 0, 128)]; % [0,0,1] refers to blue; 525 and 140 is calculated from colorbar limits [-49,59], to make 0 be white
    redColorMap = [linspace(0, 1, 128), linspace(1, 0, 128)];% [1,0,0] refers to red
    whiteColorMap = [linspace(0, 1, 128), ones(1,128)]; % [1,1,1] refers to white
    colorMap = [blueColorMap; redColorMap; whiteColorMap]'; % creat a 3-D matrix to install RGB values
    colormap(ax,colorMap); % change default color to blue-->white-->red
    hold on;
    
    ax.CLim=[-49 49];
    B = colorbar('southoutside'); % put colorbar at the bottom
    set(B, 'Position', [.58 .1 .38 .03]); %[left,bottom,width,height]
    B.Label.String = 'mm/decade'; % name colorbar title
    B.Label.FontSize = 8;
    
    export_fig Figures/annual_mean_G525 -eps -transparent -depsc % use outside function: export_fig

   % saveas(gcf,'Figures/test4.epsc');

    if(sav==1)
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
        saveas(gcf,'Figures/annual_mean_trend_G.jpg'); % need to manually save .jpg by replacing with same name
    end
end

                % 2.GHCN-D:plot 1979-2014 mean extreme prec
if plot_extr_mean_G == 1
    % All commentted acripts below is to produce this .mat. if saved out, load it directly
%    load '/adata/Data/Observations/GHCN-D/station_year_extr_NE.mat'; % 36 yr*525 st
    % calculate mean values of annual extreme prec, achieve a 36*525 .mat
%    load '/adata/Data/Observations/GHCN-D/station_data_525_adjust.mat'; % 41638*120
    % sort the sequence and find the 1% heaviest prec threshold for every station
    station_year_extr_NE = [];
    threshold = [];
    for i=1:525
        threshold_temp = sort(station_data_525_adjust(:,i+4),'descend');
        threshold_temp = threshold_temp(~isnan(threshold_temp)); % exclude NaNs
        %r = round(0.01*size(find(threshold_temp),1)); % find non-zero elements,get number of 1% heaviest daily events
        count = 0;
        for j=1:size(threshold_temp,1)
            if 0.1*threshold_temp(j)<0.254 % <0.01 inch is trace precip
                break
            end
            count = count+1;
        end
        r = fix(0.01*count);
        threshold(i) = threshold_temp(r+1); % find 1%+1 heaviest prec value
    end
  
    % calculate staion extreme prec
    for i=1:525
        for j=1:nyr_G % 36
            station_year_extr_temp=[];
            extr_temp = 0;
            yr_daynum_start = day_mon_G((j-1)*12+1) - day_mon_G(1) + 1; % start date number
            k = day_mon_G(j*12+1) - day_mon_G((j-1)*12+1);
            yr_daynum_end = yr_daynum_start + k -1;   % end date number
            station_year_extr_temp = sort(station_data_525_adjust(yr_daynum_start:yr_daynum_end,i+4),'descend');
            x = ~isnan(station_year_extr_temp);
            y = find(x==1);
            if sum(y)==0
                station_year_extr_NE(j,i) = NaN; % set 100% missing value years to NaN
            else
                for k=1:size(station_year_extr_temp) % 365/366,calculate by year by station
                    if station_year_extr_temp(k)>threshold(i)
                       extr_temp = extr_temp + station_year_extr_temp(k); % larger than threshold,sum up 
                    end
                end
                station_year_extr_NE(j,i) = extr_temp; % 36*525
            end
        end
    end

    % grid station extreme prec to 1*1 degree
    lon0 = prec_refvec(3,1); % upper left corner longitude -82.6563
    lat0 = prec_refvec(3,2); % upper left corner latitude 47.5313
    grid_lat = fix(length(lat_NE)*prec_refvec(2,1))+1; % NE span 11 degree
    grid_lon = fix(length(lon_NE)*prec_refvec(2,1))+1; % NE span 16 degree
    prcp_year_extr_NE = NaN(grid_lat,grid_lon,nyr_G); % grid station data to 11*16    
    for i=1:grid_lat
        for j=1:grid_lon
            count = 0;
            index_grid = [];
            for k=1:length(geo) % 525
                if geo(1,k)>lat0-i & geo(1,k)<=lat0-i+1 & ...
                        geo(2,k)<lon0+j & geo(2,k)>=lon0+j-1
                    count = count+1;
                    index_grid(count) = k;
                end
            end
            
            if length(index_grid)~=0
                prcp_year_extr_NE(i,j,:) = squeeze(nanmean(station_year_extr_NE(:,index_grid),2));
            else
                prcp_year_extr_NE(i,j,:) = NaN;
            end
        end
    end
STOP    
%      % JMW Begin
%  start_year = 1958;
%  plot([start_year:2012],nanmean(station_year_extr_NE(start_year-1900:112,:),2)*.1,'r')
%  hold on
%  plot([start_year:2012],squeeze(nanmean(nanmean(prcp_year_extr_NE(:,:,start_year-1900:112),1),2))*.1,'b')
%  stat_lin = polyfit([start_year:2012],nanmean(station_year_extr_NE(start_year-1900:112,:),2)'*.1,1);
%  grid_lin = polyfit([start_year:2012],squeeze(nanmean(nanmean(prcp_year_extr_NE(:,:,start_year-1900:112),1),2))'*.1,1);
%  plot([start_year:2012],stat_lin(1)*[start_year:2012]+stat_lin(2),'r')
%  plot([start_year:2012],grid_lin(1)*[start_year:2012]+grid_lin(2),'b')
%  legend('Station','Grid')
%  
%  mean_grid_del = (grid_lin(1)*[2012]+grid_lin(2) - (grid_lin(1)*[start_year]+grid_lin(2)))/mean([grid_lin(1)*[start_year]+grid_lin(2) grid_lin(1)*[2012]+grid_lin(2)])
%  mean_stat_del = (stat_lin(1)*[2012]+stat_lin(2) - (stat_lin(1)*[start_year]+stat_lin(2)))/mean([stat_lin(1)*[start_year]+stat_lin(2) stat_lin(1)*[2012]+stat_lin(2)])
%  year_grid_del = (grid_lin(1)*[2012]+grid_lin(2) - (grid_lin(1)*[start_year]+grid_lin(2)))/(grid_lin(1)*[start_year]+grid_lin(2))
%  year_stat_del = (stat_lin(1)*[2012]+stat_lin(2) - (stat_lin(1)*[start_year]+stat_lin(2)))/(stat_lin(1)*[start_year]+stat_lin(2))
%  % JMW End
% STOP

    % calculate regional trend
    %   M-K test and Thei-sen slope
    datain = [linspace(1979,2014,36)',0.1*squeeze(nanmean(nanmean(prcp_year_extr_NE,2)))];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
    slope_MK = 10*sen
    p_value_MK = sig
    H_MK = h % =1 means significant
    B2 = linspace(1979,2014,36)';
    C2 = 0.1*squeeze(nanmean(nanmean(prcp_year_extr_NE(:,:,:),2)));
    X2 = 1979:1:2014; % X data range
    [a,b] = rlinfit(B2,C2) %[slope,intercept] from Theil-sen
    %   linear trend estimate
    disp('    Hereafter is linear trend estimate');
    mdl = fitlm(linspace(1979,2014,36)',0.1*squeeze(nanmean(nanmean(prcp_year_extr_NE,2)))); % provide linear trend and p-value
%     figure
%     plotResiduals(mdl)
    [H1,p] = lillietest(mdl.Residuals.Raw) % H1=0 means a normal distribution at default alpha level 0.05
    h1 = lbqtest(mdl.Residuals.Raw) % h1=0 means residuals are independent/non-autocorrelated
    slope = 10*table2array(mdl.Coefficients(2,1)) % convert to mm/decade
    p_value = table2array(mdl.Coefficients(2,4)) % for linear regression model
STOP
    axes(ha(1)); % use subplot i by naming axes
    scatter(geo(2,:),geo(1,:),17,0.1.*nanmean(station_year_extr_NE),'filled','MarkerEdgeColor','k','LineWidth',0.1);
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    ax = gca; % get the handle for the current axes
    title(ax,'a) GHCN-D 1979-2014 mean','FontSize',8,'FontWeight','normal');
    set(ax,'FontSize',8);
    xlim([-84 -66]);
    ylim([36 48]);
    set(gca,'Xtick',linspace(-84,-66,7));
    set(gca,'Ytick',linspace(36,48,5));
    colormap(ax,flipud(jet)); % change default color to jet:blue-->red
 
    ax.CLim=[45 155];
    A = colorbar('southoutside'); % put colorbar at the bottom
    set(A, 'Position', [.085 .1 .38 .03]); %[left,bottom,width,height]
    A.Label.String = 'mm/yr'; % name colorbar title
    A.Label.FontSize = 8;
%     if(sav==1)
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
%         saveas(gcf,'Figures/extr_mean_trend_G.jpg'); % need to manually save .jpg by replacing with same name
%     end
%     figure 
%     hold on
end

             % 2'.GHCN-D:plot trend in 1979-2014 mean extreme prec
if plot_extr_trend_G == 1
    % calculate and plot linear trend
    extr_regress_coeff_temp = [];
    extr_regress_coeff = [];
    for x=1:525        
        ind = 1:length(station_year_extr_NE(:,x));
        k = ~isnan(station_year_extr_NE(:,x)); % find index with non-NaN values
%         extr_regress_coeff_temp = polyfit(ind(k)',0.1*station_year_extr_NE(k,x),1);
%         extr_regress_coeff(x) = extr_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
        datain = [];
        datain = [ind(k)',0.1*station_year_extr_NE(k,x)];
        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
        extr_regress_coeff(x) = sen; % Thei-sen slope
    end
    
    extr_decade_trend = 10.*extr_regress_coeff;  % changing yearly trends to be mm/decade
    
    % to compute regional trend:10*polyfit(linspace(1979,1914,36),0.1*nanmean(station_year_extr_NE,2)',1)
    axes(ha(2));
    scatter(geo(2,:),geo(1,:),17,extr_decade_trend,'filled','MarkerEdgeColor','k','LineWidth',0.1); % multiply by 0.1 to convert to unit 1mm
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    ax = gca; % get the handle for the current axes
    title(ax,'b) GHCN-D 1979-2014 trend','FontSize',8,'FontWeight','normal');
    set(ax,'FontSize',8);
    xlim([-84 -66]);
    ylim([36 48]);
    set(gca,'Xtick',linspace(-84,-66,7));
    set(gca,'Ytick',linspace(36,48,5));

    ax.CLim=[-6 10]; % set the min and max limits for the colorbar
    blueColorMap = [ones(1,96), linspace(1, 0, 160)]; % [0,0,1] refers to blue; 525 and 140 is calculated from colorbar limits [-49,59], to make 0 be white
    redColorMap = [linspace(0, 1, 96), linspace(1, 0, 160)];% [1,0,0] refers to red
    whiteColorMap = [linspace(0, 1, 96), ones(1,160)]; % [1,1,1] refers to white
    colorMap = [blueColorMap; redColorMap; whiteColorMap]'; % creat a 3-D matrix to install RGB values
    colormap(ax,colorMap); % change default color to blue-->white-->red
    hold on;
    
    B = colorbar('southoutside'); % put colorbar at the bottom
    set(B, 'Position', [.58 .1 .38 .03]); %[left,bottom,width,height]
    B.Label.String = 'mm/decade'; % name colorbar title
    B.Label.FontSize = 8;
    
    export_fig Figures/extr_mean_G -eps -transparent -depsc % use outside function: export_fig
 %    saveas(gcf,'Figures/test4.epsc');
        
    if(sav==1)
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
        saveas(gcf,'Figures/extr_mean_trend_G.jpg'); % need to manually save .jpg by replacing with same name
    end
end

              % 3. GHCN-D:plot 1979-2014 seasonal total prec
if plot_season_mean_G == 1
    for i=1:36
        for j=1:525
          spring_start = day_mon_G(i*12-9)-datenum_start_G;
          spring_end = day_mon_G(i*12-6)-datenum_start_G-1;
          spring_days = spring_end-spring_start+1;
          station_year_spring(i,j) = spring_days*nanmean(station_data_525_adjust(spring_start:spring_end,j+4),1);

          summer_start = day_mon_G(i*12-6)-datenum_start_G;
          summer_end = day_mon_G(i*12-3)-datenum_start_G-1;
          summer_days = summer_end-summer_start+1;
          station_year_summer(i,j) = summer_days*nanmean(station_data_525_adjust(summer_start:summer_end,j+4),1);

          autumn_start = day_mon_G(i*12-3)-datenum_start_G;
          autumn_end = day_mon_G(i*12)-datenum_start_G-1;
          autumn_days = autumn_end-autumn_start+1;
          station_year_autumn(i,j) = autumn_days*nanmean(station_data_525_adjust(autumn_start:autumn_end,j+4),1);

          if i==36
              station_year_winter(i,j) = NaN;
          else
              winter_start = day_mon_G(i*12)-datenum_start_G;
              winter_end = day_mon_G(i*12+3)-datenum_start_G-1;
              winter_days = winter_end-winter_start+1;
              station_year_winter(i,j) = winter_days*nanmean(station_data_525_adjust(winter_start:winter_end,j+4),1);
          end
        end
    end
    
        % grid station extreme prec to 1*1 degree
    lon0 = prec_refvec(3,1); % upper left corner longitude -82.6563
    lat0 = prec_refvec(3,2); % upper left corner latitude 47.5313
    grid_lat = fix(length(lat_NE)*prec_refvec(2,1))+1; % NE span 11 degree
    grid_lon = fix(length(lon_NE)*prec_refvec(2,1))+1; % NE span 16 degree
    prcp_season_NE = NaN(grid_lat,grid_lon,nyr_G); % grid station data to 11*16    
    for i=1:grid_lat
        for j=1:grid_lon
            count = 0;
            index_grid = [];
            for k=1:length(geo) % 525
                if geo(1,k)>lat0-i & geo(1,k)<=lat0-i+1 & ...
                        geo(2,k)<lon0+j & geo(2,k)>=lon0+j-1
                    count = count+1;
                    index_grid(count) = k;
                end
            end
            
            if length(index_grid)~=0
                prcp_season_NE(i,j,:) = squeeze(nanmean(station_year_winter(:,index_grid),2));
            else
                prcp_season_NE(i,j,:) = NaN;
            end
        end
    end
    
    % calculate regional trend
    %   M-K test and Thei-sen slope
    datain = [linspace(1979,2014,36)',0.1*squeeze(nanmean(nanmean(prcp_season_NE(:,:,1:36),2)))];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
    slope_MK = 10*sen
    p_value_MK = sig
    stop
    H_MK = h % =1 means significant
    %   linear trend estimate
    disp('    Hereafter is linear trend estimate');
    mdl = fitlm(linspace(1979,2014,36)',0.1*squeeze(nanmean(nanmean(prcp_season_NE(:,:,1:36),2)))); % provide linear trend and p-value
%     figure
%     plotResiduals(mdl)
    [H1,p] = lillietest(mdl.Residuals.Raw) % H1=0 means a normal distribution at default alpha level 0.05
    h1 = lbqtest(mdl.Residuals.Raw) % h1=0 means residuals are independent/non-autocorrelated
    slope = 10*table2array(mdl.Coefficients(2,1)) % convert to mm/decade
    p_value = table2array(mdl.Coefficients(2,4)) % for linear regression model

    B2 = linspace(1979,2013,35)';
    C2 = 0.1*squeeze(nanmean(nanmean(prcp_season_NE(:,:,1:35),2)));
    my_poly2=polyfit(B2,C2,1)
    X2 = 1979 % X data range
    Y2 = polyval(my_poly2,X2)
    PE=10*my_poly2(1,1)/Y2
    
    STOP
    
    scatter(geo(2,:),geo(1,:),100,0.1.*nanmean(station_year_summer),'filled','MarkerEdgeColor','k');
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    colorbar;
    title('GHCN-D 1979-2014 annual mean summer precipitation (mm)');
    ax = gca; % get the handle for the current axes
    set(ax,'FontSize',28);
    % ax.CLim=[-49 59]; % set the min and max limits for the colorbar
    colormap(flipud(jet)); % change default color to jet:blue-->red
    
    hold on
    figure
end

           % 3'.GHCN-D:plot trend in 1979-2014 seasonal PRCP (mm/decade)
if plot_season_trend_G==1
    % calculate and plot linear trend
    prcp_season_regress_coeff_temp = [];
    prcp_season_regress_coeff = [];
    for x=1:525  % if compare full period,use commentted two lines to replace the line 'extr_regress_coeff_temp...'
        ind = 1:length(station_year_winter(:,x));
        k = ~isnan(station_year_winter(:,x)); % find index with non-NaN values
        prcp_season_regress_coeff_temp = polyfit(ind(k)',0.1*station_year_winter(k,x),1);
        prcp_season_regress_coeff(x) = prcp_season_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
    end
    
    prcp_season_decade_trend = 10.*prcp_season_regress_coeff;  % changing yearly trends to be mm/decade
    
    scatter(geo(2,:),geo(1,:),100,prcp_season_decade_trend,'filled','MarkerEdgeColor','k', 'LineWidth',1);
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    colorbar;
    title('GHCN-D trends in 1979-2014 Winter precipitation (mm/decade)');
    ax = gca; % get the handle for the current axes
    set(ax,'FontSize',28);
    colormap(flipud(jet));
%     ax.CLim=[-49 59]; % set the min and max limits for the colorbar
%     blueColorMap = [ones(1,525), linspace(1, 0, 140)]; % [0,0,1] refers to blue; 525 and 140 is calculated from colorbar limits [-49,59], to make 0 be white
%     redColorMap = [linspace(0, 1, 525), linspace(1, 0, 140)];% [1,0,0] refers to red
%     whiteColorMap = [linspace(0, 1, 525), ones(1,140)]; % [1,1,1] refers to white
%     colorMap = [blueColorMap; redColorMap; whiteColorMap]'; % creat a 3-D matrix to install RGB values
%     colormap(colorMap); % change default color to blue-->white-->red
%         
%     if(sav==1)
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
%         saveas(gcf,'Figures/annual_mean_trend_G.jpg'); % need to manually save .jpg by replacing with same name
%     end
end


                % 4.GHCN-D:plot 1979-2014 seasonal extreme prec
if plot_season_extr_mean_G == 1
    % extract daily data by season for every station
    season = 'winter';
    station_season_data = []; % store seasonal daily data
    for i=1:36
        station_season_data_temp = []; % store one year-season only
        extr_temp = 0;
        if strcmp(season,'spring')
            season_start_G = day_mon_G(i*12-9)-datenum_start_G;
            season_end_G = day_mon_G(i*12-6)-datenum_start_G-1;
        elseif strcmp(season,'summer')
            season_start_G = day_mon_G(i*12-6)-datenum_start_G;
            season_end_G = day_mon_G(i*12-3)-datenum_start_G-1;
        elseif strcmp(season,'autumn')
            season_start_G = day_mon_G(i*12-3)-datenum_start_G;
            season_end_G = day_mon_G(i*12)-datenum_start_G-1;
        elseif strcmp(season,'winter')
            if i==36
                break
            else
                season_start_G = day_mon_G(i*12)-datenum_start_G;
                season_end_G = day_mon_G(i*12+3)-datenum_start_G-1;
            end
        end
        
        station_season_data_temp = station_data_525_adjust(season_start_G:season_end_G,5:529);
        station_season_data = [station_season_data;station_season_data_temp]; %
    end
    clear season_start_G season_end_G
    
    % sort the sequence and find the 1% heaviest prec threshold for every station        
    threshold = [];
    for i=1:525
        threshold_temp = sort(station_season_data(:,i),'descend');
        threshold_temp = threshold_temp(~isnan(threshold_temp)); % exclude NaNs
        %r = round(0.01*size(find(threshold_temp),1)); % find non-zero elements,get number of 1% heaviest daily events
        count = 0;
        for j=1:size(threshold_temp,1)
            if 0.1*threshold_temp(j)<0.254 % <0.01 inch is trace precip
                break  % we used 'sort' above, so exit the inner layer for loop
            end
            count = count+1;
        end
        r = fix(0.01*count);
        threshold(i) = threshold_temp(r+1); % find 1%+1 heaviest prec value
    end    

    % calculate seasonal extreme prec by stations
    for i=1:36
        for j=1:525
            station_season_extr_temp=[];
            extr_temp = 0;   
            if strcmp(season,'spring')
                season_start_G = day_mon_G(i*12-9)-datenum_start_G;
                season_end_G = day_mon_G(i*12-6)-datenum_start_G-1;
            elseif strcmp(season,'summer')
                season_start_G = day_mon_G(i*12-6)-datenum_start_G;
                season_end_G = day_mon_G(i*12-3)-datenum_start_G-1;
            elseif strcmp(season,'autumn')
                season_start_G = day_mon_G(i*12-3)-datenum_start_G;
                season_end_G = day_mon_G(i*12)-datenum_start_G-1;
            elseif strcmp(season,'winter')
                if i==36
                    station_season_extr_NE(i,j) = NaN;
                else
                    season_start_G = day_mon_G(i*12)-datenum_start_G;
                    season_end_G = day_mon_G(i*12+3)-datenum_start_G-1;
                end
            end 
            
            station_season_extr_temp = sort(station_data_525_adjust(season_start_G:season_end_G,j+4),'descend');
            x = ~isnan(station_season_extr_temp);
            y = find(x==1);
            if sum(y)==0
                station_season_extr_NE(i,j) = NaN; % set 100% missing value years to NaN
            else
                for k=1:season_end_G-season_start_G+1 % 91/92,calculate by year by station
                    if station_season_extr_temp(k)>threshold(j)
                       extr_temp = extr_temp + station_season_extr_temp(k); % larger than threshold,sum up 
                    end
                end
                station_season_extr_NE(i,j) = extr_temp; % 97*176
            end
            if strcmp(season,'winter')
                if i==36
                    station_season_extr_NE(i,j) = NaN;
                end
            end
        end
    end
    
        % grid station extreme prec to 1*1 degree
    lon0 = prec_refvec(3,1); % upper left corner longitude -82.6563
    lat0 = prec_refvec(3,2); % upper left corner latitude 47.5313
    grid_lat = fix(length(lat_NE)*prec_refvec(2,1))+1; % NE span 11 degree
    grid_lon = fix(length(lon_NE)*prec_refvec(2,1))+1; % NE span 16 degree
    prcp_season_extr_NE = NaN(grid_lat,grid_lon,nyr_G); % grid station data to 11*16    
    for i=1:grid_lat
        for j=1:grid_lon
            count = 0;
            index_grid = [];
            for k=1:length(geo) % 116
                if geo(1,k)>lat0-i & geo(1,k)<=lat0-i+1 & ...
                        geo(2,k)<lon0+j & geo(2,k)>=lon0+j-1
                    count = count+1;
                    index_grid(count) = k;
                end
            end
            
            if length(index_grid)~=0
                prcp_season_extr_NE(i,j,:) = squeeze(nanmean(station_season_extr_NE(:,index_grid),2));
            else
                prcp_season_extr_NE(i,j,:) = NaN;
            end
        end
    end
    
    % calculate regional trend
    %   M-K test and Thei-sen slope
    datain = [linspace(1979,2014,36)',0.1*squeeze(nanmean(nanmean(prcp_season_extr_NE(:,:,1:36),2)))];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
    slope_MK = 10*sen
    p_value_MK = sig
    H_MK = h % =1 means significant
    %   linear trend estimate
    disp('    Hereafter is linear trend estimate');
    mdl = fitlm(linspace(1979,2014,36)',0.1*squeeze(nanmean(nanmean(prcp_season_extr_NE(:,:,1:36),2)))); % provide linear trend and p-value
%     figure
%     plotResiduals(mdl)
    [H1,p] = lillietest(mdl.Residuals.Raw) % H1=0 means a normal distribution at default alpha level 0.05
    h1 = lbqtest(mdl.Residuals.Raw) % h1=0 means residuals are independent/non-autocorrelated
    slope = 10*table2array(mdl.Coefficients(2,1)) % convert to mm/decade
    p_value = table2array(mdl.Coefficients(2,4)) % for linear regression model
    
    threshold_mean = 0.1*sum(threshold)/size(threshold,2)
    season_extr_mean = 0.1*nanmean(nanmean(nanmean(prcp_season_extr_NE(:,:,1:36),2)))
    slope_MK
    H_MK
        % calculate change in %/decade
    B2 = linspace(1979,2013,35)';
    C2 = 0.1*squeeze(nanmean(nanmean(prcp_season_extr_NE(:,:,1:35),2)));
    X2 = 1979:1:2013; % X data range
    [a,b] = rlinfit(B2,C2) %[slope,intercept] from Theil-sen
    change = slope_MK/(1979*a+b)
    STOP
    
    scatter(geo(2,:),geo(1,:),20,0.1.*nanmean(station_season_extr_NE(1:113,:)),...
        'filled','MarkerEdgeColor','k', 'LineWidth',0.5); % if winter, eliminate 36th year
    geoshow(lat_NE_map,lon_NE_map,'Color','k','LineWidth',0.5);
    colorbar;
    title('GHCN-D 1979-2014 spring mean extreme precipitation (mm)');
    ax = gca; % get the handle for the current axes
    %set(ax,'FontSize',28);
    % ax.CLim=[-49 59]; % set the min and max limits for the colorbar
    colormap(flipud(jet)); % change default color to jet:blue-->red
    
%     hold on
%     figure
    
    if(sav==1)
        %set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
        %set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
        %saveas(gcf,'Figures/annual_mean_G.jpg'); % need to manually save .jpg by replacing with same name        
        %H=fig('units','inches','width',7,'height',2,'font','Helvetica','fontsize',16);
        saveas(gcf,'Figures/test.epsc'); % need to manually save .jpg by replacing with same name
        %set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20]);
        %         print('-dpdf', 'Figures/annual_mean_G.pdf');
    end
end

           % 4'.GHCN-D:plot trend in 1979-2014 seasonal extr PRCP (mm/decade)
if plot_season_extr_trend_G==1
    % calculate and plot linear trend
    prcp_season_extr_regress_coeff_temp = [];
    prcp_season_extr_regress_coeff = [];
    for x=1:525  % if compare full period,use commentted two lines to replace the line 'extr_regress_coeff_temp...'
        ind = 1:length(station_season_extr_NE(:,x));
        k = ~isnan(station_season_extr_NE(:,x)); % find index with non-NaN values
        prcp_season_extr_regress_coeff_temp = polyfit(ind(k)',0.1*station_season_extr_NE(k,x),1);
        prcp_season_extr_regress_coeff(x) = prcp_season_extr_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
    end
    
    prcp_season_extr_decade_trend = 10.*prcp_season_extr_regress_coeff;  % changing yearly trends to be mm/decade
    
    scatter(geo(2,:),geo(1,:),100,prcp_season_extr_decade_trend,'filled','MarkerEdgeColor','k', 'LineWidth',1);
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    colorbar;
    title('GHCN-D trends in 1979-2014 winter extreme precipitation (mm/decade)');
    ax = gca; % get the handle for the current axes
    set(ax,'FontSize',28);
    colormap(flipud(jet));
%     ax.CLim=[-49 59]; % set the min and max limits for the colorbar
%     blueColorMap = [ones(1,525), linspace(1, 0, 140)]; % [0,0,1] refers to blue; 525 and 140 is calculated from colorbar limits [-49,59], to make 0 be white
%     redColorMap = [linspace(0, 1, 525), linspace(1, 0, 140)];% [1,0,0] refers to red
%     whiteColorMap = [linspace(0, 1, 525), ones(1,140)]; % [1,1,1] refers to white
%     colorMap = [blueColorMap; redColorMap; whiteColorMap]'; % creat a 3-D matrix to install RGB values
%     colormap(colorMap); % change default color to blue-->white-->red
%         
%     if(sav==1)
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
%         saveas(gcf,'Figures/annual_mean_trend_G.jpg'); % need to manually save .jpg by replacing with same name
%     end
end