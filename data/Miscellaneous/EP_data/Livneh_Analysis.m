% using definition for extreme prec by each year(sorting heaviest 4 days),
% decade is defined as 1-10, just like what NCA did

        % Task list:
% 0.Reproduce fig 2.12 in NCA
% 0'.Reproduce fig 2.17 in NCA (line 72)
% 1.Livneh+GHCND:plot 1915-2011 annual mean prec (line 125)
% 1'.Livneh+GHCND:calculate and plot trend in 1915-2011 annual mean prec
% 2.Livneh+GHCND:plot 1915-2011 mean extreme prec (line 190)
% 2'.Livneh+GHCND:plot trend in 1915-2011 mean extreme prec (line 222)
% 3.Livneh+GHCND+NARR:plot bar charts for decadel changes in annual mean prec (line 265)
% 3'.Livneh+GHCND+NARR:plot bar charts for decadel changes in extreme prec

clear all
close all

plot_annual_change_Livneh = 0; % 0.reproduce fig.2.12 in NCA
plot_extr_change_Livneh = 0; % 0.reproduce fig.2.17 in NCA
plot_annual_mean_GL = 1; % control running, 0 means off,1 is on
plot_annual_trend_GL = 0;  
plot_extr_mean_GL = 0; 
plot_extr_trend_GL = 0;
plot_mean_decadel_change_GLN = 0;
plot_extr_decadel_change_GLN = 0;
plot_season_mean_GL = 0;
plot_season_trend_GL = 0;
plot_season_extr_mean_GL = 0;
plot_season_extr_trend_GL = 0;

sav=0;
path(path,'/ahome/hhuang/Matlab/Scripts'); % store extra functions downloaded from File Exchange
load('/ahome/hhuang/Matlab/NE_PRECIP/Proc_Data/Livneh_Prec_NE.mat'); 

%  load GHCN-D data, and latitude/longitude
load '/adata/Data/Observations/GHCN-D/stations_NE176_pro.mat';  % 6*176
load '/adata/Data/Observations/GHCN-D/station_data_yr_176.mat'; % 97*176
geo = [];
geo(1,:) = cell2mat(stations_NE176_pro(2,:)); % convert cell array element to numeric number
geo(2,:) = cell2mat(stations_NE176_pro(3,:));

hFig = figure(1); % get the name of entire figure
set(hFig, 'Position', [200 200 720 600]); % [left,bottom,width,height]:define the size of the entire figure. Very critical,since deciding the proportion of X/Y and then gap between subplots
ax0 = gca; % get the handle of entire figure
set(ax0,'visible','off'); % set the whole background axes to be invisible
text(-0.15,0.45,'Latitude(°)','FontSize',8,'FontWeight','normal','Rotation',90); % add Y label to entire figure,[left,bottom] distance to the entire axis point 0
text(0.45,-0.02,'Longitude(°)','FontSize',8,'FontWeight','normal'); % add X label to entire figure,[left,bottom] distance to the entire axis point 0

%ha = tight_subplot(3,2,[.06 .04],[.12 .04],[.06 .01]); % outside funtion:tight_subplot. limit the gaps and margins
ha = tight_subplot(2,2,[.06 .06],[.12 .04],[.06 .01]); % outside funtion:tight_subplot. limit the gaps and margins

            % 0.Reproduce fig 2.12 in NCA
if plot_annual_change_Livneh==1
    prec_baseline = mean(prec_annual_NE(:,:,1:46),3); % 1915-1960 average
    prec_current = mean(prec_annual_NE(:,:,77:97),3); % 1991-2011 average
    prec_mean_change = 100.*(prec_current-prec_baseline)./prec_baseline;

    geoshow(prec_mean_change(:,:),prec_refvec,'DisplayType','surface',...
        'ZData',zeros(size(prec_mean_change(:,:))),'CData',prec_mean_change(:,:));
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    colorbar
    title('1991-2011 mean precipitation change (%, relative to 1915-1960 annual mean)'); % Fig 2.12
    set(gca,'FontSize',18)
    ax = gca;
    ax.YLabel.String = 'Latitude';
    ax.XLabel.String = 'Longitude';

    hold on
    %figure

    prec_baseline_NE = mean(nanmean(prec_baseline)); % average all grids ignoring NaN
    prec_decade_mean_NE = [];
    prec_decade_mean_NE(1) = mean(nanmean(mean(prec_annual_NE(:,:,1:6),3))); % 1920 is 6
    prec_decade_mean_NE(10) = mean(nanmean(mean(prec_annual_NE(:,:,87:97),3)));% 2001 is 87
    for i=2:9
        prec_decade_mean_NE(i) = mean(nanmean(mean(prec_annual_NE(:,:,i*10-13:i*10-4),3)));%if decade 0-9,should be -14, -5
    end
    prec_decade_extr_change_NE = 100.*(prec_decade_mean_NE-prec_baseline_NE)./prec_baseline_NE;
    bar(1910:10:2000,prec_decade_extr_change_NE);
    set(gca,'FontSize',18);
    title('1991-2011 precipitation decadal change (%, relative to 1915-1960 annual mean)');
    xlabel('Decade');
    ylabel('Change (%)');

    hold on
    figure
end

               % 0'.Reproduce fig 2.17 in NCA
if plot_extr_change_Livneh==1
   % sort the sequence and select the 1% heaviest prec each decade
   prec_year_extr_NE = [];
   prec_year_extr_cell = [];
   for j=1:97
       prec_year_extr_temp =[];
       yr_daynum_start = day_mon((j-1)*12+1) - day_mon(1) + 1; % start date number
       k = day_mon(j*12+1) - day_mon((j-1)*12+1);
       yr_daynum_end = yr_daynum_start + k -1;   % end date number
       prec_year_extr_temp = sort(prec_daily_NE(:,:,yr_daynum_start:yr_daynum_end),3,'descend');
       prec_year_extr_NE(:,:,j) = mean(prec_year_extr_temp(:,:,1:4),3);
   end
   prec_extr_mean = squeeze(mean(prec_year_extr_NE,3));

   % calculate baseline value, get decadal changes
   prec_baseline_extr_cell = mean(prec_year_extr_NE(:,:,1:46),3); % 1915-1960 annual mean value for each cell
   prec_current_extr_cell = mean(prec_year_extr_NE(:,:,77:97),3); % 1991-2011 cell values
   prec_extr_change_scale = 100.*(prec_current_extr_cell-prec_baseline_extr_cell)./prec_baseline_extr_cell;

   geoshow(prec_extr_change_scale(:,:),prec_refvec,'DisplayType','surface',...
       'ZData',zeros(size(prec_extr_change_scale(:,:))),'CData',prec_extr_change_scale(:,:));
   geoshow(lat_NE_map,lon_NE_map,'Color','k');
   colorbar
   title('1991-2011 extreme precipitation change (%, relative to 1915-1960 annual mean)');
   set(gca,'FontSize',18);
   ylabel('Latitude');
   xlabel('Longitude');

   hold on
   figure

   % plot 2.17 bar chart
   prec_baseline_extr_NE = mean(prec_year_extr_NE(1:46)); % 1915-1960 annual mean value for whole cells
   prec_decade_extr_NE = [];
   prec_decade_extr_NE(1) = mean(prec_year_extr_NE(1:6)); % 1915-1919 should be 1:5
   prec_decade_extr_NE(10) = mean(prec_year_extr_NE(87:97)); % 2000-2011 should be 86:97
   for m=2:9
       prec_decade_extr_NE(m) = mean(prec_year_extr_NE(10*(m-1)-3:10*(m-1)+6)); %if decade 0-9,should be -4,+5
   end

   prec_decade_extr_change_NE = 100.*(prec_decade_extr_NE-prec_baseline_extr_NE)./prec_baseline_extr_NE;
   bar(1910:10:2000,prec_decade_extr_change_NE);
   title('1991-2011 extreme precipitation decadal change (%, relative to 1915-1960 annual mean)');
   set(gca,'FontSize',18);
   xlabel('Decade');
   ylabel('Change (%)');

   hold on
   figure
end


              % 1.Livneh+GHCND:plot 1915-2011 annual mean prec
if plot_annual_mean_GL==1
  % calculate NE regional trend
%   10.*polyfit(linspace(1915,2011,97)',squeeze(squeeze(nanmean(nanmean(prec_annual_NE,2)))),1)
%   [H,p_value]=Mann_Kendall(squeeze(squeeze(nanmean(nanmean(prec_annual_NE,2)))),0.05)
    figure 
    B1 = linspace(1915,2011,97)';
    C1 = squeeze(nanmean(nanmean(prec_annual_NE(:,:,:),2)));
    ipt = findchangepts(C1,'MaxNumChanges',1)
    findchangepts(C1,'MaxNumChanges',1)
    figure
    plot(B1,C1,'Color',[0.15 0.15 0.15],'Marker','o','MarkerSize',3);
    stop
    % calculate regional trend
    %   M-K test and Thei-sen slope
    datain = [linspace(1915,2011,97)',squeeze(nanmean(nanmean(prec_annual_NE(:,:,1:97),2)))];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
    slope_MK = 10*sen
    p_value_MK = sig
    H_MK = h
%     %   linear trend estimate
%     disp('    Hereafter is linear trend estimate');
%     mdl = fitlm(linspace(1979,2011,33)',squeeze(nanmean(nanmean(prec_annual_NE(:,:,65:97),2)))); % provide linear trend and p-value
%     figure
%     plotResiduals(mdl)
%     [H1,p] = lillietest(mdl.Residuals.Raw) % H1=0 means a normal distribution at default alpha level 0.05
%     h1 = lbqtest(mdl.Residuals.Raw) % h1=0 means residuals are independent/non-autocorrelated
%     slope = 10*table2array(mdl.Coefficients(2,1)) % convert to mm/decade
%     p_value = table2array(mdl.Coefficients(2,4)) % for linear regression model
% 
%     STOP
    
  prec_annual_NE_scale = mean(prec_annual_NE,3);
  axes(ha(1));
  geoshow(prec_annual_NE_scale,prec_refvec,'DisplayType','surface',...
      'ZData',zeros(size(prec_annual_NE_scale)),'CData',prec_annual_NE_scale);
  shading flat;
  hold on

  % Plot GHCN-D 1915-2011 annual mean daily PRCP (mm/day)
  prcp_1511 = nanmean(station_data_yr(:,:),1);
  scatter(geo(2,:),geo(1,:),13,0.1*prcp_1511(1,:),'filled','MarkerEdgeColor','k','LineWidth',0.1);
  geoshow(lat_NE_map,lon_NE_map,'Color','k');
  
  ax = gca; % get the handle for the current axes
  title(ax,'a) LI2013 & GHCN-D 1915-2011 mean','FontSize',8,'FontWeight','normal');
  set(ax,'FontSize',8);
  xlim([-84 -66]);
  ylim([36 48]);
  set(gca,'Xtick',linspace(-84,-66,7));
  set(gca,'Ytick',linspace(36,48,5));
  ax.CLim=[750 1450];
  colormap(ax,flipud(jet)); % change default color to jet:blue-->red
  
  if(sav==1)
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
        saveas(gcf,'Figures/annual_mean_GL.jpg'); % need to manually save .jpg by replacing with same name
  end
  %figure
  hold on
end

            % 1'.Livneh+GHCND:calculate and plot trend in 1915-2011 annual mean prec
if plot_annual_trend_GL==1
    prec_extr_regress_coeff_temp = [];
    prec_extr_regress_coeff = [];
    for x=1:size(prec_annual_NE,1)
        x
        for y=1:size(prec_annual_NE,2)
            if prec_annual_NE(x,y,size(prec_annual_NE,3)) == NaN
                prec_extr_regress_coeff(x,y) = NaN;
            else
                prec_extr_regress_coeff_temp = polyfit([1:1:size(prec_annual_NE,3)],squeeze(prec_annual_NE(x,y,:))',1);
                prec_extr_regress_coeff(x,y) = prec_extr_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
            end
        end
    end
    prec_decade_trend(:,:) = 10.*prec_extr_regress_coeff;  % changing yearly trends to be mm/decade

    % plot trend of GHCND, analyze 1915-2011 trend only
    for x=1:176
        ind = 1:length(station_data_yr(:,x)); % analyze 1979-2014 trend only
        k = ~isnan(station_data_yr(:,x)); % find index with non-NaN values
        %pref_regress_coeff_temp = polyfit(ind(k)',36.5246.*station_data_yr(k,x+4),1);
        y = find(k==1);
        prcp_regress_coeff_temp_1511 = polyfit(ind(k)',0.1.*station_data_yr(y,x),1);
        prcp_regress_coeff_1511(x) = prcp_regress_coeff_temp_1511(1); % polyfit produces two parameters, here we need coffieient, not intercept
        
        mdl = fitlm(ind(k)',0.1*station_data_yr(y,x)); % provide linear trend and p-value
        prcp_regress_coeff(x,1) = 10*table2array(mdl.Coefficients(2,1)) % convert to mm/decade
        prcp_regress_coeff(x,2) = table2array(mdl.Coefficients(2,4)) % p-value, for linear regression model       
    end
%     prcp_decade_trend_1511 = 10.*prcp_regress_coeff_1511;  % changing yearly trends to be mm/decade
    
    axes(ha(2));
    geoshow(prec_decade_trend(:,:),prec_refvec,'DisplayType','surface',...
        'ZData',zeros(size(prec_decade_trend(:,:))),'CData',prec_decade_trend(:,:));
    shading flat;
    hold on
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    
        % plot stations with significant and insignificant trends in different markers     
    pp=0;   
    qq=0;
    for x=1:176
        if prcp_regress_coeff(x,2)<0.05 % if significant trend
            pp=pp+1;
            prcp_regress_coeff1(pp)=prcp_regress_coeff(x,1);
            geo1(2,pp)=geo(2,x);
            geo1(1,pp)=geo(1,x);
        else                            % if insignificant trend
            qq=qq+1;
            prcp_regress_coeff2(qq)=prcp_regress_coeff(x,1);
            geo2(2,qq)=geo(2,x);
            geo2(1,qq)=geo(1,x);
        end
    end
    
    scatter(geo1(2,:),geo1(1,:),17,prcp_regress_coeff1,'square','filled','MarkerEdgeColor','k','LineWidth',0.1);
    hold on;
    scatter(geo2(2,:),geo2(1,:),9,prcp_regress_coeff2,'diamond','filled','MarkerEdgeColor','k','LineWidth',0.1);

%     scatter(geo(2,:),geo(1,:),13,prcp_decade_trend_1511,'filled','MarkerEdgeColor','k','LineWidth',0.1);
    ax = gca; % get the handle for the current axes
    title(ax,'b) LI2013 & GHCN-D 1915-2011 trend','FontSize',8,'FontWeight','normal');
    set(ax,'FontSize',8);
    xlim([-84 -66]);
    ylim([36 48]);
    set(gca,'Xtick',linspace(-84,-66,7));
    set(gca,'Ytick',linspace(36,48,5));
    ax.CLim=[-35 85];
    
    blueColorMap = [ones(1,75), linspace(1, 0, 181)]; % [0,0,1] refers to blue; 116 and 140 is calculated from colorbar limits [-49,59], to make 0 be white
    redColorMap = [linspace(0, 1, 75), linspace(1, 0, 181)];% [1,0,0] refers to red
    whiteColorMap = [linspace(0, 1, 75), ones(1,181)]; % [1,1,1] refers to white
    colorMap = [blueColorMap; redColorMap; whiteColorMap]'; % creat a 3-D matrix to install RGB values
    colormap(ax,colorMap);
    hold on;
    
    clear station_data_yr prcp_regress_coeff geo1 geo2 prcp_regress_coeff1 prcp_regress_coeff2
    
    % saveas(gcf,'Figures/test2.epsc')
    if(sav==1)
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
        saveas(gcf,'Figures/annual_mean_trend_GL.jpg'); % need to manually save .jpg by replacing with same name
    end
end

                % 2.Livneh+GHCND:plot 1915-2011 mean extreme prec
if plot_extr_mean_GL==1
    % find 1% threshold all daily events,sort the sequence 
    load '/adata/Data/Observations/GHCN-D/station_year_extr_NE176.mat'; % 97 yr*176 st
% if prec_year_extr_NE.mat saved out, comment following paragraph
    load('/ahome/hhuang/Matlab/NE_PRECIP/Proc_Data/prec_year_extr_NE.mat');    
%     prec_year_extr_NE = []; 
%     threshold_L = [];
%     for i=1:size(prec_daily_NE,1)
%         for j=1:size(prec_daily_NE,2)
%             threshold_temp_L = squeeze(sort(prec_daily_NE(i,j,:),'descend'));
%             if isnan(threshold_temp_L(1))
%                 threshold_L(i,j) = NaN;
%             else
%                 count = 0;
%                 for k=1:size(threshold_temp_L)
%                     if threshold_temp_L(k)<0.254 % <0.01 inch is trace precip
%                         break
%                     end
%                     count = count+1;
%                 end
%                 r = fix(0.01*count);
%                 threshold_L(i,j) = threshold_temp_L(r+1); % find 1%+1 heaviest prec value
%             end
%         end
%     end
% 
%     for i=1:size(prec_daily_NE,1)
%         i
%         for j=1:size(prec_daily_NE,2)                                                         
%             for k=1:97
%                 prec_year_extr_temp = [];
%                 extr_temp_L = 0;
%                 yr_daynum_start = day_mon((k-1)*12+1) - day_mon(1) + 1; % start date number
%                 num = day_mon(k*12+1) - day_mon((k-1)*12+1);
%                 yr_daynum_end = yr_daynum_start + num -1;   % end date number
%                 prec_year_extr_temp = squeeze(sort(prec_daily_NE(i,j,yr_daynum_start:yr_daynum_end),3,'descend'));
%                 if isnan(prec_year_extr_temp(1))
%                     prec_year_extr_NE(i,j,k) = NaN;
%                 else
%                     for d=1:size(prec_year_extr_temp) % 365/366,calculate by year by station
%                         if prec_year_extr_temp(d)>threshold_L(i,j)
%                             extr_temp_L = extr_temp_L + prec_year_extr_temp(d); % larger than threshold,sum up
%                         end
%                     end
%                     prec_year_extr_NE(i,j,k) = extr_temp_L;
%                 end
%             end
%         end
%     end    
%     STOP

    % if define the wettest 1% events(4 days) each year as extr prec, use following codes
%     for j=1:97
%         prec_year_extr_temp =[];
%         yr_daynum_start = day_mon((j-1)*12+1) - day_mon(1) + 1; % start date number
%         k = day_mon(j*12+1) - day_mon((j-1)*12+1);
%         yr_daynum_end = yr_daynum_start + k -1;   % end date number
%         prec_year_extr_temp = sort(prec_daily_NE(:,:,yr_daynum_start:yr_daynum_end),3,'descend');
%         prec_year_extr_NE(:,:,j) = mean(prec_year_extr_temp(:,:,1:4),3);
%     end

%     figure 
%     B1 = linspace(1915,2011,97)';
%     C1 = squeeze(nanmean(nanmean(prec_year_extr_NE(:,:,:),2)));
%     ipt = findchangepts(C1,'MaxNumChanges',1)
%     findchangepts(C1,'MaxNumChanges',1)
%     figure
%     plot(B1,C1,'Color',[0.15 0.15 0.15],'Marker','o','MarkerSize',3);
%     stop
    
    prec_extr_mean = squeeze(mean(prec_year_extr_NE,3));
    
    axes(ha(1));
    geoshow(prec_extr_mean,prec_refvec,'DisplayType','surface',...
        'ZData',zeros(size(prec_extr_mean)),'CData',prec_extr_mean);
    shading flat;
    hold on
    scatter(geo(2,:),geo(1,:),13,0.1.*nanmean(station_year_extr_NE(:,:)),'filled','MarkerEdgeColor','k','LineWidth',0.1);
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    ax = gca; % get the handle for the current axes
    title(ax,'a) LI2013 & GHCN-D 1915-2011 mean','FontSize',8,'FontWeight','normal');
    set(ax,'FontSize',8);
    xlim([-84 -66]);
    ylim([36 48]);
    set(gca,'Xtick',linspace(-84,-66,7));
    set(gca,'Ytick',linspace(36,48,5));
    ax.CLim=[45 125];
    colormap(ax,flipud(jet)); % change default color to jet:blue-->red
    
    if(sav==1)
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
        saveas(gcf,'Figures/extr_mean_GL.jpg'); % need to manually save .jpg by replacing with same name
    end

end


           % 2'.Livneh+GHCND:plot trend in 1915-2011 mean extreme prec
if plot_extr_trend_GL==1
%     % calculate regional trend
%     %   M-K test and Thei-sen slope
%     datain = [linspace(1979,2011,33)',squeeze(nanmean(nanmean(prec_year_extr_NE(:,:,65:97),2)))]; 
%     [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
%     slope_MK = 10*sen
%     p_value_MK = sig
%     H_MK = h
%     %   linear trend estimate
%     disp('    Hereafter is linear trend estimate');
%     mdl = fitlm(linspace(1979,2011,33)',squeeze(prec_year_extr_NE(41,162,65:97))); % provide linear trend and p-value
%     figure
%     plotResiduals(mdl)
%     [H1,p] = lillietest(mdl.Residuals.Raw) % H1=0 means a normal distribution at default alpha level 0.05
%     h1 = lbqtest(mdl.Residuals.Raw) % h1=0 means residuals are independent/non-autocorrelated
%     slope = 10*table2array(mdl.Coefficients(2,1)) % convert to mm/decade
%     p_value = table2array(mdl.Coefficients(2,4)) % for linear regression model
% 
%     STOP
    %%
    % calculate grid trend
    prec_extr_regress_coeff_temp = [];
    prec_extr_regress_coeff = [];
    for x=1:size(prec_year_extr_NE,1)
        x
        for y=1:size(prec_year_extr_NE,2)
            if isnan(prec_year_extr_NE(x,y,size(prec_year_extr_NE,3)))
                prec_extr_regress_coeff(x,y) = NaN;
            else
                datain = [];
                datain = [linspace(1915,2011,97)',squeeze(prec_year_extr_NE(x,y,65:97))];
                [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
                prec_extr_regress_coeff(x,y) = sen;
%                 prec_regress_coeff_temp = polyfit([1:1:size(prec_year_extr_NE(:,:,65:97),3)],squeeze(prec_year_extr_NE(x,y,65:97))',1);
%                 prec_regress_coeff(x,y) = prec_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
            end
        end
    end
    prec_extr_decade_trend(:,:) = 10.*prec_extr_regress_coeff;  % changing yearly trends to be mm/decade
    %%
    % calculate trend in GHCND 1915-2011 annual extreme precipitation
    prcp_extr_regress_coeff_temp = [];
    prcp_extr_regress_coeff = [];
    for x=1:176
        ind = 1:length(station_year_extr_NE(:,x));
        k = ~isnan(station_year_extr_NE(:,x)); % find index with non-NaN values
        y = find(k==1);
        datain = [ind(k)',0.1*station_year_extr_NE(y,x)]; % ind(k)'=ind(y)',y is used to find non-NaN extr
        [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
        prcp_extr_regress_coeff(x,1) = 10*sen; % Thei-sen slope
        prcp_extr_regress_coeff(x,2) = sig;
%         prcp_extr_regress_coeff_temp = polyfit(ind(k)',0.1.*station_year_extr_NE(14+y,x),1);
%         prcp_extr_regress_coeff(x) = prcp_extr_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
    end
%     prcp_extr_decade_trend = 10.*prcp_extr_regress_coeff;  % changing yearly trends to be mm/decade
    
    axes(ha(2));
    figure(2), clf
    geoshow(prec_extr_decade_trend(:,:),prec_refvec,'DisplayType','surface',...
        'ZData',zeros(size(prec_extr_decade_trend(:,:))),'CData',prec_extr_decade_trend(:,:));
    shading flat;
    hold on
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
%     scatter(geo(2,:),geo(1,:),13,prcp_extr_decade_trend,'filled','MarkerEdgeColor','k','LineWidth',0.1);
        % plot stations with significant and insignificant trends in different markers     
    pp=0;   
    qq=0;
    for x=1:176
        if prcp_extr_regress_coeff(x,2)<0.05 % if significant trend
            pp=pp+1;
            prcp_extr_regress_coeff1(pp)=prcp_extr_regress_coeff(x,1);
            geo1(2,pp)=geo(2,x);
            geo1(1,pp)=geo(1,x);
        else                            % if insignificant trend
            qq=qq+1;
            prcp_extr_regress_coeff2(qq)=prcp_extr_regress_coeff(x,1);
            geo2(2,qq)=geo(2,x);
            geo2(1,qq)=geo(1,x);
        end
    end
       
    scatter(geo1(2,:),geo1(1,:),17,prcp_extr_regress_coeff1,'square','filled','MarkerEdgeColor','k','LineWidth',0.1);
    hold on;
    scatter(geo2(2,:),geo2(1,:),9,prcp_extr_regress_coeff2,'diamond','filled','MarkerEdgeColor','k','LineWidth',0.1);

    ax = gca; % get the handle for the current axes
    title(ax,'b) LI2013 & GHCN-D 1915-2011 trend','FontSize',8,'FontWeight','normal');
    set(ax,'FontSize',8);
    xlim([-84 -66]);
    ylim([36 48]);
    set(gca,'Xtick',linspace(-84,-66,7));
    set(gca,'Ytick',linspace(36,48,5));
    ax.CLim=[-19 29];
    blueColorMap = [ones(1,101), linspace(1, 0, 155)]; % [0,0,1] refers to blue; 116 and 140 is calculated from colorbar limits [-49,59], to make 0 be white
    redColorMap = [linspace(0, 1, 101), linspace(1, 0, 155)];% [1,0,0] refers to red
    whiteColorMap = [linspace(0, 1, 101), ones(1,155)]; % [1,1,1] refers to white
    colorMap = [blueColorMap; redColorMap; whiteColorMap]'; % creat a 3-D matrix to install RGB values
    colormap(ax,colorMap);
    hold on;
    
    clear station_year_extr_NE geo1 geo2 prcp_extr_regress_coeff1 prcp_extr_regress_coeff2
    if(sav==1)
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
        saveas(gcf,'Figures/extr_mean_trend_GL.jpg'); % need to manually save .jpg by replacing with same name
    end
end

              % 3.plot bar charts for decadel changes in mean prec
if plot_mean_decadel_change_GLN==1
    mean_baseline_NE = 36.525*nanmean(nanmean(station_data_yr(1:60,5:150))); % GHCN-D 1901-1960 annual mean prec as baseline
    % compute GHCN-D annual mean prec in every decade 1901-2010,2011-2014
    for i=1:11
        prcp_decade_mean_NE(i) = 36.525.*nanmean(nanmean(station_data_yr((i-1)*10+1:i*10,5:150)));
    end
    prcp_decade_mean_NE(12) = 36.525.*nanmean(nanmean(station_data_yr(111:114,5:150))); % 2011-2014
    
    % compute Livneh annual mean prec in every decade 1921-2010,9 decades
    prec_decade_mean_NE = zeros(1,12);
    for i=3:11
        prec_decade_mean_NE(i) = nanmean(nanmean(mean(prec_annual_NE(:,:,(i-2)*10-3:(i-2)*10+6),3)));
    end
    
    % compute NARR annual mean prec in every decade 1981-2010,2011-2014
    load('/ahome/hhuang/Matlab/NE_PRECIP/Proc_Data/NARR_Prec_NE.mat');
    yr_start_L = 1979;
    yr_end_L = 2014;
    nmon_L = (yr_end_L-yr_start_L+1)*12;
    day_mon_L = datenum(yr_start_L*ones(nmon_L+1,1),[1:nmon_L+1]',1);
    tday_L = [day_mon_L(1):day_mon_L(1)+length(apcp_daily_NE)-1]';
    for mon = 1:nmon_L
        indx_N = find(tday_L>=day_mon_L(mon) & tday_L<day_mon_L(mon+1));
        apcp_mon_NE(:,:,mon) = sum(apcp_daily_NE(:,:,indx_N),3); 
    end
    % reshape to a 4th dimension that's 12 months long 45 x 69 x 12 x 36
    apcp_annual_temp = reshape(apcp_mon_NE,[165,250,12,36]);
    apcp_annual_NE = squeeze(sum(apcp_annual_temp,3)); 
    apcp_decade_mean_NE = zeros(1,12);
    for i=9:11
        apcp_decade_mean_NE(i) = nanmean(nanmean(nanmean(apcp_annual_NE(:,:,(i-8)*10-7:(i-8)*10+2),3)));
    end
    apcp_decade_mean_NE(12) = nanmean(nanmean(nanmean(apcp_annual_NE(:,:,33:36),3))); % 2011-2014

    % calculate decadel changes in mean prec
    prcp_decade_mean_change_NE = 100.*(prcp_decade_mean_NE-...
        mean_baseline_NE)./mean_baseline_NE;
    
    prec_decade_mean_change_NE = 100.*(prec_decade_mean_NE-...
        mean_baseline_NE)./mean_baseline_NE;
    prec_decade_mean_change_NE(1) = 0;
    prec_decade_mean_change_NE(2) = 0;
    prec_decade_mean_change_NE(12) = 0;
    
    apcp_decade_mean_change_NE = 100.*(apcp_decade_mean_NE-...
        mean_baseline_NE)./mean_baseline_NE;
    for i=1:8
        apcp_decade_mean_change_NE(i)=0;
    end
    
    decade_mean_change_NE_combined = [prcp_decade_mean_change_NE(:),...
        prec_decade_mean_change_NE(:),apcp_decade_mean_change_NE(:)];
    bar([1900:10:2010],decade_mean_change_NE_combined,'grouped');
    legend('GHCN-D','Livneh','NARR');
    set(gca,'FontSize',28);
    title('1901-2014 annual mean precipitation decadal change (%, relative to GHCN-D 1901-1960 annual mean)');
    xlabel('Decade');
    ylabel('Change (%)');
    
    hold on
    figure
end

              % 3'.plot bar charts for decadel changes in extreme prec
if plot_extr_decadel_change_GLN==1
    % compute GHCN-D extreme prec in every decade 1901-2010,2011-2014
    load '/adata/Data/Observations/GHCN-D/station_year_extr_NE.mat'; % 114 yr*176 st
    extr_baseline_NE = 0.1*nanmean(nanmean(station_year_extr_NE(1:60,:)));% GHCN-D 1901-1960 annual mean prec as baseline
    for i=1:11
        prcp_extr_decade_mean_NE(i) = 0.1*nanmean(nanmean(station_year_extr_NE((i-1)*10+1:i*10,:)));
    end
    prcp_extr_decade_mean_NE(12) = 0.1*nanmean(nanmean(station_year_extr_NE(111:114,:))); % 2011-2014
    
    % compute Livneh extr prec in every decade 1921-2010,9 decades
    load('/ahome/hhuang/Matlab/NE_PRECIP/Proc_Data/prec_year_extr_NE.mat'); 
%     prec_year_extr_NE = [];
%     for j=1:97
%         prec_year_extr_temp =[];
%         yr_daynum_start = day_mon((j-1)*12+1) - day_mon(1) + 1; % start date number
%         k = day_mon(j*12+1) - day_mon((j-1)*12+1);
%         yr_daynum_end = yr_daynum_start + k -1;   % end date number
%         prec_year_extr_temp = sort(prec_daily_NE(:,:,yr_daynum_start:yr_daynum_end),3,'descend');
%         prec_year_extr_NE(:,:,j) = mean(prec_year_extr_temp(:,:,1:4),3);
%     end   
    prec_extr_decade_mean_NE = zeros(1,12);
    for i=3:11
        prec_extr_decade_mean_NE(i) = nanmean(nanmean(mean(prec_year_extr_NE(:,:,(i-2)*10-3:(i-2)*10+6),3)));
    end  
    
    % compute NARR extreme prec in every decade 1981-2010,2011-2014
    load('/ahome/hhuang/Matlab/NE_PRECIP/Proc_Data/apcp_year_extr_NE.mat');
%     yr_start_N = 1979;
%     yr_end_N = 2014;
%     nmon_N = (yr_end_N-yr_start_N+1)*12;
%     day_mon_N = datenum(yr_start_N*ones(nmon_N+1,1),[1:nmon_N+1]',1);
%     apcp_year_extr_NE = [];
%     for j=1:36
%         apcp_year_extr_temp =[];
%         yr_daynum_start_N = day_mon_N((j-1)*12+1) - day_mon_N(1) + 1; % start date number
%         k_N = day_mon_N(j*12+1) - day_mon_N((j-1)*12+1);
%         yr_daynum_end_N = yr_daynum_start_N + k_N -1;   % end date number
%         apcp_year_extr_temp = sort(apcp_daily_NE(:,:,yr_daynum_start_N:yr_daynum_end_N),3,'descend');
%         apcp_year_extr_NE(:,:,j) = mean(apcp_year_extr_temp(:,:,1:4),3); % 45*69*36 matrix
%     end
    apcp_extr_decade_mean_NE = zeros(1,12);
    for i=9:11
        apcp_extr_decade_mean_NE(i) = nanmean(nanmean(nanmean(apcp_year_extr_NE(:,:,(i-8)*10-7:(i-8)*10+2),3)));
    end
    apcp_extr_decade_mean_NE(12) = nanmean(nanmean(nanmean(apcp_year_extr_NE(:,:,33:36),3))); % 2011-2014

    % calculate decadel changes in extreme prec
    prcp_extr_decade_mean_change_NE = 100.*(prcp_extr_decade_mean_NE-...
        extr_baseline_NE)./extr_baseline_NE;
    
    prec_extr_decade_mean_change_NE = 100.*(prec_extr_decade_mean_NE-...
        extr_baseline_NE)./extr_baseline_NE;
    prec_extr_decade_mean_change_NE(1) = 0;
    prec_extr_decade_mean_change_NE(2) = 0;
    prec_extr_decade_mean_change_NE(12) = 0;
    
    apcp_extr_decade_mean_change_NE = 100.*(apcp_extr_decade_mean_NE-...
        extr_baseline_NE)./extr_baseline_NE;
    for i=1:8
        apcp_extr_decade_mean_change_NE(i)=0;
    end
    
    decade_extr_change_NE_combined = [prcp_extr_decade_mean_change_NE(:),...
        prec_extr_decade_mean_change_NE(:),apcp_extr_decade_mean_change_NE(:)];
    bar([1900:10:2010],decade_extr_change_NE_combined,'grouped');
    legend('GHCN-D','Livneh','NARR');
    set(gca,'FontSize',28);
    title('1901-2014 extreme precipitation decadal change (%, relative to GHCN-D 1901-1960 annual mean)');
    xlabel('Decade');
    ylabel('Change (%)');
end   
    

            % 4. Livneh & GHCN-D:plot 1915-2011 seasonal total prec
if plot_season_mean_GL == 1
    datenum_start = day_mon(1);
    for i=1:97
        i
        for j=1:size(prec_daily_NE,1)
            for k=1:size(prec_daily_NE,2)
                spring_start = day_mon(i*12-9)-datenum_start;
                spring_end = day_mon(i*12-6)-datenum_start-1;
                prec_year_spring(j,k,i) = sum(prec_daily_NE(j,k,spring_start:spring_end),3);

                summer_start = day_mon(i*12-6)-datenum_start;
                summer_end = day_mon(i*12-3)-datenum_start-1;
                prec_year_summer(j,k,i) = sum(prec_daily_NE(j,k,summer_start:summer_end),3);

                autumn_start = day_mon(i*12-3)-datenum_start;
                autumn_end = day_mon(i*12)-datenum_start-1;
                prec_year_autumn(j,k,i) = sum(prec_daily_NE(j,k,autumn_start:autumn_end),3);

                if i==97
                    prec_year_winter(j,k,i) = NaN;
                else
                    winter_start = day_mon(i*12)-datenum_start;
                    winter_end = day_mon(i*12+3)-datenum_start-1;
                    prec_year_winter(j,k,i) = sum(prec_daily_NE(j,k,winter_start:winter_end),3);
                end
            end
        end
    end

%    % calculate NE regional trend
%    10.*polyfit(linspace(1915,2011,97)',squeeze(squeeze(nanmean(nanmean(prec_year_spring,2)))),1)
%    [H,p_value]=Mann_Kendall(squeeze(squeeze(nanmean(nanmean(prec_year_spring,2)))),0.05)

   % calculate regional trend
    %   M-K test and Thei-sen slope
    datain = [linspace(1915,2011,97)',squeeze(nanmean(nanmean(prec_year_spring(:,:,1:97),2)))];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
    slope_MK = 10*sen
    p_value_MK = sig
    stop
    H_MK = h
    %   linear trend estimate
    disp('    Hereafter is linear trend estimate');
    mdl = fitlm(linspace(1979,2011,33)',squeeze(nanmean(nanmean(prec_year_spring(:,:,65:97),2)))) % provide linear trend and p-value
    figure
    plotResiduals(mdl)
    [H1,p] = lillietest(mdl.Residuals.Raw) % H1=0 means a normal distribution at default alpha level 0.05
    h1 = lbqtest(mdl.Residuals.Raw) % h1=0 means residuals are independent/non-autocorrelated
    slope = 10*table2array(mdl.Coefficients(2,1)) % convert to mm/decade
    p_value = table2array(mdl.Coefficients(2,4)) % for linear regression model

    STOP
    
    load '/adata/Data/Observations/GHCN-D/station_year_spring.mat';  % 114*176
    load '/adata/Data/Observations/GHCN-D/station_year_summer.mat';  % 114*176
    load '/adata/Data/Observations/GHCN-D/station_year_autumn.mat';  % 114*176
    load '/adata/Data/Observations/GHCN-D/station_year_winter.mat';  % 114*176
    prec_summer_mean = squeeze(nanmean(prec_year_summer,3));

%     geoshow(prec_summer_mean,prec_refvec,'DisplayType','surface',...
%         'ZData',zeros(size(prec_summer_mean)),'CData',prec_summer_mean);
%     shading flat;
%     hold on
%     scatter(geo(2,:),geo(1,:),100,0.1.*nanmean(station_year_summer(15:111,:)),'filled','MarkerEdgeColor','k');
%     geoshow(lat_NE_map,lon_NE_map,'Color','k');
%     colorbar;
%     title('Livneh & GHCN-D 1915-2011 annual mean summer precipitation (mm)');
%     ax = gca; % get the handle for the current axes
%     set(ax,'FontSize',28);
%     % ax.CLim=[-49 59]; % set the min and max limits for the colorbar
%     colormap(flipud(jet)); % change default color to jet:blue-->red
end

             % 4'. Plot Livneh+GHCN-D trend in 1915-2011 seasonal prec
if plot_season_trend_GL==1    
    % calculate linear trend in Livneh seasonal precipitation
    prec_season_regress_coeff_temp = [];
    prec_season_regress_coeff = [];
    for x=1:size(prec_year_spring,1)
        for y=1:size(prec_year_spring,2)
            if isnan(prec_year_spring(x,y,1))
                prec_season_regress_coeff(x,y) = NaN;
            else
                prec_season_regress_coeff_temp = polyfit([1:1:size(prec_year_spring,3)],squeeze(prec_year_spring(x,y,:))',1);
                prec_season_regress_coeff(x,y) = prec_season_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
            end
        end
    end 
    prec_season_decade_trend(:,:) = 10.*prec_season_regress_coeff;  % changing yearly trends to be mm/decade

    % calculate trend in GHCND annual mean precipitation,analyze 1979-2014 trend only
    for x=1:176
        ind = 1:length(station_year_spring(15:111,x)); % analyze 1979-2014 trend only
        k = ~isnan(station_year_spring(15:111,x)); % find index with non-NaN values
        y = find(k==1);
        prcp_season_regress_coeff_temp = polyfit(ind(k)',0.1*station_year_spring(14+y,x),1);
        prcp_season_regress_coeff(x) = prcp_season_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
    end
    prcp_season_decade_trend = 10.*prcp_season_regress_coeff;  % changing yearly trends to be mm/decade
    
    % plot trend of Livneh, and then GHCND
    geoshow(prec_season_decade_trend(:,:),prec_refvec,'DisplayType','surface',...
        'ZData',zeros(size(prec_season_decade_trend(:,:))),'CData',prec_season_decade_trend(:,:));
    shading flat;
    hold on
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    scatter(geo(2,:),geo(1,:),100,prcp_season_decade_trend,'filled','MarkerEdgeColor','k');
    colorbar
    title('Livneh & GHCND trends in 1915-2011 Spring Precipitation (mm/decade)');
    ax = gca; % get the handle for the current axes
    set(ax,'FontSize',28);
    %ax.CLim=[-49 59]; % set the min and max limits for the colorbar
    colormap(flipud(jet)); % change default color to jet:blue-->red
    
%     if(sav==1)
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
%         saveas(gcf,'Figures/annual_mean_trend_GN.jpg'); % need to manually save .jpg by replacing with same name
%     end
end

                % 5.GHCN-D:plot 1901-2014 seasonal extreme prec
if plot_season_extr_mean_GL == 1
    season = 'summer'; 
    yr_start_L = 1915;
    yr_end_L = 2011;
    nmon_L = (yr_end_L-yr_start_L+1)*12;
    day_mon_L = datenum(yr_start_L*ones(nmon_L+1,1),[1:nmon_L+1]',1);
    %tday_L = [day_mon_L(1):day_mon_L(1)+length(apcp_daily_NE)-1]';
    datenum_start_L = 699440;
    
    yr_start_G = 1901; % GHCN-D
    yr_end_G = 2014;
    nyr_G = yr_end_G - yr_start_G +1;
    nmon_G = (yr_end_G-yr_start_G+1)*12;
    datenum_start_G = 694327;
    day_mon_G = datenum(yr_start_G*ones(nmon_G+1,1),[1:nmon_G+1]',1);
    
    % extract daily data by season for every cell
    for i=1:size(prec_daily_NE,1)
        %i
        for j=1:size(prec_daily_NE,2)
            prec_season_data_temp = [];
            prec_season_data = [];% store seasonal daily data
            if isnan(prec_daily_NE(i,j,1))
                threshold_L(i,j) = NaN;
            else
                for k=1:97 % if winter, use 96
                    if strcmp(season,'spring')
                        season_start_L = day_mon_L(k*12-9)-datenum_start_L;
                        season_end_L = day_mon_L(k*12-6)-datenum_start_L-1;
                    elseif strcmp(season,'summer')
                        season_start_L = day_mon_L(k*12-6)-datenum_start_L;
                        season_end_L = day_mon_L(k*12-3)-datenum_start_L-1;
                    elseif strcmp(season,'autumn')
                        season_start_L = day_mon_L(k*12-3)-datenum_start_L;
                        season_end_L = day_mon_L(k*12)-datenum_start_L-1;
                    elseif strcmp(season,'winter')
                        if k==97
                            prec_season_extr_NE(i,j,k) = NaN;
                        else
                            season_start_L = day_mon_L(k*12)-datenum_start_L;
                            season_end_L = day_mon_L(k*12+3)-datenum_start_L-1;
                        end
                    end
                    prec_season_data_temp = squeeze(squeeze(prec_daily_NE(i,j,season_start_L:season_end_L)));
                    prec_season_data = [prec_season_data;prec_season_data_temp]; %
                end
                threshold_temp_L = squeeze(sort(prec_season_data,'descend'));
                if isnan(threshold_temp_L(1))
                    threshold_L(i,j) = NaN;
                else
                    count = 0;
                    for k=1:size(threshold_temp_L)
                        if threshold_temp_L(k)<0.254 % <0.01 inch is trace precip
                            break
                        end
                        count = count+1;
                    end
                    r = fix(0.01*count);
                    threshold_L(i,j) = threshold_temp_L(r+1); % find 1%+1 heaviest prec value
                end                
            end            
        end
    end
    clear season_start_L season_end_L
        
%     % calculate NARR extreme threshold and then extract seasonal extreme prec
%     prec_season_extr_NE = [];   
%     threshold_L = [];
%     for i=1:size(prec_daily_NE,1)
%         for j=1:size(prec_daily_NE,2)
%             threshold_temp_L = squeeze(sort(prec_season_data(i,j,:),'descend'));
%             if isnan(threshold_temp_L(1))
%                 threshold_L(i,j) = NaN;
%             else
%                 count = 0;
%                 for k=1:size(threshold_temp_L)
%                     if threshold_temp_L(k)<0.254 % <0.01 inch is trace precip
%                         break
%                     end
%                     count = count+1;
%                 end
%                 r = fix(0.01*count);
%                 threshold_L(i,j) = threshold_temp_L(r+1); % find 1%+1 heaviest prec value
%             end
%         end
%     end

    for i=1:size(prec_daily_NE,1)
        %i
        for j=1:size(prec_daily_NE,2)                                                         
            for k=1:97 % if winter, use 96
                prec_season_extr_temp = [];
                extr_temp_L = 0;
                if strcmp(season,'spring')
                    season_start_L = day_mon_L(k*12-9)-datenum_start_L;
                    season_end_L = day_mon_L(k*12-6)-datenum_start_L-1;
                elseif strcmp(season,'summer')
                    season_start_L = day_mon_L(k*12-6)-datenum_start_L;
                    season_end_L = day_mon_L(k*12-3)-datenum_start_L-1;
                elseif strcmp(season,'autumn')
                    season_start_L = day_mon_L(k*12-3)-datenum_start_L;
                    season_end_L = day_mon_L(k*12)-datenum_start_L-1;
                elseif strcmp(season,'winter')
                    if k==97
                        prec_season_extr_NE(i,j,k) = NaN;
                    else
                        season_start_L = day_mon_L(k*12)-datenum_start_L;
                        season_end_L = day_mon_L(k*12+3)-datenum_start_L-1;
                    end
                end
            
                prec_season_extr_temp = squeeze(sort(prec_daily_NE(i,j,season_start_L:season_end_L),3,'descend'));
                if isnan(prec_season_extr_temp(1))
                    prec_season_extr_NE(i,j,k) = NaN;
                else
                    for d=1:size(prec_season_extr_temp) % 365/366,calculate by year by station
                        if prec_season_extr_temp(d)>threshold_L(i,j)
                            extr_temp_L = extr_temp_L + prec_season_extr_temp(d); % larger than threshold,sum up
                        end
                    end
                    prec_season_extr_NE(i,j,k) = extr_temp_L;
                end
                if strcmp(season,'winter')
                    if k==97
                        prec_season_extr_NE(i,j,k) = NaN;
                    end
                end
            end
        end
    end    
    prec_season_extr_mean = squeeze(mean(prec_season_extr_NE,3)); % 165*250 matrix
    
    % calculate regional trend
    %   M-K test and Thei-sen slope
    datain = [linspace(1915,2011,97)',squeeze(nanmean(nanmean(prec_season_extr_NE(:,:,1:97),2)))];
    [taub tau h sig Z S sigma sen n senplot CIlower CIupper D Dall C3]=ktaub(datain,0.05,0);
    slope_MK = 10*sen
    p_value_MK = sig
    H_MK = h
    mean = nanmean(nanmean(nanmean(prec_season_extr_NE,2)))
    threshold_mean = nanmean(nanmean(threshold_L))
    figure
    STOP
%     %   linear trend estimate
%     disp('    Hereafter is linear trend estimate');
%     mdl = fitlm(linspace(1979,2011,33)',squeeze(nanmean(nanmean(prec_season_extr_NE(:,:,65:97),2)))); % provide linear trend and p-value
%     figure
%     plotResiduals(mdl)
%     [H1,p] = lillietest(mdl.Residuals.Raw) % H1=0 means a normal distribution at default alpha level 0.05
%     h1 = lbqtest(mdl.Residuals.Raw) % h1=0 means residuals are independent/non-autocorrelated
%     slope = 10*table2array(mdl.Coefficients(2,1)) % convert to mm/decade
%     p_value = table2array(mdl.Coefficients(2,4)) % for linear regression model

    STOP
    
    % calculate GHCN-D extreme threshold and then extract seasonal extreme prec
    load '/adata/Data/Observations/GHCN-D/station_data_adjust.mat'; % 41638*150
    % sort the sequence and find the 1% heaviest prec threshold for every station
    station_season_extr_NE = [];
    threshold = [];
    r = round(0.01*size(station_data_adjust,1)); % get 416 heaviest daily events
    for j=1:176
        threshold_temp = sort(station_data_adjust(:,j+4),'descend');
        threshold_temp = threshold_temp(~isnan(threshold_temp)); % exclude NaNs
        threshold(j) = threshold_temp(r+1); % find 417th heaviest prec value
    end
    
    for i=1:114 % if winter, use 113
        i
        for j=1:176
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
                if i==114
                    station_season_extr_NE(i,j) = NaN;
                else
                    season_start_G = day_mon_G(i*12)-datenum_start_G;
                    season_end_G = day_mon_G(i*12+3)-datenum_start_G-1;
                end
            end 
            
            station_season_extr_temp = sort(station_data_adjust(season_start_G:season_end_G,j+4),'descend');
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
                station_season_extr_NE(i,j) = extr_temp; % 114*176
            end
            if strcmp(season,'winter')
                if i==114
                    station_season_extr_NE(i,j) = NaN;
                end
            end
        end
    end    
    
    geoshow(prec_season_extr_mean,prec_refvec,'DisplayType','surface',...
        'ZData',zeros(size(prec_season_extr_mean)),'CData',prec_season_extr_mean);
    shading flat;
    hold on    
    scatter(geo(2,:),geo(1,:),100,0.1.*nanmean(station_season_extr_NE(15:110,:)),'filled','MarkerEdgeColor','k');
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    colorbar;
    title('Livneh & GHCN-D 1915-2011 spring extreme precipitation (mm)');
    ax = gca; % get the handle for the current axes
    set(ax,'FontSize',28);
    % ax.CLim=[-49 59]; % set the min and max limits for the colorbar
    colormap(flipud(jet)); % change default color to jet:blue-->red
    
    figure
    hold on
end

           % 5'.GHCN-D:plot trend in 1901-2014 seasonal extr PRCP (mm/decade)
if plot_season_extr_trend_GL==1
    % calculate linear trend in NARR seasonal extreme precipitation
    prec_season_extr_regress_coeff_temp = [];
    prec_season_extr_regress_coeff = [];
    for x=1:size(prec_season_extr_NE,1)
        for y=1:size(prec_season_extr_NE,2)
            if isnan(prec_season_extr_NE(x,y,1))
                prec_season_extr_regress_coeff(x,y) = NaN;
            else
                prec_season_extr_regress_coeff_temp = polyfit([1:1:size(prec_season_extr_NE,3)],squeeze(prec_season_extr_NE(x,y,:))',1);
                prec_season_extr_regress_coeff(x,y) = prec_season_extr_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
            end
        end
    end 
    prec_season_extr_decade_trend(:,:) = 10.*prec_season_extr_regress_coeff;  % changing yearly trends to be mm/decade

    % calculate and plot GHCN-D trend
    prcp_season_extr_regress_coeff_temp = [];
    prcp_season_extr_regress_coeff = [];
    for x=1:176  % if compare full period,use commentted two lines to replace the line 'extr_regress_coeff_temp...'
        ind = 1:length(station_season_extr_NE(15:110,x));
        k = ~isnan(station_season_extr_NE(15:110,x)); % find index with non-NaN values
        y = find(k==1);
        prcp_season_extr_regress_coeff_temp = polyfit(ind(k)',0.1*station_season_extr_NE(14+y,x),1);
        prcp_season_extr_regress_coeff(x) = prcp_season_extr_regress_coeff_temp(1); % polyfit produces two parameters, here we need coffieient, not intercept
    end
    prcp_season_extr_decade_trend = 10.*prcp_season_extr_regress_coeff;  % changing yearly trends to be mm/decade
    
    geoshow(prec_season_extr_decade_trend,prec_refvec,'DisplayType','surface',...
        'ZData',zeros(size(prec_season_extr_decade_trend)),'CData',prec_season_extr_decade_trend);
    shading flat;
    hold on
    scatter(geo(2,:),geo(1,:),100,prcp_season_extr_decade_trend,'filled','MarkerEdgeColor','k', 'LineWidth',1);
    geoshow(lat_NE_map,lon_NE_map,'Color','k');
    colorbar;
    title('Livneh & GHCN-D trends in 1915-2011 spring extreme precipitation (mm/decade)');
    ax = gca; % get the handle for the current axes
    set(ax,'FontSize',28);
    colormap(flipud(jet));
%     ax.CLim=[-49 59]; % set the min and max limits for the colorbar
%     blueColorMap = [ones(1,116), linspace(1, 0, 140)]; % [0,0,1] refers to blue; 116 and 140 is calculated from colorbar limits [-49,59], to make 0 be white
%     redColorMap = [linspace(0, 1, 116), linspace(1, 0, 140)];% [1,0,0] refers to red
%     whiteColorMap = [linspace(0, 1, 116), ones(1,140)]; % [1,1,1] refers to white
%     colorMap = [blueColorMap; redColorMap; whiteColorMap]'; % creat a 3-D matrix to install RGB values
%     colormap(colorMap); % change default color to blue-->white-->red
%         
%     if(sav==1)
%         set(gcf, 'Position', get(0,'Screensize')); % Maximize figure. 
%         saveas(gcf,'Figures/annual_mean_trend_G.jpg'); % need to manually save .jpg by replacing with same name
%     end
end