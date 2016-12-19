%% plotting a figure "BF extracellular data_1"
close all
clear all

cd 'D:\Inbal M.Sc\Data PhD\Figures\PhD Proposal Figures\Figure BF extracellular_1'

%opening saved figures:
optopatcher = imread('Optopatcher','JPG');
optopatcher_pos = [0.02 , 0.74 , 0.24 , 0.20];
optopatcher_pos_top = optopatcher_pos(1,2)+optopatcher_pos(1,4);
% 
slice = imread('ChAT-ChR2_27_slide9_slice_003','tiff');
slice_pos = [0.28 , 0.74 , 0.24 , 0.20];
slice_pos_top = slice_pos(1,2)+slice_pos(1,4);

mean_LFP = open('2013-02-27-002_h1+2_BF_LFP_mean.fig');
mean_LFP_ax = get(gcf, 'children');
mean_LFP_pos = [0.86, 0.90, 0.07, 0.03; ...
                0.62 , 0.77 , 0.32 , 0.16];
mean_LFP_pos_top = mean_LFP_pos(2,2)+mean_LFP_pos(2,4);

clust = open('2013-03-20_clustering.fig');  %contains 2 axes
clust_ax = get(gcf, 'children');
clust_pos = [0.36 , 0.56 , 0.24 , 0.10; ...
             0.08 , 0.56 , 0.20 , 0.10];
clust_pos_top = [clust_pos(1,2)+clust_pos(1,4);...
                 clust_pos(2,2)+clust_pos(2,4)];

raw_data_trace = open('2013-03-20_h12_trace10_15Hz.fig');
trace_ax = get(gcf, 'children');
trace_pos = [0.68 , 0.56 , 0.28 , 0.10];
trace_pos_top = trace_pos(1,2)+trace_pos(1,4);

freq_1 = open('2013-03-20_raster_PSTH_2Hz.fig');    %contains 2 axes
freq_1_ax = get(gcf, 'children');
freq_1_pos = [0.08 , 0.22 , 0.24 , 0.09; ...
              0.08 , 0.33 , 0.24 , 0.09];
freq_1_pos_top = freq_1_pos(2,2)+freq_1_pos(2,4);

freq_2 = open('2013-03-20_raster_PSTH_5Hz.fig');    %contains 2 axes
freq_2_ax = get(gcf, 'children');
freq_2_pos = [0.39 , 0.22 , 0.24 , 0.09; ...
              0.39 , 0.33 , 0.24 , 0.09];
freq_2_pos_top = freq_2_pos(2,2)+freq_2_pos(2,4);

freq_3 = open('2013-03-20_raster_PSTH_15Hz.fig');    %contains 2 axes
freq_3_ax = get(gcf, 'children');
freq_3_pos = [0.70 , 0.22 , 0.24 , 0.09; ...
              0.70 , 0.33 , 0.24 , 0.09];
freq_3_pos_top = freq_3_pos(2,2)+freq_3_pos(2,4);

%open a new figure:
Fig_BF1 = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
annotation('textbox', [optopatcher_pos(1,1) optopatcher_pos_top 0 0]+[0.02 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [slice_pos(1,1) slice_pos_top 0 0]+[-0.05 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [mean_LFP_pos(2,1) slice_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [clust_pos(2,1) clust_pos_top(2) 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [clust_pos(1,1) clust_pos_top(1) 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')       
annotation('textbox', [trace_pos(1,1) trace_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [freq_1_pos(2,1) freq_1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [freq_2_pos(2,1) freq_2_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [freq_3_pos(2,1) freq_3_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [freq_1_pos(2,1)+freq_1_pos(2,3)./2 freq_1_pos_top 0 0]+[-0.03 0.02 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', '2Hz',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')       
annotation('textbox', [freq_2_pos(1,1)+freq_2_pos(1,3)./2 freq_2_pos_top 0 0]+[-0.03 0.02 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', '5Hz',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [freq_3_pos(1,1)+freq_3_pos(1,3)./2 freq_3_pos_top 0 0]+[-0.03 0.02 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', '15Hz',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
       
       %Placing plots in the figure:
optopatcher_ax = axes('position',optopatcher_pos);
imshow(optopatcher, 'parent', optopatcher_ax) 

slice_ax = axes('position',slice_pos);
imshow(slice, 'parent', slice_ax) 

mean_LFP_ax_copy = copyobj(mean_LFP_ax,Fig_BF1); % Copy trace_ax to new fig
set(mean_LFP_ax_copy(1),'position',mean_LFP_pos(1,:)) % Set its position
set(mean_LFP_ax_copy(2),'position',mean_LFP_pos(2,:)) % Set its position

trace_ax_copy = copyobj(trace_ax,Fig_BF1); % Copy trace_ax to new fig
set(trace_ax_copy,'position',trace_pos(1,:)) % Set its position

clust_ax_copy = copyobj(clust_ax,Fig_BF1); % Copy clust_ax to new fig
set(clust_ax_copy(1),'position',clust_pos(1,:)) % Set its position
set(clust_ax_copy(2),'position',clust_pos(2,:)) % Set its position

freq_1_ax_copy = copyobj(freq_1_ax,Fig_BF1); % Copy freq_1_ax to new fig
set(freq_1_ax_copy(1),'position',freq_1_pos(1,:)) % Set its position
set(freq_1_ax_copy(2),'position',freq_1_pos(2,:)) % Set its position

freq_2_ax_copy = copyobj(freq_2_ax,Fig_BF1); % Copy freq_2_ax to new fig
set(freq_2_ax_copy(1),'position',freq_2_pos(1,:)) % Set its position
set(freq_2_ax_copy(2),'position',freq_2_pos(2,:)) % Set its position

freq_3_ax_copy = copyobj(freq_3_ax,Fig_BF1); % Copy freq_3_ax to new fig
set(freq_3_ax_copy(1),'position',freq_3_pos(1,:)) % Set its position
set(freq_3_ax_copy(2),'position',freq_3_pos(2,:)) % Set its position

% saveas(7,'BF1.eps');
print -depsc2 'Fig Extracellular BF_1'


% %Specifying the axis:
% xlimits = [0.7 4.3];
% xticks = 1 : 4 ;
% ylimits = [-28 2];
% yticks = [-28 0];
% 
% set( gca, 'xlim', xlimits, 'xtick', xticks, 'ylim', ylimits, 'ytick',...
% [ylimits(1) 0 ylimits(2)], 'ticklength', [0.030 0.030], 'box', 'off' );
% % Set the limits and ticks you defined earlier