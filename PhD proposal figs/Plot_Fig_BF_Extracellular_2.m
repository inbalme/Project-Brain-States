%% plotting a figure "BF extracellular data_2"

close all
clear all

cd 'D:\Inbal M.Sc\Data PhD\Figures\PhD Proposal Figures\Figure BF extracellular_2'

%opening saved figures:
freq_1_zoom1 = open('2013-03-20_raster_PSTH_2Hz_zoom-in_pulse_1.fig');    
freq_1_zoom1_ax = get(gcf, 'children');
freq_1_zoom1_pos = [0.08 , 0.86 , 0.24 , 0.09];
freq_1_zoom1_pos_top = freq_1_zoom1_pos(1,2)+freq_1_zoom1_pos(1,4);

freq_1_zoom2 = open('2013-03-20_raster_PSTH_2Hz_zoom-in_pulse_2.fig');   
freq_1_zoom2_ax = get(gcf, 'children');
freq_1_zoom2_pos = [0.08 , 0.68 , 0.24 , 0.09];
freq_1_zoom2_pos_top = freq_1_zoom2_pos(1,2)+freq_1_zoom2_pos(1,4);


freq_2_zoom1 = open('2013-03-20_raster_PSTH_5Hz_zoom-in_pulse_1.fig');    
freq_2_zoom1_ax = get(gcf, 'children');
freq_2_zoom1_pos = [0.39 , 0.86 , 0.24 , 0.09];
freq_2_zoom1_pos_top = freq_2_zoom1_pos(1,2)+freq_2_zoom1_pos(1,4);

freq_2_zoom2 = open('2013-03-20_raster_PSTH_5Hz_zoom-in_pulse_5.fig');   
freq_2_zoom2_ax = get(gcf, 'children');
freq_2_zoom2_pos = [0.39 , 0.68 , 0.24 , 0.09];
freq_2_zoom2_pos_top = freq_2_zoom2_pos(1,2)+freq_2_zoom2_pos(1,4);
          
freq_3_zoom1 = open('2013-03-20_raster_PSTH_15Hz_zoom-in_pulse_1.fig');    
freq_3_zoom1_ax = get(gcf, 'children');
freq_3_zoom1_pos = [0.70 , 0.86 , 0.24 , 0.09];
freq_3_zoom1_pos_top = freq_3_zoom1_pos(1,2)+freq_3_zoom1_pos(1,4);

freq_3_zoom2 = open('2013-03-20_raster_PSTH_15Hz_zoom-in_pulse_5.fig');   
freq_3_zoom2_ax = get(gcf, 'children');
freq_3_zoom2_pos = [0.70 , 0.68 , 0.24 , 0.09];
freq_3_zoom2_pos_top = freq_3_zoom2_pos(1,2)+freq_3_zoom2_pos(1,4);

latency_2Hz = open('Mean_Latency_3_cells_2Hz_no_xlabel_y.fig');    
latency_2Hz_ax = get(gcf, 'children');
latency_2Hz_pos = [0.08 , 0.46 , 0.24 , 0.13];
latency_2Hz_pos_top = latency_2Hz_pos(1,2)+latency_2Hz_pos(1,4);

latency_5Hz = open('Mean_Latency_3_cells_5Hz_no_xlabel_y.fig');    
latency_5Hz_ax = get(gcf, 'children');
latency_5Hz_pos = [0.39 , 0.46 , 0.24 , 0.13];
latency_5Hz_pos_top = latency_5Hz_pos(1,2)+latency_5Hz_pos(1,4);

latency_15Hz = open('Mean_Latency_1_cell_15Hz_no_xlabel_y.fig');    
latency_15Hz_ax = get(gcf, 'children');
latency_15Hz_pos = [0.70 , 0.46 , 0.24 , 0.13];
latency_15Hz_pos_top = latency_15Hz_pos(1,2)+latency_15Hz_pos(1,4);

variance_2Hz = open('Variance_3_cells_2Hz_y.fig');    
variance_2Hz_ax = get(gcf, 'children');
variance_2Hz_pos = [0.08 , 0.26 , 0.24 , 0.13];
variance_2Hz_pos_top = variance_2Hz_pos(1,2)+variance_2Hz_pos(1,4);

variance_5Hz = open('Variance_3_cells_5Hz_y.fig');    
variance_5Hz_ax = get(gcf, 'children');
variance_5Hz_pos = [0.39 , 0.26 , 0.24 , 0.13];
variance_5Hz_pos_top = variance_5Hz_pos(1,2)+variance_5Hz_pos(1,4);

variance_15Hz = open('Variance_1_cell_15Hz_y.fig');    
variance_15Hz_ax = get(gcf, 'children');
variance_15Hz_pos = [0.70 , 0.26 , 0.24 , 0.13];
variance_15Hz_pos_top = variance_15Hz_pos(1,2)+variance_15Hz_pos(1,4);

%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
annotation('textbox', [freq_1_zoom1_pos(1,1) freq_1_zoom1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [freq_1_zoom2_pos(1,1) freq_1_zoom2_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')    
annotation('textbox', [freq_2_zoom1_pos(1,1) freq_2_zoom1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')  
annotation('textbox', [freq_2_zoom2_pos(1,1) freq_2_zoom2_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')           
annotation('textbox', [freq_3_zoom1_pos(1,1) freq_3_zoom1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')    
annotation('textbox', [freq_3_zoom2_pos(1,1) freq_3_zoom2_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')    
annotation('textbox', [latency_2Hz_pos(1,1) latency_2Hz_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')  
annotation('textbox', [latency_5Hz_pos(1,1) latency_5Hz_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')  
annotation('textbox', [latency_15Hz_pos(1,1) latency_15Hz_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')  
annotation('textbox', [variance_2Hz_pos(1,1) variance_2Hz_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')  
annotation('textbox', [variance_5Hz_pos(1,1) variance_5Hz_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'K',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')  
annotation('textbox', [variance_15Hz_pos(1,1) variance_15Hz_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'L',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')         
annotation('textbox', [freq_1_zoom1_pos(1,1)+freq_1_zoom1_pos(1,3)./2 freq_1_zoom1_pos_top 0 0]+[-0.03 0.02 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', '2Hz',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [freq_2_zoom1_pos(1,1)+freq_2_zoom1_pos(1,3)./2 freq_2_zoom1_pos_top 0 0]+[-0.03 0.02 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', '5Hz',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [freq_3_zoom1_pos(1,1)+freq_3_zoom1_pos(1,3)./2 freq_3_zoom1_pos_top 0 0]+[-0.03 0.02 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', '15Hz',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
       
%Placing plots in the figure:
freq_1_zoom1_ax_copy = copyobj(freq_1_zoom1_ax,F); % Copy freq_1_zoom1_ax to new fig
set(freq_1_zoom1_ax_copy(1),'position',freq_1_zoom1_pos(1,:)) % Set its position

freq_1_zoom2_ax_copy = copyobj(freq_1_zoom2_ax,F); % Copy freq_1_zoom2_ax to new fig
set(freq_1_zoom2_ax_copy(1),'position',freq_1_zoom2_pos(1,:)) % Set its position

freq_2_zoom1_ax_copy = copyobj(freq_2_zoom1_ax,F); % Copy freq_2_zoom1_ax to new fig
set(freq_2_zoom1_ax_copy(1),'position',freq_2_zoom1_pos(1,:)) % Set its position

freq_2_zoom2_ax_copy = copyobj(freq_2_zoom2_ax,F); % Copy freq_2_zoom2_ax to new fig
set(freq_2_zoom2_ax_copy(1),'position',freq_2_zoom2_pos(1,:)) % Set its position

freq_3_zoom1_ax_copy = copyobj(freq_3_zoom1_ax,F); % Copy freq_3_zoom1_ax to new fig
set(freq_3_zoom1_ax_copy(1),'position',freq_3_zoom1_pos(1,:)) % Set its position

freq_3_zoom2_ax_copy = copyobj(freq_3_zoom2_ax,F); % Copy freq_3_zoom2_ax to new fig
set(freq_3_zoom2_ax_copy(1),'position',freq_3_zoom2_pos(1,:)) % Set its position

latency_2Hz_ax_copy = copyobj(latency_2Hz_ax,F); % Copy latency_2Hz_ax to new fig
set(latency_2Hz_ax_copy(1),'position',latency_2Hz_pos(1,:)) % Set its position

latency_5Hz_ax_copy = copyobj(latency_5Hz_ax,F); % Copy latency_5Hz_ax to new fig
set(latency_5Hz_ax_copy(1),'position',latency_5Hz_pos(1,:)) % Set its position

latency_15Hz_ax_copy = copyobj(latency_15Hz_ax,F); % Copy latency_15Hz_ax to new fig
set(latency_15Hz_ax_copy(1),'position',latency_15Hz_pos(1,:)) % Set its position

variance_2Hz_ax_copy = copyobj(variance_2Hz_ax,F); % Copy variance_2Hz_ax to new fig
set(variance_2Hz_ax_copy(1),'position',variance_2Hz_pos(1,:)) % Set its position

variance_5Hz_ax_copy = copyobj(variance_5Hz_ax,F); % Copy variance_5Hz_ax to new fig
set(variance_5Hz_ax_copy(1),'position',variance_5Hz_pos(1,:)) % Set its position

variance_15Hz_ax_copy = copyobj(variance_15Hz_ax,F); % Copy variance_15Hz_ax to new fig
set(variance_15Hz_ax_copy(1),'position',variance_15Hz_pos(1,:)) % Set its position

print -depsc2 'Fig Extracellular BF_2'
