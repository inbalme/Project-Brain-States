%% Plot_Fig_Intracellular

close all
clear all

cd 'D:\Inbal M.Sc\Data PhD\Figures\PhD Proposal Figures\Figure intracellular'

Light_Off_t1 = open('2013-01-09-001_h2_t1.fig');    
Light_Off_t1_ax = get(gcf, 'children');
Light_Off_t1_pos = [0.08 , 0.84 , 0.40 , 0.09];
Light_Off_t1_pos_top = Light_Off_t1_pos(1,2)+Light_Off_t1_pos(1,4);

Light_Off_t2 = open('2013-01-09-001_h2_t2.fig');    
Light_Off_t2_ax = get(gcf, 'children');
Light_Off_t2_pos = [0.08 , 0.70 , 0.40 , 0.09];
Light_Off_t2_pos_top = Light_Off_t2_pos(1,2)+Light_Off_t2_pos(1,4);

Light_Off_t3 = open('2013-01-09-001_h2_t3.fig');    
Light_Off_t3_ax = get(gcf, 'children');
Light_Off_t3_pos = [0.08 , 0.56 , 0.40 , 0.09];
Light_Off_t3_pos_top = Light_Off_t3_pos(1,2)+Light_Off_t3_pos(1,4);

Light_Off_t4 = open('2013-01-09-001_h2_t4.fig');    
Light_Off_t4_ax = get(gcf, 'children');
Light_Off_t4_pos = [0.08 , 0.42 , 0.40 , 0.09];
Light_Off_t4_pos_top = Light_Off_t4_pos(1,2)+Light_Off_t4_pos(1,4);

Light_On_t21 = open('2013-01-09-001_h2_t21.fig');    
Light_On_t21_ax = get(gcf, 'children');
Light_On_t21_pos = [0.56 , 0.84 , 0.40 , 0.09];
Light_On_t21_pos_top = Light_On_t21_pos(1,2)+Light_On_t21_pos(1,4);

Light_On_t22 = open('2013-01-09-001_h2_t22.fig');    
Light_On_t22_ax = get(gcf, 'children');
Light_On_t22_pos = [0.56 , 0.70 , 0.40 , 0.09];
Light_On_t22_pos_top = Light_On_t22_pos(1,2)+Light_On_t22_pos(1,4);

Light_On_t23 = open('2013-01-09-001_h2_t23.fig');    
Light_On_t23_ax = get(gcf, 'children');
Light_On_t23_pos = [0.56 , 0.56 , 0.40 , 0.09];
Light_On_t23_pos_top = Light_On_t23_pos(1,2)+Light_On_t23_pos(1,4);

Light_On_t24 = open('2013-01-09-001_h2_t24.fig');    
Light_On_t24_ax = get(gcf, 'children');
Light_On_t24_pos = [0.56, 0.42 , 0.40 , 0.09];
Light_On_t24_pos_top = Light_On_t24_pos(1,2)+Light_On_t24_pos(1,4);


%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
     
annotation('textbox', [Light_Off_t1_pos(1,1) Light_Off_t1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [Light_On_t21_pos(1,1) Light_On_t21_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')       
       
%Placing plots in the figure:
Light_Off_t1_ax_copy = copyobj(Light_Off_t1_ax,F); % Copy Light_Off_t1_ax to new fig
set(Light_Off_t1_ax_copy(1),'position',Light_Off_t1_pos(1,:)) % Set its position 

Light_Off_t2_ax_copy = copyobj(Light_Off_t2_ax,F); % Copy Light_Off_t2_ax to new fig
set(Light_Off_t2_ax_copy(1),'position',Light_Off_t2_pos(1,:)) % Set its position  

Light_Off_t3_ax_copy = copyobj(Light_Off_t3_ax,F); % Copy Light_Off_t3_ax to new fig
set(Light_Off_t3_ax_copy(1),'position',Light_Off_t3_pos(1,:)) % Set its position  

Light_Off_t4_ax_copy = copyobj(Light_Off_t4_ax,F); % Copy Light_Off_t4_ax to new fig
set(Light_Off_t4_ax_copy(1),'position',Light_Off_t4_pos(1,:)) % Set its position  

Light_On_t21_ax_copy = copyobj(Light_On_t21_ax,F); % Copy Light_On_t21_ax to new fig
set(Light_On_t21_ax_copy(1),'position',Light_On_t21_pos(1,:)) % Set its position 

Light_On_t22_ax_copy = copyobj(Light_On_t22_ax,F); % Copy Light_On_t22_ax to new fig
set(Light_On_t22_ax_copy(1),'position',Light_On_t22_pos(1,:)) % Set its position 

Light_On_t23_ax_copy = copyobj(Light_On_t23_ax,F); % Copy Light_On_t23_ax to new fig
set(Light_On_t23_ax_copy(1),'position',Light_On_t23_pos(1,:)) % Set its position 

Light_On_t24_ax_copy = copyobj(Light_On_t24_ax,F); % Copy Light_On_t24_ax to new fig
set(Light_On_t24_ax_copy(1),'position',Light_On_t24_pos(1,:)) % Set its position 

 print -depsc2 'Fig Intracellular'
      