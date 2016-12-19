%% Plot LFP in the cortex in ChAT mice

close all
clear all

cd 'D:\Inbal M.Sc\Data PhD\Figures\PhD Proposal Figures\Figure LFP'

BC_LFP_trace17 = open('2013-02-13-002_h3_BC_LFP_trace17.fig');    
BC_LFP_trace17_ax = get(gcf, 'children');
BC_LFP_trace17_pos = [0.08 , 0.84 , 0.40 , 0.09];
BC_LFP_trace17_pos_top = BC_LFP_trace17_pos(1,2)+BC_LFP_trace17_pos(1,4);

BC_LFP_trace17_filt = open('2013-02-13-002_h3_BC_LFP_trace17_BPF05-5.fig');    
BC_LFP_trace17_filt_ax = get(gcf, 'children');
BC_LFP_trace17_filt_pos = [0.08 , 0.70 , 0.40 , 0.09];
BC_LFP_trace17_filt_pos_top = BC_LFP_trace17_filt_pos(1,2)+BC_LFP_trace17_filt_pos(1,4);

BC_LFP_mean = open('2013-02-13-002_h3_BC_LFP_mean+std.fig');    
BC_LFP_mean_ax = get(gcf, 'children');
BC_LFP_mean_pos = [0.08 , 0.56 , 0.40 , 0.09];
BC_LFP_mean_pos_top = BC_LFP_mean_pos(1,2)+BC_LFP_mean_pos(1,4);

BC_LFP_mean_filt = open('2013-02-13-002_h3_BC_LFP_mean+std_BPF05-5.fig');    
BC_LFP_mean_filt_ax = get(gcf, 'children');
BC_LFP_mean_filt_pos = [0.08 , 0.42 , 0.40 , 0.09];
BC_LFP_mean_filt_pos_top = BC_LFP_mean_filt_pos(1,2)+BC_LFP_mean_filt_pos(1,4);

M1_LFP_trace17 = open('2013-03-20-001_h6_M1_LFP_trace17.fig');    
M1_LFP_trace17_ax = get(gcf, 'children');
M1_LFP_trace17_pos = [0.56 , 0.84 , 0.40 , 0.09];
M1_LFP_trace17_pos_top = M1_LFP_trace17_pos(1,2)+M1_LFP_trace17_pos(1,4);

M1_LFP_trace17_filt = open('2013-03-20-001_h6_M1_LFP_trace17_BPF05-5.fig');    
M1_LFP_trace17_filt_ax = get(gcf, 'children');
M1_LFP_trace17_filt_pos = [0.56 , 0.70 , 0.40 , 0.09];
M1_LFP_trace17_filt_pos_top = M1_LFP_trace17_filt_pos(1,2)+M1_LFP_trace17_filt_pos(1,4);

M1_LFP_mean = open('2013-03-20-001_h6_M1_LFP_mean+std.fig');    
M1_LFP_mean_ax = get(gcf, 'children');
M1_LFP_mean_pos = [0.56 , 0.56 , 0.40 , 0.09];
M1_LFP_mean_pos_top = M1_LFP_mean_pos(1,2)+M1_LFP_mean_pos(1,4);

M1_LFP_mean_filt = open('2013-03-20-001_h6_M1_LFP_mean+std_BPF05-5.fig');    
M1_LFP_mean_filt_ax = get(gcf, 'children');
M1_LFP_mean_filt_pos = [0.56 , 0.42 , 0.40 , 0.09];
M1_LFP_mean_filt_pos_top = M1_LFP_mean_filt_pos(1,2)+M1_LFP_mean_filt_pos(1,4);

%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
     
annotation('textbox', [BC_LFP_trace17_pos(1,1) BC_LFP_trace17_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [BC_LFP_trace17_filt_pos(1,1) BC_LFP_trace17_filt_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [BC_LFP_mean_pos(1,1) BC_LFP_mean_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [BC_LFP_mean_filt_pos(1,1) BC_LFP_mean_filt_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')         
annotation('textbox', [M1_LFP_trace17_pos(1,1) M1_LFP_trace17_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [M1_LFP_trace17_filt_pos(1,1) M1_LFP_trace17_filt_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')       
annotation('textbox', [M1_LFP_mean_pos(1,1) M1_LFP_mean_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [M1_LFP_mean_filt_pos(1,1) M1_LFP_mean_filt_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [BC_LFP_trace17_pos(1,1)+BC_LFP_trace17_pos(1,3)./2 BC_LFP_trace17_pos_top 0 0]+[-0.08 0.02 0.20 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Barrel Cortex',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [M1_LFP_trace17_pos(1,1)+M1_LFP_trace17_pos(1,3)./2 M1_LFP_trace17_pos_top 0 0]+[-0.08 0.02 0.20 0.04],...
           'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Motor Cortex',...
           'FontName','helvetica', 'fontsize', 12, 'fontweight', 'bold')
       
%Placing plots in the figure:
BC_LFP_mean_ax_copy = copyobj(BC_LFP_mean_ax,F); % Copy BC_LFP_mean_ax to new fig
set(BC_LFP_mean_ax_copy(1),'position',BC_LFP_mean_pos(1,:)) % Set its position  

BC_LFP_mean_filt_ax_copy = copyobj(BC_LFP_mean_filt_ax,F); % Copy BC_LFP_mean_filt_ax to new fig
set(BC_LFP_mean_filt_ax_copy(1),'position',BC_LFP_mean_filt_pos(1,:)) % Set its position  

BC_LFP_trace17_ax_copy = copyobj(BC_LFP_trace17_ax,F); % Copy BC_LFP_trace17_ax to new fig
set(BC_LFP_trace17_ax_copy(1),'position',BC_LFP_trace17_pos(1,:)) % Set its position  

BC_LFP_trace17_filt_ax_copy = copyobj(BC_LFP_trace17_filt_ax,F); % Copy BC_LFP_trace17_filt_ax to new fig
set(BC_LFP_trace17_filt_ax_copy(1),'position',BC_LFP_trace17_filt_pos(1,:)) % Set its position  

M1_LFP_mean_ax_copy = copyobj(M1_LFP_mean_ax,F); % Copy M1_LFP_mean_ax to new fig
set(M1_LFP_mean_ax_copy(1),'position',M1_LFP_mean_pos(1,:)) % Set its position  

M1_LFP_mean_filt_ax_copy = copyobj(M1_LFP_mean_filt_ax,F); % Copy M1_LFP_mean_filt_ax to new fig
set(M1_LFP_mean_filt_ax_copy(1),'position',M1_LFP_mean_filt_pos(1,:)) % Set its position  

M1_LFP_trace17_ax_copy = copyobj(M1_LFP_trace17_ax,F); % Copy M1_LFP_trace17_ax to new fig
set(M1_LFP_trace17_ax_copy(1),'position',M1_LFP_trace17_pos(1,:)) % Set its position  

M1_LFP_trace17_filt_ax_copy = copyobj(M1_LFP_trace17_filt_ax,F); % Copy M1_LFP_trace17_filt_ax to new fig
set(M1_LFP_trace17_filt_ax_copy(1),'position',M1_LFP_trace17_filt_pos(1,:)) % Set its position  

 print -depsc2 'Fig LFP Cortex'