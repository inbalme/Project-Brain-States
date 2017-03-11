%% plotting a figure 

close all
clear all
save_flag=1;
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP filtered 49-51Hz Presentation'
cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm-LFP correlations\LFP_50Hz+BP0.1-200 Vm_50Hz';
%opening saved figures:
%%
% cross-correlations

%population parameters:
%actual_data
cc_paired_plot = open('Vm-LFPcc_spont+evoked_max-peak_paired_population.fig');    
cc_paired_plot_ax = get(gcf, 'children');
%shuffled data
cc_paired_plot_shuff = open('Vm-LFPcc_spont+evoked_max-peak_paired_population_shuff.fig');    
cc_paired_plot_shuff_ax = get(gcf, 'children');

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 8]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
 
%% Positions:

cc_paired_plot_shuff_pos(1,:) = [0.1 , 0.1 , 0.35 , 0.75]; 
cc_paired_plot_pos(1,:) =  [cc_paired_plot_shuff_pos(1,1)+cc_paired_plot_shuff_pos(1,3)+0.15 , cc_paired_plot_shuff_pos(1,2) ,cc_paired_plot_shuff_pos(1,3) , cc_paired_plot_shuff_pos(1,4)];
cc_paired_plot_shuff_pos_top = cc_paired_plot_shuff_pos(1,2)+cc_paired_plot_shuff_pos(1,4);
cc_paired_plot_pos_top = cc_paired_plot_pos(1,2)+cc_paired_plot_pos(1,4);

%%
%Placing plots in the figure:

%cc population paired plots
cc_paired_plot_ax_copy = copyobj(cc_paired_plot_ax,F); % copy axes to new fig
set(cc_paired_plot_ax_copy,'position',cc_paired_plot_pos(1,:))
   set(F, 'currentaxes', cc_paired_plot_ax_copy); yl=ylabel('');
cc_paired_plot_ax_copy.FontSize=14;
cc_paired_plot_ax_copy.XTickLabel={'Off', 'On','Off', 'On'}; 

cc_paired_plot_shuff_ax_copy = copyobj(cc_paired_plot_shuff_ax,F); % copy axes to new fig
set(cc_paired_plot_shuff_ax_copy,'position',cc_paired_plot_shuff_pos(1,:))
cc_paired_plot_shuff_ax_copy.FontSize=14;
cc_paired_plot_shuff_ax_copy.XTickLabel={'Off', 'On','Off', 'On'}; 
    %   'fontname', 'arial','fontsize',13,'linewidth',1.5,'box','off'); %, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks, 'xticklabel',x6ticklab);
%  set(F, 'currentaxes', Response_modulation_ax_copy); t=title(''); yl=ylabel('Spikes/Stim. train','fontsize',13);


%   annotation:
annotation('textbox', [cc_paired_plot_shuff_pos(1,1) cc_paired_plot_shuff_pos_top 0 0]+[0.01 0.08 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Shift-Predictor Correlations', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
annotation('textbox', [cc_paired_plot_pos(1,1) cc_paired_plot_pos_top 0 0]+[0.01 0.08 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Shift-Corrected Correlations', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
%  
 a_pos1=[-0.03 -0.02 0.04 0.04];
%  a_pos2=[-0.03 0 0.04 0.04];
 a_pos3=[-0.04 -0.01 0.04 0.04];
 a_pos4=[-0.06 0.04 0.04 0.04];

  annotation('textbox', [cc_paired_plot_shuff_pos(1,1) cc_paired_plot_shuff_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
  annotation('textbox', [cc_paired_plot_pos(1,1) cc_paired_plot_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')


%  annotation('textbox', [cc_paired_plot_pos(1,1) cc_paired_plot_pos_top 0 0]+[0.02 0.0 0.2 0.04],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Noise Correlations', 'FontName','arial', 'fontsize', 12,'color', [0 0 153]/256)
%   annotation('textbox', [cc_paired_plot_shuff_pos(1,1) cc_paired_plot_shuff_pos_top 0 0]+[0.02 0.0 0.2 0.04],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Signal correlations', 'FontName','arial', 'fontsize', 12,'color', [0 0 153]/256)
%% 
 
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures'
if save_flag==1
    filename='Fig 5b ChAT ccVmLFP';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end
