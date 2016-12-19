%% plotting a figure "VmLFPcc" for file 44

close all
clear all

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP filtered 49-51Hz'

%opening saved figures:
%spontaneous:
spont_trace_Off_f44 = open('Vm-LFP_spont_ES_Off_file_44.fig');    
spont_trace_Off_f44_ax = get(gcf, 'children');
spont_trace_Off_f44_pos = [0.08 , 0.76 , 0.31 , 0.15];
spont_trace_Off_f44_pos_top = spont_trace_Off_f44_pos(1,2)+spont_trace_Off_f44_pos(1,4);

spont_trace_On_f44 = open('Vm-LFP_spont_ES_On_file_44.fig');    
spont_trace_On_f44_ax = get(gcf, 'children');
spont_trace_On_f44_pos = [0.43 , 0.76 , 0.31 , 0.15];
spont_trace_On_f44_pos_top = spont_trace_On_f44_pos(1,2)+spont_trace_On_f44_pos(1,4);

spont_VmLFPcc_f44 = open('Vm-LFPcc_spont_file_44.fig');    
spont_VmLFPcc_f44_ax = get(gcf, 'children');
spont_VmLFPcc_f44_pos = [0.84 , 0.76 , 0.12 , 0.15];
spont_VmLFPcc_f44_pos_top = spont_VmLFPcc_f44_pos(1,2)+spont_VmLFPcc_f44_pos(1,4);
%evoked
evoked_trace_Off_f44 = open('Vm-LFP_evoked_ES_Off_file_44.fig');    
evoked_trace_Off_f44_ax = get(gcf, 'children');
evoked_trace_Off_f44_pos = [0.08 , 0.5 , 0.31 , 0.15];
evoked_trace_Off_f44_pos_top = evoked_trace_Off_f44_pos(1,2)+evoked_trace_Off_f44_pos(1,4);

evoked_trace_On_f44 = open('Vm-LFP_evoked_ES_On_file_44.fig');    
evoked_trace_On_f44_ax = get(gcf, 'children');
evoked_trace_On_f44_pos = [0.43 , 0.5 , 0.31 , 0.15];
evoked_trace_On_f44_pos_top = evoked_trace_On_f44_pos(1,2)+evoked_trace_On_f44_pos(1,4);

evoked_VmLFPcc_f44 = open('Vm-LFPcc_evoked_file_44.fig');    
evoked_VmLFPcc_f44_ax = get(gcf, 'children');
evoked_VmLFPcc_f44_pos = [0.84 , 0.5 , 0.12 , 0.15];
evoked_VmLFPcc_f44_pos_top = evoked_VmLFPcc_f44_pos(1,2)+evoked_VmLFPcc_f44_pos(1,4);
%population plots
population_spont = open('Vm-LFPcc_spont_zero-lag_paired_population.fig');    
population_spont_ax = get(gcf, 'children');
population_spont_pos = [0.23 , 0.1 , 0.15 , 0.2];
population_spont_pos_top = population_spont_pos(1,2)+population_spont_pos(1,4);

population_evoked = open('Vm-LFPcc_evoked_zero-lag_paired_population.fig');    
population_evoked_ax = get(gcf, 'children');
population_evoked_pos = [0.63, 0.1 , 0.15 , 0.2];
population_evoked_pos_top = population_evoked_pos(1,2)+population_evoked_pos(1,4);


%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
% annotation('textbox', [freq_1_zoom1_pos(1,1) freq_1_zoom1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
%            'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A',...
%            'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
% annotation('textbox', [freq_1_zoom2_pos(1,1) freq_1_zoom2_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
%            'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B',...
%            'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')    
% annotation('textbox', [freq_2_zoom1_pos(1,1) freq_2_zoom1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
%            'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C',...
%            'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')  
%% xlimits, y limits, ticks etc.
 x1limits = [0.5 2.5];
x1ticks = [0.5, 1, 1.5, 2, 2.5];
x1ticklab=[0, 0.5, 1, 1.5, 2];
y1limits = [-10 20];
y1ticks = [-10 0 10 20];   
         x2limits = [5 7];
        x2ticks = [5, 5.5, 6, 6.5, 7];
        x2ticklab=[0, 0.5, 1, 1.5, 2];
        y2limits = [-10 20];
        y2ticks = [-10 0 10 20];
        y2ticklab=[];
             x3limits = [-0.2 0.2];
            x3ticks = [-0.2, 0, 0.2];
            x3ticklab=[-0.2, 0, 0.2];
            y3limits = [-0.2 0.3];
            y3ticks = [-0.2 0 0.3];
 
x4limits = [6 7];
x4ticks = [6, 6.5, 7];
x4ticklab=[0, 0.5, 1];
y4limits = [-10 20];
y4ticks = [-10 0 10 20];   
        x5limits = [6 7];
        x5ticks = [6, 6.5, 7];
        x5ticklab=[0, 0.5, 1];
        y5limits = [-10 20];
        y5ticks = [-10 0 10 20];
        y5ticklab=[];
             x6limits = [-0.2 0.2];
            x6ticks = [-0.2, 0, 0.2];
            x6ticklab=[-0.2, 0, 0.2];
            y6limits = [-0.2 0.3];
            y6ticks = [-0.2 0 0.3];
            %%
%Placing plots in the figure:
spont_trace_Off_f44_ax_copy = copyobj(spont_trace_Off_f44_ax,F); % copy axes to new fig
set(spont_trace_Off_f44_ax_copy,'position',spont_trace_Off_f44_pos(1,:),...
       'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y1limits,'ytick', y1ticks,'xlim',x1limits,'xtick',x1ticks,...
       'xticklabel',x1ticklab,'ticklength',[0.01 0.025]);
set(F, 'currentaxes', spont_trace_Off_f44_ax_copy); t=title('NB-'); xl=xlabel('');
   l=legend('Vm','LFP'); set(l,'box','off'),  set(l,'FontSize', 8), set(l,'linewidth',1.5); set(l,'position',[0.35 0.89 0.03 0.03]); 
    
spont_trace_On_f44_ax_copy = copyobj(spont_trace_On_f44_ax,F); % copy axes to new fig
set(spont_trace_On_f44_ax_copy,'position',spont_trace_On_f44_pos(1,:),...
       'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y2limits,'ytick', y2ticks, 'yticklabel',y2ticklab,'xlim',x2limits,'xtick',x2ticks,...
       'xticklabel',x2ticklab,'ticklength',[0.01 0.025]);
  set(F, 'currentaxes', spont_trace_On_f44_ax_copy); t=title('NB+'); yl=ylabel(''); xl=xlabel('');
   l=legend('Vm','LFP'); set(l,'box','off'), set(l,'FontSize', 8) ; set(l,'linewidth',1.5);set(l,'position',[0.7 0.89 0.03 0.03]);

spont_VmLFPcc_f44_ax_copy = copyobj(spont_VmLFPcc_f44_ax,F); % copy axes to new fig
set(spont_VmLFPcc_f44_ax_copy,'position',spont_VmLFPcc_f44_pos(1,:),...
       'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y3limits,'ytick', y3ticks,'xlim',x3limits,'xtick',x3ticks,...
       'xticklabel',x3ticklab,'ticklength',[0.01 0.025]);
 set(F, 'currentaxes', spont_VmLFPcc_f44_ax_copy); t=title(''); xl=xlabel('');
   l=legend('NB-','NB+'); set(l,'box','off') ;set(l,'FontSize', 8);set(l,'linewidth',1.5);set(l,'position',[0.93 0.89 0.03 0.03]);

   %evoked
   evoked_trace_Off_f44_ax_copy = copyobj(evoked_trace_Off_f44_ax,F); % copy axes to new fig
set(evoked_trace_Off_f44_ax_copy,'position',evoked_trace_Off_f44_pos(1,:),...
       'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y4limits,'ytick', y4ticks,'xlim',x4limits,'xtick',x4ticks,...
       'xticklabel',x4ticklab,'ticklength',[0.01 0.025]);
set(F, 'currentaxes', evoked_trace_Off_f44_ax_copy); t=title('');
   l=legend('Vm','LFP'); set(l,'box','off');  set(l,'FontSize', 8);set(l,'linewidth',1.5); set(l,'position',[0.35 0.63 0.03 0.03]);
    
evoked_trace_On_f44_ax_copy = copyobj(evoked_trace_On_f44_ax,F); % copy axes to new fig
set(evoked_trace_On_f44_ax_copy,'position',evoked_trace_On_f44_pos(1,:),...
       'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y5limits,'ytick', y5ticks, 'yticklabel',y5ticklab,'xlim',x5limits,'xtick',x5ticks,...
       'xticklabel',x5ticklab,'ticklength',[0.01 0.025]);
  set(F, 'currentaxes', evoked_trace_On_f44_ax_copy); t=title(''); yl=ylabel(''); %delete(t); delete(yl);
   l=legend('Vm','LFP'); set(l,'box','off'), set(l,'FontSize', 8); set(l,'linewidth',1.5);set(l,'position',[0.7 0.63 0.03 0.03]);

evoked_VmLFPcc_f44_ax_copy = copyobj(evoked_VmLFPcc_f44_ax,F); % copy axes to new fig
set(evoked_VmLFPcc_f44_ax_copy,'position',evoked_VmLFPcc_f44_pos(1,:),...
       'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks,...
       'xticklabel',x6ticklab,'ticklength',[0.01 0.025]);
 set(F, 'currentaxes', evoked_VmLFPcc_f44_ax_copy); t=title('');
   l=legend('NB-','NB+'); set(l,'box','off'); set(l,'FontSize', 8);set(l,'linewidth',1.5); set(l,'position',[0.93 0.63 0.03 0.03]);
%population
population_spont_ax_copy = copyobj(population_spont_ax,F); % copy axes to new fig
set(population_spont_ax_copy,'position',population_spont_pos(1,:),...
       'fontname', 'arial','fontsize',13,'linewidth',1.5,'box','off');
 set(F, 'currentaxes', population_spont_ax_copy); t=title('Spontaneous','fontsize',16);ylabel('Correlation','fontsize',16);
 
population_evoked_ax_copy = copyobj(population_evoked_ax,F); % copy axes to new fig
set(population_evoked_ax_copy,'position',population_evoked_pos(1,:),...
       'fontname', 'arial','fontsize',13,'linewidth',1.5,'box','off');
 set(F, 'currentaxes', population_evoked_ax_copy); t=title('Evoked','fontsize',16);ylabel('Correlation','fontsize',16);

      

filename='Fig_VmLFPcc_file44';
saveas(F,'Fig_VmLFPcc_file44.fig'); 
print(F,filename,'-dpng','-r600','-opengl') 
% print(F, '-depsc2', filename); 


