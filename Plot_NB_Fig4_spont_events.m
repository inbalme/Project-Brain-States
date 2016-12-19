%% plotting a figure 

close all
clear all


%opening saved figures:
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'

spont_Vm_5prcentile = open('Vm_5prcentile_ongoing.fig');    
spont_Vm_5prcentile_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz\single trials\ongoing'

spont_event_freq = open('Spontaneous event frequency.fig');    
spont_event_freq_ax = get(gcf, 'children');

spont_event_amp = open('Spontaneous event amplitude.fig');    
spont_event_amp_ax = get(gcf, 'children');

spont_event_halfwidth = open('Spontaneous event half-width.fig');    
spont_event_halfwidth_ax = get(gcf, 'children');

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 15 12]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
 
%% Positions:
spont_Vm_5prcentile_pos(1,:) = [0.15 , 0.55 , 0.35 , 0.35];
spont_event_freq_pos(1,:) = [0.15 , 0.6 , 0.35 , 0.35]; spont_event_freq_pos(1,2)=spont_Vm_5prcentile_pos(1,2)-spont_Vm_5prcentile_pos(1,4)-0.1;

h_dist1=spont_Vm_5prcentile_pos(1,1)+spont_Vm_5prcentile_pos(1,3)+0.12;
spont_event_amp_pos(1,:) = [h_dist1 , spont_Vm_5prcentile_pos(1,2) ,  spont_Vm_5prcentile_pos(1,3) ,  spont_Vm_5prcentile_pos(1,4)];
spont_event_halfwidth_pos(1,:) = [h_dist1, spont_event_freq_pos(1,2) ,  spont_event_freq_pos(1,3) ,  spont_event_freq_pos(1,4)];

%top positions:
spont_Vm_5prcentile_pos_top = spont_Vm_5prcentile_pos(1,2)+spont_Vm_5prcentile_pos(1,4);
spont_event_freq_pos_top = spont_event_freq_pos(1,2)+spont_event_freq_pos(1,4);
spont_event_amp_pos_top = spont_event_amp_pos(1,2)+spont_event_amp_pos(1,4);
spont_event_halfwidth_pos_top = spont_event_halfwidth_pos(1,2)+spont_event_halfwidth_pos(1,4);


%%
%Placing plots in the figure:
spont_Vm_5prcentile_ax_copy = copyobj(spont_Vm_5prcentile_ax,F); % copy axes to new fig
set(spont_Vm_5prcentile_ax_copy,'position',spont_Vm_5prcentile_pos(1,:))
set(F, 'currentaxes', spont_Vm_5prcentile_ax_copy);  xl=xlabel('');  tl=title('');
yl=ylabel('5 percentile [mV]');
spont_Vm_5prcentile_ax_copy.FontSize=12;
spont_Vm_5prcentile_ax_copy.XTickLabel=[];

spont_event_freq_ax_copy = copyobj(spont_event_freq_ax,F); % copy axes to new fig
set(spont_event_freq_ax_copy,'position',spont_event_freq_pos(1,:))
set(F, 'currentaxes', spont_event_freq_ax_copy);  xl=xlabel('');  tl=title('');
yl=ylabel('frequency [Hz]');
spont_event_freq_ax_copy.FontSize=12;

  spont_event_amp_ax_copy = copyobj(spont_event_amp_ax,F); % copy axes to new fig
set(spont_event_amp_ax_copy,'position',spont_event_amp_pos(1,:))
set(F, 'currentaxes', spont_event_amp_ax_copy);  xl=xlabel('');  tl=title(''); 
yl=ylabel('Amplitude [mV]');
spont_event_amp_ax_copy.FontSize=12;
spont_event_amp_ax_copy.XTickLabel=[];

spont_event_halfwidth_ax_copy = copyobj(spont_event_halfwidth_ax,F); % copy axes to new fig
set(spont_event_halfwidth_ax_copy,'position',spont_event_halfwidth_pos(1,:))
set(F, 'currentaxes', spont_event_halfwidth_ax_copy);  xl=xlabel('');  tl=title('');
yl=ylabel('Half-Width [mS]');
spont_event_halfwidth_ax_copy.FontSize=12;

 a_pos1=[-0.06 0.05 0.04 0.04];
 
annotation('textbox', [spont_Vm_5prcentile_pos(1,1),spont_Vm_5prcentile_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_event_freq_pos(1,1) spont_event_freq_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_event_amp_pos(1,1) spont_event_amp_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_event_halfwidth_pos(1,1) spont_event_halfwidth_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')


cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures'
filename='Fig 4 spont events';
saveas(F,'Fig 4 spont events.fig'); 
print(F,filename,'-dpng','-r600','-opengl') 
print(F, '-depsc2', filename);
