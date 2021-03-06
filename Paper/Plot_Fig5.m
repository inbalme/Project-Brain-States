%% plotting a figure 
close all
clear all
save_flag=0;
 no_numbering_flag=0;
 ax_fontsize=10;
 barWidth=1.5;
letter_fontsize=12; 
%opening saved figures:
%Short traces:
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Zoom-in Trace Presentation\style2'

%evoked
evoked_trace_Off_f46_zoom = open('f46_traces_x2_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10)
linexdata=get(b1(1),'xdata');
set(b1(1),'xdata',linexdata-70/1000,'linewidth',barWidth);
delete(b1(2,1)); 
delete(b2(2,1));
textxdata=get(b2(1),'position');
textxdata(1)=textxdata(1)-100/1000;
set(b2(1),'position',textxdata);
rec=findobj(gca,'type','rectangle'); 
delete(rec);
evoked_trace_Off_f46_zoom_ax = get(gcf, 'children');

evoked_trace_On_f46_zoom = open('f46_traces_x3_v2.fig');  
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10)
linexdata=get(b1(1),'xdata');
set(b1(1),'xdata',linexdata-70/1000,'linewidth',barWidth);
delete(b1(2,1)); 
delete(b2(2,1));
textxdata=get(b2(1),'position');
textxdata(1)=textxdata(1)-100/1000;
set(b2(1),'position',textxdata);
rec=findobj(gca,'type','rectangle'); 
delete(rec);
evoked_trace_On_f46_zoom_ax = get(gcf, 'children');

evoked_mean_f46_zoom = open('f46_mean_x2+3_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
set(b2(2),'string','1 s');
linexdata=get(b1(1),'xdata');
set(b1(1),'xdata',linexdata-70/1000,'linewidth',barWidth);
lineydata=get(b1(2),'ydata');
set(b1(2),'ydata',lineydata-0.5,'linewidth',barWidth);
textxdata=get(b2(1),'position');
textxdata(1)=textxdata(1)-100/1000;
set(b2(1),'position',textxdata);
rec=findobj(gca,'type','rectangle'); 
delete(rec);
evoked_mean_f46_zoom_ax = get(gcf, 'children');

evoked_std_f46_zoom = open('f46_std_x2+3_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10)
linexdata=get(b1(1),'xdata');
set(b1(1),'xdata',linexdata-70/1000,'linewidth',barWidth);
delete(b1(2,1)); 
delete(b2(2,1));
textxdata=get(b2(1),'position');
textxdata(1)=textxdata(1)-100/1000;
set(b2(1),'position',textxdata);
rec=findobj(gca,'type','rectangle'); 
delete(rec);
evoked_std_f46_zoom_ax = get(gcf, 'children');

f46_WS = open('f46_WS alone.fig');    
set(gca,'visible','off');
f46_WS_ax = get(gcf, 'children');

%ChAT:

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Zoom-in Trace Presentation\style2'

%evoked
evoked_trace_Off_f80_zoom = open('f80_traces_x2_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10)
linexdata=get(b1(1),'xdata');
set(b1(1),'xdata',linexdata-70/1000,'linewidth',barWidth);
delete(b1(2,1)); 
delete(b2(2,1));
textxdata=get(b2(1),'position');
textxdata(1)=textxdata(1)-100/1000;
set(b2(1),'position',textxdata);
rec=findobj(gca,'type','rectangle'); 
delete(rec);
evoked_trace_Off_f80_zoom_ax = get(gcf, 'children');

evoked_trace_On_f80_zoom = open('f80_traces_x3_v2.fig');  
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
linexdata=get(b1(1),'xdata');
set(b1(1),'xdata',linexdata-70/1000,'linewidth',barWidth);
delete(b1(2,1)); 
delete(b2(2,1));
textxdata=get(b2(1),'position');
textxdata(1)=textxdata(1)-100/1000;
set(b2(1),'position',textxdata);
rec=findobj(gca,'type','rectangle'); 
delete(rec);
evoked_trace_On_f80_zoom_ax = get(gcf, 'children');

evoked_mean_f80_zoom = open('f80_mean_x2+3_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
set(b2(2),'string','1 s');
linexdata=get(b1(1),'xdata');
set(b1(1),'xdata',linexdata-70/1000,'linewidth',barWidth);
lineydata=get(b1(2),'ydata');
set(b1(2),'ydata',lineydata-0.2,'linewidth',barWidth);
textxdata=get(b2(1),'position');
textxdata(1)=textxdata(1)-100/1000;
set(b2(1),'position',textxdata);
rec=findobj(gca,'type','rectangle'); 
delete(rec);
evoked_mean_f80_zoom_ax = get(gcf, 'children');

evoked_std_f80_zoom = open('f80_std_x2+3_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
linexdata=get(b1(1),'xdata');
set(b1(1),'xdata',linexdata-70/1000,'linewidth',barWidth);
delete(b1(2,1)); 
delete(b2(2,1));
textxdata=get(b2(1),'position');
textxdata(1)=textxdata(1)-100/1000;
set(b2(1),'position',textxdata);
rec=findobj(gca,'type','rectangle'); 
delete(rec);
evoked_std_f80_zoom_ax = get(gcf, 'children');

f80_WS = open('f80_WS alone.fig');    
set(gca,'visible','off');
f80_WS_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
NB_illustration = imread('NBES_Schematic_Illustration_one-electrode_tight','TIF');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
ChAT_illustration = imread('ChAT_Schematic_Illustration_one-electrode_tight','TIF');
%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',11);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 24]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

%% xlimits, y limits, ticks etc.
 
%% Positions:
    
v_dist=0.06;
evoked_trace_Off_f46_zoom_pos(1,:) = [0.03 , 0.6 ,  0.45 ,  0.21];
evoked_trace_On_f46_zoom_pos(1,:) = [evoked_trace_Off_f46_zoom_pos(1,1), 0.5, evoked_trace_Off_f46_zoom_pos(1,3),evoked_trace_Off_f46_zoom_pos(1,4)]; 
evoked_trace_On_f46_zoom_pos(1,2)=evoked_trace_Off_f46_zoom_pos(1,2)-evoked_trace_On_f46_zoom_pos(1,4)-0.03;
evoked_mean_f46_zoom_pos(1,:) = [evoked_trace_Off_f46_zoom_pos(1,1), evoked_trace_On_f46_zoom_pos(1,2)-0.11-v_dist, evoked_trace_Off_f46_zoom_pos(1,3), 0.11];
evoked_std_f46_zoom_pos(1,:) = [evoked_trace_Off_f46_zoom_pos(1,1), evoked_mean_f46_zoom_pos(1,2)-evoked_mean_f46_zoom_pos(1,4)-v_dist, evoked_trace_Off_f46_zoom_pos(1,3), evoked_mean_f46_zoom_pos(1,4)];
f46_WS_pos(1,:) = [evoked_trace_Off_f46_zoom_pos(1,1),evoked_trace_Off_f46_zoom_pos(1,2)+evoked_trace_Off_f46_zoom_pos(1,4)+0.01,evoked_trace_Off_f46_zoom_pos(1,3),0.03];
%ChAT
h_dist1=evoked_trace_Off_f46_zoom_pos(1,1)+evoked_trace_Off_f46_zoom_pos(1,3)+0.06;

evoked_trace_Off_f80_zoom_pos(1,:) = [h_dist1, evoked_trace_Off_f46_zoom_pos(1,2) ,  evoked_trace_Off_f46_zoom_pos(1,3) ,  evoked_trace_Off_f46_zoom_pos(1,4)];
evoked_trace_On_f80_zoom_pos(1,:) = [h_dist1, evoked_trace_On_f46_zoom_pos(1,2) ,  evoked_trace_Off_f80_zoom_pos(1,3) ,  evoked_trace_Off_f80_zoom_pos(1,4)+0.01];
evoked_mean_f80_zoom_pos(1,:) = [h_dist1, evoked_mean_f46_zoom_pos(1,2) ,   evoked_trace_Off_f80_zoom_pos(1,3) ,   evoked_mean_f46_zoom_pos(1,4)+0.01];
evoked_std_f80_zoom_pos(1,:) = [h_dist1, evoked_std_f46_zoom_pos(1,2) ,  evoked_trace_Off_f80_zoom_pos(1,3) ,  evoked_std_f46_zoom_pos(1,4)+0.01];
f80_WS_pos(1,:) = [evoked_trace_Off_f80_zoom_pos(1,1),evoked_trace_Off_f80_zoom_pos(1,2)+evoked_trace_Off_f80_zoom_pos(1,4)+0.01,evoked_trace_Off_f80_zoom_pos(1,3),0.03];

%position of top of each panel
evoked_trace_Off_f46_zoom_pos_top = evoked_trace_Off_f46_zoom_pos(1,2)+evoked_trace_Off_f46_zoom_pos(1,4);
evoked_trace_On_f46_zoom_pos_top = evoked_trace_On_f46_zoom_pos(1,2)+evoked_trace_On_f46_zoom_pos(1,4);
evoked_mean_f46_zoom_pos_top = evoked_mean_f46_zoom_pos(1,2)+evoked_mean_f46_zoom_pos(1,4);
evoked_std_f46_zoom_pos_top = evoked_std_f46_zoom_pos(1,2)+evoked_std_f46_zoom_pos(1,4);
evoked_trace_Off_f80_zoom_pos_top = evoked_trace_Off_f80_zoom_pos(1,2)+evoked_trace_Off_f80_zoom_pos(1,4);
evoked_trace_On_f80_zoom_pos_top = evoked_trace_On_f80_zoom_pos(1,2)+evoked_trace_On_f80_zoom_pos(1,4);
evoked_mean_f80_zoom_pos_top = evoked_mean_f80_zoom_pos(1,2)+evoked_mean_f80_zoom_pos(1,4);
evoked_std_f80_zoom_pos_top = evoked_std_f80_zoom_pos(1,2)+evoked_std_f80_zoom_pos(1,4);
f46_WS_pos_top = f46_WS_pos(1,2)+f46_WS_pos(1,4);
f80_WS_pos_top = f80_WS_pos(1,2)+f80_WS_pos(1,4);

NB_illustration_pos(1,:)=[evoked_trace_Off_f46_zoom_pos(1,1)+0.12,evoked_trace_Off_f46_zoom_pos_top+0.07, 0.11, 0.12];
ChAT_illustration_pos(1,:)=[evoked_trace_Off_f80_zoom_pos(1,1)+0.14,evoked_trace_Off_f80_zoom_pos_top+0.07, NB_illustration_pos(1,3), NB_illustration_pos(1,4)];

NB_illustration_pos_top=NB_illustration_pos(1,2)+NB_illustration_pos(1,4);
ChAT_illustration_pos_top=ChAT_illustration_pos(1,2)+ChAT_illustration_pos(1,4);
%%
%Placing plots in the figure:

NB_illustration_ax = axes('position',NB_illustration_pos);
imshow(NB_illustration, 'parent', NB_illustration_ax) 

ChAT_illustration_ax = axes('position',ChAT_illustration_pos);
imshow(ChAT_illustration, 'parent', ChAT_illustration_ax) 

evoked_trace_Off_f46_zoom_ax_copy = copyobj(evoked_trace_Off_f46_zoom_ax,F); % copy axes to new fig
set(evoked_trace_Off_f46_zoom_ax_copy,'position',evoked_trace_Off_f46_zoom_pos(1,:))

f46_WS_ax_copy = copyobj(f46_WS_ax,F); % copy axes to new fig
set(f46_WS_ax_copy,'position',f46_WS_pos(1,:))

evoked_trace_On_f46_zoom_ax_copy = copyobj(evoked_trace_On_f46_zoom_ax,F); % copy axes to new fig
set(evoked_trace_On_f46_zoom_ax_copy,'position',evoked_trace_On_f46_zoom_pos(1,:))

evoked_mean_f46_zoom_ax_copy = copyobj(evoked_mean_f46_zoom_ax,F); % copy axes to new fig
set(evoked_mean_f46_zoom_ax_copy,'position',evoked_mean_f46_zoom_pos(1,:))

evoked_std_f46_zoom_ax_copy = copyobj(evoked_std_f46_zoom_ax,F); % copy axes to new fig
set(evoked_std_f46_zoom_ax_copy,'position',evoked_std_f46_zoom_pos(1,:))

%position legend:
set(evoked_mean_f46_zoom_ax_copy(1),'position',[0.055 evoked_trace_Off_f46_zoom_pos_top+0.06 0.08 0.005])
evoked_mean_f46_zoom_ax_copy(1).FontSize=9;
% set(evoked_mean_f46_zoom_ax_copy(1),'position',[0.9 spont_std_f46_pos_top+0.01 0.08 0.05])

%ChAT

evoked_trace_Off_f80_zoom_ax_copy = copyobj(evoked_trace_Off_f80_zoom_ax,F); % copy axes to new fig
set(evoked_trace_Off_f80_zoom_ax_copy,'position',evoked_trace_Off_f80_zoom_pos(1,:))

f80_WS_ax_copy = copyobj(f80_WS_ax,F); % copy axes to new fig
set(f80_WS_ax_copy,'position',f80_WS_pos(1,:))

evoked_trace_On_f80_zoom_ax_copy = copyobj(evoked_trace_On_f80_zoom_ax,F); % copy axes to new fig
set(evoked_trace_On_f80_zoom_ax_copy,'position',evoked_trace_On_f80_zoom_pos(1,:))

evoked_mean_f80_zoom_ax_copy = copyobj(evoked_mean_f80_zoom_ax,F); % copy axes to new fig
set(evoked_mean_f80_zoom_ax_copy,'position',evoked_mean_f80_zoom_pos(1,:))

evoked_std_f80_zoom_ax_copy = copyobj(evoked_std_f80_zoom_ax,F); % copy axes to new fig
set(evoked_std_f80_zoom_ax_copy,'position',evoked_std_f80_zoom_pos(1,:))

%position legend:
set(evoked_mean_f80_zoom_ax_copy(1),'position',[0.84 evoked_trace_Off_f46_zoom_pos_top+0.06 0.08 0.005])
evoked_mean_f80_zoom_ax_copy(1).FontSize=9;

 a_pos1=[-0.03 -0.01 0.04 0.04];
 a_pos4=[0.4 0 0.16 0.04]; %[0 -0.01 0.16 0.04];
 a_pos2=[-0.03 -0.02 0.04 0.04]; %[-0.02 0.01 0.04 0.04]
 a_pos3=[0.03 -0.01 0.04 0.04];

 annotation('textbox', [0.27, NB_illustration_pos_top-0.05, 0.6 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Sensory evoked synaptic responses', 'FontName','arial', 'fontsize', 12); %, 'fontweight', 'bold'
 
  annotation('textbox', [f46_WS_pos(1,1)+0.14 f46_WS_pos_top+0.02 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'WS', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
 annotation('textbox', [f80_WS_pos(1,1)+0.14 f80_WS_pos_top+0.02 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'WS', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
 annotation('textbox', [evoked_trace_Off_f46_zoom_pos(1,1)-0.025 evoked_trace_Off_f46_zoom_pos_top+0.04 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Single trials', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
annotation('textbox', [evoked_mean_f46_zoom_pos(1,1) evoked_mean_f46_zoom_pos_top+0.01 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [evoked_std_f46_zoom_pos(1,1)+0.005, evoked_std_f46_zoom_pos_top+0.01 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 12)

 if no_numbering_flag==0;

annotation('textbox', [evoked_trace_Off_f46_zoom_pos(1,1)-0.03 evoked_trace_Off_f46_zoom_pos_top+0.04 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [evoked_trace_On_f46_zoom_pos(1,1)-0.03 evoked_trace_On_f46_zoom_pos_top+0.03 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [evoked_mean_f46_zoom_pos(1,1)-0.03 evoked_mean_f46_zoom_pos_top+0.03 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [evoked_std_f46_zoom_pos(1,1)-0.03 evoked_std_f46_zoom_pos_top+0.03 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [evoked_trace_Off_f80_zoom_pos(1,1)-0.03 evoked_trace_Off_f80_zoom_pos_top+0.04 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [evoked_trace_On_f80_zoom_pos(1,1)-0.03 evoked_trace_On_f46_zoom_pos_top+0.03 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [evoked_mean_f80_zoom_pos(1,1)-0.03 evoked_mean_f46_zoom_pos_top+0.03 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [evoked_std_f80_zoom_pos(1,1)-0.03 evoked_std_f46_zoom_pos_top+0.03 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
 end
  
      cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'
 
  if save_flag==1
    filename='Fig5_f46_f80';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print -depsc2 Fig5_f46_f80
%     print(F, '-depsc2', filename);
end
