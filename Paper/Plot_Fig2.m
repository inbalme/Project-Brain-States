%% plotting a figure 
close all
clear all
save_flag=1;
 no_numbering_flag=1;
%opening saved figures:
%Long traces NB:
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Long Trace Presentation\style2'
%spontaneous:
spont_trace_f46 = open('f46_traces_x1_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
spont_trace_f46_ax = get(gcf, 'children');

spont_mean_f46 = open('f46_mean_x1_v2.fig');    
spont_mean_f46_ax = get(gcf, 'children');

spont_std_f46 = open('f46_std_x1_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
spont_std_f46_ax = get(gcf, 'children');

%% long traces ChAT

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Long Trace Presentation\style2'
%spontaneous:
spont_trace_f80 = open('f80_traces_x1_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
spont_trace_f80_ax = get(gcf, 'children');

spont_mean_f80 = open('f80_mean_x1_v2.fig');    
spont_mean_f80_ax = get(gcf, 'children');

spont_std_f80 = open('f80_std_x1_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
spont_std_f80_ax = get(gcf, 'children');


cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous';

VmM_NB = open('Spontaneous VmM.fig');   
VmM_NB_ax = get(gcf, 'children');

 
VmSTD_NB = open('Spontaneous VmSTD.fig');   
VmSTD_NB_ax = get(gcf, 'children');  

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous';
VmM_ChAT = open('Spontaneous VmM.fig');    
VmM_ChAT_ax = get(gcf, 'children');

VmSTD_ChAT = open('Spontaneous VmSTD.fig');    
VmSTD_ChAT_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'

ChAT_illustration = imread('ChAT_Schematic_Illustration_tight','TIF');

ChAT_histology = imread('ChAT_Histology_tight','TIF');


%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',11);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 20]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
% annotation('textbox', [freq_1_zoom1_pos(1,1) freq_1_zoom1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
%            'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A',...
%            'FontName','arial', 'fontsize', 11, 'fontweight', 'bold')
% annotation('textbox', [freq_1_zoom2_pos(1,1) freq_1_zoom2_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
%            'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B',...
%            'FontName','arial', 'fontsize', 11, 'fontweight', 'bold')    
% annotation('textbox', [freq_2_zoom1_pos(1,1) freq_2_zoom1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
%            'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C',...
%            'FontName','arial', 'fontsize', 11, 'fontweight', 'bold')  
%% xlimits, y limits, ticks etc.
 
%% Positions:
    
v_dist=0.01;
spont_trace_f46_pos(1,:) = [0.04 , 0.73 ,  0.6 ,  0.2];
spont_mean_f46_pos(1,:) = [spont_trace_f46_pos(1,1), spont_trace_f46_pos(1,2)-0.10-v_dist, spont_trace_f46_pos(1,3), 0.10];
spont_std_f46_pos(1,:) = [spont_trace_f46_pos(1,1), spont_mean_f46_pos(1,2)-spont_mean_f46_pos(1,4)-v_dist, spont_trace_f46_pos(1,3), spont_mean_f46_pos(1,4)];

spont_trace_f80_pos(1,:) =[spont_trace_f46_pos(1,1) , 0.5 ,  spont_trace_f46_pos(1,3) ,  spont_trace_f46_pos(1,4)]; spont_trace_f80_pos(1,2)=spont_std_f46_pos(1,2)-spont_trace_f80_pos(1,4)-0.05;
spont_mean_f80_pos(1,:) =[spont_trace_f80_pos(1,1), spont_trace_f80_pos(1,2)-0.10-v_dist, spont_trace_f80_pos(1,3), 0.10];
spont_std_f80_pos(1,:) =[spont_trace_f80_pos(1,1), spont_mean_f80_pos(1,2)-spont_mean_f80_pos(1,4)-v_dist, spont_trace_f80_pos(1,3), spont_mean_f80_pos(1,4)];

h_dist1=spont_trace_f46_pos(1,1)+spont_trace_f46_pos(1,3)+0.07;

v_dist2=0.05;
VmM_NB_pos(1,:)=[h_dist1+0.03, spont_trace_f46_pos(1,2)+0.02, 0.07, 0.13]; 
VmSTD_NB_pos(1,:)= VmM_NB_pos; VmSTD_NB_pos(1,1)=VmM_NB_pos(1,1)+VmM_NB_pos(1,3)+0.06;
VmM_ChAT_pos(1,:)=[h_dist1+0.03, spont_std_f80_pos(1,2)+0.02, VmM_NB_pos(1,3), VmM_NB_pos(1,4)]; 
VmSTD_ChAT_pos(1,:)= VmM_ChAT_pos; VmSTD_ChAT_pos(1,1)=VmM_ChAT_pos(1,1)+VmM_ChAT_pos(1,3)+0.06;

ChAT_illustration_pos(1,:)=[h_dist1+0.03,VmM_NB_pos(1,2)-0.2, 0.15, 0.15];
ChAT_histology_pos(1,:)=[h_dist1-0.04,ChAT_illustration_pos(1,2)-0.21, 0.3, 0.2];

%position of top of each panel
spont_trace_f46_pos_top = spont_trace_f46_pos(1,2)+spont_trace_f46_pos(1,4);
spont_mean_f46_pos_top = spont_mean_f46_pos(1,2)+spont_mean_f46_pos(1,4);
spont_std_f46_pos_top = spont_std_f46_pos(1,2)+spont_std_f46_pos(1,4);
spont_trace_f80_pos_top = spont_trace_f80_pos(1,2)+spont_trace_f80_pos(1,4);
spont_mean_f80_pos_top = spont_mean_f80_pos(1,2)+spont_mean_f80_pos(1,4);
spont_std_f80_pos_top = spont_std_f80_pos(1,2)+spont_std_f80_pos(1,4);
VmM_NB_pos_top = VmM_NB_pos(1,2)+VmM_NB_pos(1,4);
VmSTD_NB_pos_top = VmSTD_NB_pos(1,2)+VmSTD_NB_pos(1,4);
VmM_ChAT_pos_top = VmM_ChAT_pos(1,2)+VmM_ChAT_pos(1,4);
VmSTD_ChAT_pos_top = VmSTD_ChAT_pos(1,2)+VmSTD_ChAT_pos(1,4);
%%
%Placing plots in the figure:
%Cell 44 - spont

spont_trace_f46_ax_copy = copyobj(spont_trace_f46_ax,F); % copy axes to new fig
set(spont_trace_f46_ax_copy,'position',spont_trace_f46_pos(1,:))
   %    'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:));
spont_mean_f46_ax_copy = copyobj(spont_mean_f46_ax,F); % copy axes to new fig
set(spont_mean_f46_ax_copy(2),'position',spont_mean_f46_pos(1,:))

spont_std_f46_ax_copy = copyobj(spont_std_f46_ax,F); % copy axes to new fig
set(spont_std_f46_ax_copy,'position',spont_std_f46_pos(1,:))

spont_trace_f80_ax_copy = copyobj(spont_trace_f80_ax,F); % copy axes to new fig
set(spont_trace_f80_ax_copy,'position',spont_trace_f80_pos(1,:))

spont_mean_f80_ax_copy = copyobj(spont_mean_f80_ax,F); % copy axes to new fig
set(spont_mean_f80_ax_copy(2),'position',spont_mean_f80_pos(1,:))

spont_std_f80_ax_copy = copyobj(spont_std_f80_ax,F); % copy axes to new fig
set(spont_std_f80_ax_copy,'position',spont_std_f80_pos(1,:))


%population panels:
VmM_NB_ax_copy = copyobj(VmM_NB_ax,F); % copy axes to new fig
set(VmM_NB_ax_copy(1),'position',VmM_NB_pos(1,:))
set(F, 'currentaxes', VmM_NB_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
VmM_NB_ax_copy(1).FontSize=11;
% VmM_NB_ax_copy(1).XTickLabel=[];

VmSTD_NB_ax_copy = copyobj(VmSTD_NB_ax(1),F); % copy axes to new fig
set(VmSTD_NB_ax_copy,'position',VmSTD_NB_pos(1,:))
set(F, 'currentaxes', VmSTD_NB_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
VmSTD_NB_ax_copy.FontSize=11;
% VmSTD_NB_ax_copy.XTickLabel=[];

VmM_ChAT_ax_copy = copyobj(VmM_ChAT_ax,F); % copy axes to new fig
set(VmM_ChAT_ax_copy(1),'position',VmM_ChAT_pos(1,:))
set(F, 'currentaxes', VmM_ChAT_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
VmM_ChAT_ax_copy(1).FontSize=11;
VmM_ChAT_ax_copy(1).XTickLabel={'Off','On'};

VmSTD_ChAT_ax_copy = copyobj(VmSTD_ChAT_ax(1),F); % copy axes to new fig
set(VmSTD_ChAT_ax_copy,'position',VmSTD_ChAT_pos(1,:))
set(F, 'currentaxes', VmSTD_ChAT_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
VmSTD_ChAT_ax_copy.FontSize=11;
VmSTD_ChAT_ax_copy.XTickLabel={'Off','On'};

ChAT_illustration_ax = axes('position',ChAT_illustration_pos);
imshow(ChAT_illustration, 'parent', ChAT_illustration_ax) 

ChAT_histology_ax = axes('position',ChAT_histology_pos);
imshow(ChAT_histology, 'parent', ChAT_histology_ax) 

%position legend:
set(spont_mean_f46_ax_copy(1),'position',[0.5 spont_trace_f46_pos_top+0.02 0.08 0.005])
spont_mean_f46_ax_copy(1).FontSize=10;

set(spont_mean_f80_ax_copy(1),'position',[0.5 spont_trace_f80_pos_top 0.08 0.005])
spont_mean_f80_ax_copy(1).FontSize=10;
%   annotation:

annotation('textbox', [spont_trace_f46_pos(1,1), spont_trace_f46_pos_top 0 0]+[0 0.03 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous Activity', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
 
 a_pos1=[-0.03 -0.01 0.04 0.04];
 a_pos4=[0 -0.01 0.16 0.04];
 a_pos2=[-0.02 -0.01 0.04 0.04];
 a_pos3=[-0.035 0 0.04 0.04];
 
 annotation('textbox', [spont_trace_f46_pos(1,1) spont_trace_f46_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Vm Traces', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
annotation('textbox', [spont_mean_f46_pos(1,1) spont_mean_f46_pos_top 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [spont_std_f46_pos(1,1) spont_std_f46_pos_top 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [spont_trace_f80_pos(1,1) spont_trace_f80_pos_top-0.018 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Vm Traces', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
annotation('textbox', [spont_mean_f80_pos(1,1) spont_mean_f80_pos_top-0.01 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [spont_std_f80_pos(1,1) spont_std_f80_pos_top-0.01 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 12)
    
annotation('textbox', [VmM_NB_pos(1,1)-0.02 VmM_NB_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [VmSTD_NB_pos(1,1)-0.02 VmSTD_NB_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [VmM_ChAT_pos(1,1)-0.02 VmM_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [VmSTD_ChAT_pos(1,1)-0.02 VmSTD_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 12)

 if no_numbering_flag==0;
annotation('textbox', [spont_trace_f46_pos(1,1) spont_trace_f46_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [spont_mean_f46_pos(1,1) spont_mean_f46_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [spont_std_f46_pos(1,1) spont_std_f46_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [VmM_NB_pos(1,1) VmM_NB_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [VmSTD_NB_pos(1,1) VmM_NB_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 end
 
      cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'
  
if save_flag==1
    filename='Fig 2 intracellular traces spont_f46_f80';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end
