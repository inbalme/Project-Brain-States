%% plotting a figure 
close all
clear all
save_flag=1;
 no_numbering_flag=1;
%opening saved figures:
%Long traces:
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Long Trace Presentation'
%spontaneous:
spont_trace_f16 = open('f16_traces_x1_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
spont_trace_f16_ax = get(gcf, 'children');

spont_mean_f16 = open('f16_mean_x1_v2.fig');    
spont_mean_f16_ax = get(gcf, 'children');

spont_std_f16 = open('f16_std_x1_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
spont_std_f16_ax = get(gcf, 'children');

%%
%Short traces:
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Zoom-in Trace Presentation'

%evoked
evoked_trace_Off_f16_zoom = open('f16_traces_x2_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
evoked_trace_Off_f16_zoom_ax = get(gcf, 'children');

evoked_trace_On_f16_zoom = open('f16_traces_x3_v2.fig');  
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
evoked_trace_On_f16_zoom_ax = get(gcf, 'children');

evoked_mean_f16_zoom = open('f16_mean_x2+3_v2.fig');    
evoked_mean_f16_zoom_ax = get(gcf, 'children');

evoked_std_f16_zoom = open('f16_std_x2+3_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
evoked_std_f16_zoom_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\LFP_50Hz Vm_50Hz'

VmM = open('Vm M_Before During After Sensory stim.fig');    
VmM_ax = get(gcf, 'children');

VmSTD = open('Vm STD_Before During After Sensory stim.fig');    
VmSTD_ax = get(gcf, 'children');

SNR = open('SNR1.fig');    
SNR_ax = get(gcf, 'children');

Amplitude_Signal = open('Amplitude_Signal.fig');    
Amplitude_Signal_ax = get(gcf, 'children');

Amplitude_Noise = open('Amplitude_Noise1.fig');    
Amplitude_Noise_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Evoked';
% Peak_amp = open('Train Amplitude Local_v2.fig');    
Peak_amp = open('Evoked Amplitude.fig');    
Peak_amp_ax = get(gcf, 'children');

Hist_amp = open('Amp_Diff_Hist.fig'); 
color2=[51,102,0]./256; %green
color4=[204, 102,0]./256; %light brown
c2=findall(gcf,'type','bar');
        set(c2,'facecolor',color4);
Hist_amp_ax = get(gcf, 'children');


%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',11);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 23 29]); %[left, bottom, width, height] 
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
    
v_dist=0.04;
spont_trace_f16_pos(1,:) = [0.04 , 0.72 ,  0.44 ,  0.20];
spont_mean_f16_pos(1,:) = [spont_trace_f16_pos(1,1), spont_trace_f16_pos(1,2)-0.10-v_dist, spont_trace_f16_pos(1,3), 0.10];
spont_std_f16_pos(1,:) = [spont_trace_f16_pos(1,1), spont_mean_f16_pos(1,2)-spont_mean_f16_pos(1,4)-v_dist, spont_trace_f16_pos(1,3), spont_mean_f16_pos(1,4)];

h_dist1=spont_trace_f16_pos(1,1)+spont_trace_f16_pos(1,3)+0.04;

evoked_trace_Off_f16_zoom_pos(1,:) = [h_dist1, spont_trace_f16_pos(1,2)+0.1 ,  0.48 ,  0.12];
evoked_trace_On_f16_zoom_pos(1,:) = [h_dist1, evoked_trace_Off_f16_zoom_pos(1,2)-evoked_trace_Off_f16_zoom_pos(1,4)-0.02 ,  evoked_trace_Off_f16_zoom_pos(1,3) ,  evoked_trace_Off_f16_zoom_pos(1,4)];
evoked_mean_f16_zoom_pos(1,:) = [h_dist1, evoked_trace_On_f16_zoom_pos(1,2)-0.10-0.03 ,   evoked_trace_Off_f16_zoom_pos(1,3) ,   0.10];
evoked_std_f16_zoom_pos(1,:) = [h_dist1, evoked_mean_f16_zoom_pos(1,2)-evoked_mean_f16_zoom_pos(1,4)-0.03 ,  evoked_trace_Off_f16_zoom_pos(1,3) ,  evoked_mean_f16_zoom_pos(1,4)];

v_dist2=0.05;
VmM_pos(1,:)=[0.08, 0.5, 0.14, 0.10]; VmM_pos(1,2)=evoked_std_f16_zoom_pos(1,2)-VmM_pos(1,4)-v_dist2;
VmSTD_pos(1,:)= VmM_pos; VmSTD_pos(1,1)=VmM_pos(1,1)+VmM_pos(1,3)+0.06;
SNR_pos(1,:)= VmSTD_pos; SNR_pos(1,1)=VmSTD_pos(1,1)+VmSTD_pos(1,3)+0.08; SNR_pos(1,3)=0.12;
Amplitude_Signal_pos(1,:)=SNR_pos; Amplitude_Signal_pos(1,1)=SNR_pos(1,1)+SNR_pos(1,3)+0.06;
Amplitude_Noise_pos(1,:)=Amplitude_Signal_pos; Amplitude_Noise_pos(1,1)=Amplitude_Signal_pos(1,1)+Amplitude_Signal_pos(1,3)+0.06;

Peak_amp_pos(1,:)=[VmM_pos(1,1)+0.02,0.5,0.3,0.12];  Peak_amp_pos(1,2)=VmM_pos(1,2)-Peak_amp_pos(1,4)-0.08;
Hist_amp_pos(1,:)=Peak_amp_pos(1,:); Hist_amp_pos(1,1)=Peak_amp_pos(1,1)+Peak_amp_pos(1,3)+0.15;
%position of top of each panel
spont_trace_f16_pos_top = spont_trace_f16_pos(1,2)+spont_trace_f16_pos(1,4);
spont_mean_f16_pos_top = spont_mean_f16_pos(1,2)+spont_mean_f16_pos(1,4);
spont_std_f16_pos_top = spont_std_f16_pos(1,2)+spont_std_f16_pos(1,4);
evoked_trace_Off_f16_zoom_pos_top = evoked_trace_Off_f16_zoom_pos(1,2)+evoked_trace_Off_f16_zoom_pos(1,4);
evoked_trace_On_f16_zoom_pos_top = evoked_trace_On_f16_zoom_pos(1,2)+evoked_trace_On_f16_zoom_pos(1,4);
evoked_mean_f16_zoom_pos_top = evoked_mean_f16_zoom_pos(1,2)+evoked_mean_f16_zoom_pos(1,4);
evoked_std_f16_zoom_pos_top = evoked_std_f16_zoom_pos(1,2)+evoked_std_f16_zoom_pos(1,4);
VmM_pos_top = VmM_pos(1,2)+VmM_pos(1,4);
VmSTD_pos_top = VmSTD_pos(1,2)+VmSTD_pos(1,4);
Peak_amp_pos_top = Peak_amp_pos(1,2)+Peak_amp_pos(1,4);
Hist_amp_pos_top = Hist_amp_pos(1,2)+Hist_amp_pos(1,4);
SNR_pos_top = SNR_pos(1,2)+SNR_pos(1,4);
Amplitude_Signal_pos_top = Amplitude_Signal_pos(1,2)+Amplitude_Signal_pos(1,4);
Amplitude_Noise_pos_top = Amplitude_Noise_pos(1,2)+Amplitude_Noise_pos(1,4);
%%
%Placing plots in the figure:
%Cell 16 - spont

spont_trace_f16_ax_copy = copyobj(spont_trace_f16_ax,F); % copy axes to new fig
set(spont_trace_f16_ax_copy,'position',spont_trace_f16_pos(1,:))
   %    'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:));
spont_mean_f16_ax_copy = copyobj(spont_mean_f16_ax,F); % copy axes to new fig
set(spont_mean_f16_ax_copy,'position',spont_mean_f16_pos(1,:))

spont_std_f16_ax_copy = copyobj(spont_std_f16_ax,F); % copy axes to new fig
set(spont_std_f16_ax_copy,'position',spont_std_f16_pos(1,:))

%% zoom-in traces - evoked
evoked_trace_Off_f16_zoom_ax_copy = copyobj(evoked_trace_Off_f16_zoom_ax,F); % copy axes to new fig
set(evoked_trace_Off_f16_zoom_ax_copy,'position',evoked_trace_Off_f16_zoom_pos(1,:))

evoked_trace_On_f16_zoom_ax_copy = copyobj(evoked_trace_On_f16_zoom_ax,F); % copy axes to new fig
set(evoked_trace_On_f16_zoom_ax_copy,'position',evoked_trace_On_f16_zoom_pos(1,:))

evoked_mean_f16_zoom_ax_copy = copyobj(evoked_mean_f16_zoom_ax,F); % copy axes to new fig
set(evoked_mean_f16_zoom_ax_copy,'position',evoked_mean_f16_zoom_pos(1,:))

evoked_std_f16_zoom_ax_copy = copyobj(evoked_std_f16_zoom_ax,F); % copy axes to new fig
set(evoked_std_f16_zoom_ax_copy,'position',evoked_std_f16_zoom_pos(1,:))

%population panels:
VmM_ax_copy = copyobj(VmM_ax,F); % copy axes to new fig
set(VmM_ax_copy(2),'position',VmM_pos(1,:))
set(F, 'currentaxes', VmM_ax_copy(2));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel('sensory stimulation','fontsize',10); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
VmM_ax_copy(2).FontSize=11;
VmM_ax_copy(2).XTickLabel=[];VmM_ax_copy(2).XTickLabel(1,:)='Bef.';VmM_ax_copy(2).XTickLabel(2,:)='Dur.';VmM_ax_copy(2).XTickLabel(3,:)='Aft.';

VmSTD_ax_copy = copyobj(VmSTD_ax(2),F); % copy axes to new fig
set(VmSTD_ax_copy,'position',VmSTD_pos(1,:))
set(F, 'currentaxes', VmSTD_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('sensory stimulation','fontsize',10);%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
VmSTD_ax_copy.FontSize=11;
VmSTD_ax_copy.XTickLabel=[];VmSTD_ax_copy.XTickLabel(1,:)='Bef.';VmSTD_ax_copy.XTickLabel(2,:)='Dur.';VmSTD_ax_copy.XTickLabel(3,:)='Aft.';


Peak_amp_ax_copy = copyobj(Peak_amp_ax,F); % copy axes to new fig
set(Peak_amp_ax_copy,'position',Peak_amp_pos(1,:))
set(F, 'currentaxes', Peak_amp_ax_copy); tl=title(''); yl=ylabel('Amp. [mV]');
Peak_amp_ax_copy.FontSize=11;
Peak_amp_ax_copy.XTickLabel(11,:)=[]; Peak_amp_ax_copy.XTickLabel(11,1)='R';

Hist_amp_ax_copy = copyobj(Hist_amp_ax,F); % copy axes to new fig
set(Hist_amp_ax_copy,'position',Hist_amp_pos(1,:))
set(F, 'currentaxes', Hist_amp_ax_copy); tl=title('');
% Hist_amp_ax_copy.FontSize=12;

SNR_ax_copy = copyobj(SNR_ax,F); % copy axes to new fig
set(SNR_ax_copy,'position',SNR_pos(1,:))
set(F, 'currentaxes', SNR_ax_copy);  tl=title('');  yl=ylabel('A.U.','fontsize',11); %tl=title('SNR','fontsize',11,'fontweight','normal');
SNR_ax_copy.FontSize=11;
% SNR_ax_copy.XTickLabel=[];


Amplitude_Signal_ax_copy = copyobj(Amplitude_Signal_ax,F); % copy axes to new fig
set(Amplitude_Signal_ax_copy,'position',Amplitude_Signal_pos(1,:))
set(F, 'currentaxes', Amplitude_Signal_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11); %tl=title('Signal Amp.','fontsize',11,'fontweight','normal');
Amplitude_Signal_ax_copy.FontSize=11;
% Amplitude_Signal_ax_copy.XTickLabel=[];

Amplitude_Noise_ax_copy = copyobj(Amplitude_Noise_ax,F); % copy axes to new fig
set(Amplitude_Noise_ax_copy,'position',Amplitude_Noise_pos(1,:))
set(F, 'currentaxes', Amplitude_Noise_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11);  %tl=title('Noise Amp.','fontsize',11,'fontweight','normal');
Amplitude_Noise_ax_copy.FontSize=11;

%position legend:
set(VmM_ax_copy(1),'position',[0.9 evoked_trace_Off_f16_zoom_pos_top+0.03 0.08 0.005])
VmM_ax_copy(1).FontSize=10;
% set(evoked_mean_f16_zoom_ax_copy(1),'position',[0.9 spont_std_f16_pos_top+0.01 0.08 0.05])
%   annotation:
annotation('textbox', [evoked_trace_Off_f16_zoom_pos(1,1) evoked_trace_Off_f16_zoom_pos_top 0 0]+[0.03 0.01 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Sensory-evoked Responses', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
annotation('textbox', [spont_trace_f16_pos(1,1), evoked_trace_Off_f16_zoom_pos_top 0 0]+[0.06 0.01 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous Activity', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
 
 a_pos1=[-0.03 -0.01 0.04 0.04];
 a_pos4=[0 -0.01 0.16 0.04];
 a_pos2=[-0.02 -0.01 0.04 0.04];
 a_pos3=[-0.035 0 0.04 0.04];
 
 annotation('textbox', [spont_trace_f16_pos(1,1) spont_trace_f16_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Vm Traces', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
annotation('textbox', [spont_mean_f16_pos(1,1) spont_mean_f16_pos_top 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [spont_std_f16_pos(1,1) spont_std_f16_pos_top 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [VmM_pos(1,1)+0.02 VmM_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [VmSTD_pos(1,1)+0.025 VmSTD_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [SNR_pos(1,1)+0.03 SNR_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'SNR', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [Amplitude_Signal_pos(1,1) Amplitude_Signal_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Signal Amp.', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [Amplitude_Noise_pos(1,1) Amplitude_Noise_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Noise Amp.', 'FontName','arial', 'fontsize', 12)

 if no_numbering_flag==0;
annotation('textbox', [spont_trace_f16_pos(1,1) spont_trace_f16_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [spont_mean_f16_pos(1,1) spont_mean_f16_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [spont_std_f16_pos(1,1) spont_std_f16_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')

annotation('textbox', [evoked_trace_Off_f16_zoom_pos(1,1) evoked_trace_Off_f16_zoom_pos_top 0 0]+a_pos2,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [evoked_trace_On_f16_zoom_pos(1,1) evoked_trace_On_f16_zoom_pos_top 0 0]+a_pos2,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [evoked_mean_f16_zoom_pos(1,1) evoked_mean_f16_zoom_pos_top 0 0]+a_pos2,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [evoked_std_f16_zoom_pos(1,1) evoked_std_f16_zoom_pos_top 0 0]+a_pos2,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [VmM_pos(1,1) VmM_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [VmSTD_pos(1,1) VmM_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')

annotation('textbox', [SNR_pos(1,1) VmM_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [Amplitude_Signal_pos(1,1) VmM_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'K', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [Amplitude_Noise_pos(1,1) VmM_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'L', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [Peak_amp_pos(1,1)-0.06 Peak_amp_pos_top 0.04 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'M', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [Hist_amp_pos(1,1)-0.06 Peak_amp_pos_top 0.04 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'N', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 end
  if no_numbering_flag==1;
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\No Numbering'
  else
      cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures'
  end
if save_flag==1
    filename='Fig 3 intracellular traces+SNR_f16_v4';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end
