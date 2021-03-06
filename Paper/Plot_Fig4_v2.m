%% plotting a figure 
close all
clear all
save_flag=0;
 no_numbering_flag=1;
%opening saved figures:
%Short traces:
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Zoom-in Trace Presentation'

%evoked
evoked_trace_Off_f46_zoom = open('f46_traces_x2_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
evoked_trace_Off_f46_zoom_ax = get(gcf, 'children');

evoked_trace_On_f46_zoom = open('f46_traces_x3_v2.fig');  
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
evoked_trace_On_f46_zoom_ax = get(gcf, 'children');

evoked_mean_f46_zoom = open('f46_mean_x2+3_v2.fig');    
evoked_mean_f46_zoom_ax = get(gcf, 'children');

evoked_std_f46_zoom = open('f46_std_x2+3_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
evoked_std_f46_zoom_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\LFP_50Hz Vm_50Hz'

VmM = open('Vm M_During Sensory stim.fig');    
VmM_ax = get(gcf, 'children');

VmSTD = open('Vm STD_During Sensory stim.fig');    
VmSTD_ax = get(gcf, 'children');

SNR = open('SNR1.fig');    
SNR_ax = get(gcf, 'children');

Amplitude_Signal = open('Amplitude_Signal.fig');    
Amplitude_Signal_ax = get(gcf, 'children');

Amplitude_Noise = open('Amplitude_Noise1.fig');    
Amplitude_Noise_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Evoked';
Peak_amp = open('Evoked Amplitude Normalized.fig');    
Peak_amp_ax = get(gcf, 'children');

%ChAT:

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Zoom-in Trace Presentation'

%evoked
evoked_trace_Off_f80_zoom = open('f80_traces_x2_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
evoked_trace_Off_f80_zoom_ax = get(gcf, 'children');

evoked_trace_On_f80_zoom = open('f80_traces_x3_v2.fig');  
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
evoked_trace_On_f80_zoom_ax = get(gcf, 'children');

evoked_mean_f80_zoom = open('f80_mean_x2+3_v2.fig');    
evoked_mean_f80_zoom_ax = get(gcf, 'children');

evoked_std_f80_zoom = open('f80_std_x2+3_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
delete(b1(2,1)); 
delete(b2(2,1));
evoked_std_f80_zoom_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary\LFP_50Hz Vm_50Hz'

VmM_ChAT = open('Vm M_During Sensory stim.fig');    
VmM_ChAT_ax = get(gcf, 'children');

VmSTD_ChAT = open('Vm STD_During Sensory stim.fig');    
VmSTD_ChAT_ax = get(gcf, 'children');

SNR_ChAT = open('SNR1.fig');    
SNR_ChAT_ax = get(gcf, 'children');

Amplitude_Signal_ChAT = open('Amplitude_Signal.fig');    
Amplitude_Signal_ChAT_ax = get(gcf, 'children');

Amplitude_Noise_ChAT = open('Amplitude_Noise1.fig');    
Amplitude_Noise_ChAT_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Evoked';
Peak_amp_ChAT = open('Evoked Amplitude Normalized.fig');    
Peak_amp_ChAT_ax = get(gcf, 'children');

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

%% xlimits, y limits, ticks etc.
 
%% Positions:
    
v_dist=0.04;
evoked_trace_Off_f46_zoom_pos(1,:) = [0.02 , 0.72 ,  0.34 ,  0.20];
evoked_trace_On_f46_zoom_pos(1,:) = [evoked_trace_Off_f46_zoom_pos(1,1), 0.5, evoked_trace_Off_f46_zoom_pos(1,3),evoked_trace_Off_f46_zoom_pos(1,4)]; 
        evoked_trace_On_f46_zoom_pos(1,2)=evoked_trace_Off_f46_zoom_pos(1,2)-evoked_trace_On_f46_zoom_pos(1,4)-v_dist;
evoked_mean_f46_zoom_pos(1,:) = [evoked_trace_Off_f46_zoom_pos(1,1), evoked_trace_On_f46_zoom_pos(1,2)-0.10-v_dist, evoked_trace_Off_f46_zoom_pos(1,3), 0.10];
evoked_std_f46_zoom_pos(1,:) = [evoked_trace_Off_f46_zoom_pos(1,1), evoked_mean_f46_zoom_pos(1,2)-evoked_mean_f46_zoom_pos(1,4)-v_dist, evoked_trace_Off_f46_zoom_pos(1,3), evoked_mean_f46_zoom_pos(1,4)];

% Peak_amp_pos(1,:)=[0.05, 0.5, 0.12, 0.10];  Peak_amp_pos(1,2)=evoked_std_f46_zoom_pos(1,2)-Peak_amp_pos(1,4)-v_dist2;
VmM_pos(1,:)=[0.07, 0.5, 0.06, 0.11]; VmM_pos(1,2)=evoked_std_f46_zoom_pos(1,2)-VmM_pos(1,4)-0.05;
VmSTD_pos(1,:)= VmM_pos; VmSTD_pos(1,1)=VmM_pos(1,1)+VmM_pos(1,3)+0.05;
Peak_amp_pos(1,:)=[0.5,VmM_pos(1,2)+0.03, 0.15, VmM_pos(1,4)-0.03];  Peak_amp_pos(1,1)=VmSTD_pos(1,1)+VmSTD_pos(1,3)+0.08;

h_dist2=evoked_trace_Off_f46_zoom_pos(1,1)+evoked_trace_Off_f46_zoom_pos(1,3)+0.04;
SNR_pos(1,:)= [h_dist2, evoked_trace_Off_f46_zoom_pos(1,2)+0.03,0.08, 0.12]; 
Amplitude_Signal_pos(1,:)=SNR_pos; Amplitude_Signal_pos(1,2)=SNR_pos(1,2)-Amplitude_Signal_pos(1,4)-0.1;
Amplitude_Noise_pos(1,:)=Amplitude_Signal_pos; Amplitude_Noise_pos(1,2)=Amplitude_Signal_pos(1,2)-Amplitude_Noise_pos(1,4)-0.1;

%ChAT
h_dist1=SNR_pos(1,1)+SNR_pos(1,3)+0.03;

evoked_trace_Off_f80_zoom_pos(1,:) = [h_dist1, evoked_trace_Off_f46_zoom_pos(1,2) ,  evoked_trace_Off_f46_zoom_pos(1,3) ,  evoked_trace_Off_f46_zoom_pos(1,4)];
evoked_trace_On_f80_zoom_pos(1,:) = [h_dist1, evoked_trace_On_f46_zoom_pos(1,2) ,  evoked_trace_Off_f80_zoom_pos(1,3) ,  evoked_trace_Off_f80_zoom_pos(1,4)+0.01];
evoked_mean_f80_zoom_pos(1,:) = [h_dist1, evoked_mean_f46_zoom_pos(1,2) ,   evoked_trace_Off_f80_zoom_pos(1,3) ,   evoked_mean_f46_zoom_pos(1,4)+0.01];
evoked_std_f80_zoom_pos(1,:) = [h_dist1, evoked_std_f46_zoom_pos(1,2) ,  evoked_trace_Off_f80_zoom_pos(1,3) ,  evoked_std_f46_zoom_pos(1,4)+0.01];
% evoked_trace_On_f80_zoom_pos(1,:) = [h_dist1, evoked_trace_Off_f80_zoom_pos(1,2)-evoked_trace_Off_f80_zoom_pos(1,4)-v_dist ,  evoked_trace_Off_f80_zoom_pos(1,3) ,  evoked_trace_Off_f80_zoom_pos(1,4)+0.01];
% evoked_mean_f80_zoom_pos(1,:) = [h_dist1, evoked_trace_On_f80_zoom_pos(1,2)-0.10-v_dist ,   evoked_trace_Off_f80_zoom_pos(1,3) ,   evoked_mean_f46_zoom_pos(1,4)+0.01];
% evoked_std_f80_zoom_pos(1,:) = [h_dist1, evoked_mean_f80_zoom_pos(1,2)-evoked_mean_f80_zoom_pos(1,4)-v_dist ,  evoked_trace_Off_f80_zoom_pos(1,3) ,  evoked_std_f46_zoom_pos(1,4)+0.01];

VmM_ChAT_pos(1,:)=[evoked_trace_Off_f80_zoom_pos(1,1)+0.05, VmM_pos(1,2), VmM_pos(1,3), VmM_pos(1,4)]; %VmM_ChAT_pos(1,2)=evoked_std_f80_zoom_pos(1,2)-VmM_ChAT_pos(1,4)-0.05;
VmSTD_ChAT_pos(1,:)= VmM_ChAT_pos; VmSTD_ChAT_pos(1,1)=VmM_ChAT_pos(1,1)+VmM_ChAT_pos(1,3)+0.05;
Peak_amp_ChAT_pos(1,:)=[0.5,Peak_amp_pos(1,2), Peak_amp_pos(1,3), Peak_amp_pos(1,4)];  Peak_amp_ChAT_pos(1,1)=VmSTD_ChAT_pos(1,1)+VmSTD_ChAT_pos(1,3)+0.08;

h_dist2=evoked_trace_Off_f80_zoom_pos(1,1)+evoked_trace_Off_f80_zoom_pos(1,3)+0.04;
SNR_ChAT_pos(1,:)= [h_dist2, SNR_pos(1,2),SNR_pos(1,3), SNR_pos(1,4)]; 
Amplitude_Signal_ChAT_pos(1,:)=SNR_ChAT_pos; Amplitude_Signal_ChAT_pos(1,2)=Amplitude_Signal_pos(1,2);
Amplitude_Noise_ChAT_pos(1,:)=Amplitude_Signal_ChAT_pos; Amplitude_Noise_ChAT_pos(1,2)=Amplitude_Noise_pos(1,2);

%position of top of each panel
evoked_trace_Off_f46_zoom_pos_top = evoked_trace_Off_f46_zoom_pos(1,2)+evoked_trace_Off_f46_zoom_pos(1,4);
evoked_trace_On_f46_zoom_pos_top = evoked_trace_On_f46_zoom_pos(1,2)+evoked_trace_On_f46_zoom_pos(1,4);
evoked_mean_f46_zoom_pos_top = evoked_mean_f46_zoom_pos(1,2)+evoked_mean_f46_zoom_pos(1,4);
evoked_std_f46_zoom_pos_top = evoked_std_f46_zoom_pos(1,2)+evoked_std_f46_zoom_pos(1,4);
evoked_trace_Off_f80_zoom_pos_top = evoked_trace_Off_f80_zoom_pos(1,2)+evoked_trace_Off_f80_zoom_pos(1,4);
evoked_trace_On_f80_zoom_pos_top = evoked_trace_On_f80_zoom_pos(1,2)+evoked_trace_On_f80_zoom_pos(1,4);
evoked_mean_f80_zoom_pos_top = evoked_mean_f80_zoom_pos(1,2)+evoked_mean_f80_zoom_pos(1,4);
evoked_std_f80_zoom_pos_top = evoked_std_f80_zoom_pos(1,2)+evoked_std_f80_zoom_pos(1,4);
VmM_pos_top = VmM_pos(1,2)+VmM_pos(1,4);
VmSTD_pos_top = VmSTD_pos(1,2)+VmSTD_pos(1,4);
Peak_amp_pos_top = Peak_amp_pos(1,2)+Peak_amp_pos(1,4);
SNR_pos_top = SNR_pos(1,2)+SNR_pos(1,4);
Amplitude_Signal_pos_top = Amplitude_Signal_pos(1,2)+Amplitude_Signal_pos(1,4);
Amplitude_Noise_pos_top = Amplitude_Noise_pos(1,2)+Amplitude_Noise_pos(1,4);
VmM_ChAT_pos_top = VmM_ChAT_pos(1,2)+VmM_ChAT_pos(1,4);
VmSTD_ChAT_pos_top = VmSTD_ChAT_pos(1,2)+VmSTD_ChAT_pos(1,4);
Peak_amp_ChAT_pos_top = Peak_amp_ChAT_pos(1,2)+Peak_amp_ChAT_pos(1,4);
SNR_ChAT_pos_top = SNR_ChAT_pos(1,2)+SNR_ChAT_pos(1,4);
Amplitude_Signal_ChAT_pos_top = Amplitude_Signal_ChAT_pos(1,2)+Amplitude_Signal_ChAT_pos(1,4);
Amplitude_Noise_ChAT_pos_top = Amplitude_Noise_ChAT_pos(1,2)+Amplitude_Noise_ChAT_pos(1,4);
%%
%Placing plots in the figure:

evoked_trace_Off_f46_zoom_ax_copy = copyobj(evoked_trace_Off_f46_zoom_ax,F); % copy axes to new fig
set(evoked_trace_Off_f46_zoom_ax_copy,'position',evoked_trace_Off_f46_zoom_pos(1,:))

evoked_trace_On_f46_zoom_ax_copy = copyobj(evoked_trace_On_f46_zoom_ax,F); % copy axes to new fig
set(evoked_trace_On_f46_zoom_ax_copy,'position',evoked_trace_On_f46_zoom_pos(1,:))

evoked_mean_f46_zoom_ax_copy = copyobj(evoked_mean_f46_zoom_ax,F); % copy axes to new fig
set(evoked_mean_f46_zoom_ax_copy,'position',evoked_mean_f46_zoom_pos(1,:))

evoked_std_f46_zoom_ax_copy = copyobj(evoked_std_f46_zoom_ax,F); % copy axes to new fig
set(evoked_std_f46_zoom_ax_copy,'position',evoked_std_f46_zoom_pos(1,:))

%population panels:
VmM_ax_copy = copyobj(VmM_ax,F); % copy axes to new fig
set(VmM_ax_copy(1),'position',VmM_pos(1,:))
set(F, 'currentaxes', VmM_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
VmM_ax_copy(1).FontSize=11;
% VmM_ax_copy(1).XTickLabel=[];
% delete(VmM_ax_copy(1));

VmSTD_ax_copy = copyobj(VmSTD_ax(1),F); % copy axes to new fig
set(VmSTD_ax_copy,'position',VmSTD_pos(1,:))
set(F, 'currentaxes', VmSTD_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
VmSTD_ax_copy.FontSize=11;
% VmSTD_ax_copy.XTickLabel=[];

Peak_amp_ax_copy = copyobj(Peak_amp_ax,F); % copy axes to new fig
set(Peak_amp_ax_copy,'position',Peak_amp_pos(1,:))
set(F, 'currentaxes', Peak_amp_ax_copy); tl=title(''); yl=ylabel('Amp. Norm.');
Peak_amp_ax_copy.FontSize=11;
Peak_amp_ax_copy.XTickLabel(11,:)=[]; Peak_amp_ax_copy.XTickLabel(11,1)='R';

SNR_ax_copy = copyobj(SNR_ax,F); % copy axes to new fig
set(SNR_ax_copy,'position',SNR_pos(1,:))
set(F, 'currentaxes', SNR_ax_copy);  tl=title('');  yl=ylabel('A.U.','fontsize',11); %tl=title('SNR','fontsize',11,'fontweight','normal');
SNR_ax_copy.FontSize=11;
SNR_ax_copy.XTickLabel=[];

Amplitude_Signal_ax_copy = copyobj(Amplitude_Signal_ax,F); % copy axes to new fig
set(Amplitude_Signal_ax_copy,'position',Amplitude_Signal_pos(1,:))
set(F, 'currentaxes', Amplitude_Signal_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11); %tl=title('Signal Amp.','fontsize',11,'fontweight','normal');
Amplitude_Signal_ax_copy.FontSize=11;
Amplitude_Signal_ax_copy.XTickLabel=[];

Amplitude_Noise_ax_copy = copyobj(Amplitude_Noise_ax,F); % copy axes to new fig
set(Amplitude_Noise_ax_copy,'position',Amplitude_Noise_pos(1,:))
set(F, 'currentaxes', Amplitude_Noise_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11);  %tl=title('Noise Amp.','fontsize',11,'fontweight','normal');
Amplitude_Noise_ax_copy.FontSize=11;

%position legend:
set(evoked_mean_f46_zoom_ax_copy(1),'position',[0.32 evoked_trace_Off_f46_zoom_pos_top+0.03 0.08 0.005])
evoked_mean_f46_zoom_ax_copy(1).FontSize=10;
% set(evoked_mean_f46_zoom_ax_copy(1),'position',[0.9 spont_std_f46_pos_top+0.01 0.08 0.05])


%ChAT

evoked_trace_Off_f80_zoom_ax_copy = copyobj(evoked_trace_Off_f80_zoom_ax,F); % copy axes to new fig
set(evoked_trace_Off_f80_zoom_ax_copy,'position',evoked_trace_Off_f80_zoom_pos(1,:))

evoked_trace_On_f80_zoom_ax_copy = copyobj(evoked_trace_On_f80_zoom_ax,F); % copy axes to new fig
set(evoked_trace_On_f80_zoom_ax_copy,'position',evoked_trace_On_f80_zoom_pos(1,:))

evoked_mean_f80_zoom_ax_copy = copyobj(evoked_mean_f80_zoom_ax,F); % copy axes to new fig
set(evoked_mean_f80_zoom_ax_copy,'position',evoked_mean_f80_zoom_pos(1,:))

evoked_std_f80_zoom_ax_copy = copyobj(evoked_std_f80_zoom_ax,F); % copy axes to new fig
set(evoked_std_f80_zoom_ax_copy,'position',evoked_std_f80_zoom_pos(1,:))

VmM_ChAT_ax_copy = copyobj(VmM_ChAT_ax,F); % copy axes to new fig
set(VmM_ChAT_ax_copy(1),'position',VmM_ChAT_pos(1,:))
set(F, 'currentaxes', VmM_ChAT_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
VmM_ChAT_ax_copy(1).FontSize=11;
VmM_ChAT_ax_copy(1).XTickLabel={'Off','On'};
% delete(VmM_ChAT_ax_copy(1));

VmSTD_ChAT_ax_copy = copyobj(VmSTD_ChAT_ax(1),F); % copy axes to new fig
set(VmSTD_ChAT_ax_copy,'position',VmSTD_ChAT_pos(1,:))
set(F, 'currentaxes', VmSTD_ChAT_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
VmSTD_ChAT_ax_copy.FontSize=11;
VmSTD_ChAT_ax_copy.XTickLabel={'Off','On'};

Peak_amp_ChAT_ax_copy = copyobj(Peak_amp_ChAT_ax,F); % copy axes to new fig
set(Peak_amp_ChAT_ax_copy,'position',Peak_amp_ChAT_pos(1,:))
set(F, 'currentaxes', Peak_amp_ChAT_ax_copy); tl=title(''); yl=ylabel('Amp. Norm.');
Peak_amp_ChAT_ax_copy.FontSize=11;
Peak_amp_ChAT_ax_copy.XTickLabel(11,:)=[]; Peak_amp_ax_copy.XTickLabel(11,1)='R';

SNR_ChAT_ax_copy = copyobj(SNR_ChAT_ax,F); % copy axes to new fig
set(SNR_ChAT_ax_copy,'position',SNR_ChAT_pos(1,:))
set(F, 'currentaxes', SNR_ChAT_ax_copy);  tl=title('');  yl=ylabel('A.U.','fontsize',11); %tl=title('SNR','fontsize',11,'fontweight','normal');
SNR_ChAT_ax_copy.FontSize=11;
SNR_ChAT_ax_copy.XTickLabel=[];

Amplitude_Signal_ChAT_ax_copy = copyobj(Amplitude_Signal_ChAT_ax,F); % copy axes to new fig
set(Amplitude_Signal_ChAT_ax_copy,'position',Amplitude_Signal_ChAT_pos(1,:))
set(F, 'currentaxes', Amplitude_Signal_ChAT_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11); %tl=title('Signal Amp.','fontsize',11,'fontweight','normal');
Amplitude_Signal_ChAT_ax_copy.FontSize=11;
Amplitude_Signal_ChAT_ax_copy.XTickLabel=[];

Amplitude_Noise_ChAT_ax_copy = copyobj(Amplitude_Noise_ChAT_ax,F); % copy axes to new fig
set(Amplitude_Noise_ChAT_ax_copy,'position',Amplitude_Noise_ChAT_pos(1,:))
set(F, 'currentaxes', Amplitude_Noise_ChAT_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11);  %tl=title('Noise Amp.','fontsize',11,'fontweight','normal');
Amplitude_Noise_ChAT_ax_copy.FontSize=11;
Amplitude_Noise_ChAT_ax_copy.XTickLabel={'Off','On'};

%position legend:
set(evoked_mean_f80_zoom_ax_copy(1),'position',[0.83 evoked_trace_Off_f80_zoom_pos_top+0.03 0.08 0.005])
evoked_mean_f80_zoom_ax_copy(1).FontSize=10;

 a_pos1=[-0.03 -0.01 0.04 0.04];
 a_pos4=[0 -0.01 0.16 0.04];
 a_pos2=[-0.02 -0.01 0.04 0.04];
 a_pos3=[-0.035 0 0.04 0.04];
 
 annotation('textbox', [evoked_trace_Off_f46_zoom_pos(1,1) evoked_trace_Off_f46_zoom_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Vm Traces', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
annotation('textbox', [evoked_mean_f46_zoom_pos(1,1) evoked_mean_f46_zoom_pos_top 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [evoked_std_f46_zoom_pos(1,1), evoked_std_f46_zoom_pos_top 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [VmM_pos(1,1)-0.005 VmM_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [VmSTD_pos(1,1)-0.005 VmSTD_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [SNR_pos(1,1)+0.01 SNR_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'SNR', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [Amplitude_Signal_pos(1,1)-0.01 Amplitude_Signal_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Signal Amp.', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [Amplitude_Noise_pos(1,1)-0.01 Amplitude_Noise_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Noise Amp.', 'FontName','arial', 'fontsize', 12)
%ChAT:
annotation('textbox', [VmM_ChAT_pos(1,1)-0.005 VmM_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [VmSTD_ChAT_pos(1,1)-0.005 VmSTD_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [SNR_ChAT_pos(1,1)+0.01 SNR_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'SNR', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [Amplitude_Signal_ChAT_pos(1,1)-0.01 Amplitude_Signal_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Signal Amp.', 'FontName','arial', 'fontsize', 12)
annotation('textbox', [Amplitude_Noise_ChAT_pos(1,1)-0.01 Amplitude_Noise_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Noise Amp.', 'FontName','arial', 'fontsize', 12)

 if no_numbering_flag==0;
annotation('textbox', [spont_trace_f46_pos(1,1) spont_trace_f46_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [spont_mean_f46_pos(1,1) spont_mean_f46_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [spont_std_f46_pos(1,1) spont_std_f46_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')

annotation('textbox', [evoked_trace_Off_f46_zoom_pos(1,1) evoked_trace_Off_f46_zoom_pos_top 0 0]+a_pos2,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [evoked_trace_On_f46_zoom_pos(1,1) evoked_trace_On_f46_zoom_pos_top 0 0]+a_pos2,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [evoked_mean_f46_zoom_pos(1,1) evoked_mean_f46_zoom_pos_top 0 0]+a_pos2,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [evoked_std_f46_zoom_pos(1,1) evoked_std_f46_zoom_pos_top 0 0]+a_pos2,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
% annotation('textbox', [VmM_pos(1,1) VmM_pos_top 0 0]+a_pos3,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
% annotation('textbox', [VmSTD_pos(1,1) VmM_pos_top 0 0]+a_pos3,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')

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
  
      cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'
 
  if save_flag==1
    filename='Fig 4_f46_f80';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end
