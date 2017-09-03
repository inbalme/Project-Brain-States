%% plotting a figure 
close all
clear all
save_flag=1;
 no_numbering_flag=0;
 abslen=0.05;
 ax_fontsize=10;
 color1=[255, 102,102]./256; %pink
 color2=[102, 172,255]./256; %light blue
%opening saved figures:

%evoked population

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\LFP_50Hz Vm_50Hz'

VmM = open('Vm M_During Sensory stim.fig');   
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(10),'color',color1); uistack(l1(10),'top')
VmM_ax = get(gcf, 'children');

VmSTD = open('Vm STD_During Sensory stim.fig'); 
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(10),'color',color1); uistack(l1(10),'top')
VmSTD_ax = get(gcf, 'children');

SNR = open('SNR1.fig');    
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(10),'color',color1); uistack(l1(10),'top')
SNR_ax = get(gcf, 'children');

Amplitude_Signal = open('Amplitude_Signal.fig');    
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(10),'color',color1); uistack(l1(10),'top')
Amplitude_Signal_ax = get(gcf, 'children');

Amplitude_Noise = open('Amplitude_Noise1.fig');    
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(10),'color',color1); uistack(l1(10),'top')
Amplitude_Noise_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Evoked';
Peak_amp = open('Evoked Amplitude Normalized.fig');    
Peak_amp_ax = get(gcf, 'children');

%ChAT:

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary\LFP_50Hz Vm_50Hz'

VmM_ChAT = open('Vm M_During Sensory stim.fig');    
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(9),'color',color2); uistack(l1(9),'top')
VmM_ChAT_ax = get(gcf, 'children');

VmSTD_ChAT = open('Vm STD_During Sensory stim.fig');    
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(9),'color',color2); uistack(l1(9),'top')
VmSTD_ChAT_ax = get(gcf, 'children');

SNR_ChAT = open('SNR1.fig');    
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(9),'color',color2); uistack(l1(9),'top')
SNR_ChAT_ax = get(gcf, 'children');

Amplitude_Signal_ChAT = open('Amplitude_Signal.fig');    
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(9),'color',color2); uistack(l1(9),'top')
Amplitude_Signal_ChAT_ax = get(gcf, 'children');

Amplitude_Noise_ChAT = open('Amplitude_Noise1.fig');    
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(9),'color',color2); uistack(l1(9),'top')
Amplitude_Noise_ChAT_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Evoked';
Peak_amp_ChAT = open('Evoked Amplitude Normalized.fig');    
Peak_amp_ChAT_ax = get(gcf, 'children');

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
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 16]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

%% xlimits, y limits, ticks etc.
 
%% Positions:
    
v_dist=0.12;
VmM_pos(1,:)=[0.12, 0.68, 0.1, 0.12]; 
VmSTD_pos(1,:)= VmM_pos; VmSTD_pos(1,1)=VmM_pos(1,1)+VmM_pos(1,3)+0.07;
Peak_amp_pos(1,:)=[VmM_pos(1,1),0.5, 0.3, 0.22 ]; Peak_amp_pos(1,2)=VmM_pos(1,2)-Peak_amp_pos(1,4)-v_dist;
% Peak_amp_pos(1,:)=[0.5,VmM_pos(1,2)+0.03, 0.15, VmM_pos(1,4)-0.03];  Peak_amp_pos(1,1)=VmSTD_pos(1,1)+VmSTD_pos(1,3)+0.08;

SNR_pos(1,:)= [VmM_pos(1,1), 0.5,0.08, 0.12]; SNR_pos(1,2)=Peak_amp_pos(1,2)-SNR_pos(1,4)-v_dist-0.01;
Amplitude_Signal_pos(1,:)=SNR_pos; Amplitude_Signal_pos(1,1)=SNR_pos(1,1)+SNR_pos(1,3)+0.05;
Amplitude_Noise_pos(1,:)=Amplitude_Signal_pos; Amplitude_Noise_pos(1,1)=Amplitude_Signal_pos(1,1)+Amplitude_Signal_pos(1,3)+0.05;

%ChAT
h_dist1=VmSTD_pos(1,1)+VmSTD_pos(1,3)+0.1;

VmM_ChAT_pos(1,:)=VmM_pos; VmM_ChAT_pos(1,1)=0.57;
VmSTD_ChAT_pos(1,:)= VmSTD_pos; VmSTD_ChAT_pos(1,1)=VmM_ChAT_pos(1,1)+VmM_ChAT_pos(1,3)+0.07;
Peak_amp_ChAT_pos(1,:)=Peak_amp_pos; Peak_amp_ChAT_pos(1,1)=VmM_ChAT_pos(1,1);

SNR_ChAT_pos(1,:)= [VmM_ChAT_pos(1,1), SNR_pos(1,2),0.08, 0.12]; 
Amplitude_Signal_ChAT_pos(1,:)=SNR_ChAT_pos; Amplitude_Signal_ChAT_pos(1,1)=SNR_ChAT_pos(1,1)+SNR_ChAT_pos(1,3)+0.05;
Amplitude_Noise_ChAT_pos(1,:)=Amplitude_Signal_ChAT_pos; Amplitude_Noise_ChAT_pos(1,1)=Amplitude_Signal_ChAT_pos(1,1)+Amplitude_Signal_ChAT_pos(1,3)+0.05;

%position of top of each panel

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

NB_illustration_pos(1,:)=[VmM_pos(1,1)+0.02,VmM_pos_top+0.07, 0.11, 0.12];
ChAT_illustration_pos(1,:)=[VmM_ChAT_pos(1,1)+0.12,VmM_ChAT_pos_top+0.07, NB_illustration_pos(1,3), NB_illustration_pos(1,4)];

NB_illustration_pos_top=NB_illustration_pos(1,2)+NB_illustration_pos(1,4);
ChAT_illustration_pos_top=ChAT_illustration_pos(1,2)+ChAT_illustration_pos(1,4);
%%
%Placing plots in the figure:

NB_illustration_ax = axes('position',NB_illustration_pos);
imshow(NB_illustration, 'parent', NB_illustration_ax) 

ChAT_illustration_ax = axes('position',ChAT_illustration_pos);
imshow(ChAT_illustration, 'parent', ChAT_illustration_ax) 

%population panels:
VmM_ax_copy = copyobj(VmM_ax,F); % copy axes to new fig
set(VmM_ax_copy(1),'position',VmM_pos(1,:))
set(F, 'currentaxes', VmM_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
VmM_ax_copy(1).FontSize=ax_fontsize;
 set(gca,'tickdir','out','ylim',[-80,-20],'ytick',[-80,-60,-40,-20])
ticklen=fn_get_abs_ticklength(gca, abslen);
% VmM_ax_copy(1).XTickLabel=[];
% delete(VmM_ax_copy(1));

VmSTD_ax_copy = copyobj(VmSTD_ax(1),F); % copy axes to new fig
set(VmSTD_ax_copy,'position',VmSTD_pos(1,:))
set(F, 'currentaxes', VmSTD_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
VmSTD_ax_copy.FontSize=ax_fontsize;
ticklen=fn_get_abs_ticklength(gca, abslen);
% VmSTD_ax_copy.XTickLabel=[];

Peak_amp_ax_copy = copyobj(Peak_amp_ax,F); % copy axes to new fig
set(Peak_amp_ax_copy,'position',Peak_amp_pos(1,:))
set(F, 'currentaxes', Peak_amp_ax_copy); tl=title(''); yl=ylabel('Amp. Norm.');
Peak_amp_ax_copy.FontSize=ax_fontsize;
Peak_amp_ax_copy.XTickLabel(11,:)=[]; Peak_amp_ax_copy.XTickLabel(11,1)='R';
 set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);

SNR_ax_copy = copyobj(SNR_ax,F); % copy axes to new fig
set(SNR_ax_copy,'position',SNR_pos(1,:))
set(F, 'currentaxes', SNR_ax_copy);  tl=title('');  yl=ylabel('A.U.','fontsize',11); %tl=title('SNR','fontsize',11,'fontweight','normal');
SNR_ax_copy.FontSize=ax_fontsize;
 set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);
% SNR_ax_copy.XTickLabel=[];

Amplitude_Signal_ax_copy = copyobj(Amplitude_Signal_ax,F); % copy axes to new fig
set(Amplitude_Signal_ax_copy,'position',Amplitude_Signal_pos(1,:))
set(F, 'currentaxes', Amplitude_Signal_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11); %tl=title('Signal Amp.','fontsize',11,'fontweight','normal');
Amplitude_Signal_ax_copy.FontSize=ax_fontsize;
 set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);
% Amplitude_Signal_ax_copy.XTickLabel=[];

Amplitude_Noise_ax_copy = copyobj(Amplitude_Noise_ax,F); % copy axes to new fig
set(Amplitude_Noise_ax_copy,'position',Amplitude_Noise_pos(1,:))
set(F, 'currentaxes', Amplitude_Noise_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11);  %tl=title('Noise Amp.','fontsize',11,'fontweight','normal');
Amplitude_Noise_ax_copy.FontSize=ax_fontsize;
 set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);

%position legend:
% set(evoked_mean_f46_zoom_ax_copy(1),'position',[0.32 evoked_trace_Off_f46_zoom_pos_top+0.03 0.08 0.005])
% evoked_mean_f46_zoom_ax_copy(1).FontSize=10;
% set(evoked_mean_f46_zoom_ax_copy(1),'position',[0.9 spont_std_f46_pos_top+0.01 0.08 0.05])

%ChAT

VmM_ChAT_ax_copy = copyobj(VmM_ChAT_ax,F); % copy axes to new fig
set(VmM_ChAT_ax_copy(1),'position',VmM_ChAT_pos(1,:))
set(F, 'currentaxes', VmM_ChAT_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
VmM_ChAT_ax_copy(1).FontSize=ax_fontsize;
VmM_ChAT_ax_copy(1).XTickLabel={'Off','On'};
 set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);
% delete(VmM_ChAT_ax_copy(1));

VmSTD_ChAT_ax_copy = copyobj(VmSTD_ChAT_ax(1),F); % copy axes to new fig
set(VmSTD_ChAT_ax_copy,'position',VmSTD_ChAT_pos(1,:))
set(F, 'currentaxes', VmSTD_ChAT_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
VmSTD_ChAT_ax_copy.FontSize=ax_fontsize;
VmSTD_ChAT_ax_copy.XTickLabel={'Off','On'};
  set(gca,'tickdir','out','ytick',[0,3,6])
ticklen=fn_get_abs_ticklength(gca, abslen);

Peak_amp_ChAT_ax_copy = copyobj(Peak_amp_ChAT_ax,F); % copy axes to new fig
set(Peak_amp_ChAT_ax_copy,'position',Peak_amp_ChAT_pos(1,:))
set(F, 'currentaxes', Peak_amp_ChAT_ax_copy); tl=title(''); yl=ylabel('Amp. Norm.');
Peak_amp_ChAT_ax_copy.FontSize=ax_fontsize;
Peak_amp_ChAT_ax_copy.XTickLabel(11,:)=[]; Peak_amp_ChAT_ax_copy.XTickLabel(11,1)='R';
 set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);

SNR_ChAT_ax_copy = copyobj(SNR_ChAT_ax,F); % copy axes to new fig
set(SNR_ChAT_ax_copy,'position',SNR_ChAT_pos(1,:))
set(F, 'currentaxes', SNR_ChAT_ax_copy);  tl=title('');  yl=ylabel('A.U.','fontsize',11); %tl=title('SNR','fontsize',11,'fontweight','normal');
SNR_ChAT_ax_copy.FontSize=ax_fontsize;
SNR_ChAT_ax_copy.XTickLabel={'Off','On'};
 set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);

Amplitude_Signal_ChAT_ax_copy = copyobj(Amplitude_Signal_ChAT_ax,F); % copy axes to new fig
set(Amplitude_Signal_ChAT_ax_copy,'position',Amplitude_Signal_ChAT_pos(1,:))
set(F, 'currentaxes', Amplitude_Signal_ChAT_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11); %tl=title('Signal Amp.','fontsize',11,'fontweight','normal');
Amplitude_Signal_ChAT_ax_copy.FontSize=ax_fontsize;
Amplitude_Signal_ChAT_ax_copy.XTickLabel={'Off','On'};
 set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);

Amplitude_Noise_ChAT_ax_copy = copyobj(Amplitude_Noise_ChAT_ax,F); % copy axes to new fig
set(Amplitude_Noise_ChAT_ax_copy,'position',Amplitude_Noise_ChAT_pos(1,:))
set(F, 'currentaxes', Amplitude_Noise_ChAT_ax_copy); tl=title('');  yl=ylabel('mV','fontsize',11);  %tl=title('Noise Amp.','fontsize',11,'fontweight','normal');
Amplitude_Noise_ChAT_ax_copy.FontSize=ax_fontsize;
Amplitude_Noise_ChAT_ax_copy.XTickLabel={'Off','On'};
 set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);

%Changing the fontsize of xlabel and ylabel of all axes in the figure:
 haxes=get(gcf,'children');
 haxes1=get(haxes,'type');
 for i=1:size(haxes1,1)
    findax(i)=strcmp({'axes'},haxes1{i});
 end
 haxlabels=get(haxes(findax),{'XLabel' 'YLabel'});
 for i=1:numel(haxlabels)
    set(haxlabels{i},'fontsize',ax_fontsize);
 end
 
 %position legend:
% set(evoked_mean_f80_zoom_ax_copy(1),'position',[0.83 evoked_trace_Off_f80_zoom_pos_top+0.03 0.08 0.005])
% evoked_mean_f80_zoom_ax_copy(1).FontSize=10;

 a_pos1=[-0.03 -0.01 0.04 0.04];
 a_pos4=[0 -0.01 0.16 0.04];
 a_pos2=[-0.02 -0.01 0.04 0.04];
 a_pos3=[-0.035 0.015 0.04 0.04];

 annotation('textbox', [0.255, NB_illustration_pos_top-0.05, 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Sensory evoked synaptic responses', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 
 annotation('textbox', [VmM_pos(1,1)-0.005 VmM_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [VmSTD_pos(1,1)-0.005 VmSTD_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [Peak_amp_pos(1,1)+0.05, Peak_amp_pos_top+0.005, 0.2, 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Response Amplitude', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [SNR_pos(1,1)+0.01 SNR_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'SNR', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [Amplitude_Signal_pos(1,1)-0.01 Amplitude_Signal_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Signal Amp.', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [Amplitude_Noise_pos(1,1)-0.01 Amplitude_Noise_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Noise Amp.', 'FontName','arial', 'fontsize', ax_fontsize)
%ChAT:
annotation('textbox', [VmM_ChAT_pos(1,1)-0.005 VmM_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [VmSTD_ChAT_pos(1,1)-0.005 VmSTD_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [Peak_amp_ChAT_pos(1,1)+0.05, Peak_amp_ChAT_pos_top+0.005, 0.2, 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Response Amplitude', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [SNR_ChAT_pos(1,1)+0.01 SNR_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'SNR', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [Amplitude_Signal_ChAT_pos(1,1)-0.01 Amplitude_Signal_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Signal Amp.', 'FontName','arial', 'fontsize', ax_fontsize)
annotation('textbox', [Amplitude_Noise_ChAT_pos(1,1)-0.01 Amplitude_Noise_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Noise Amp.', 'FontName','arial', 'fontsize', ax_fontsize)

 if no_numbering_flag==0;
annotation('textbox', [VmM_pos(1,1) VmM_pos_top 0 0]+a_pos3,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [VmSTD_pos(1,1) VmSTD_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [Peak_amp_pos(1,1) Peak_amp_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [SNR_pos(1,1) SNR_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [Amplitude_Signal_pos(1,1) Amplitude_Signal_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [Amplitude_Noise_pos(1,1) Amplitude_Noise_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')

annotation('textbox', [VmM_ChAT_pos(1,1) VmM_ChAT_pos_top 0 0]+a_pos3,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [VmSTD_ChAT_pos(1,1) VmSTD_ChAT_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [Peak_amp_ChAT_pos(1,1) Peak_amp_ChAT_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [SNR_ChAT_pos(1,1) SNR_ChAT_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [Amplitude_Signal_ChAT_pos(1,1) Amplitude_Signal_ChAT_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'K', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [Amplitude_Noise_ChAT_pos(1,1) Amplitude_Noise_ChAT_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'L', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 end
  
      cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'
 
  if save_flag==1
    filename='Fig 5_f46_f80';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end
