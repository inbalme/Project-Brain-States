%% plotting a figure 

%% new version

close all
clear all
save_flag=0;
no_numbering_flag=0;
abslen=0.05; %[in cm]
ax_fontsize=10;
pairPlotLineWidth=1.2;
letter_fontsize=12; 
 barWidth=1.5;
color1=[255, 102,102]./256; %pink
color2=[102, 172,255]./256; %light blue
color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [102, 172,255]./256; [30,75,14]/256; [112,172,90]/256];  

%% opening saved figures:

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
ChAT_illustration = imread('ChAT_Schematic_Illustration_one-electrode_tight','TIF');
ChAT_histology = imread('ChAT_Histology_tight','PNG');

% long traces ChAT
cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Long Trace Presentation\style2'
%spontaneous:
spont_trace_f121 = open('f121_VmWhisk_traces_x1.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
textposition=get(b2(12),'position');
set(b2(12),'position', [textposition(1)-0.5,textposition(2)+6.5,textposition(3)]);
linexdata=get(b1(6),'xdata');
lineydata=get(b1(6),'ydata');
set(b1(6),'xdata',linexdata-200/1000,'linewidth',barWidth,'ydata',lineydata+5);
delete(b1(4)); 
delete(b2(8));
spont_trace_f121_ax = get(gcf, 'children');

spont_mean_f121 = open('f121_mean_x1_v2.fig');   
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
set(b2(1),'string','5 mV');
set(b2(2),'string','1 s');
textposition=get(b2(1),'position');
set(b2(1),'position', [textposition(1)-0.6,textposition(2)+0.5,textposition(3)]);
b1=findall(gcf,'type','line');
linexdata=get(b1(1),'xdata');
lineydata=get(b1(1),'ydata');
set(b1(1),'xdata',linexdata-200/1000,'linewidth',barWidth,'ydata',lineydata+0.5);
lineydata=get(b1(2),'ydata');
set(b1(2),'ydata',lineydata-0.2,'linewidth',barWidth);
rec1=findall(gcf,'type','rectangle');
delete(rec1)
spont_mean_f121_ax = get(gcf, 'children');

spont_std_f121 = open('f121_std_x1_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
set(b2(1),'string','2 mV');
textposition=get(b2(1),'position');
set(b2(1),'position', [textposition(1)-0.65,textposition(2)+1,textposition(3)]);
linexdata=get(b1(1),'xdata');
lineydata=get(b1(1),'ydata');
set(b1(1),'xdata',linexdata-200/1000,'linewidth',barWidth,'ydata',lineydata+1);
delete(b1(2,1)); 
delete(b2(2,1));
rec1=findall(gcf,'type','rectangle');
delete(rec1)
spont_std_f121_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary_awake\LFP_50Hz Vm_50Hz';
VmM_ChAT = open('Vm M_Before Sensory stim.fig');    
         l1=findall(gcf,'color',[0.7 0.7 0.7]);
         set(l1,'linewidth',pairPlotLineWidth);
         set(l1(6),'color',color2); uistack(l1(6),'top'); %l1(1)for f129,  l1(6) for f121
VmM_ChAT_ax = get(gcf, 'children');

VmSTD_ChAT = open('Vm STD_Before Sensory stim.fig');    
         l1=findall(gcf,'color',[0.7 0.7 0.7]);
         set(l1,'linewidth',pairPlotLineWidth);
         set(l1(6),'color',color2); uistack(l1(6),'top')
VmSTD_ChAT_ax = get(gcf, 'children');

 cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm Histograms and Powerspec_awake'

        spont_Vm_hist_ChAT = open('f121_Vm_histogram_ongoing.fig');
        %placing legend:
hVm=get(gca,'Children');
[l2,OBJH2,OUTH2,OUTM2] = legend([hVm(1) hVm(2)],{'Light Off','Light On'},'fontsize',9, 'box', 'off'); % returns a handle LEGH to the legend axes; a vector OBJH containing handles for the text, lines, and patches in the legend; a vector OUTH of handles to thelines and patches in the plot; and a cell array OUTM containingthe text in the legend.
lpatch2=findobj(OBJH2,'type','patch');
lpatch2(1).Vertices(2:3,2)=lpatch2(1).Vertices(1,2)+0.15;
lpatch2(2).Vertices(2:3,2)=lpatch2(2).Vertices(1,2)+0.15;
lpatch2(1).Vertices(:,2)=lpatch2(1).Vertices(:,2)+0.06;
lpatch2(2).Vertices(:,2)=lpatch2(2).Vertices(:,2)+0.08;
lpatch2(1).Vertices([1,2,5],1)=lpatch2(1).Vertices(3,1)-0.2;
lpatch2(2).Vertices([1,2,5],1)=lpatch2(1).Vertices(3,1)-0.2;
lpatch2(1).FaceColor=color_table(1,:);
lpatch2(1).FaceAlpha=0.7;
lpatch2(2).FaceColor=color_table(2,:);
lpatch2(2).FaceAlpha=0.3;
% saveas(gcf,'f121_Vm_histogram_ongoing','fig'); 
spont_Vm_hist_ChAT_ax = get(gcf, 'children');

        spont_Vm_5prcentile_ChAT = open('Vm_5prcentile_ongoing.fig'); 
            l1=findall(gcf,'color',[0.7 0.7 0.7]);
            set(l1,'linewidth',pairPlotLineWidth);
            set(l1(6),'color',color2); uistack(l1(6),'top') 
        spont_Vm_5prcentile_ChAT_ax = get(gcf, 'children');
 
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis_awake\Spontaneous'
    
       spont_event_freq_ChAT = open('Spontaneous event frequency.fig');  
     set(gca,'ylim',[0 10])
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1,'linewidth',pairPlotLineWidth);
        set(l1(6),'color',color2); uistack(l1(6),'top')
    spont_event_freq_ChAT_ax = get(gcf, 'children');

    spont_event_amp_ChAT = open('Spontaneous event amplitude.fig');  
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1,'linewidth',pairPlotLineWidth);
        set(l1(6),'color',color2); uistack(l1(6),'top')
    spont_event_amp_ChAT_ax = get(gcf, 'children');
% cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Rin'
% 
%     spont_Rin_ChAT = open('Rin.fig');    
%     spont_Rin_ChAT_ax = get(gcf, 'children');

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',11);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 18]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
% set(gcf, 'PaperPositionMode', 'auto');

%% Positions: (reduce 20% of hight)
 
ChAT_illustration_pos(1,:)=[0.11,0.79, 0.24, 0.2];
ChAT_histology_pos(1,:)=[ChAT_illustration_pos(1,1)+ChAT_illustration_pos(1,3)+0.05,ChAT_illustration_pos(1,2), 0.6, ChAT_illustration_pos(1,4)];

v_dist=0.01;
spont_trace_f121_pos(7,:) = [0.1 , 0.7 ,  0.6 ,  0.06]; 
spont_trace_f121_pos(6,:) = spont_trace_f121_pos(7,:)-[0 , spont_trace_f121_pos(7,4) , 0 , 0];
spont_trace_f121_pos(5,:) = spont_trace_f121_pos(6,:)-[0 , spont_trace_f121_pos(6,4)-0.02 , 0 , 0];
spont_trace_f121_pos(4,:) = spont_trace_f121_pos(5,:)-[0 , spont_trace_f121_pos(5,4) , 0 , 0];
spont_trace_f121_pos(3,:) = spont_trace_f121_pos(4,:)-[0 , spont_trace_f121_pos(4,4)-0.02 , 0 , 0];
spont_trace_f121_pos(1,:) = spont_trace_f121_pos(3,:)-[0 , spont_trace_f121_pos(3,4) , 0 , 0];
spont_trace_f121_pos(2,:) =spont_trace_f121_pos(1,:);

spont_mean_f121_pos(1,:) = [spont_trace_f121_pos(1,1), spont_trace_f121_pos(1,2)-0.09-v_dist, spont_trace_f121_pos(1,3), 0.09];
spont_std_f121_pos(1,:) = [spont_trace_f121_pos(1,1), spont_mean_f121_pos(1,2)-spont_mean_f121_pos(1,4)-v_dist-0.02, spont_trace_f121_pos(1,3), spont_mean_f121_pos(1,4)];

spont_Vm_hist_ChAT_pos(1,:)=[0.15 , 0.5 ,  0.27 ,  0.09]; spont_Vm_hist_ChAT_pos(1,2)=spont_std_f121_pos(1,2)-spont_Vm_hist_ChAT_pos(1,4)-0.06;
h_dist1=spont_Vm_hist_ChAT_pos(1,1)+spont_Vm_hist_ChAT_pos(1,3)+0.16;
spont_Vm_5prcentile_ChAT_pos(1,:) = [h_dist1 , spont_Vm_hist_ChAT_pos(1,2)-0.02 ,  0.12 ,  0.14];
spont_Rin_ChAT_pos(1,:) = [h_dist1+spont_Vm_5prcentile_ChAT_pos(1,3)+0.12, spont_Vm_5prcentile_ChAT_pos(1,2) ,  spont_Vm_5prcentile_ChAT_pos(1,3) ,  spont_Vm_5prcentile_ChAT_pos(1,4)];

v_dist2=0.1;
VmM_ChAT_pos(1,:)=[spont_Rin_ChAT_pos(1,1), spont_trace_f121_pos(1,2)+0.04, 0.12, 0.14]; 
VmSTD_ChAT_pos(1,:)= VmM_ChAT_pos; VmSTD_ChAT_pos(1,2)=VmM_ChAT_pos(1,2)-VmM_ChAT_pos(1,4)-0.1;

spont_event_freq_ChAT_pos(1,:) =[spont_amp_hist_ChAT_pos(1,1),spont_event_freq_pos(1,2) , spont_event_freq_pos(1,3),spont_event_freq_pos(1,4)]; 
spont_event_amp_ChAT_pos(1,:)=[spont_event_freq_ChAT_pos(1,1)+spont_event_freq_ChAT_pos(1,3)+0.09,spont_event_freq_ChAT_pos(1,2),spont_event_freq_ChAT_pos(1,3),spont_event_freq_ChAT_pos(1,4)];

%position of top of each panel
spont_trace_f121_pos_top = spont_trace_f121_pos(7,2)+spont_trace_f121_pos(7,4);
spont_mean_f121_pos_top = spont_mean_f121_pos(1,2)+spont_mean_f121_pos(1,4);
spont_std_f121_pos_top = spont_std_f121_pos(1,2)+spont_std_f121_pos(1,4);
VmM_ChAT_pos_top = VmM_ChAT_pos(1,2)+VmM_ChAT_pos(1,4);
VmSTD_ChAT_pos_top = VmSTD_ChAT_pos(1,2)+VmSTD_ChAT_pos(1,4);
spont_Vm_hist_ChAT_pos_top=spont_Vm_hist_ChAT_pos(1,2)+spont_Vm_hist_ChAT_pos(1,4);
spont_Vm_5prcentile_ChAT_pos_top=spont_Vm_5prcentile_ChAT_pos(1,2)+spont_Vm_5prcentile_ChAT_pos(1,4);
spont_event_freq_ChAT_pos_top = spont_event_freq_ChAT_pos(1,2)+spont_event_freq_ChAT_pos(1,4);
spont_event_amp_ChAT_pos_top = spont_event_amp_ChAT_pos(1,2)+spont_event_amp_ChAT_pos(1,4);
% spont_Rin_ChAT_pos_top=spont_Rin_ChAT_pos(1,2)+spont_Rin_ChAT_pos(1,4);

ChAT_illustration_pos_top=ChAT_illustration_pos(1,2)+ChAT_illustration_pos(1,4);
ChAT_histology_pos_top=ChAT_histology_pos(1,2)+ChAT_histology_pos(1,4);

%%
%Placing plots in the figure:
%Cell 80 - spont

    spont_trace_f121_ax_copy = copyobj(spont_trace_f121_ax,F); % copy axes to new fig
for i=1:length(spont_trace_f121_ax)
    set(spont_trace_f121_ax_copy(i),'position',spont_trace_f121_pos(i,:))
end

% spont_trace_f121_ax_copy = copyobj(spont_trace_f121_ax,F); % copy axes to new fig
% set(spont_trace_f121_ax_copy,'position',spont_trace_f121_pos(1,:))

spont_mean_f121_ax_copy = copyobj(spont_mean_f121_ax,F); % copy axes to new fig
set(spont_mean_f121_ax_copy(2),'position',spont_mean_f121_pos(1,:))

spont_std_f121_ax_copy = copyobj(spont_std_f121_ax,F); % copy axes to new fig
set(spont_std_f121_ax_copy,'position',spont_std_f121_pos(1,:))

% position legend:
set(spont_trace_f121_ax_copy(2),'position',[0.65 spont_trace_f121_pos_top 0.08 0.005])
spont_trace_f121_ax_copy(2).FontSize=9;
% delete(spont_trace_f121_ax_copy(2));
%position legend:
% set(spont_mean_f121_ax_copy(1),'position',[0.65 spont_trace_f121_pos_top 0.08 0.005])
% spont_mean_f121_ax_copy(1).FontSize=9;
delete(spont_mean_f121_ax_copy(1));
%population panels:
VmM_ChAT_ax_copy = copyobj(VmM_ChAT_ax,F); % copy axes to new fig
set(VmM_ChAT_ax_copy(1),'position',VmM_ChAT_pos(1,:))
set(F, 'currentaxes', VmM_ChAT_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
VmM_ChAT_ax_copy(1).FontSize=ax_fontsize;
VmM_ChAT_ax_copy(1).XTickLabel={'Off','On'};
set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);

VmSTD_ChAT_ax_copy = copyobj(VmSTD_ChAT_ax(1),F); % copy axes to new fig
set(VmSTD_ChAT_ax_copy,'position',VmSTD_ChAT_pos(1,:))
set(F, 'currentaxes', VmSTD_ChAT_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
VmSTD_ChAT_ax_copy.FontSize=ax_fontsize;
VmSTD_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);

% histogram+5th percentile+Rin

spont_Vm_hist_ChAT_ax_copy = copyobj(spont_Vm_hist_ChAT_ax,F); % copy axes to new fig
set(spont_Vm_hist_ChAT_ax_copy(2),'position',spont_Vm_hist_ChAT_pos(1,:))
set(F, 'currentaxes', spont_Vm_hist_ChAT_ax_copy(2));  tl=title('');
yl=ylabel({'Prob.'},'fontsize',ax_fontsize);
xl=xlabel('Vm (mV)','fontsize',ax_fontsize);
spont_Vm_hist_ChAT_ax_copy(2).FontSize=ax_fontsize;
set(gca,'tickdir','out');
ticklen=fn_get_abs_ticklength(gca, abslen);
% spont_Vm_hist_ax_copy.XTickLabel=[];

%placing legend:
hVm=get(gca,'Children');
[l2,OBJH2,OUTH2,OUTM2] = legend([hVm(1) hVm(2)],{'Light Off','Light On'}, 'position',[spont_Vm_hist_ChAT_pos(1,1)+spont_Vm_hist_ChAT_pos(1,3)-0.14,spont_Vm_hist_ChAT_pos_top-0.02, 0.2,0.02],'fontsize',9, 'box', 'off'); % returns a handle LEGH to the legend axes; a vector OBJH containing handles for the text, lines, and patches in the legend; a vector OUTH of handles to thelines and patches in the plot; and a cell array OUTM containingthe text in the legend.
lpatch2=findobj(OBJH2,'type','patch');
lpatch2(1).Vertices(2:3,2)=lpatch2(1).Vertices(1,2)+0.15;
lpatch2(2).Vertices(2:3,2)=lpatch2(2).Vertices(1,2)+0.15;
lpatch2(1).Vertices(:,2)=lpatch2(1).Vertices(:,2)+0.06;
lpatch2(2).Vertices(:,2)=lpatch2(2).Vertices(:,2)+0.08;
lpatch2(1).Vertices([1,2,5],1)=lpatch2(1).Vertices(3,1)-0.2;
lpatch2(2).Vertices([1,2,5],1)=lpatch2(1).Vertices(3,1)-0.2;
lpatch2(1).FaceColor=color_table(1,:);
lpatch2(1).FaceAlpha=0.7;
lpatch2(2).FaceColor=color_table(2,:);
lpatch2(2).FaceAlpha=0.3;

spont_Vm_5prcentile_ChAT_ax_copy = copyobj(spont_Vm_5prcentile_ChAT_ax,F); % copy axes to new fig
set(spont_Vm_5prcentile_ChAT_ax_copy,'position',spont_Vm_5prcentile_ChAT_pos(1,:))
set(F, 'currentaxes', spont_Vm_5prcentile_ChAT_ax_copy);   tl=title('');
% yl=ylabel({'Lower 5th'; 'perc. (mV)'},'fontsize',ax_fontsize);
spont_Vm_5prcentile_ChAT_ax_copy.FontSize=ax_fontsize;
spont_Vm_5prcentile_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out');
ticklen=fn_get_abs_ticklength(gca, abslen);

spont_event_freq_ChAT_ax_copy = copyobj(spont_event_freq_ChAT_ax,F); % copy axes to new fig
set(spont_event_freq_ChAT_ax_copy,'position',spont_event_freq_ChAT_pos(1,:))
set(F, 'currentaxes', spont_event_freq_ChAT_ax_copy);  xl=xlabel('');  tl=title('');
yl=ylabel('Freq. (Hz)');
spont_event_freq_ChAT_ax_copy.FontSize=ax_fontsize;
spont_event_freq_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out');
ticklen=fn_get_abs_ticklength(gca, abslen);

  spont_event_amp_ChAT_ax_copy = copyobj(spont_event_amp_ChAT_ax,F); % copy axes to new fig
set(spont_event_amp_ChAT_ax_copy,'position',spont_event_amp_ChAT_pos(1,:))
set(F, 'currentaxes', spont_event_amp_ChAT_ax_copy);  xl=xlabel('');  tl=title(''); 
yl=ylabel('Amp. (mV)');
spont_event_amp_ChAT_ax_copy.FontSize=ax_fontsize;
spont_event_amp_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out');
ticklen=fn_get_abs_ticklength(gca, abslen);

% spont_Rin_ChAT_ax_copy = copyobj(spont_Rin_ChAT_ax,F); % copy axes to new fig
% set(spont_Rin_ChAT_ax_copy,'position',spont_Rin_ChAT_pos(1,:))
% set(F, 'currentaxes', spont_Rin_ChAT_ax_copy);   tl=title(''); 
% yl=ylabel('Rin (MOhm)','fontsize',ax_fontsize);
% spont_Rin_ChAT_ax_copy.FontSize=ax_fontsize;
% spont_Rin_ChAT_ax_copy.XTickLabel={'Off','On'};
% set(gca,'tickdir','out');
% ticklen=fn_get_abs_ticklength(gca, abslen);


% ChAT_illustration_ax = axes('position',ChAT_illustration_pos);
% imshow(ChAT_illustration, 'parent', ChAT_illustration_ax) 

% ChAT_histology_ax = axes('position',ChAT_histology_pos);
% imshow(ChAT_histology, 'parent', ChAT_histology_ax) 

%   annotation:

annotation('textbox', [spont_trace_f121_pos(7,1)+0.05, spont_trace_f121_pos_top 0 0]+[0 0.015 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous Activity', 'FontName','arial', 'fontsize', 12) %, 'fontweight', 'bold'
 
 a_pos1=[-0.03 -0.01 0.04 0.04];
 a_pos4=[0 -0.01 0.16 0.04];
 a_pos2=[-0.02 -0.01 0.04 0.04];
 a_pos3=[-0.035 0 0.04 0.04];

% annotation('textbox', [spont_trace_f121_pos(1,1)+0.15 spont_trace_f121_pos_top+0.04 0.4 0],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Optogenetic stimulation', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
annotation('textbox', [spont_trace_f121_pos(7,1) spont_trace_f121_pos_top-0.016 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Vm Traces', 'FontName','arial', 'fontsize', 11) %'fontweight', 'bold'
annotation('textbox', [spont_mean_f121_pos(1,1) spont_mean_f121_pos_top-0.0085 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 11)
annotation('textbox', [spont_std_f121_pos(1,1) spont_std_f121_pos_top-0.009 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 11)
    
annotation('textbox', [VmM_ChAT_pos(1,1)+0.01 VmM_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 10)
annotation('textbox', [VmSTD_ChAT_pos(1,1)+0.01, VmSTD_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 10)

 if no_numbering_flag==0;

annotation('textbox', [spont_trace_f121_pos(1,1)-0.05 ChAT_illustration_pos_top+0.01 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
%  annotation('textbox', [ChAT_histology_pos(1,1)+0.05 ChAT_histology_pos_top-0.02 0 0],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_trace_f121_pos(1,1)-0.05 spont_trace_f121_pos_top+0.03 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [spont_mean_f121_pos(1,1)-0.05 spont_mean_f121_pos_top+0.03 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [spont_std_f121_pos(1,1)-0.05 spont_std_f121_pos_top+0.03 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [VmM_ChAT_pos(1,1)-0.08 VmM_ChAT_pos_top+0.06 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [VmM_ChAT_pos(1,1)-0.08 VmSTD_ChAT_pos_top+0.06 0 0],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [spont_std_f121_pos(1,1)-0.05 spont_Vm_hist_ChAT_pos_top+0.05 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
annotation('textbox', [spont_Vm_5prcentile_ChAT_pos(1,1)-0.09 spont_Vm_hist_ChAT_pos_top+0.05 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', letter_fontsize, 'fontweight', 'bold')
% annotation('textbox', [spont_Rin_ChAT_pos(1,1) spont_Vm_hist_ChAT_pos_top 0 0]+a_pos1,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 end

      cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'
  
if save_flag==1
    filename='Fig9_intracellular_traces_spont_f121';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print -depsc2 Fig9_intracellular_traces_spont_f121
%     print(F, '-depsc2', filename);
end

%     savefig('Fig 2b intracellular traces spont_f121_v2.fig')
