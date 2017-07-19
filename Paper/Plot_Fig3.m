%% plotting a figure 

close all
clear all
save_flag=0;
no_numbering_flag=1;
dotsize=5; %markersize of dots for detected events
%opening saved figures:

        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm Histograms and Powerspec\Spontaneous histograms'
  
    spont_Vm_hist = open('f46_Vm_histogram_ongoing.fig');
   
    spont_Vm_hist_ax = get(gcf, 'children'); 
%     spont_Vm_median = open('f46_median_ci_ongoing.fig');
%     hm=get(gca,'Children');
%     set(hm,'MarkerSize',2)
%     spont_Vm_median_ax = get(gcf, 'children'); 
    spont_Vm_5prcentile = open('Vm_5prcentile_ongoing.fig');    
    spont_Vm_5prcentile_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous\detection f46'
         color_onset=[255,178,102]./256; color_top=[255,153,16]./256;
        detection_example1= open('detection example - ongoing NB-, trace 3.fig');    %trace 5
        a=findall(gca,'type','line'); delete(a(2));
        b=findall(gca,'type','text'); delete(b(2));  set(b(1),'rotation',0)
        c=findall(gca,'type','scatter');
        set(c(1),'markerfacecolor',color_top,'sizedata',dotsize);
        set(c(2),'markerfacecolor',color_onset,'sizedata',dotsize);
        detection_example1_ax = get(gcf, 'children');
        detection_example2= open('detection example - ongoing NB+, trace 6.fig'); 
        a=findall(gca,'type','line'); delete(a(1:2));
        b=findall(gca,'type','text'); delete(b(1:2));
        c=findall(gca,'type','scatter');
       set(c(1),'markerfacecolor',color_top,'sizedata',dotsize);
        set(c(2),'markerfacecolor',color_onset,'sizedata',dotsize);
        detection_example2_ax = get(gcf, 'children');
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous\detection f44'
        detection_example3= open('detection example - ongoing NB-, trace 7.fig');    %1
        a=findall(gca,'type','line'); %delete(a(2));
        b=findall(gca,'type','text'); set(b(1),'rotation',0)%delete(b(2));
        c=findall(gca,'type','scatter');
        set(c(1),'markerfacecolor',color_top,'sizedata',dotsize);
        set(c(2),'markerfacecolor',color_onset,'sizedata',dotsize);
        detection_example3_ax = get(gcf, 'children');
        detection_example4= open('detection example - ongoing NB+, trace 4.fig');    
        a=findall(gca,'type','line'); delete(a(1:2));
        b=findall(gca,'type','text'); delete(b(1:2));
        c=findall(gca,'type','scatter');
       set(c(1),'markerfacecolor',color_top,'sizedata',dotsize);
        set(c(2),'markerfacecolor',color_onset,'sizedata',dotsize);
        detection_example4_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous'
%         eventoverlay = open('f46_events overlay.fig');    
%         eventoverlay_ax = get(gcf, 'children');
    spont_amp_hist = open('Event_Amp_Hist.fig');
    spont_amp_hist_ax = get(gcf, 'children');
    spont_amp_median = open('Event_Amp_Median_ci.fig');
    hm=get(gca,'Children');
    set(hm,'MarkerSize',2)
    spont_amp_median_ax = get(gcf, 'children');
    
    spont_event_freq = open('Spontaneous event frequency.fig');    
    spont_event_freq_ax = get(gcf, 'children');

    spont_event_amp = open('Spontaneous event amplitude.fig');    
    spont_event_amp_ax = get(gcf, 'children');
    
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Rin'
    spont_Rin = open('Rin.fig');    
    spont_Rin_ax = get(gcf, 'children');

    %% ChAT:
    cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm Histograms and Powerspec\Spontaneous histograms'

        spont_Vm_hist_ChAT = open('f80_Vm_histogram_ongoing.fig');
        spont_Vm_hist_ChAT_ax = get(gcf, 'children');
%         spont_Vm_median_ChAT = open('f80_median_ci_ongoing.fig');
%         hm=get(gca,'Children');
%         set(hm,'MarkerSize',2)
%         spont_Vm_median_ChAT_ax = get(gcf, 'children'); 
        spont_Vm_5prcentile_ChAT = open('Vm_5prcentile_ongoing.fig');    
        spont_Vm_5prcentile_ChAT_ax = get(gcf, 'children');
cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous\detection f80'
color_onset=[255,178,102]./256; color_top=[255,153,16]./256;
        detection_example1_ChAT= open('detection example - ongoing NB-, trace 5.fig');    
        a=findall(gca,'type','line'); delete(a(2));
        b=findall(gca,'type','text'); delete(b(2));
        c=findall(gca,'type','scatter');
        set(c(1),'sizedata',dotsize);
        set(c(2),'sizedata',dotsize);
        detection_example1_ChAT_ax = get(gcf, 'children');
        detection_example2_ChAT= open('detection example - ongoing NB+, trace 6.fig'); 
        a=findall(gca,'type','line'); delete(a(1:2));
        b=findall(gca,'type','text'); delete(b(1:2));
        c=findall(gca,'type','scatter');
        set(c(1),'sizedata',dotsize);
        set(c(2),'sizedata',dotsize);
        detection_example2_ChAT_ax = get(gcf, 'children');
         cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous\detection f74'
        detection_example3_ChAT= open('detection example - ongoing NB-, trace 3.fig');    %8,9
        a=findall(gca,'type','line'); %delete(a(2));
        b=findall(gca,'type','text'); %delete(b(2));
        c=findall(gca,'type','scatter');
        set(c(1),'sizedata',dotsize);
        set(c(2),'sizedata',dotsize);
        detection_example3_ChAT_ax = get(gcf, 'children');
        detection_example4_ChAT= open('detection example - ongoing NB+, trace 10.fig');    
        a=findall(gca,'type','line'); delete(a(1:2));
        b=findall(gca,'type','text'); delete(b(1:2));
        c=findall(gca,'type','scatter');
        set(c(1),'sizedata',dotsize);
        set(c(2),'sizedata',dotsize);
        detection_example4_ChAT_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous'
%         eventoverlay_ChAT = open('f80_events overlay.fig');    
%         eventoverlay_ChAT_ax = get(gcf, 'children');
   
    spont_amp_hist_ChAT = open('Event_Amp_Hist.fig');
    spont_amp_hist_ChAT_ax = get(gcf, 'children');
    spont_amp_median_ChAT = open('Event_Amp_Median_ci.fig');
    hm=get(gca,'Children');
    set(hm,'MarkerSize',2)
    spont_amp_median_ChAT_ax = get(gcf, 'children');
    
    spont_event_freq_ChAT = open('Spontaneous event frequency.fig');    
    spont_event_freq_ChAT_ax = get(gcf, 'children');

    spont_event_amp_ChAT = open('Spontaneous event amplitude.fig');    
    spont_event_amp_ChAT_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Rin'

    spont_Rin_ChAT = open('Rin.fig');    
    spont_Rin_ChAT_ax = get(gcf, 'children');

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 25]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
 
%% Positions:
%NB
spont_Vm_hist_pos(1,:)=[0.1, 0.88, 0.35, 0.08];
% spont_Vm_median_pos(1,:)=[spont_Vm_hist_pos(1,1)+0.22, spont_Vm_hist_pos(1,2)+0.02, 0.08, 0.04];
h_dist1=spont_Vm_hist_pos(1,1)+spont_Vm_hist_pos(1,3)+0.18;
spont_Vm_5prcentile_pos(1,:) = [h_dist1 , spont_Vm_hist_pos(1,2) ,  0.1 ,  0.08];
spont_Rin_pos(1,:) = [h_dist1+spont_Vm_5prcentile_pos(1,3)+0.12, spont_Vm_5prcentile_pos(1,2) ,  spont_Vm_5prcentile_pos(1,3) ,  spont_Vm_5prcentile_pos(1,4)];

detection_example1_pos(1,:) = [0.05 , 0.5 , 0.45 , 0.12]; detection_example1_pos(1,2)=spont_Vm_hist_pos(1,2)-detection_example1_pos(1,4)-0.05;
hdist3=detection_example1_pos(1,1)+detection_example1_pos(1,3)+0.05;
detection_example2_pos(1,:) = [hdist3 , detection_example1_pos(1,2) , detection_example1_pos(1,3) , detection_example1_pos(1,4)]; %detection_example2_pos(1,2)=detection_example1_pos(1,2)-detection_example1_pos(1,4)+0.05;
detection_example3_pos(1,:) = [detection_example1_pos(1,1) , 0.5 , detection_example1_pos(1,3) , detection_example1_pos(1,4)]; detection_example3_pos(1,2)=detection_example1_pos(1,2)-detection_example3_pos(1,4)+0.07;
detection_example4_pos(1,:) = [hdist3 , detection_example3_pos(1,2) , detection_example1_pos(1,3) , detection_example1_pos(1,4)]; %detection_example4_pos(1,2)=detection_example3_pos(1,2)-detection_example3_pos(1,4)+0.05;

spont_amp_hist_pos(1,:) = [spont_Vm_hist_pos(1,1), 0.5 ,  spont_Vm_hist_pos(1,3) ,  0.08]; spont_amp_hist_pos(1,2)=detection_example4_pos(1,2)-spont_amp_hist_pos(1,4)-0.01;
spont_amp_median_pos(1,:) = [spont_amp_hist_pos(1,1)+0.24, spont_amp_hist_pos(1,2)+0.04,0.05,0.05];
h_dist2=spont_amp_hist_pos(1,1)+spont_amp_hist_pos(1,3)+0.15;
spont_event_freq_pos(1,:) =[spont_Vm_5prcentile_pos(1,1),spont_amp_hist_pos(1,2)-0.02, spont_Vm_5prcentile_pos(1,3),spont_Vm_5prcentile_pos(1,4)];
spont_event_amp_pos(1,:)=[spont_event_freq_pos(1,1)+spont_event_freq_pos(1,3)+0.1,spont_event_freq_pos(1,2),spont_event_freq_pos(1,3),spont_event_freq_pos(1,4)];

%ChAT:
spont_Vm_hist_ChAT_pos(1,:)=[0.1, 0.4, 0.35, 0.08];
% spont_Vm_median_ChAT_pos(1,:) = [spont_Vm_hist_ChAT_pos(1,1)+0.22, spont_Vm_hist_ChAT_pos(1,2)+0.02,0.08,0.04];
% h_dist1=spont_Vm_hist_pos(1,1)+spont_Vm_hist_pos(1,3)+0.14;
spont_Vm_5prcentile_ChAT_pos(1,:) = [h_dist1 , spont_Vm_hist_ChAT_pos(1,2) ,  0.1 ,  0.08];
spont_Rin_ChAT_pos(1,:) = [h_dist1+spont_Vm_5prcentile_ChAT_pos(1,3)+0.12, spont_Vm_5prcentile_ChAT_pos(1,2) ,  spont_Vm_5prcentile_ChAT_pos(1,3) ,  spont_Vm_5prcentile_ChAT_pos(1,4)];

detection_example1_ChAT_pos(1,:) = [0.05 , 0.5 , 0.45 , 0.22]; detection_example1_ChAT_pos(1,2)=spont_Vm_hist_ChAT_pos(1,2)-detection_example1_ChAT_pos(1,4)+0.05;
hdist3=detection_example1_ChAT_pos(1,1)+detection_example1_ChAT_pos(1,3)+0.05;
detection_example2_ChAT_pos(1,:) = [hdist3 , detection_example1_ChAT_pos(1,2) , detection_example1_ChAT_pos(1,3) , detection_example1_ChAT_pos(1,4)]; %detection_example2_pos(1,2)=detection_example1_pos(1,2)-detection_example1_pos(1,4)+0.05;
detection_example3_ChAT_pos(1,:) = [detection_example1_ChAT_pos(1,1) , 0.5 , detection_example1_ChAT_pos(1,3) , detection_example1_ChAT_pos(1,4)]; detection_example3_ChAT_pos(1,2)=detection_example1_ChAT_pos(1,2)-detection_example3_ChAT_pos(1,4)+0.15;
detection_example4_ChAT_pos(1,:) = [hdist3 , detection_example3_ChAT_pos(1,2) , detection_example1_ChAT_pos(1,3) , detection_example1_ChAT_pos(1,4)]; %detection_example4_pos(1,2)=detection_example3_pos(1,2)-detection_example3_pos(1,4)+0.05;

spont_amp_hist_ChAT_pos(1,:) = [spont_Vm_hist_ChAT_pos(1,1), 0.5 ,  spont_Vm_hist_ChAT_pos(1,3) ,  0.08]; spont_amp_hist_ChAT_pos(1,2)=detection_example4_ChAT_pos(1,2)-spont_amp_hist_ChAT_pos(1,4)-0.01;
spont_amp_median_ChAT_pos(1,:) = [spont_amp_hist_ChAT_pos(1,1)+0.24, spont_amp_hist_ChAT_pos(1,2)+0.04,0.05,0.05];
h_dist2=spont_amp_hist_pos(1,1)+spont_amp_hist_pos(1,3)+0.15;
spont_event_freq_ChAT_pos(1,:) =[spont_Vm_5prcentile_ChAT_pos(1,1),spont_amp_hist_ChAT_pos(1,2)-0.02, spont_Vm_5prcentile_ChAT_pos(1,3),spont_Vm_5prcentile_ChAT_pos(1,4)];
spont_event_amp_ChAT_pos(1,:)=[spont_event_freq_ChAT_pos(1,1)+spont_event_freq_ChAT_pos(1,3)+0.1,spont_event_freq_ChAT_pos(1,2),spont_event_freq_ChAT_pos(1,3),spont_event_freq_ChAT_pos(1,4)];

%top positions:
spont_Vm_hist_pos_top = spont_Vm_hist_pos(1,2)+spont_Vm_hist_pos(1,4);
spont_Vm_5prcentile_pos_top = spont_Vm_5prcentile_pos(1,2)+spont_Vm_5prcentile_pos(1,4);
spont_event_freq_pos_top = spont_event_freq_pos(1,2)+spont_event_freq_pos(1,4);
spont_event_amp_pos_top = spont_event_amp_pos(1,2)+spont_event_amp_pos(1,4);
spont_Rin_pos_top = spont_Rin_pos(1,2)+spont_Rin_pos(1,4);
detection_example1_pos_top = detection_example1_pos(1,2)+detection_example1_pos(1,4);
detection_example2_pos_top = detection_example2_pos(1,2)+detection_example2_pos(1,4);
detection_example3_pos_top = detection_example3_pos(1,2)+detection_example3_pos(1,4);
detection_example4_pos_top = detection_example4_pos(1,2)+detection_example4_pos(1,4);
spont_amp_hist_pos_top = spont_amp_hist_pos(1,2)+spont_amp_hist_pos(1,4);
%ChAT:
spont_Vm_hist_ChAT_pos_top = spont_Vm_hist_ChAT_pos(1,2)+spont_Vm_hist_ChAT_pos(1,4);
spont_Vm_5prcentile_ChAT_pos_top = spont_Vm_5prcentile_ChAT_pos(1,2)+spont_Vm_5prcentile_ChAT_pos(1,4);
spont_event_freq_ChAT_pos_top = spont_event_freq_ChAT_pos(1,2)+spont_event_freq_ChAT_pos(1,4);
spont_event_amp_ChAT_pos_top = spont_event_amp_ChAT_pos(1,2)+spont_event_amp_ChAT_pos(1,4);
spont_Rin_ChAT_pos_top = spont_Rin_ChAT_pos(1,2)+spont_Rin_ChAT_pos(1,4);
detection_example1_ChAT_pos_top = detection_example1_ChAT_pos(1,2)+detection_example1_ChAT_pos(1,4);
detection_example2_ChAT_pos_top = detection_example2_ChAT_pos(1,2)+detection_example2_ChAT_pos(1,4);
detection_example3_ChAT_pos_top = detection_example3_ChAT_pos(1,2)+detection_example3_ChAT_pos(1,4);
detection_example4_ChAT_pos_top = detection_example4_ChAT_pos(1,2)+detection_example4_ChAT_pos(1,4);
spont_amp_hist_ChAT_pos_top = spont_amp_hist_ChAT_pos(1,2)+spont_amp_hist_ChAT_pos(1,4);

% spont_event_halfwidth_pos_top = spont_event_halfwidth_pos(1,2)+spont_event_halfwidth_pos(1,4);
% spont_event_risetime_pos_top = spont_event_risetime_pos(1,2)+spont_event_risetime_pos(1,4);

%%
%Placing plots in the figure:
ax_fontsize=11;
color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];        
spont_Vm_hist_ax_copy = copyobj(spont_Vm_hist_ax,F); % copy axes to new fig
set(spont_Vm_hist_ax_copy(2),'position',spont_Vm_hist_pos(1,:))
set(F, 'currentaxes', spont_Vm_hist_ax_copy(2));  tl=title('');
spont_Vm_hist_ax_copy(2).FontSize=ax_fontsize;
set(gca,'tickdir','out');
% spont_Vm_hist_ax_copy.XTickLabel=[];

%placing legend:
 hVm=get(gca,'Children');
[l1,OBJH1,OUTH1,OUTM1] = legend([hVm(1) hVm(2)],{'NB-','NB+'}, 'position',[spont_Vm_hist_pos(1,1)+spont_Vm_hist_pos(1,3)-0.13,spont_Vm_hist_pos_top-0.02, 0.2,0.02],'fontsize',10, 'box', 'off'); % returns a handle LEGH to the legend axes; a vector OBJH containing handles for the text, lines, and patches in the legend; a vector OUTH of handles to thelines and patches in the plot; and a cell array OUTM containingthe text in the legend.
lpatch1=findobj(OBJH1,'type','patch');
lpatch1(1).Vertices(2:3,2)=lpatch1(1).Vertices(1,2)+0.15;
lpatch1(2).Vertices(2:3,2)=lpatch1(2).Vertices(1,2)+0.15;
lpatch1(1).Vertices(:,2)=lpatch1(1).Vertices(:,2)+0.06;
lpatch1(2).Vertices(:,2)=lpatch1(2).Vertices(:,2)+0.08;
lpatch1(1).Vertices([1,2,5],1)=lpatch1(1).Vertices(3,1)-0.2;
lpatch1(2).Vertices([1,2,5],1)=lpatch1(1).Vertices(3,1)-0.2;
lpatch1(1).FaceColor=color_table(1,:);
lpatch1(1).FaceAlpha=0.7;
lpatch1(2).FaceColor=color_table(2,:);
lpatch1(2).FaceAlpha=0.3;

% spont_Vm_median_ax_copy = copyobj(spont_Vm_median_ax,F); % copy axes to new fig
% set(spont_Vm_median_ax_copy,'position',spont_Vm_median_pos(1,:))
% set(F, 'currentaxes', spont_Vm_median_ax_copy);  tl=title('');
% set(gca,'xlim',[0.9 2.1]);
% spont_Vm_median_ax_copy.FontSize=10;

spont_Vm_5prcentile_ax_copy = copyobj(spont_Vm_5prcentile_ax,F); % copy axes to new fig
set(spont_Vm_5prcentile_ax_copy,'position',spont_Vm_5prcentile_pos(1,:))
set(F, 'currentaxes', spont_Vm_5prcentile_ax_copy);   tl=title('');
yl=ylabel({'Lower 5th'; 'perc. [mV]'});
spont_Vm_5prcentile_ax_copy.FontSize=ax_fontsize;
set(gca,'tickdir','out');

spont_amp_hist_ax_copy = copyobj(spont_amp_hist_ax(2),F); % copy axes to new fig
set(spont_amp_hist_ax_copy,'position',spont_amp_hist_pos(1,:),'tickdir','out');
spont_amp_hist_ax_copy.FontSize=ax_fontsize;

spont_amp_median_ax_copy = copyobj(spont_amp_median_ax,F); % copy axes to new fig
set(spont_amp_median_ax_copy,'position',spont_amp_median_pos(1,:))
set(F, 'currentaxes', spont_amp_median_ax_copy);  tl=title('');
set(gca,'xlim',[0.7 2.3],'tickdir','out');
spont_amp_median_ax_copy.FontSize=8;

spont_event_freq_ax_copy = copyobj(spont_event_freq_ax,F); % copy axes to new fig
set(spont_event_freq_ax_copy,'position',spont_event_freq_pos(1,:))
set(F, 'currentaxes', spont_event_freq_ax_copy);  xl=xlabel('');  tl=title('');
yl=ylabel('Freq. [Hz]');
spont_event_freq_ax_copy.FontSize=ax_fontsize;
set(gca,'tickdir','out');
% spont_event_freq_ax_copy.XTickLabel=[];

  spont_event_amp_ax_copy = copyobj(spont_event_amp_ax,F); % copy axes to new fig
set(spont_event_amp_ax_copy,'position',spont_event_amp_pos(1,:))
set(F, 'currentaxes', spont_event_amp_ax_copy);  xl=xlabel('');  tl=title(''); 
yl=ylabel('Amp. [mV]');
spont_event_amp_ax_copy.FontSize=ax_fontsize;
set(gca,'tickdir','out');
% spont_event_amp_ax_copy.XTickLabel=[];

spont_Rin_ax_copy = copyobj(spont_Rin_ax,F); % copy axes to new fig
set(spont_Rin_ax_copy,'position',spont_Rin_pos(1,:))
set(F, 'currentaxes', spont_Rin_ax_copy);   tl=title(''); 
% yl=ylabel('Rin [MOhm]');
spont_Rin_ax_copy.FontSize=ax_fontsize;
set(gca,'tickdir','out');

detection_example1_ax_copy = copyobj(detection_example1_ax,F); % copy axes to new fig
set(detection_example1_ax_copy,'position',detection_example1_pos(1,:))

detection_example2_ax_copy = copyobj(detection_example2_ax,F); % copy axes to new fig
set(detection_example2_ax_copy,'position',detection_example2_pos(1,:))

detection_example3_ax_copy = copyobj(detection_example3_ax,F); % copy axes to new fig
set(detection_example3_ax_copy,'position',detection_example3_pos(1,:))

detection_example4_ax_copy = copyobj(detection_example4_ax,F); % copy axes to new fig
set(detection_example4_ax_copy,'position',detection_example4_pos(1,:))

set(findall(gcf,'type','scatter'),'sizedata',20);

%ChAT:
color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [102, 172,255]./256; [30,75,14]/256; [112,172,90]/256];  
spont_Vm_hist_ChAT_ax_copy = copyobj(spont_Vm_hist_ChAT_ax,F); % copy axes to new fig
set(spont_Vm_hist_ChAT_ax_copy(2),'position',spont_Vm_hist_ChAT_pos(1,:))
set(F, 'currentaxes', spont_Vm_hist_ChAT_ax_copy(2));  tl=title('');
spont_Vm_hist_ChAT_ax_copy(2).FontSize=ax_fontsize;
set(gca,'tickdir','out');
% spont_Vm_hist_ax_copy.XTickLabel=[];

%placing legend:

hVm=get(gca,'Children');
[l2,OBJH2,OUTH2,OUTM2] = legend([hVm(1) hVm(2)],{'Light Off','Light On'}, 'position',[spont_Vm_hist_ChAT_pos(1,1)+spont_Vm_hist_ChAT_pos(1,3)-0.14,spont_Vm_hist_ChAT_pos_top-0.02, 0.2,0.02],'fontsize',10, 'box', 'off'); % returns a handle LEGH to the legend axes; a vector OBJH containing handles for the text, lines, and patches in the legend; a vector OUTH of handles to thelines and patches in the plot; and a cell array OUTM containingthe text in the legend.
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

% spont_Vm_median_ChAT_ax_copy = copyobj(spont_Vm_median_ChAT_ax,F); % copy axes to new fig
% set(spont_Vm_median_ChAT_ax_copy,'position',spont_Vm_median_ChAT_pos(1,:))
% set(F, 'currentaxes', spont_Vm_median_ChAT_ax_copy);  tl=title('');
% set(gca,'xlim',[0.9 2.1]);
% spont_Vm_median_ChAT_ax_copy.FontSize=10;

spont_Vm_5prcentile_ChAT_ax_copy = copyobj(spont_Vm_5prcentile_ChAT_ax,F); % copy axes to new fig
set(spont_Vm_5prcentile_ChAT_ax_copy,'position',spont_Vm_5prcentile_ChAT_pos(1,:))
set(F, 'currentaxes', spont_Vm_5prcentile_ChAT_ax_copy);   tl=title('');
yl=ylabel({'Lower 5th'; 'perc. [mV]'});
spont_Vm_5prcentile_ChAT_ax_copy.FontSize=ax_fontsize;
spont_Vm_5prcentile_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out');

spont_amp_hist_ChAT_ax_copy = copyobj(spont_amp_hist_ChAT_ax(2),F); % copy axes to new fig
set(spont_amp_hist_ChAT_ax_copy,'position',spont_amp_hist_ChAT_pos(1,:))
spont_amp_hist_ChAT_ax_copy.FontSize=ax_fontsize;

spont_amp_median_ChAT_ax_copy = copyobj(spont_amp_median_ChAT_ax,F); % copy axes to new fig
set(spont_amp_median_ChAT_ax_copy,'position',spont_amp_median_ChAT_pos(1,:))
set(F, 'currentaxes', spont_amp_median_ChAT_ax_copy);  tl=title('');
set(gca,'xlim',[0.7 2.3],'xticklabel',{'Off','On'},'tickdir','out');
spont_amp_median_ChAT_ax_copy.FontSize=8;

spont_event_freq_ChAT_ax_copy = copyobj(spont_event_freq_ChAT_ax,F); % copy axes to new fig
set(spont_event_freq_ChAT_ax_copy,'position',spont_event_freq_ChAT_pos(1,:))
set(F, 'currentaxes', spont_event_freq_ChAT_ax_copy);  xl=xlabel('');  tl=title('');
yl=ylabel('Freq. [Hz]');
spont_event_freq_ChAT_ax_copy.FontSize=ax_fontsize;
spont_event_freq_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out');

  spont_event_amp_ChAT_ax_copy = copyobj(spont_event_amp_ChAT_ax,F); % copy axes to new fig
set(spont_event_amp_ChAT_ax_copy,'position',spont_event_amp_ChAT_pos(1,:))
set(F, 'currentaxes', spont_event_amp_ChAT_ax_copy);  xl=xlabel('');  tl=title(''); 
yl=ylabel('Amp. [mV]');
spont_event_amp_ChAT_ax_copy.FontSize=ax_fontsize;
spont_event_amp_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out');

spont_Rin_ChAT_ax_copy = copyobj(spont_Rin_ChAT_ax,F); % copy axes to new fig
set(spont_Rin_ChAT_ax_copy,'position',spont_Rin_ChAT_pos(1,:))
set(F, 'currentaxes', spont_Rin_ChAT_ax_copy);   tl=title(''); 
% yl=ylabel('Rin [MOhm]');
spont_Rin_ChAT_ax_copy.FontSize=ax_fontsize;
spont_Rin_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out');

detection_example1_ChAT_ax_copy = copyobj(detection_example1_ChAT_ax,F); % copy axes to new fig
set(detection_example1_ChAT_ax_copy,'position',detection_example1_ChAT_pos(1,:))

detection_example2_ChAT_ax_copy = copyobj(detection_example2_ChAT_ax,F); % copy axes to new fig
set(detection_example2_ChAT_ax_copy,'position',detection_example2_ChAT_pos(1,:))

detection_example3_ChAT_ax_copy = copyobj(detection_example3_ChAT_ax,F); % copy axes to new fig
set(detection_example3_ChAT_ax_copy,'position',detection_example3_ChAT_pos(1,:))

detection_example4_ChAT_ax_copy = copyobj(detection_example4_ChAT_ax,F); % copy axes to new fig
set(detection_example4_ChAT_ax_copy,'position',detection_example4_ChAT_pos(1,:))

set(findall(gcf,'type','scatter'),'sizedata',20);

color_text=[102,51,0]./256;
 a_pos1=[-0.06 0 0.04 0.04];
%  annotation('textbox', [spont_Vm_hist_pos(1,1),spont_Vm_hist_pos_top+0.03 0.4 0.05],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous activity', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
annotation('textbox', [detection_example1_pos(1,1)+0.09,detection_example1_pos_top-0.05, 0.1 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'NB-', 'FontName','arial', 'fontsize', 12, 'color','k');
 annotation('textbox', [detection_example2_pos(1,1)+0.09,detection_example1_pos_top-0.05, 0.1 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'NB+', 'FontName','arial', 'fontsize', 12, 'color','k');
 
 annotation('textbox', [detection_example1_pos(1,1)-0.05,detection_example1_pos_top-0.08, 0.1 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Cell#1', 'FontName','arial', 'fontsize', 12, 'color','k');
 annotation('textbox', [detection_example3_pos(1,1)-0.05,detection_example3_pos_top-0.08, 0.1 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Cell#2', 'FontName','arial', 'fontsize', 12, 'color','k');
 
 annotation('textbox', [detection_example1_ChAT_pos(1,1)+0.05,detection_example1_ChAT_pos_top-0.14, 0.2 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Light Off', 'FontName','arial', 'fontsize', 12, 'color','k');
 annotation('textbox', [detection_example2_ChAT_pos(1,1)+0.05,detection_example1_ChAT_pos_top-0.14, 0.2 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Light On', 'FontName','arial', 'fontsize', 12, 'color','k');

 annotation('textbox', [detection_example1_ChAT_pos(1,1)-0.05,detection_example1_ChAT_pos_top-0.15, 0.1 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Cell#3', 'FontName','arial', 'fontsize', 12, 'color','k');
 annotation('textbox', [detection_example3_ChAT_pos(1,1)-0.05,detection_example3_ChAT_pos_top-0.18, 0.1 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Cell#4', 'FontName','arial', 'fontsize', 12, 'color','k');

 an1pos=[spont_Vm_hist_pos(1,1),spont_Vm_hist_pos_top 0.04 0.04]+[-0.06, 0, 0, 0];
if no_numbering_flag==0;
annotation('textbox',an1pos ,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold');
annotation('textbox', [spont_Vm_5prcentile_pos(1,1)-0.06, an1pos(1,2) 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [spont_Rin_pos(1,1)-0.06, an1pos(1,2) 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [ an1pos(1,1), detection_example1_pos_top-0.04 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [detection_example2_pos(1,1)-0.03 detection_example1_pos_top-0.04 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
  annotation('textbox', [ an1pos(1,1), detection_example3_pos_top-0.12 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [detection_example4_pos(1,1)-0.03 detection_example3_pos_top-0.12 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [an1pos(1,1), spont_event_freq_pos_top+0.02 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_event_freq_pos(1,1)-0.06 spont_event_freq_pos_top+0.02 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_event_amp_pos(1,1)-0.06 spont_event_amp_pos_top+0.02 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
end 
%  annotation('textbox', [spont_event_halfwidth_pos(1,1) spont_event_halfwidth_pos_top 0 0]+a_pos1,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
%   annotation('textbox', [spont_event_risetime_pos(1,1) spont_event_risetime_pos_top 0 0]+a_pos1,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
close(1:20)
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'

if save_flag==1   
    filename='Fig 3_f46_f80_v2';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end