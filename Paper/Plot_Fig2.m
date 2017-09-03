%% plotting a figure 

%% new version

close all
clear all
save_flag=0;
no_numbering_flag=0;
abslen=0.05; %[in cm]
ax_fontsize=10;
color1=[255, 102,102]./256; %pink
color2=[102, 172,255]./256; %light blue
 
%opening saved figures:
%Long traces NB:
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Long Trace Presentation\style2'
%spontaneous:
spont_trace_f46 = open('f46_traces_x1_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
delete(b1(2,1)); 
delete(b2(2,1));
spont_trace_f46_ax = get(gcf, 'children');

spont_mean_f46 = open('f46_mean_x1_v2.fig');    
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
spont_mean_f46_ax = get(gcf, 'children');

spont_std_f46 = open('f46_std_x1_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
delete(b1(2,1)); 
delete(b2(2,1));
spont_std_f46_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm Histograms and Powerspec\Spontaneous histograms'
  
    spont_Vm_hist = open('f46_Vm_histogram_ongoing.fig');
    spont_Vm_hist_ax = get(gcf, 'children'); 

    spont_Vm_5prcentile = open('Vm_5prcentile_ongoing.fig');   
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(10),'color',color1); uistack(l1(10),'top')
    spont_Vm_5prcentile_ax = get(gcf, 'children');
    
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Rin'
    spont_Rin = open('Rin.fig');    
    spont_Rin_ax = get(gcf, 'children');

%% long traces ChAT

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Long Trace Presentation\style2'
%spontaneous:
spont_trace_f80 = open('f80_traces_x1_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
delete(b1(2,1)); 
delete(b2(2,1));
spont_trace_f80_ax = get(gcf, 'children');

spont_mean_f80 = open('f80_mean_x1_v2.fig');   
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
spont_mean_f80_ax = get(gcf, 'children');

spont_std_f80 = open('f80_std_x1_mean-subt_v2.fig');    
b1=findall(gcf,'type','line');
b2=findall(gcf,'type','text');
set(b2,'fontsize',10);
delete(b1(2,1)); 
delete(b2(2,1));
spont_std_f80_ax = get(gcf, 'children');


cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\LFP_50Hz Vm_50Hz';

VmM_NB = open('Vm M_Before Sensory stim.fig');   
            l1=findall(gcf,'color',[0.7 0.7 0.7]);
            set(l1(10),'color',color1); uistack(l1(10),'top')
VmM_NB_ax = get(gcf, 'children');

 
VmSTD_NB = open('Vm STD_Before Sensory stim.fig');   
        l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(10),'color',color1); uistack(l1(10),'top')
VmSTD_NB_ax = get(gcf, 'children');  

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary\LFP_50Hz Vm_50Hz';
VmM_ChAT = open('Vm M_Before Sensory stim.fig');    
         l1=findall(gcf,'color',[0.7 0.7 0.7]);
         set(l1(9),'color',color2); uistack(l1(9),'top')
VmM_ChAT_ax = get(gcf, 'children');

VmSTD_ChAT = open('Vm STD_Before Sensory stim.fig');    
         l1=findall(gcf,'color',[0.7 0.7 0.7]);
         set(l1(9),'color',color2); uistack(l1(9),'top')
VmSTD_ChAT_ax = get(gcf, 'children');

 cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm Histograms and Powerspec\Spontaneous histograms'

        spont_Vm_hist_ChAT = open('f80_Vm_histogram_ongoing.fig');
        spont_Vm_hist_ChAT_ax = get(gcf, 'children');

        spont_Vm_5prcentile_ChAT = open('Vm_5prcentile_ongoing.fig'); 
            l1=findall(gcf,'color',[0.7 0.7 0.7]);
            set(l1(9),'color',color2); uistack(l1(9),'top')
        spont_Vm_5prcentile_ChAT_ax = get(gcf, 'children');
        
cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Rin'

    spont_Rin_ChAT = open('Rin.fig');    
    spont_Rin_ChAT_ax = get(gcf, 'children');


cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
NB_illustration = imread('NBES_Schematic_Illustration_one-electrode_tight','TIF');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
ChAT_illustration = imread('ChAT_Schematic_Illustration_one-electrode_tight','TIF');
% ChAT_histology = imread('ChAT_Histology_tight','TIF');


%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',11);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 27]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

%% Positions: (reduce 20% of hight)
    
v_dist=0.01;
spont_trace_f46_pos(1,:) = [0.15 , 0.79 ,  0.6 ,  0.14];
spont_mean_f46_pos(1,:) = [spont_trace_f46_pos(1,1), spont_trace_f46_pos(1,2)-0.06-v_dist, spont_trace_f46_pos(1,3), 0.06];
spont_std_f46_pos(1,:) = [spont_trace_f46_pos(1,1), spont_mean_f46_pos(1,2)-spont_mean_f46_pos(1,4)-v_dist, spont_trace_f46_pos(1,3), spont_mean_f46_pos(1,4)];

%NB
spont_Vm_hist_pos(1,:)=[0.2 , 0.5 ,  0.27 ,  0.06]; spont_Vm_hist_pos(1,2)=spont_std_f46_pos(1,2)-spont_Vm_hist_pos(1,4)-0.04;
h_dist1=spont_Vm_hist_pos(1,1)+spont_Vm_hist_pos(1,3)+0.16;
spont_Vm_5prcentile_pos(1,:) = [h_dist1 , spont_Vm_hist_pos(1,2)-0.02 ,  0.1 ,  0.09];
spont_Rin_pos(1,:) = [h_dist1+spont_Vm_5prcentile_pos(1,3)+0.12, spont_Vm_5prcentile_pos(1,2) ,  spont_Vm_5prcentile_pos(1,3) ,  spont_Vm_5prcentile_pos(1,4)];

spont_trace_f80_pos(1,:) =[spont_trace_f46_pos(1,1) , 0.5 ,  spont_trace_f46_pos(1,3) ,  spont_trace_f46_pos(1,4)]; spont_trace_f80_pos(1,2)=spont_Vm_hist_pos(1,2)-spont_trace_f80_pos(1,4)-0.11;
spont_mean_f80_pos(1,:) =[spont_trace_f80_pos(1,1), spont_trace_f80_pos(1,2)-0.06-v_dist, spont_trace_f80_pos(1,3), 0.06];
spont_std_f80_pos(1,:) =[spont_trace_f80_pos(1,1), spont_mean_f80_pos(1,2)-spont_mean_f80_pos(1,4)-v_dist, spont_trace_f80_pos(1,3), spont_mean_f80_pos(1,4)];

spont_Vm_hist_ChAT_pos(1,:)=[spont_Vm_hist_pos(1,1), 0.5, 0.27, 0.06]; spont_Vm_hist_ChAT_pos(1,2)=spont_std_f80_pos(1,2)-spont_Vm_hist_ChAT_pos(1,4)-0.04;
spont_Vm_5prcentile_ChAT_pos(1,:) = [h_dist1 , spont_Vm_hist_ChAT_pos(1,2)-0.02 ,  0.1 ,  0.09];
spont_Rin_ChAT_pos(1,:) = [h_dist1+spont_Vm_5prcentile_ChAT_pos(1,3)+0.12, spont_Vm_5prcentile_ChAT_pos(1,2) ,  spont_Vm_5prcentile_ChAT_pos(1,3) ,  spont_Vm_5prcentile_ChAT_pos(1,4)];

h_dist1=spont_trace_f46_pos(1,1)+spont_trace_f46_pos(1,3)+0.07;

v_dist2=0.05;
VmM_NB_pos(1,:)=[h_dist1+0.03, spont_trace_f46_pos(1,2)+0.04, 0.1, 0.09]; 
VmSTD_NB_pos(1,:)= VmM_NB_pos; VmSTD_NB_pos(1,2)=VmM_NB_pos(1,2)-VmM_NB_pos(1,4)-0.07;
VmM_ChAT_pos(1,:)=[h_dist1+0.03, spont_trace_f80_pos(1,2)+0.04, VmM_NB_pos(1,3), VmM_NB_pos(1,4)]; 
VmSTD_ChAT_pos(1,:)= VmM_ChAT_pos; VmSTD_ChAT_pos(1,2)=VmM_ChAT_pos(1,2)-VmM_ChAT_pos(1,4)-0.07;

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
spont_Vm_hist_pos_top=spont_Vm_hist_pos(1,2)+spont_Vm_hist_pos(1,4);
spont_Vm_5prcentile_pos_top=spont_Vm_5prcentile_pos(1,2)+spont_Vm_5prcentile_pos(1,4);
spont_Rin_pos_top=spont_Rin_pos(1,2)+spont_Rin_pos(1,4);
spont_Vm_hist_ChAT_pos_top=spont_Vm_hist_ChAT_pos(1,2)+spont_Vm_hist_ChAT_pos(1,4);
spont_Vm_5prcentile_ChAT_pos_top=spont_Vm_5prcentile_ChAT_pos(1,2)+spont_Vm_5prcentile_ChAT_pos(1,4);
spont_Rin_ChAT_pos_top=spont_Rin_ChAT_pos(1,2)+spont_Rin_ChAT_pos(1,4);

NB_illustration_pos(1,:)=[0.02,spont_trace_f46_pos_top-0.05, 0.11, 0.1];
ChAT_illustration_pos(1,:)=[NB_illustration_pos(1,1),spont_trace_f80_pos_top-0.05, NB_illustration_pos(1,3), NB_illustration_pos(1,4)];

NB_illustration_pos_top=NB_illustration_pos(1,2)+NB_illustration_pos(1,4);
ChAT_illustration_pos_top=ChAT_illustration_pos(1,2)+ChAT_illustration_pos(1,4);

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

%position legend:
set(spont_mean_f46_ax_copy(1),'position',[0.65 spont_trace_f46_pos_top+0.02 0.08 0.005])
spont_mean_f46_ax_copy(1).FontSize=9;

set(spont_mean_f80_ax_copy(1),'position',[0.65 spont_trace_f80_pos_top 0.08 0.005])
spont_mean_f80_ax_copy(1).FontSize=9;

%population panels:
VmM_NB_ax_copy = copyobj(VmM_NB_ax,F); % copy axes to new fig
set(VmM_NB_ax_copy(1),'position',VmM_NB_pos(1,:))
set(F, 'currentaxes', VmM_NB_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
VmM_NB_ax_copy(1).FontSize=ax_fontsize;
set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);
% VmM_NB_ax_copy(1).XTickLabel=[];

VmSTD_NB_ax_copy = copyobj(VmSTD_NB_ax(1),F); % copy axes to new fig
set(VmSTD_NB_ax_copy,'position',VmSTD_NB_pos(1,:))
set(F, 'currentaxes', VmSTD_NB_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
VmSTD_NB_ax_copy.FontSize=ax_fontsize;
set(gca,'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);
% VmSTD_NB_ax_copy.XTickLabel=[];

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
color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];        
spont_Vm_hist_ax_copy = copyobj(spont_Vm_hist_ax,F); % copy axes to new fig
set(spont_Vm_hist_ax_copy(2),'position',spont_Vm_hist_pos(1,:))
set(F, 'currentaxes', spont_Vm_hist_ax_copy(2));  tl=title('');
yl=ylabel({'Prob.'},'fontsize',ax_fontsize);
xl=xlabel('Vm (mV)','fontsize',ax_fontsize);
spont_Vm_hist_ax_copy(2).FontSize=ax_fontsize;
set(gca,'tickdir','out');
ticklen=fn_get_abs_ticklength(gca, abslen);

% spont_Vm_hist_ax_copy.XTickLabel=[];

%placing legend:
 hVm=get(gca,'Children');
[l1,OBJH1,OUTH1,OUTM1] = legend([hVm(1) hVm(2)],{'NB-','NB+'}, 'position',[spont_Vm_hist_pos(1,1)+spont_Vm_hist_pos(1,3)-0.13,spont_Vm_hist_pos_top-0.02, 0.2,0.02],'fontsize',9, 'box', 'off'); % returns a handle LEGH to the legend axes; a vector OBJH containing handles for the text, lines, and patches in the legend; a vector OUTH of handles to thelines and patches in the plot; and a cell array OUTM containingthe text in the legend.
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

spont_Vm_5prcentile_ax_copy = copyobj(spont_Vm_5prcentile_ax,F); % copy axes to new fig
set(spont_Vm_5prcentile_ax_copy,'position',spont_Vm_5prcentile_pos(1,:))
set(F, 'currentaxes', spont_Vm_5prcentile_ax_copy);   tl=title('');
yl=ylabel({'Lower 5th'; 'perc. (mV)'},'fontsize',ax_fontsize);
spont_Vm_5prcentile_ax_copy.FontSize=ax_fontsize;
set(gca,'tickdir','out');
ticklen=fn_get_abs_ticklength(gca, abslen);

spont_Rin_ax_copy = copyobj(spont_Rin_ax,F); % copy axes to new fig
set(spont_Rin_ax_copy,'position',spont_Rin_pos(1,:))
set(F, 'currentaxes', spont_Rin_ax_copy);   tl=title(''); 
yl=ylabel('Rin (MOhm)','fontsize',ax_fontsize);
spont_Rin_ax_copy.FontSize=ax_fontsize;
set(gca,'tickdir','out');
ticklen=fn_get_abs_ticklength(gca, abslen);

color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [102, 172,255]./256; [30,75,14]/256; [112,172,90]/256];  
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
yl=ylabel({'Lower 5th'; 'perc. (mV)'},'fontsize',ax_fontsize);
spont_Vm_5prcentile_ChAT_ax_copy.FontSize=ax_fontsize;
spont_Vm_5prcentile_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out');
ticklen=fn_get_abs_ticklength(gca, abslen);

spont_Rin_ChAT_ax_copy = copyobj(spont_Rin_ChAT_ax,F); % copy axes to new fig
set(spont_Rin_ChAT_ax_copy,'position',spont_Rin_ChAT_pos(1,:))
set(F, 'currentaxes', spont_Rin_ChAT_ax_copy);   tl=title(''); 
yl=ylabel('Rin (MOhm)','fontsize',ax_fontsize);
spont_Rin_ChAT_ax_copy.FontSize=ax_fontsize;
spont_Rin_ChAT_ax_copy.XTickLabel={'Off','On'};
set(gca,'tickdir','out');
ticklen=fn_get_abs_ticklength(gca, abslen);


NB_illustration_ax = axes('position',NB_illustration_pos);
imshow(NB_illustration, 'parent', NB_illustration_ax) 

ChAT_illustration_ax = axes('position',ChAT_illustration_pos);
imshow(ChAT_illustration, 'parent', ChAT_illustration_ax) 

%   annotation:

annotation('textbox', [spont_trace_f46_pos(1,1), spont_trace_f46_pos_top 0 0]+[0 0.015 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous Activity', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold') %
 
 a_pos1=[-0.03 -0.01 0.04 0.04];
 a_pos4=[0 -0.01 0.16 0.04];
 a_pos2=[-0.02 -0.01 0.04 0.04];
 a_pos3=[-0.035 0 0.04 0.04];
 annotation('textbox', [spont_trace_f46_pos(1,1)+0.14 spont_trace_f46_pos_top+0.04 0.4 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'NB electrical stimulation', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
 annotation('textbox', [spont_trace_f46_pos(1,1) spont_trace_f46_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Vm Traces', 'FontName','arial', 'fontsize', 11) %'fontweight', 'bold'
annotation('textbox', [spont_mean_f46_pos(1,1) spont_mean_f46_pos_top 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 11)
annotation('textbox', [spont_std_f46_pos(1,1) spont_std_f46_pos_top 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 11)
annotation('textbox', [spont_trace_f80_pos(1,1)+0.15 spont_trace_f80_pos_top+0.04 0.4 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Optogenetic stimulation', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
annotation('textbox', [spont_trace_f80_pos(1,1) spont_trace_f80_pos_top-0.016 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Vm Traces', 'FontName','arial', 'fontsize', 11) %'fontweight', 'bold'
annotation('textbox', [spont_mean_f80_pos(1,1) spont_mean_f80_pos_top-0.0085 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 11)
annotation('textbox', [spont_std_f80_pos(1,1) spont_std_f80_pos_top-0.009 0 0]+a_pos4,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 11)
    
annotation('textbox', [VmM_NB_pos(1,1) VmM_NB_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 10)
annotation('textbox', [VmSTD_NB_pos(1,1) VmSTD_NB_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 10)
annotation('textbox', [VmM_ChAT_pos(1,1) VmM_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 10)
annotation('textbox', [VmSTD_ChAT_pos(1,1) VmSTD_ChAT_pos_top+0.005 0.16 0.04],...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 10)

 if no_numbering_flag==0;
 annotation('textbox', [NB_illustration_pos(1,1)-0.01 NB_illustration_pos_top-0.02 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold') 
annotation('textbox', [spont_trace_f46_pos(1,1) spont_trace_f46_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_mean_f46_pos(1,1) spont_mean_f46_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_std_f46_pos(1,1) spont_std_f46_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [VmM_NB_pos(1,1) VmM_NB_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [VmSTD_NB_pos(1,1) VmSTD_NB_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_Vm_hist_pos(1,1) spont_Vm_hist_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_Vm_5prcentile_pos(1,1) spont_Vm_5prcentile_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_Rin_pos(1,1) spont_Rin_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')

annotation('textbox', [ChAT_illustration_pos(1,1)-0.01 ChAT_illustration_pos_top-0.02 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_trace_f80_pos(1,1) spont_trace_f80_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'K', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_mean_f80_pos(1,1) spont_mean_f80_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'L', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_std_f80_pos(1,1) spont_std_f80_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'M', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [VmM_ChAT_pos(1,1) VmM_ChAT_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'N', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [VmSTD_ChAT_pos(1,1) VmSTD_ChAT_pos_top 0 0]+a_pos3,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'O', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_Vm_hist_ChAT_pos(1,1) spont_Vm_hist_ChAT_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'P', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_Vm_5prcentile_ChAT_pos(1,1) spont_Vm_5prcentile_ChAT_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Q', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
annotation('textbox', [spont_Rin_ChAT_pos(1,1) spont_Rin_ChAT_pos_top 0 0]+a_pos1,...
    'FitHeightToText', 'on', 'edgecolor', 'none','string', 'R', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 end

      cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'
  
if save_flag==1
    filename='Fig 2 intracellular traces spont_f46_f80_v4';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end

%% old version:
% close all
% clear all
% save_flag=0;
%  no_numbering_flag=0;
%  abslen=0.05; %[in cm]
%  ax_fontsize=10;
%  color1=[255, 102,102]./256; %pink
% color2=[102, 172,255]./256; %light blue
%  
% %opening saved figures:
% %Long traces NB:
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Long Trace Presentation\style2'
% %spontaneous:
% spont_trace_f46 = open('f46_traces_x1_v2.fig');    
% b1=findall(gcf,'type','line');
% b2=findall(gcf,'type','text');
% delete(b1(2,1)); 
% delete(b2(2,1));
% spont_trace_f46_ax = get(gcf, 'children');
% 
% spont_mean_f46 = open('f46_mean_x1_v2.fig');    
% spont_mean_f46_ax = get(gcf, 'children');
% 
% spont_std_f46 = open('f46_std_x1_mean-subt_v2.fig');    
% b1=findall(gcf,'type','line');
% b2=findall(gcf,'type','text');
% delete(b1(2,1)); 
% delete(b2(2,1));
% spont_std_f46_ax = get(gcf, 'children');
% 
% %% long traces ChAT
% 
% cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Long Trace Presentation\style2'
% %spontaneous:
% spont_trace_f80 = open('f80_traces_x1_v2.fig');    
% b1=findall(gcf,'type','line');
% b2=findall(gcf,'type','text');
% delete(b1(2,1)); 
% delete(b2(2,1));
% spont_trace_f80_ax = get(gcf, 'children');
% 
% spont_mean_f80 = open('f80_mean_x1_v2.fig');    
% spont_mean_f80_ax = get(gcf, 'children');
% 
% spont_std_f80 = open('f80_std_x1_mean-subt_v2.fig');    
% b1=findall(gcf,'type','line');
% b2=findall(gcf,'type','text');
% delete(b1(2,1)); 
% delete(b2(2,1));
% spont_std_f80_ax = get(gcf, 'children');
% 
% 
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\LFP_50Hz Vm_50Hz';
% 
% VmM_NB = open('Vm M_Before Sensory stim.fig');   
%             l1=findall(gcf,'color',[0.7 0.7 0.7]);
%             set(l1(10),'color',color1); uistack(l1(10),'top')
% VmM_NB_ax = get(gcf, 'children');
% 
%  
% VmSTD_NB = open('Vm STD_Before Sensory stim.fig');   
%         l1=findall(gcf,'color',[0.7 0.7 0.7]);
%         set(l1(10),'color',color1); uistack(l1(10),'top')
% VmSTD_NB_ax = get(gcf, 'children');  
% 
% cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary\LFP_50Hz Vm_50Hz';
% VmM_ChAT = open('Vm M_Before Sensory stim.fig');    
%          l1=findall(gcf,'color',[0.7 0.7 0.7]);
%          set(l1(9),'color',color2); uistack(l1(9),'top')
% VmM_ChAT_ax = get(gcf, 'children');
% 
% VmSTD_ChAT = open('Vm STD_Before Sensory stim.fig');    
%          l1=findall(gcf,'color',[0.7 0.7 0.7]);
%          set(l1(9),'color',color2); uistack(l1(9),'top')
% VmSTD_ChAT_ax = get(gcf, 'children');
% 
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
% NB_illustration = imread('NBES_Schematic_Illustration_one-electrode_tight','TIF');
% 
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
% ChAT_illustration = imread('ChAT_Schematic_Illustration_one-electrode_tight','TIF');
% % ChAT_histology = imread('ChAT_Histology_tight','TIF');
% 
% 
% %%
% %open a new figure:
% F = figure;
% set(gcf,'color','w');
% set(gcf,'DefaultAxesFontSize',11);
% set(gcf,'DefaultAxesFontName','arial');
% set(gcf, 'PaperType', 'A4');
% set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 20]); %[left, bottom, width, height] 
% set(gcf,'PaperOrientation','portrait');
% set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
% % annotation('textbox', [freq_1_zoom1_pos(1,1) freq_1_zoom1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
% %            'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A',...
% %            'FontName','arial', 'fontsize', 11, 'fontweight', 'bold')
% % annotation('textbox', [freq_1_zoom2_pos(1,1) freq_1_zoom2_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
% %            'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B',...
% %            'FontName','arial', 'fontsize', 11, 'fontweight', 'bold')    
% % annotation('textbox', [freq_2_zoom1_pos(1,1) freq_2_zoom1_pos_top 0 0]+[-0.06 0.01 0.04 0.04],...
% %            'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C',...
% %            'FontName','arial', 'fontsize', 11, 'fontweight', 'bold')  
% %% xlimits, y limits, ticks etc.
%  
% %% Positions:
%     
% v_dist=0.01;
% spont_trace_f46_pos(1,:) = [0.15 , 0.73 ,  0.6 ,  0.2];
% spont_mean_f46_pos(1,:) = [spont_trace_f46_pos(1,1), spont_trace_f46_pos(1,2)-0.10-v_dist, spont_trace_f46_pos(1,3), 0.10];
% spont_std_f46_pos(1,:) = [spont_trace_f46_pos(1,1), spont_mean_f46_pos(1,2)-spont_mean_f46_pos(1,4)-v_dist, spont_trace_f46_pos(1,3), spont_mean_f46_pos(1,4)];
% 
% spont_trace_f80_pos(1,:) =[spont_trace_f46_pos(1,1) , 0.5 ,  spont_trace_f46_pos(1,3) ,  spont_trace_f46_pos(1,4)]; spont_trace_f80_pos(1,2)=spont_std_f46_pos(1,2)-spont_trace_f80_pos(1,4)-0.07;
% spont_mean_f80_pos(1,:) =[spont_trace_f80_pos(1,1), spont_trace_f80_pos(1,2)-0.10-v_dist, spont_trace_f80_pos(1,3), 0.10];
% spont_std_f80_pos(1,:) =[spont_trace_f80_pos(1,1), spont_mean_f80_pos(1,2)-spont_mean_f80_pos(1,4)-v_dist, spont_trace_f80_pos(1,3), spont_mean_f80_pos(1,4)];
% 
% h_dist1=spont_trace_f46_pos(1,1)+spont_trace_f46_pos(1,3)+0.07;
% 
% v_dist2=0.05;
% VmM_NB_pos(1,:)=[h_dist1+0.03, spont_trace_f46_pos(1,2)+0.03, 0.12, 0.13]; 
% VmSTD_NB_pos(1,:)= VmM_NB_pos; VmSTD_NB_pos(1,2)=VmM_NB_pos(1,2)-VmM_NB_pos(1,4)-0.1;
% VmM_ChAT_pos(1,:)=[h_dist1+0.03, spont_trace_f80_pos(1,2)+0.03, VmM_NB_pos(1,3), VmM_NB_pos(1,4)]; 
% VmSTD_ChAT_pos(1,:)= VmM_ChAT_pos; VmSTD_ChAT_pos(1,2)=VmM_ChAT_pos(1,2)-VmM_ChAT_pos(1,4)-0.1;
% 
% %position of top of each panel
% spont_trace_f46_pos_top = spont_trace_f46_pos(1,2)+spont_trace_f46_pos(1,4);
% spont_mean_f46_pos_top = spont_mean_f46_pos(1,2)+spont_mean_f46_pos(1,4);
% spont_std_f46_pos_top = spont_std_f46_pos(1,2)+spont_std_f46_pos(1,4);
% spont_trace_f80_pos_top = spont_trace_f80_pos(1,2)+spont_trace_f80_pos(1,4);
% spont_mean_f80_pos_top = spont_mean_f80_pos(1,2)+spont_mean_f80_pos(1,4);
% spont_std_f80_pos_top = spont_std_f80_pos(1,2)+spont_std_f80_pos(1,4);
% VmM_NB_pos_top = VmM_NB_pos(1,2)+VmM_NB_pos(1,4);
% VmSTD_NB_pos_top = VmSTD_NB_pos(1,2)+VmSTD_NB_pos(1,4);
% VmM_ChAT_pos_top = VmM_ChAT_pos(1,2)+VmM_ChAT_pos(1,4);
% VmSTD_ChAT_pos_top = VmSTD_ChAT_pos(1,2)+VmSTD_ChAT_pos(1,4);
% 
% NB_illustration_pos(1,:)=[0.02,spont_trace_f46_pos_top-0.05, 0.11, 0.12];
% ChAT_illustration_pos(1,:)=[NB_illustration_pos(1,1),spont_trace_f80_pos_top-0.05, NB_illustration_pos(1,3), NB_illustration_pos(1,4)];
% 
% NB_illustration_pos_top=NB_illustration_pos(1,2)+NB_illustration_pos(1,4);
% ChAT_illustration_pos_top=ChAT_illustration_pos(1,2)+ChAT_illustration_pos(1,4);
% 
% %%
% %Placing plots in the figure:
% %Cell 44 - spont
% 
% spont_trace_f46_ax_copy = copyobj(spont_trace_f46_ax,F); % copy axes to new fig
% set(spont_trace_f46_ax_copy,'position',spont_trace_f46_pos(1,:))
%    %    'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:));
% spont_mean_f46_ax_copy = copyobj(spont_mean_f46_ax,F); % copy axes to new fig
% set(spont_mean_f46_ax_copy(2),'position',spont_mean_f46_pos(1,:))
% 
% spont_std_f46_ax_copy = copyobj(spont_std_f46_ax,F); % copy axes to new fig
% set(spont_std_f46_ax_copy,'position',spont_std_f46_pos(1,:))
% 
% spont_trace_f80_ax_copy = copyobj(spont_trace_f80_ax,F); % copy axes to new fig
% set(spont_trace_f80_ax_copy,'position',spont_trace_f80_pos(1,:))
% 
% spont_mean_f80_ax_copy = copyobj(spont_mean_f80_ax,F); % copy axes to new fig
% set(spont_mean_f80_ax_copy(2),'position',spont_mean_f80_pos(1,:))
% 
% spont_std_f80_ax_copy = copyobj(spont_std_f80_ax,F); % copy axes to new fig
% set(spont_std_f80_ax_copy,'position',spont_std_f80_pos(1,:))
% 
% 
% %population panels:
% VmM_NB_ax_copy = copyobj(VmM_NB_ax,F); % copy axes to new fig
% set(VmM_NB_ax_copy(1),'position',VmM_NB_pos(1,:))
% set(F, 'currentaxes', VmM_NB_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
% VmM_NB_ax_copy(1).FontSize=ax_fontsize;
% set(gca,'tickdir','out')
% ticklen=fn_get_abs_ticklength(gca, abslen);
% % VmM_NB_ax_copy(1).XTickLabel=[];
% 
% VmSTD_NB_ax_copy = copyobj(VmSTD_NB_ax(1),F); % copy axes to new fig
% set(VmSTD_NB_ax_copy,'position',VmSTD_NB_pos(1,:))
% set(F, 'currentaxes', VmSTD_NB_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
% VmSTD_NB_ax_copy.FontSize=ax_fontsize;
% set(gca,'tickdir','out')
% ticklen=fn_get_abs_ticklength(gca, abslen);
% % VmSTD_NB_ax_copy.XTickLabel=[];
% 
% VmM_ChAT_ax_copy = copyobj(VmM_ChAT_ax,F); % copy axes to new fig
% set(VmM_ChAT_ax_copy(1),'position',VmM_ChAT_pos(1,:))
% set(F, 'currentaxes', VmM_ChAT_ax_copy(1));  tl=title(''); yl=ylabel('mV','fontsize',11); xl=xlabel(''); %tl=title('Mean Vm','fontsize',11,'fontweight','normal');
% VmM_ChAT_ax_copy(1).FontSize=ax_fontsize;
% VmM_ChAT_ax_copy(1).XTickLabel={'Off','On'};
% set(gca,'tickdir','out')
% ticklen=fn_get_abs_ticklength(gca, abslen);
% 
% VmSTD_ChAT_ax_copy = copyobj(VmSTD_ChAT_ax(1),F); % copy axes to new fig
% set(VmSTD_ChAT_ax_copy,'position',VmSTD_ChAT_pos(1,:))
% set(F, 'currentaxes', VmSTD_ChAT_ax_copy);  tl=title('');yl=ylabel('mV','fontsize',11); xl=xlabel('');%tl=title('STD Vm','fontsize',11,'fontweight','normal'); 
% VmSTD_ChAT_ax_copy.FontSize=ax_fontsize;
% VmSTD_ChAT_ax_copy.XTickLabel={'Off','On'};
% set(gca,'tickdir','out')
% ticklen=fn_get_abs_ticklength(gca, abslen);
% 
% NB_illustration_ax = axes('position',NB_illustration_pos);
% imshow(NB_illustration, 'parent', NB_illustration_ax) 
% 
% ChAT_illustration_ax = axes('position',ChAT_illustration_pos);
% imshow(ChAT_illustration, 'parent', ChAT_illustration_ax) 
% 
% 
% %position legend:
% set(spont_mean_f46_ax_copy(1),'position',[0.65 spont_trace_f46_pos_top+0.02 0.08 0.005])
% spont_mean_f46_ax_copy(1).FontSize=9;
% 
% set(spont_mean_f80_ax_copy(1),'position',[0.65 spont_trace_f80_pos_top 0.08 0.005])
% spont_mean_f80_ax_copy(1).FontSize=9;
% %   annotation:
% 
% annotation('textbox', [spont_trace_f46_pos(1,1), spont_trace_f46_pos_top 0 0]+[0 0.03 0.5 0.05],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous Activity', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold') %
%  
%  a_pos1=[-0.03 -0.01 0.04 0.04];
%  a_pos4=[0 -0.01 0.16 0.04];
%  a_pos2=[-0.02 -0.01 0.04 0.04];
%  a_pos3=[-0.035 0 0.04 0.04];
%  annotation('textbox', [spont_trace_f46_pos(1,1)+0.14 spont_trace_f46_pos_top+0.05 0.4 0],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'NB electrical stimulation', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
%  annotation('textbox', [spont_trace_f46_pos(1,1) spont_trace_f46_pos_top 0 0]+a_pos4,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Vm Traces', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
% annotation('textbox', [spont_mean_f46_pos(1,1) spont_mean_f46_pos_top 0 0]+a_pos4,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 12)
% annotation('textbox', [spont_std_f46_pos(1,1) spont_std_f46_pos_top 0 0]+a_pos4,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 12)
% annotation('textbox', [spont_trace_f80_pos(1,1)+0.15 spont_trace_f80_pos_top+0.05 0.4 0],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Optogenetic stimulation', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
% annotation('textbox', [spont_trace_f80_pos(1,1) spont_trace_f80_pos_top-0.018 0 0]+a_pos4,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Vm Traces', 'FontName','arial', 'fontsize', 12) %'fontweight', 'bold'
% annotation('textbox', [spont_mean_f80_pos(1,1) spont_mean_f80_pos_top-0.01 0 0]+a_pos4,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean', 'FontName','arial', 'fontsize', 12)
% annotation('textbox', [spont_std_f80_pos(1,1) spont_std_f80_pos_top-0.01 0 0]+a_pos4,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD', 'FontName','arial', 'fontsize', 12)
%     
% annotation('textbox', [VmM_NB_pos(1,1) VmM_NB_pos_top+0.005 0.16 0.04],...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 11)
% annotation('textbox', [VmSTD_NB_pos(1,1) VmSTD_NB_pos_top+0.005 0.16 0.04],...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 11)
% annotation('textbox', [VmM_ChAT_pos(1,1) VmM_ChAT_pos_top+0.005 0.16 0.04],...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Mean Vm', 'FontName','arial', 'fontsize', 11)
% annotation('textbox', [VmSTD_ChAT_pos(1,1) VmSTD_ChAT_pos_top+0.005 0.16 0.04],...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'STD Vm', 'FontName','arial', 'fontsize', 11)
% 
%  if no_numbering_flag==0;
%  annotation('textbox', [NB_illustration_pos(1,1)-0.01 NB_illustration_pos_top-0.02 0 0],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold') 
% annotation('textbox', [spont_trace_f46_pos(1,1) spont_trace_f46_pos_top 0 0]+a_pos1,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [spont_mean_f46_pos(1,1) spont_mean_f46_pos_top 0 0]+a_pos1,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [spont_std_f46_pos(1,1) spont_std_f46_pos_top 0 0]+a_pos1,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [VmM_NB_pos(1,1) VmM_NB_pos_top 0 0]+a_pos3,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [VmSTD_NB_pos(1,1) VmSTD_NB_pos_top 0 0]+a_pos3,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [ChAT_illustration_pos(1,1)-0.01 ChAT_illustration_pos_top-0.02 0 0],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [spont_trace_f80_pos(1,1) spont_trace_f80_pos_top 0 0]+a_pos1,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [spont_mean_f80_pos(1,1) spont_mean_f80_pos_top 0 0]+a_pos1,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [spont_std_f80_pos(1,1) spont_std_f80_pos_top 0 0]+a_pos1,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [VmM_ChAT_pos(1,1) VmM_ChAT_pos_top 0 0]+a_pos3,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'K', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
% annotation('textbox', [VmSTD_ChAT_pos(1,1) VmSTD_ChAT_pos_top 0 0]+a_pos3,...
%     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'L', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
%  end
%  
%       cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'
%   
% if save_flag==1
%     filename='Fig 2 intracellular traces spont_f46_f80_v2';
%     saveas(F,filename,'fig'); 
%     print(F,filename,'-dpng','-r600','-opengl') 
%     print(F, '-depsc2', filename);
% end
