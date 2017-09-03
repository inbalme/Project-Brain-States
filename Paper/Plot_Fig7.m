%% plotting a figure 

close all
clear all
save_flag=0;
no_numbering_flag=0;
abslen=0.05; %[in cm]
ax_fontsize=10;
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP filtered 49-51Hz Presentation'
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP_50Hz+BP0.1-200 Vm_50Hz';
cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm-LFP correlations\LFP_50Hz+BP1-200 Vm_50Hz';
%opening saved figures:
% [0 0 0; [30,75,14]/256; [136 137 138]/256; [112,172,90]/256; [216 22 22]/256; [255 153 153]/256];
color1=[30,75,14]/256;
color2=[112,172,90]/256;
color3=[160,160,160]./256; %grey; 
color4=[102, 172,255]./256; %light blue

% color3=[255,153,16]./256; %orange
% color4=[255,178,102]./256; %light orange
% color3=[102,51,0]./256; %brown; 204 102
% color4=[153,0,0]./256; %brown-red
% color3=[153, 76,0]./256; %brown
% color4=[204, 102,0]./256; %light brown

%spontaneous:
spont_trace_Off_f80 = open('Vm-LFP_spont_stim1_Off_f80_t2  4  5.fig');    
spont_trace_Off_f80_ax = get(gcf, 'children');
        c1=findall(gcf,'color',color1);
        set(c1,'color',color3);
spont_trace_On_f80 = open('Vm-LFP_spont_stim1_On_f80_t2  4  5.fig');    
spont_trace_On_f80_ax = get(gcf, 'children');
        c2=findall(gcf,'color',color2);
        set(c2,'color',color4);

%evoked
evoked_trace_Off_f80 = open('Vm-LFP_evoked_stim1_Off_f80_t2  4  5.fig');    
evoked_trace_Off_f80_ax = get(gcf, 'children');
        c1=findall(gcf,'color',color1);
        set(c1,'color',color3);

evoked_trace_On_f80 = open('Vm-LFP_evoked_stim1_On_f80_t2  4  5.fig');    
evoked_trace_On_f80_ax = get(gcf, 'children');
        c2=findall(gcf,'color',color2);
        set(c2,'color',color4);

%%
% cross-correlations
%spontaneous:
spont_cc_Off_f80 = open('Vm-LFPcc_spont_f80actual_data.fig');    
spont_cc_Off_f80_ax = get(gcf, 'children');

spont_cc_On_f80 = open('Vm-LFPcc_spont_f80shuffled_data.fig');    
spont_cc_On_f80_ax = get(gcf, 'children');

%evoked
evoked_cc_Off_f80 = open('Vm-LFPcc_evoked_f80actual_data.fig');    
evoked_cc_Off_f80_ax = get(gcf, 'children');

evoked_cc_On_f80 = open('Vm-LFPcc_evoked_f80shuffled_data.fig');    
evoked_cc_On_f80_ax = get(gcf, 'children');

%population parameters:
%spontaneous shuff+real
% cc_paired_plot = open('Vm-LFPcc_spont+evoked_max-peak_paired_population.fig');    
% cc_paired_plot = open('Vm-LFPcc_spont+evoked_lag0_paired_population.fig');    
cc_paired_plot = open('Vm-LFPcc_evoked_lag0_paired_population_real-shuff.fig');
l1=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l1(7),'color',color4); uistack(l1(7),'top');
        set(l1(16),'color',color4); uistack(l1(16),'top');
        
cc_paired_plot_ax = get(gcf, 'children');
%evoked shuff+real
% cc_paired_plot_shuff = open('Vm-LFPcc_spont+evoked_max-peak_paired_population_shuff.fig');    
% cc_paired_plot_shuff = open('Vm-LFPcc_spont+evoked_lag0_paired_population_shuff.fig');    
cc_paired_plot_shuff = open('Vm-LFPcc_spont_lag0_paired_population_real-shuff.fig'); 
l2=findall(gcf,'color',[0.7 0.7 0.7]);
        set(l2(7),'color',color4); uistack(l2(7),'top')
        set(l2(16),'color',color4); uistack(l2(16),'top')
cc_paired_plot_shuff_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
ChAT_illustration = imread('ChAT_Schematic_Illustration_two-electrodes_tight','TIF');
%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 20]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
 
%% Positions:
spont_trace_Off_f80_pos(1,:) = [0.05 , 0.87 , 0.43 , 0.06];
for i=1:(length(spont_trace_Off_f80_ax)-1)/2
spont_trace_Off_f80_pos(2*i,:) = spont_trace_Off_f80_pos(2*i-1,:)-[0 , spont_trace_Off_f80_pos(2*i-1,4)-0.02 , 0 , 0];
spont_trace_Off_f80_pos(2*i+1,:) = spont_trace_Off_f80_pos(2*i,:)-[0 , spont_trace_Off_f80_pos(2*i,4) , 0 , 0];
end

spont_trace_On_f80_pos(1,:) = spont_trace_Off_f80_pos(1,:); spont_trace_On_f80_pos(1,2)=spont_trace_Off_f80_pos(6,2)-spont_trace_Off_f80_pos(1,4)-0.03;
for i=1:(length(spont_trace_On_f80_ax)-1)/2
spont_trace_On_f80_pos(2*i,:) = spont_trace_On_f80_pos(2*i-1,:)-[0 , spont_trace_On_f80_pos(2*i-1,4)-0.06 , 0 , 0];
spont_trace_On_f80_pos(2*i+1,:) = spont_trace_On_f80_pos(2*i,:)-[0 , spont_trace_On_f80_pos(2*i,4)+0.015 , 0 , 0];
end


h_dist1=spont_trace_Off_f80_pos(1,1)+spont_trace_Off_f80_pos(1,3)+0.05;
evoked_trace_Off_f80_pos(1,:) = [h_dist1 , spont_trace_Off_f80_pos(1,2) ,  spont_trace_Off_f80_pos(1,3) ,  spont_trace_Off_f80_pos(1,4)];
for i=1:(length(evoked_trace_Off_f80_ax)-1)/2
evoked_trace_Off_f80_pos(2*i,:) = evoked_trace_Off_f80_pos(2*i-1,:)-[0 , evoked_trace_Off_f80_pos(2*i-1,4)-0.02 , 0 , 0];
evoked_trace_Off_f80_pos(2*i+1,:) = evoked_trace_Off_f80_pos(2*i,:)-[0 , evoked_trace_Off_f80_pos(2*i,4) , 0 , 0];
end

evoked_trace_On_f80_pos(1,:) = [h_dist1, spont_trace_On_f80_pos(1,2) ,  evoked_trace_Off_f80_pos(1,3) ,  spont_trace_Off_f80_pos(1,4)];
for i=1:(length(evoked_trace_On_f80_ax)-1)/2
evoked_trace_On_f80_pos(2*i,:) = evoked_trace_On_f80_pos(2*i-1,:)-[0 , evoked_trace_On_f80_pos(2*i-1,4)-0.06 , 0 , 0];
evoked_trace_On_f80_pos(2*i+1,:) = evoked_trace_On_f80_pos(2*i,:)-[0 , evoked_trace_On_f80_pos(2*i,4)+0.015 , 0 , 0];
end

spont_cc_On_f80_pos(1,:) = [0.105 , spont_trace_On_f80_pos(6,2)-0.14 , 0.13 , 0.08]; %shuffled
h_dist2=spont_cc_On_f80_pos(1,1)+spont_cc_On_f80_pos(1,3)+0.06;
spont_cc_Off_f80_pos(1,:) = [h_dist2 ,spont_cc_On_f80_pos(1,2) , spont_cc_On_f80_pos(1,3) , spont_cc_On_f80_pos(1,4)]; 
h_dist3=spont_cc_Off_f80_pos(1,1)+spont_cc_Off_f80_pos(1,3)+0.165;
evoked_cc_On_f80_pos(1,:) = [h_dist3 , spont_cc_On_f80_pos(1,2) ,  spont_cc_On_f80_pos(1,3) ,  spont_cc_On_f80_pos(1,4)]; %shuffled
h_dist4=evoked_cc_On_f80_pos(1,1)+evoked_cc_On_f80_pos(1,3)+0.06;
evoked_cc_Off_f80_pos(1,:) = [h_dist4 , spont_cc_Off_f80_pos(1,2) , evoked_cc_On_f80_pos(1,3) ,  spont_cc_Off_f80_pos(1,4)];

cc_paired_plot_shuff_pos(1,:) = [0.12 , evoked_cc_Off_f80_pos(1,2)-0.22 , 0.27 , 0.12]; 
cc_paired_plot_pos(1,:) =  [cc_paired_plot_shuff_pos(1,1)+cc_paired_plot_shuff_pos(1,3)+0.2 , cc_paired_plot_shuff_pos(1,2) , cc_paired_plot_shuff_pos(1,3) , cc_paired_plot_shuff_pos(1,4)];

%top positions
spont_trace_Off_f80_pos_top = spont_trace_Off_f80_pos(1,2)+spont_trace_Off_f80_pos(1,4);
spont_trace_On_f80_pos_top = spont_trace_On_f80_pos(1,2)+spont_trace_On_f80_pos(1,4);
evoked_trace_Off_f80_pos_top = evoked_trace_Off_f80_pos(1,2)+evoked_trace_Off_f80_pos(1,4);
evoked_trace_On_f80_pos_top = evoked_trace_On_f80_pos(1,2)+evoked_trace_On_f80_pos(1,4);

spont_cc_Off_f80_pos_top = spont_cc_Off_f80_pos(1,2)+spont_cc_Off_f80_pos(1,4);
spont_cc_On_f80_pos_top = spont_cc_On_f80_pos(1,2)+spont_cc_On_f80_pos(1,4);
evoked_cc_Off_f80_pos_top = evoked_cc_Off_f80_pos(1,2)+evoked_cc_Off_f80_pos(1,4);
evoked_cc_On_f80_pos_top = evoked_cc_On_f80_pos(1,2)+evoked_cc_On_f80_pos(1,4);
cc_paired_plot_shuff_pos_top = cc_paired_plot_shuff_pos(1,2)+cc_paired_plot_shuff_pos(1,4);
cc_paired_plot_pos_top = cc_paired_plot_pos(1,2)+cc_paired_plot_pos(1,4);

ChAT_illustration_pos(1,:)=[spont_trace_Off_f80_pos(1,1)+spont_trace_Off_f80_pos(1,3)-0.04,spont_trace_Off_f80_pos_top-0.01, 0.08, 0.08];
ChAT_illustration_pos_top=ChAT_illustration_pos(1,2)+ChAT_illustration_pos(1,4);
%%
%Placing plots in the figure:
ChAT_illustration_ax = axes('position',ChAT_illustration_pos);
imshow(ChAT_illustration, 'parent', ChAT_illustration_ax) 

%Cell 46 - spont
for i=1:length(spont_trace_Off_f80_ax)-1
    spont_trace_Off_f80_ax_copy(i) = copyobj(spont_trace_Off_f80_ax(i+1),F); % copy axes to new fig
    set(spont_trace_Off_f80_ax_copy(i),'position',spont_trace_Off_f80_pos(i,:))
end
   %    'fontname', 'arial','fontsize',13,'linewidth',1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:));
for i=1:length(spont_trace_On_f80_ax)-1
    spont_trace_On_f80_ax_copy(i) = copyobj(spont_trace_On_f80_ax(i+1),F); % copy axes to new fig
    set(spont_trace_On_f80_ax_copy(i),'position',spont_trace_On_f80_pos(i,:))
end
      
%    l=legend('Vm','LFP'); set(l,'box','off'), set(l,'FontSize', 8) ; set(l,'linewidth',1.5);set(l,'position',[0.73 0.89 0.03 0.03]);

   %Cell 46 evoked
       evoked_trace_Off_f80_ax_copy = copyobj(evoked_trace_Off_f80_ax,F); % copy axes to new fig
for i=2:length(evoked_trace_Off_f80_ax)
    set(evoked_trace_Off_f80_ax_copy(i),'position',evoked_trace_Off_f80_pos(i-1,:))
    evoked_trace_Off_f80_ax_copy(i).LineWidth=1.5;
end
   
%position legend:
a=evoked_trace_Off_f80_ax_copy(1).Position;
set(evoked_trace_Off_f80_ax_copy(1),'position',[0.91 evoked_trace_Off_f80_pos_top 0.08 0.04])
evoked_trace_Off_f80_ax_copy(1).FontSize=9;
%  for ix=1:length(evoked_trace_Off_f80_ax_copy(1).String)
%       str = evoked_trace_Off_f80_ax_copy(1).String{ix};
%       h = findobj(gcf,'DisplayName',str);
%       h(1).LineWidth =1.5;
%     end
%  
 evoked_trace_On_f80_ax_copy = copyobj(evoked_trace_On_f80_ax,F); % copy axes to new fig
for i=2:length(evoked_trace_On_f80_ax)
    set(evoked_trace_On_f80_ax_copy(i),'position',evoked_trace_On_f80_pos(i-1,:))
end
%position legend
set(evoked_trace_On_f80_ax_copy(1),'position',[0.91 evoked_trace_On_f80_pos_top 0.08 0.04]) 
evoked_trace_On_f80_ax_copy(1).FontSize=9;
%  for ix=1:length(evoked_trace_On_f80_ax_copy(1).String)
%       str = evoked_trace_On_f80_ax_copy(1).String{ix};
%       h = findobj(gcf,'DisplayName',str);
%       h(1).LineWidth =1.5;
%     end
 %%
 %Spontaneous Cross-correlations file 46 actual data+shuffled data
spont_cc_Off_f80_ax_copy(1) = copyobj(spont_cc_Off_f80_ax(2),F); % copy axes to new fig
set(spont_cc_Off_f80_ax_copy(1),'position',spont_cc_Off_f80_pos(1,:), 'ylim',[-0.5,0.25])
    set(F, 'currentaxes', spont_cc_Off_f80_ax_copy(1));  yl=ylabel(''); 
spont_cc_Off_f80_ax_copy(1).FontSize=ax_fontsize;    
set(gca,'xtick',[-0.5,0,0.5],'ytick',[-0.5,0,0.2],'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);
 
spont_cc_On_f80_ax_copy(1) = copyobj(spont_cc_On_f80_ax(2),F); % copy axes to new fig
set(spont_cc_On_f80_ax_copy(1),'position',spont_cc_On_f80_pos(1,:), 'ylim',[-0.5,0.25]) %'xticklabel',[]
 set(F, 'currentaxes', spont_cc_On_f80_ax_copy(1));  %xl=xlabel('');  
spont_cc_On_f80_ax_copy(1).FontSize=ax_fontsize;
L1=findobj(gca,'type','line');
set(L1,'linestyle','-');
set(gca,'xtick',[-0.5,0,0.5],'ytick',[-0.5,0,0.2],'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);
%   set(F, 'currentaxes', spont_trace_On_f33_ax_copy); t=title('NB+'); yl=ylabel(''); xl=xlabel('');
%    l=legend('Vm','LFP'); set(l,'box','off'), set(l,'FontSize', 8) ; set(l,'linewidth',1.5);set(l,'position',[0.73 0.89 0.03 0.03]);

%Evoked Cross-correlations file 46 actual data+shuffled data
evoked_cc_Off_f80_ax_copy = copyobj(evoked_cc_Off_f80_ax,F); % copy axes to new fig
set(evoked_cc_Off_f80_ax_copy(2),'position',evoked_cc_Off_f80_pos(1,:), 'ylim',[-0.5,0.25])
set(F, 'currentaxes', evoked_cc_Off_f80_ax_copy(2));  yl=ylabel(''); 
evoked_cc_Off_f80_ax_copy(2).FontSize=ax_fontsize;
set(gca,'xtick',[-0.5,0,0.5],'ytick',[-0.5,0,0.2],'tickdir','out')
ticklen=fn_get_abs_ticklength(gca, abslen);

evoked_cc_On_f80_ax_copy = copyobj(evoked_cc_On_f80_ax,F); % copy axes to new fig
set(evoked_cc_On_f80_ax_copy(2),'position',evoked_cc_On_f80_pos(1,:))
   set(F, 'currentaxes', evoked_cc_On_f80_ax_copy(2)); 
   L2=findobj(gca,'type','line');
   set(L2,'linestyle','-');
   evoked_cc_On_f80_ax_copy(2).FontSize=ax_fontsize;
   set(gca,'xtick',[-0.5,0,0.5],'ytick',[-0.5,0,0.5],'tickdir','out')
  ticklen=fn_get_abs_ticklength(gca, abslen);
   
   %position legend shuffled
%  legend_shuff_pos(1,:)=[spont_cc_Off_f80_pos(1,1)+spont_cc_Off_f80_pos(1,3)+0.03, spont_cc_Off_f80_pos_top-0.01, 0.08, 0.05];
%  legend_shuff_pos(1,:)=[spont_cc_Off_f80_pos(1,1)+spont_cc_Off_f80_pos(1,3)-0.05, spont_cc_Off_f80_pos_top, 0.08, 0.05]; 
  delete(evoked_cc_On_f80_ax_copy(1));
% set(evoked_cc_On_f80_ax_copy(1),'position',legend_shuff_pos(1,:)) 
% evoked_cc_On_f80_ax_copy(1).FontSize=10;
% evoked_cc_On_f80_ax_copy(1).LineWidth=1.5;

%position legend not shuffled
% set(evoked_cc_Off_f80_ax_copy(1),'position',legend_shuff_pos(1,:)+[0.3,0,0,0]) 
legend_pos(1,:)=[spont_cc_Off_f80_pos(1,1)+spont_cc_Off_f80_pos(1,3)+0.005, spont_cc_Off_f80_pos_top+0.01, 0.08, 0.05]; 
set(evoked_cc_Off_f80_ax_copy(1),'position',legend_pos) 
evoked_cc_Off_f80_ax_copy(1).FontSize=9;
evoked_cc_Off_f80_ax_copy(1).LineWidth=1.5;

% evoked_cc_Off_f80_ax_copy(1).String{1}='Light Off shuff. subt.';
% evoked_cc_Off_f80_ax_copy(1).String{2}='Light On shuff. subt.';
    
%cc population paired plots
cc_paired_plot_ax_copy = copyobj(cc_paired_plot_ax,F); % copy axes to new fig
set(cc_paired_plot_ax_copy,'position',cc_paired_plot_pos(1,:))
   set(F, 'currentaxes', cc_paired_plot_ax_copy); yl=ylabel('');
    txt1=findobj(gca,'type','text'); delete(txt1(1:2));
cc_paired_plot_ax_copy.FontSize=ax_fontsize;
cc_paired_plot_ax_copy.XTickLabel={'Off','On','Off','On'};
cc_paired_plot_ax_copy.TickDir='out';
ticklen=fn_get_abs_ticklength(gca, abslen);

cc_paired_plot_shuff_ax_copy = copyobj(cc_paired_plot_shuff_ax,F); % copy axes to new fig
set(cc_paired_plot_shuff_ax_copy,'position',cc_paired_plot_shuff_pos(1,:))
 set(F, 'currentaxes', cc_paired_plot_shuff_ax_copy); 
 txt2=findobj(gca,'type','text'); delete(txt2(1:2));
cc_paired_plot_shuff_ax_copy.FontSize=ax_fontsize;
cc_paired_plot_shuff_ax_copy.XTickLabel={'Off','On','Off','On'};
cc_paired_plot_shuff_ax_copy.TickDir='out';
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

%   annotation:
annotation('textbox', [spont_trace_Off_f80_pos(1,1) spont_trace_Off_f80_pos_top 0 0]+[0.04 0.02 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous Activity', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [evoked_trace_Off_f80_pos(1,1) evoked_trace_Off_f80_pos_top 0 0]+[0.03 0.02 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Sensory-evoked Responses', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [spont_cc_On_f80_pos(1,1)+0.01 spont_cc_On_f80_pos_top+0.03 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Shuffled', 'FontName','arial', 'fontsize', 10)
 annotation('textbox', [spont_cc_Off_f80_pos(1,1)+0.01 spont_cc_Off_f80_pos_top+0.03 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Corrected', 'FontName','arial', 'fontsize', 10)
 annotation('textbox', [evoked_cc_On_f80_pos(1,1)+0.01 evoked_cc_On_f80_pos_top+0.03 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Shuffled', 'FontName','arial', 'fontsize', 10)
 annotation('textbox', [evoked_cc_Off_f80_pos(1,1)+0.01 evoked_cc_Off_f80_pos_top+0.03 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Corrected', 'FontName','arial', 'fontsize', 10)
 %population
 annotation('textbox', [cc_paired_plot_shuff_pos(1,1)+0.005 cc_paired_plot_shuff_pos_top+0.03 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Shuffled', 'FontName','arial', 'fontsize', 10)
 annotation('textbox', [cc_paired_plot_shuff_pos(1,1)+0.155 cc_paired_plot_shuff_pos_top+0.03 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Corrected', 'FontName','arial', 'fontsize', 10)
 annotation('textbox', [cc_paired_plot_pos(1,1)+0.005 cc_paired_plot_pos_top+0.03 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Shuffled', 'FontName','arial', 'fontsize', 10)
 annotation('textbox', [cc_paired_plot_pos(1,1)+0.155 cc_paired_plot_pos_top+0.03 0 0],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Corrected', 'FontName','arial', 'fontsize', 10)
 %  
 a_pos1=[-0.02 -0.02 0.04 0.04];
%  a_pos2=[-0.03 0 0.04 0.04];
 a_pos3=[-0.04 -0.01 0.04 0.04];
 a_pos4=[-0.04 0.01 0.04 0.04];
 if no_numbering_flag==0;
annotation('textbox', [spont_trace_Off_f80_pos(1,1) spont_trace_Off_f80_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 annotation('textbox', [spont_trace_On_f80_pos(1,1) spont_trace_On_f80_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 annotation('textbox', [evoked_trace_Off_f80_pos(1,1) evoked_trace_Off_f80_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 annotation('textbox', [evoked_trace_On_f80_pos(1,1) evoked_trace_On_f80_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 annotation('textbox', [spont_cc_Off_f80_pos(1,1) spont_cc_Off_f80_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 annotation('textbox', [spont_cc_On_f80_pos(1,1) spont_cc_On_f80_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 annotation('textbox', [evoked_cc_Off_f80_pos(1,1) evoked_cc_Off_f80_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 annotation('textbox', [evoked_cc_On_f80_pos(1,1) evoked_cc_On_f80_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
  annotation('textbox', [cc_paired_plot_shuff_pos(1,1) cc_paired_plot_shuff_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
  annotation('textbox', [cc_paired_plot_pos(1,1) cc_paired_plot_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 end
 annotation('textbox', [spont_trace_Off_f80_pos(1,1) spont_trace_Off_f80_pos_top 0 0]+[0 -0.01 0.25 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Light Off', 'FontName','arial', 'fontsize', 12, 'color', [0 0 0])
 annotation('textbox', [spont_trace_On_f80_pos(1,1) spont_trace_On_f80_pos_top 0 0]+[0 -0.01 0.25 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Light On', 'FontName','arial', 'fontsize', 12, 'color', [0 0 0])
%  
%  annotation('textbox', [cc_paired_plot_pos(1,1) cc_paired_plot_pos_top 0 0]+[0.02 0.0 0.2 0.04],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Noise Correlations', 'FontName','arial', 'fontsize', 12,'color', [0 0 153]/256)
%   annotation('textbox', [cc_paired_plot_shuff_pos(1,1) cc_paired_plot_shuff_pos_top 0 0]+[0.02 0.0 0.2 0.04],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Signal correlations', 'FontName','arial', 'fontsize', 12,'color', [0 0 153]/256)
%% 

    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'

if save_flag==1
    filename='Fig 7_f80';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end
