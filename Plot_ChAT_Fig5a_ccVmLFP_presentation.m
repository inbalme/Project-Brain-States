%% plotting a figure 

close all
clear all
save_flag=1;
 no_numbering_flag=1;
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP filtered 49-51Hz Presentation'
% cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm-LFP correlations\LFP_50Hz+BP0.1-200 Vm_50Hz';
% cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm-LFP correlations\LFP_50Hz+BP1-200 Vm_50Hz_7cells';

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm-LFP correlations\LFP_50Hz+BP1-200 Vm_50Hz';


%opening saved figures:

%spontaneous:
spont_trace_Off_f80 = open('Vm-LFP_spont_stim1_Off_f80_t2  3  4.fig');    
spont_trace_Off_f80_ax = get(gcf, 'children');

spont_trace_On_f80 = open('Vm-LFP_spont_stim1_On_f80_t2  3  4.fig');    
spont_trace_On_f80_ax = get(gcf, 'children');

%evoked
evoked_trace_Off_f80 = open('Vm-LFP_evoked_stim1_Off_f80_t2  4  5.fig');    
evoked_trace_Off_f80_ax = get(gcf, 'children');

evoked_trace_On_f80 = open('Vm-LFP_evoked_stim1_On_f80_t2  4  5.fig');    
evoked_trace_On_f80_ax = get(gcf, 'children');

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
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP filtered 49-51Hz'
%actual_data
% cc_paired_plot = open('Vm-LFPcc_spont+evoked_max-peak_paired_population.fig');    
cc_paired_plot = open('Vm-LFPcc_spont+evoked_lag0_paired_population.fig');    
cc_paired_plot_ax = get(gcf, 'children');
%shuffled data
% cc_paired_plot_shuff = open('Vm-LFPcc_spont+evoked_max-peak_paired_population_shuff.fig');    
cc_paired_plot_shuff = open('Vm-LFPcc_spont+evoked_lag0_paired_population_shuff.fig');    
cc_paired_plot_shuff_ax = get(gcf, 'children');

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 30 20]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
 
%% Positions:
spont_trace_Off_f80_pos(1,:) = [0.05 , 0.87 , 0.43 , 0.06];
for i=1:(length(spont_trace_Off_f80_ax)-1)/2
spont_trace_Off_f80_pos(2*i,:) = spont_trace_Off_f80_pos(2*i-1,:)-[0 , spont_trace_Off_f80_pos(2*i-1,4)-0.02 , 0 , 0];
spont_trace_Off_f80_pos(2*i+1,:) = spont_trace_Off_f80_pos(2*i,:)-[0 , spont_trace_Off_f80_pos(2*i,4)+0.03 , 0 , 0];
end

spont_trace_On_f80_pos(1,:) = spont_trace_Off_f80_pos(1,:); spont_trace_On_f80_pos(1,2)=spont_trace_Off_f80_pos(6,2)-spont_trace_Off_f80_pos(1,4)-0.04;
for i=1:(length(spont_trace_On_f80_ax)-1)/2
spont_trace_On_f80_pos(2*i,:) = spont_trace_On_f80_pos(2*i-1,:)-[0 , spont_trace_On_f80_pos(2*i-1,4)-0.04 , 0 , 0];
spont_trace_On_f80_pos(2*i+1,:) = spont_trace_On_f80_pos(2*i,:)-[0 , spont_trace_On_f80_pos(2*i,4)+0.03 , 0 , 0];
end
spont_trace_Off_f80_pos_top = spont_trace_Off_f80_pos(1,2)+spont_trace_Off_f80_pos(1,4);
spont_trace_On_f80_pos_top = spont_trace_On_f80_pos(1,2)+spont_trace_On_f80_pos(1,4);

h_dist1=spont_trace_Off_f80_pos(1,1)+spont_trace_Off_f80_pos(1,3)+0.05;
evoked_trace_Off_f80_pos(1,:) = [h_dist1 , spont_trace_Off_f80_pos(1,2) ,  spont_trace_Off_f80_pos(1,3) ,  spont_trace_Off_f80_pos(1,4)];
for i=1:(length(evoked_trace_Off_f80_ax)-1)/2
evoked_trace_Off_f80_pos(2*i,:) = evoked_trace_Off_f80_pos(2*i-1,:)-[0 , evoked_trace_Off_f80_pos(2*i-1,4)-0.02 , 0 , 0];
evoked_trace_Off_f80_pos(2*i+1,:) = evoked_trace_Off_f80_pos(2*i,:)-[0 , evoked_trace_Off_f80_pos(2*i,4)+0.03 , 0 , 0];
end
% for i=2:length(evoked_trace_Off_f80_ax)-1
% evoked_trace_Off_f80_pos(i,:) = evoked_trace_Off_f80_pos(i-1,:)-[0 , evoked_trace_Off_f80_pos(i-1,4) , 0 , 0];
% end

evoked_trace_On_f80_pos(1,:) = [h_dist1, spont_trace_On_f80_pos(1,2) ,  evoked_trace_Off_f80_pos(1,3) ,  spont_trace_Off_f80_pos(1,4)];
for i=1:(length(evoked_trace_On_f80_ax)-1)/2
evoked_trace_On_f80_pos(2*i,:) = evoked_trace_On_f80_pos(2*i-1,:)-[0 , evoked_trace_On_f80_pos(2*i-1,4)-0.04 , 0 , 0];
evoked_trace_On_f80_pos(2*i+1,:) = evoked_trace_On_f80_pos(2*i,:)-[0 , evoked_trace_On_f80_pos(2*i,4)+0.03 , 0 , 0];
end
evoked_trace_Off_f80_pos_top = evoked_trace_Off_f80_pos(1,2)+evoked_trace_Off_f80_pos(1,4);
evoked_trace_On_f80_pos_top = evoked_trace_On_f80_pos(1,2)+evoked_trace_On_f80_pos(1,4);

spont_cc_On_f80_pos(1,:) = [0.09 , spont_trace_On_f80_pos(6,2)-0.165 , 0.13 , 0.1]; %shuffled
h_dist2=spont_cc_On_f80_pos(1,1)+spont_cc_On_f80_pos(1,3)+0.05;
spont_cc_Off_f80_pos(1,:) = [h_dist2 ,spont_cc_On_f80_pos(1,2) , spont_cc_On_f80_pos(1,3) , spont_cc_On_f80_pos(1,4)]; 
h_dist3=spont_cc_Off_f80_pos(1,1)+spont_cc_Off_f80_pos(1,3)+0.19;
evoked_cc_On_f80_pos(1,:) = [h_dist3 , spont_cc_On_f80_pos(1,2) ,  spont_cc_On_f80_pos(1,3) ,  spont_cc_On_f80_pos(1,4)]; %shuffled
h_dist4=evoked_cc_On_f80_pos(1,1)+evoked_cc_On_f80_pos(1,3)+0.05;
evoked_cc_Off_f80_pos(1,:) = [h_dist4 , spont_cc_Off_f80_pos(1,2) , evoked_cc_On_f80_pos(1,3) ,  spont_cc_Off_f80_pos(1,4)];

spont_cc_Off_f80_pos_top = spont_cc_Off_f80_pos(1,2)+spont_cc_Off_f80_pos(1,4);
spont_cc_On_f80_pos_top = spont_cc_On_f80_pos(1,2)+spont_cc_On_f80_pos(1,4);
evoked_cc_Off_f80_pos_top = evoked_cc_Off_f80_pos(1,2)+evoked_cc_Off_f80_pos(1,4);
evoked_cc_On_f80_pos_top = evoked_cc_On_f80_pos(1,2)+evoked_cc_On_f80_pos(1,4);

cc_paired_plot_shuff_pos(1,:) = [h_dist3 , evoked_cc_Off_f80_pos(1,2)-0.01 , 0.15 , 0.22]; 
cc_paired_plot_shuff_pos_top = cc_paired_plot_shuff_pos(1,2)+cc_paired_plot_shuff_pos(1,4);

cc_paired_plot_pos(1,:) =  [cc_paired_plot_shuff_pos(1,1)+cc_paired_plot_shuff_pos(1,3)+0.06 , cc_paired_plot_shuff_pos(1,2) , 0.15 , cc_paired_plot_shuff_pos(1,4)];
cc_paired_plot_pos_top = cc_paired_plot_pos(1,2)+cc_paired_plot_pos(1,4);

%%
%Placing plots in the figure:
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
end
   
%position legend:
a=evoked_trace_Off_f80_ax_copy(1).Position;
set(evoked_trace_Off_f80_ax_copy(1),'position',[0.92 evoked_trace_Off_f80_pos_top 0.08 0.05])
evoked_trace_Off_f80_ax_copy(1).FontSize=11;
% for ix=1:length(evoked_trace_Off_f80_ax_copy(1).String)
%       str = evoked_trace_Off_f80_ax_copy(1).String{ix};
%       h = findobj(gcf,'DisplayName',str);
%       h(1).LineWidth =1.5;
% end
%
 evoked_trace_On_f80_ax_copy = copyobj(evoked_trace_On_f80_ax,F); % copy axes to new fig
for i=2:length(evoked_trace_On_f80_ax)
    set(evoked_trace_On_f80_ax_copy(i),'position',evoked_trace_On_f80_pos(i-1,:))
end
%position legend
set(evoked_trace_On_f80_ax_copy(1),'position',[0.92 evoked_trace_On_f80_pos_top 0.08 0.05]) 
evoked_trace_On_f80_ax_copy(1).FontSize=11;
% for ix=1:length(evoked_trace_On_f80_ax_copy(1).String)
%       str = evoked_trace_On_f80_ax_copy(1).String{ix};
%       h = findobj(gcf,'DisplayName',str);
%       h(1).LineWidth =1.5;
% end
%%
 %Spontaneous Cross-correlations file 46 actual data+shuffled data
spont_cc_Off_f80_ax_copy(1) = copyobj(spont_cc_Off_f80_ax(2),F); % copy axes to new fig
set(spont_cc_Off_f80_ax_copy(1),'position',spont_cc_Off_f80_pos(1,:))
    set(F, 'currentaxes', spont_cc_Off_f80_ax_copy(1));  yl=ylabel(''); 
spont_cc_Off_f80_ax_copy(1).FontSize=12;    
 
spont_cc_On_f80_ax_copy(1) = copyobj(spont_cc_On_f80_ax(2),F); % copy axes to new fig
set(spont_cc_On_f80_ax_copy(1),'position',spont_cc_On_f80_pos(1,:)) %,'xticklabel',[]
 set(F, 'currentaxes', spont_cc_On_f80_ax_copy(1));  %xl=xlabel('');  
spont_cc_On_f80_ax_copy(1).FontSize=12;
   
%   set(F, 'currentaxes', spont_trace_On_f33_ax_copy); t=title('NB+'); yl=ylabel(''); xl=xlabel('');
%    l=legend('Vm','LFP'); set(l,'box','off'), set(l,'FontSize', 8) ; set(l,'linewidth',1.5);set(l,'position',[0.73 0.89 0.03 0.03]);

%Evoked Cross-correlations file 46 actual data+shuffled data
evoked_cc_Off_f80_ax_copy = copyobj(evoked_cc_Off_f80_ax,F); % copy axes to new fig
set(evoked_cc_Off_f80_ax_copy(2),'position',evoked_cc_Off_f80_pos(1,:))
set(F, 'currentaxes', evoked_cc_Off_f80_ax_copy(2)); yl=ylabel(''); 
evoked_cc_Off_f80_ax_copy(2).FontSize=12;

evoked_cc_On_f80_ax_copy = copyobj(evoked_cc_On_f80_ax,F); % copy axes to new fig
set(evoked_cc_On_f80_ax_copy(2),'position',evoked_cc_On_f80_pos(1,:))
   set(F, 'currentaxes', evoked_cc_On_f80_ax_copy(2)); yl=ylabel('');  
   evoked_cc_On_f80_ax_copy(2).FontSize=12;
   
   %position legend shuffled
 legend_shuff_pos(1,:)=[spont_cc_Off_f80_pos(1,1)+spont_cc_Off_f80_pos(1,3)+0.03, spont_cc_Off_f80_pos_top-0.01, 0.08, 0.05];
 set(evoked_cc_On_f80_ax_copy(1),'position',legend_shuff_pos(1,:)) 
evoked_cc_On_f80_ax_copy(1).FontSize=11;
evoked_cc_On_f80_ax_copy(1).LineWidth=1.5;

%position legend not shuffled
set(evoked_cc_Off_f80_ax_copy(1),'position',legend_shuff_pos(1,:)-[0,0.07,0,0]) 
evoked_cc_Off_f80_ax_copy(1).FontSize=11;
evoked_cc_Off_f80_ax_copy(1).LineWidth=1.5;

%cc population paired plots
% cc_paired_plot_ax_copy = copyobj(cc_paired_plot_ax,F); % copy axes to new fig
% set(cc_paired_plot_ax_copy,'position',cc_paired_plot_pos(1,:))
% set(F, 'currentaxes', cc_paired_plot_ax_copy); yl=ylabel('');
% cc_paired_plot_ax_copy.FontSize=14;
% cc_paired_plot_ax_copy.XTickLabel={'Off', 'On','Off', 'On'}; 
% 
% cc_paired_plot_shuff_ax_copy = copyobj(cc_paired_plot_shuff_ax,F); % copy axes to new fig
% set(cc_paired_plot_shuff_ax_copy,'position',cc_paired_plot_shuff_pos(1,:))
% cc_paired_plot_shuff_ax_copy.FontSize=14;
% cc_paired_plot_shuff_ax_copy.XTickLabel={'Off', 'On','Off', 'On'}; 
    %   'fontname', 'arial','fontsize',13,'linewidth',1.5,'box','off'); %, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks, 'xticklabel',x6ticklab);
%  set(F, 'currentaxes', Response_modulation_ax_copy); t=title(''); yl=ylabel('Spikes/Stim. train','fontsize',13);

%   annotation:
annotation('textbox', [spont_trace_Off_f80_pos(1,1) spont_trace_Off_f80_pos_top 0 0]+[0.09 0 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous Activity', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
annotation('textbox', [evoked_trace_Off_f80_pos(1,1) evoked_trace_Off_f80_pos_top 0 0]+[0.09 0 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Sensory-evoked Responses', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
%  
a_pos1=[-0.03 -0.02 0.04 0.04];
%  a_pos2=[-0.03 0 0.04 0.04];
 a_pos3=[-0.04 -0.01 0.04 0.04];
 a_pos4=[-0.04 0.01 0.04 0.04];
  if no_numbering_flag==0;
annotation('textbox', [spont_trace_Off_f80_pos(1,1) spont_trace_Off_f80_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_trace_On_f80_pos(1,1) spont_trace_On_f80_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [evoked_trace_Off_f80_pos(1,1) evoked_trace_Off_f80_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [evoked_trace_On_f80_pos(1,1) evoked_trace_On_f80_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_cc_Off_f80_pos(1,1) spont_cc_Off_f80_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_cc_On_f80_pos(1,1) spont_cc_On_f80_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [evoked_cc_Off_f80_pos(1,1) evoked_cc_Off_f80_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [evoked_cc_On_f80_pos(1,1) evoked_cc_On_f80_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
%   annotation('textbox', [cc_paired_plot_shuff_pos(1,1) cc_paired_plot_shuff_pos_top 0 0]+a_pos4,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
%   annotation('textbox', [cc_paired_plot_pos(1,1) cc_paired_plot_pos_top 0 0]+a_pos4,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
  end
 annotation('textbox', [spont_trace_Off_f80_pos(1,1) spont_trace_Off_f80_pos_top 0 0]+[0 -0.01 0.1 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Light Off', 'FontName','arial', 'fontsize', 12, 'color',[0 0 0] ) %[0 0 153]/256
 annotation('textbox', [spont_trace_On_f80_pos(1,1) spont_trace_On_f80_pos_top 0 0]+[0 -0.01 0.1 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Light On', 'FontName','arial', 'fontsize', 12, 'color', [0 0 0])
%  
%  annotation('textbox', [cc_paired_plot_pos(1,1) cc_paired_plot_pos_top 0 0]+[0.02 0.0 0.2 0.04],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Noise Correlations', 'FontName','arial', 'fontsize', 12,'color', [0 0 153]/256)
%   annotation('textbox', [cc_paired_plot_shuff_pos(1,1) cc_paired_plot_shuff_pos_top 0 0]+[0.02 0.0 0.2 0.04],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Signal correlations', 'FontName','arial', 'fontsize', 12,'color', [0 0 153]/256)
%% 
 if no_numbering_flag==1;
     cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\No Numbering'
 else
     cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures'
 end
if save_flag==1
    filename='Fig 5a ChAT ccVmLFP lag0';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end
