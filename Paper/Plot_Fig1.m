%% plotting a figure 

close all
clear all
ax_fontsize=12;
%f23 - xlim=[5.5 8.5], f28 - [4 7], f55 - [4.5 7.5], f60 - [5 8];
xlim1=[4,7]; xlim2=[4.5,7.5];
save_flag=0;
no_numbering_flag=0;

%opening saved figures:
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'

NB_illustration = imread('NBES_Schematic_Illustration_tight','TIF');

NB_histology = imread('20150709_slide1(red)_sec3_whole section_2_v1','TIF');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\LFP depth'

% traces_depth3 = open('20150323-002_LFP_traces_h6_t3  6  7_4600um.fig');    
traces_depth3 = open('2016-10-20-001_LFP_traces_h28_t2  3  4  5  6_v3.fig');    
traces_depth3_ax = get(gcf, 'children');

% PSD_depth3 = open('20150323-002_LFP_PSD_1-100Hz_h6_4600um.fig');    
PSD_depth3 = open('2016-10-20-001_LFP_PSD_1-100Hz_h28_v3.fig');    
PSD_depth3_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm Histograms and Powerspec\Powerspec LFP\Not Normalized';

%non-normalized power
total_power = open('Total_Power.fig');
total_power_ax = get(gcf, 'children');

delta = open('Lowfreq_Power.fig');
delta_ax = get(gcf, 'children');

beta = open('Beta_Power.fig');
beta_ax = get(gcf, 'children');

gamma = open('Gamma_Power.fig');
gamma_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH'

%evoked exmaple1
evoked_trace_Off = open('file28 whisker evoked activity NB-.fig');    
evoked_trace_Off_ax = get(gcf, 'children');
evoked_trace_On = open('file28 whisker evoked activity NB+.fig');    
evoked_trace_On_ax = get(gcf, 'children');

%evoked exmaple2
evoked_trace_Off2 = open('file55 whisker evoked activity NB-.fig');    
evoked_trace_Off2_ax = get(gcf, 'children');
evoked_trace_On2 = open('file55 whisker evoked activity NB+.fig');    
evoked_trace_On2_ax = get(gcf, 'children');

% Modulation Index
Modulation_index = open('Modulation index.fig');    
Modulation_index_ax = get(gcf, 'children');

%SNR Index
SNR = open('SNR.fig');    
SNR_ax = get(gcf, 'children');

%Response Modulation
Response_modulation = open('Response_modulation.fig');    
Response_modulation_ax = get(gcf, 'children');

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 20]); %[left, bottom, width, height] [1.2 1.2 18 20]
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);


%% Positions:
NB_illustration_pos = [0 , 0.75 , 0.27 , 0.17];
NB_histology_pos = [0.04 , 0.52 , 0.17 , 0.2];

traces_depth3_pos(1,:) = [0.25 , 0.65 , 0.45 , 0.3];

dist1=0.12; dist2=0.06; h_dist1=traces_depth3_pos(1,1)+traces_depth3_pos(1,3)+dist1;

PSD_depth3_pos(1,:) = [h_dist1 , traces_depth3_pos(1,2)+0.13 ,  0.12,  0.15];

total_power_pos(1,:) = [0.35 , 0.5 , 0.1 , 0.1]; total_power_pos(1,2)=traces_depth3_pos(1,2)-total_power_pos(1,4)-0.03;
delta_pos(1,:) = [total_power_pos(1,1)+total_power_pos(1,3)+dist2 , total_power_pos(1,2) , total_power_pos(1,3) , total_power_pos(1,4)];
beta_pos(1,:) = [delta_pos(1,1)+delta_pos(1,3)+dist2 , total_power_pos(1,2) , total_power_pos(1,3) , total_power_pos(1,4)];
gamma_pos(1,:) = [beta_pos(1,1)+beta_pos(1,3)+dist2 , total_power_pos(1,2) , total_power_pos(1,3) , total_power_pos(1,4)];

evoked_trace_Off_pos(1,:) = [0.1 , total_power_pos(1,2)-0.15 , 0.25 , 0.09];
evoked_trace_Off_pos(2,:) = [evoked_trace_Off_pos(1,1) , 0.5 ,evoked_trace_Off_pos(1,3) ,evoked_trace_Off_pos(1,4)]; evoked_trace_Off_pos(2,2)=evoked_trace_Off_pos(1,2)-evoked_trace_Off_pos(2,4)-0.01;
evoked_trace_On_pos(1,:) = [evoked_trace_Off_pos(1,1) , 0.5 , evoked_trace_Off_pos(1,3) , evoked_trace_Off_pos(1,4)];  evoked_trace_On_pos(1,2)=evoked_trace_Off_pos(2,2)-evoked_trace_Off_pos(1,4)-0.01;
evoked_trace_On_pos(2,:) = [evoked_trace_Off_pos(1,1) , 0.5 , evoked_trace_Off_pos(1,3) , evoked_trace_Off_pos(1,4)]; evoked_trace_On_pos(2,2)=evoked_trace_On_pos(1,2)-evoked_trace_On_pos(2,4)-0.01;

evoked_trace_Off2_pos(1,:) = [evoked_trace_Off_pos(1,1)+evoked_trace_Off_pos(1,3)+0.1 , evoked_trace_Off_pos(1,2) , evoked_trace_Off_pos(1,3) , evoked_trace_Off_pos(1,4)];
evoked_trace_Off2_pos(2,:) = [evoked_trace_Off2_pos(1,1) , 0.5 ,evoked_trace_Off2_pos(1,3) ,evoked_trace_Off2_pos(1,4)]; evoked_trace_Off2_pos(2,2)=evoked_trace_Off2_pos(1,2)-evoked_trace_Off2_pos(2,4)-0.01;
evoked_trace_On2_pos(1,:) = [evoked_trace_Off2_pos(1,1) , 0.5 , evoked_trace_Off2_pos(1,3) , evoked_trace_Off2_pos(1,4)];  evoked_trace_On2_pos(1,2)=evoked_trace_Off2_pos(2,2)-evoked_trace_Off2_pos(1,4)-0.01;
evoked_trace_On2_pos(2,:) = [evoked_trace_Off2_pos(1,1) , 0.5 , evoked_trace_Off2_pos(1,3) , evoked_trace_Off2_pos(1,4)]; evoked_trace_On2_pos(2,2)=evoked_trace_On2_pos(1,2)-evoked_trace_On2_pos(2,4)-0.01;

vert_dim=0.1;
SNR_pos(1,:) = [h_dist1 , evoked_trace_On_pos(2,2), 0.13 , vert_dim];
Modulation_index_pos(1,:) = [h_dist1 , SNR_pos(1,2)+vert_dim+0.05 , 0.13 , vert_dim];
Response_modulation_pos(1,:) = [h_dist1 , Modulation_index_pos(1,2)+vert_dim+0.05 , 0.13 , vert_dim];

%top positions
PSD_depth3_pos_top = PSD_depth3_pos(1,2)+PSD_depth3_pos(1,4);
traces_depth3_pos_top = traces_depth3_pos(1,2)+traces_depth3_pos(1,4);
total_power_pos_top = total_power_pos(1,2)+total_power_pos(1,4);
delta_pos_top = delta_pos(1,2)+delta_pos(1,4);
beta_pos_top = beta_pos(1,2)+beta_pos(1,4);
gamma_pos_top = gamma_pos(1,2)+gamma_pos(1,4);
evoked_trace_Off_pos_top = evoked_trace_Off_pos(1,2)+evoked_trace_Off_pos(1,4);
evoked_trace_On_pos_top = evoked_trace_On_pos(1,2)+evoked_trace_On_pos(1,4);
Modulation_index_pos_top = Modulation_index_pos(1,2)+Modulation_index_pos(1,4);
Response_modulation_pos_top = Response_modulation_pos(1,2)+Response_modulation_pos(1,4);
SNR_pos_top = SNR_pos(1,2)+SNR_pos(1,4);
%%
%Placing plots in the figure:

NB_illustration_ax = axes('position',NB_illustration_pos);
imshow(NB_illustration, 'parent', NB_illustration_ax) 

NB_histology_ax = axes('position',NB_histology_pos);
imshow(NB_histology, 'parent', NB_histology_ax) 

% LFP traces
traces_depth3_ax_copy = copyobj(traces_depth3_ax,F); % copy axes to new fig
set(traces_depth3_ax_copy,'position',traces_depth3_pos(1,:))
     
% LFP PSD
PSD_depth3_ax_copy = copyobj(PSD_depth3_ax,F); % copy axes to new fig
set(PSD_depth3_ax_copy(2),'position',PSD_depth3_pos(1,:))
set(F, 'currentaxes', PSD_depth3_ax_copy(2)); 
PSD_depth3_ax_copy(2).FontSize=ax_fontsize;
% PSD_depth3_ax_copy(2).XTickLabel=[];
%position legend:
set(PSD_depth3_ax_copy(1),'position',[0.88 PSD_depth3_pos_top-0.025 0.08 0.05])
PSD_depth3_ax_copy(1).FontSize=ax_fontsize;
PSD_depth3_ax_copy(1).LineWidth=1.5;  

total_power_ax_copy = copyobj(total_power_ax,F); % copy axes to new fig
set(total_power_ax_copy,'position',total_power_pos(1,:))
total_power_ax_copy.FontSize=ax_fontsize;
set(F, 'currentaxes', total_power_ax_copy); 
tl=title({'Total Power'; ''},'fontweight','normal','fontsize',12);

delta_ax_copy = copyobj(delta_ax,F); % copy axes to new fig
set(delta_ax_copy,'position',delta_pos(1,:))
delta_ax_copy.FontSize=ax_fontsize;
set(F, 'currentaxes', delta_ax_copy); 
%  tl=title({'Delta Power'; ''},'fontweight','normal','fontsize',12);
 tl=title({'1-10 Hz'; ''},'fontweight','normal','fontsize',12); yl=ylabel('');
% yl=ylabel('Normalized PSD');

beta_ax_copy = copyobj(beta_ax,F); % copy axes to new fig
set(beta_ax_copy,'position',beta_pos(1,:))
beta_ax_copy.FontSize=ax_fontsize;
set(F, 'currentaxes', beta_ax_copy); 
% yl=ylabel('');  tl=title({'Beta Power'; ''},'fontweight','normal','fontsize',12);
yl=ylabel('');  tl=title({'12-25 Hz'; ''},'fontweight','normal','fontsize',12);

gamma_ax_copy = copyobj(gamma_ax,F); % copy axes to new fig
set(gamma_ax_copy,'position',gamma_pos(1,:))
gamma_ax_copy.FontSize=ax_fontsize;
set(F, 'currentaxes', gamma_ax_copy); 
% yl=ylabel('');  tl=title({'Gamma Power'; ''},'fontweight','normal','fontsize',12);
yl=ylabel('');  tl=title({'30-50 Hz'; ''},'fontweight','normal','fontsize',12);
%   annotation:
annotation('textbox', [traces_depth3_pos(1,1) traces_depth3_pos_top 0 0]+[0.12 0.01 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Barrel cortex L4 LFP traces', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')

 %PSTH
  evoked_trace_Off_ax_copy = copyobj(evoked_trace_Off_ax,F); % copy axes to new fig
set(evoked_trace_Off_ax_copy(2),'position',evoked_trace_Off_pos(1,:),'xticklabel',[], 'xlim',xlim1,...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5); % 'ylim', y1limits(1,:),'ytick', y1ticks(1,:)
   set(F, 'currentaxes', evoked_trace_Off_ax_copy(2)); t=title(''); yl=ylabel('Trial#'); xl=xlabel('');  
   set(evoked_trace_Off_ax_copy(1),'position',evoked_trace_Off_pos(2,:),'xticklabel',[],'xlim',xlim1,...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5 ); %'xlim',x1limits,'xtick',x1ticks, 'xticklabel',[], 'ylim', y2limits(2,:),'ytick', y2ticks(2,:)
  set(F, 'currentaxes', evoked_trace_Off_ax_copy(1)); t=title(''); yl=ylabel('Spikes/S'); xl=xlabel('');  
  
evoked_trace_On_ax_copy = copyobj(evoked_trace_On_ax,F); % copy axes to new fig
set(evoked_trace_On_ax_copy(2),'position',evoked_trace_On_pos(1,:),'xlim',xlim1,...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5); %,1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:)
   set(F, 'currentaxes', evoked_trace_On_ax_copy(2)); t=title(''); yl=ylabel('Trial#'); xl=xlabel('');  
    set(evoked_trace_On_ax_copy(1),'position',evoked_trace_On_pos(2,:),'xlim',xlim1,...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5); %,'xlim',x1limits,'xtick',x1ticks, 'xticklabel',x1ticklab, 'ylim', y2limits(2,:),'ytick', y2ticks(2,:)
 set(F, 'currentaxes', evoked_trace_On_ax_copy(1)); t=title(''); yl=ylabel('Spikes/S');  xl=xlabel('Time from ES offset [S]');  
 
 evoked_trace_Off2_ax_copy = copyobj(evoked_trace_Off2_ax,F); % copy axes to new fig
set(evoked_trace_Off2_ax_copy(2),'position',evoked_trace_Off2_pos(1,:),'xticklabel',[], 'xlim',xlim2,...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5); % 'ylim', y1limits(1,:),'ytick', y1ticks(1,:)
   set(F, 'currentaxes', evoked_trace_Off2_ax_copy(2)); t=title(''); yl=ylabel('Trial#'); xl=xlabel('');  
   set(evoked_trace_Off2_ax_copy(1),'position',evoked_trace_Off2_pos(2,:),'xticklabel',[],'xlim',xlim2,...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5 ); %'xlim',x1limits,'xtick',x1ticks, 'xticklabel',[], 'ylim', y2limits(2,:),'ytick', y2ticks(2,:)
  set(F, 'currentaxes', evoked_trace_Off2_ax_copy(1)); t=title(''); yl=ylabel('Spikes/S'); xl=xlabel('');  
  
evoked_trace_On2_ax_copy = copyobj(evoked_trace_On2_ax,F); % copy axes to new fig
set(evoked_trace_On2_ax_copy(2),'position',evoked_trace_On2_pos(1,:),'xlim',xlim2,...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5); %,1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:)
   set(F, 'currentaxes', evoked_trace_On2_ax_copy(2)); t=title(''); yl=ylabel('Trial#'); xl=xlabel('');  
    set(evoked_trace_On2_ax_copy(1),'position',evoked_trace_On2_pos(2,:),'xlim',xlim2,...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5); %,'xlim',x1limits,'xtick',x1ticks, 'xticklabel',x1ticklab, 'ylim', y2limits(2,:),'ytick', y2ticks(2,:)
 set(F, 'currentaxes', evoked_trace_On2_ax_copy(1)); t=title(''); yl=ylabel('Spikes/S');  xl=xlabel('Time from ES offset [S]');  

 %Response modulation
Response_modulation_ax_copy = copyobj(Response_modulation_ax,F); % copy axes to new fig
set(Response_modulation_ax_copy,'position',Response_modulation_pos(1,:),'ylim',[-5 30] ,'ytick', [0 10 20],'xticklabel',[],...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5,'box','off'); %, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks, 'xticklabel',x6ticklab);
 set(F, 'currentaxes', Response_modulation_ax_copy); t=title(''); yl=ylabel('Spikes/Train','fontsize',13); 

 %Modulation index
 Modulation_index_ax_copy = copyobj(Modulation_index_ax,F); % copy axes to new fig
set(Modulation_index_ax_copy,'position',Modulation_index_pos(1,:),'ylim',[-1 1] ,'ytick', [-1 0 1],...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5,'box','off'); %, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks, 'xticklabel',x6ticklab);
 set(F, 'currentaxes', Modulation_index_ax_copy); t=title(''); yl=ylabel('MI','fontsize',13);
 %SNR
SNR_ax_copy = copyobj(SNR_ax,F); % copy axes to new fig
set(SNR_ax_copy,'position',SNR_pos(1,:),'ylim',[-0.1 1.1] ,'ytick', [0 0.5 1],...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5,'box','off'); %, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks, 'xticklabel',x6ticklab);
 set(F, 'currentaxes', SNR_ax_copy); t=title(''); yl=ylabel('SNRI','fontsize',13); %set Ytick to be on the left side and not on the right

 %  
 a_pos1=[0 -0.02 0.04 0.04];
 a_pos2=[-0.06 0 0.04 0.04];
 a_pos3=[0.4*traces_depth3_pos(1,3), -0.025, 0.5, 0.05];
 a_pos4=[-0.03, 0.02, 0.04, 0.04];
 

%   annotation('textbox', [traces_depth3_pos(1,1) traces_depth3_pos_top 0 0]+a_pos3,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Stimulation depth 5000 \mum', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
%  
if no_numbering_flag==0

 annotation('textbox', [traces_depth3_pos(1,1) traces_depth3_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [PSD_depth3_pos(1,1) PSD_depth3_pos_top 0 0]+a_pos2,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [total_power_pos(1,1) total_power_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold') 
annotation('textbox', [delta_pos(1,1) total_power_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')  
 annotation('textbox', [beta_pos(1,1) total_power_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold') 
  annotation('textbox', [gamma_pos(1,1) total_power_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold') 
end
 
%% 
if no_numbering_flag==1
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\No Numbering'
else 
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'
end
if save_flag==1;
filename='Fig 1 LFP traces+PSD_1-100Hz+SNR_file28';
saveas(F,filename,'fig'); 
print(F,filename,'-dpng','-r600','-opengl') 
print(F, '-depsc2', filename);
end
