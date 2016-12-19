%% plotting a figure 

close all
clear all
ax_fontsize=12;
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\LFP depth'

%opening saved figures:

%spontaneous:
traces_depth1 = open('20150323-002_LFP_traces_h1_t3  6  7_2500um.fig');    
traces_depth1_ax = get(gcf, 'children');

traces_depth2 = open('20150323-002_LFP_traces_h3_t2  4  8_3500um.fig');    
traces_depth2_ax = get(gcf, 'children');

traces_depth3 = open('20150323-002_LFP_traces_h6_t3  6  7_4600um.fig');    
traces_depth3_ax = get(gcf, 'children');

%evoked
% PSD_depth1 = open('20150323-002_LFP_PSD_h1_t3  6  7_2500um.fig');    
PSD_depth1 = open('20150323-002_LFP_PSD_1-100Hz_h1_2500um.fig');   
PSD_depth1_ax = get(gcf, 'children');

% PSD_depth2 = open('20150323-002_LFP_PSD_h3_t2  4  8_3500um.fig');    
PSD_depth2 = open('20150323-002_LFP_PSD_1-100Hz_h3_3500um.fig');    
PSD_depth2_ax = get(gcf, 'children');


% PSD_depth3 = open('20150323-002_LFP_PSD_h6_t3  6  7_4600um.fig');    
PSD_depth3 = open('20150323-002_LFP_PSD_1-100Hz_h6_4600um.fig');    

PSD_depth3_ax = get(gcf, 'children');

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 21 16]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);


%% Positions:
traces_depth1_pos(1,:) = [0.05 , 0.68 , 0.58 , 0.2];
traces_depth2_pos(1,:) = [0.05 , 0.6 , traces_depth1_pos(1,3) , traces_depth1_pos(1,4)]; traces_depth2_pos(1,2)=traces_depth1_pos(1,2)-traces_depth1_pos(1,4)-0.09;
traces_depth3_pos(1,:) = [0.05 , traces_depth2_pos(1,2)-traces_depth2_pos(1,4)-0.09 , traces_depth1_pos(1,3) , traces_depth1_pos(1,4)];

traces_depth1_pos_top = traces_depth1_pos(1,2)+traces_depth1_pos(1,4);
traces_depth2_pos_top = traces_depth2_pos(1,2)+traces_depth2_pos(1,4);
traces_depth3_pos_top = traces_depth3_pos(1,2)+traces_depth3_pos(1,4);

dist1=0.12;
h_dist1=traces_depth1_pos(1,1)+traces_depth1_pos(1,3)+dist1;
PSD_depth1_pos(1,:) = [h_dist1 , traces_depth1_pos(1,2)+0.04 ,  1-(2*traces_depth1_pos(1,1)+traces_depth1_pos(1,3)+dist1) ,  traces_depth1_pos(1,4)-0.04];
PSD_depth2_pos(1,:) = [h_dist1, traces_depth2_pos(1,2)+0.04 ,  PSD_depth1_pos(1,3) ,  traces_depth2_pos(1,4)-0.04];
PSD_depth3_pos(1,:) = [h_dist1 , traces_depth3_pos(1,2)+0.04 ,  PSD_depth1_pos(1,3) ,  traces_depth3_pos(1,4)-0.04];

PSD_depth1_pos_top = PSD_depth1_pos(1,2)+PSD_depth1_pos(1,4);
PSD_depth2_pos_top = PSD_depth2_pos(1,2)+PSD_depth2_pos(1,4);
PSD_depth3_pos_top = PSD_depth3_pos(1,2)+PSD_depth3_pos(1,4);

%%
%Placing plots in the figure:
% LFP traces
traces_depth1_ax_copy = copyobj(traces_depth1_ax,F); % copy axes to new fig
set(traces_depth1_ax_copy,'position',traces_depth1_pos(1,:))

traces_depth2_ax_copy = copyobj(traces_depth2_ax,F); % copy axes to new fig
set(traces_depth2_ax_copy,'position',traces_depth2_pos(1,:))

traces_depth3_ax_copy = copyobj(traces_depth3_ax,F); % copy axes to new fig
set(traces_depth3_ax_copy,'position',traces_depth3_pos(1,:))
     
% LFP PSD
PSD_depth1_ax_copy = copyobj(PSD_depth1_ax,F); % copy axes to new fig
set(PSD_depth1_ax_copy(2),'position',PSD_depth1_pos(1,:))
set(F, 'currentaxes', PSD_depth1_ax_copy(2)); 
xl=xlabel('');  
PSD_depth1_ax_copy(2).FontSize=ax_fontsize;
PSD_depth1_ax_copy(2).XTickLabel=[];
%position legend:
set(PSD_depth1_ax_copy(1),'position',[0.88 PSD_depth1_pos_top-0.025 0.08 0.05])
PSD_depth1_ax_copy(1).FontSize=ax_fontsize;
PSD_depth1_ax_copy(1).LineWidth=1.5;  
%   
PSD_depth2_ax_copy = copyobj(PSD_depth2_ax(2),F); % copy axes to new fig
set(PSD_depth2_ax_copy,'position',PSD_depth2_pos(1,:))
set(F, 'currentaxes', PSD_depth2_ax_copy); 
xl=xlabel('');  
PSD_depth2_ax_copy.FontSize=ax_fontsize;
PSD_depth2_ax_copy.XTickLabel=[];

PSD_depth3_ax_copy = copyobj(PSD_depth3_ax(2),F); % copy axes to new fig
set(PSD_depth3_ax_copy,'position',PSD_depth3_pos(1,:))
set(F, 'currentaxes', PSD_depth3_ax_copy); 
PSD_depth3_ax_copy.FontSize=ax_fontsize;

%   annotation:
annotation('textbox', [traces_depth1_pos(1,1) traces_depth1_pos_top 0 0]+[0.17 0.07 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Barrel cortex LFP traces', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
% annotation('textbox', [PSD_depth1_pos(1,1) PSD_depth1_pos_top 0 0]+[0.12 0 0.5 0.05],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Sensory-evoked Responses', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
%  
 a_pos1=[0 -0.01 0.04 0.04];
 a_pos2=[-0.05 0 0.04 0.04];
 a_pos3=[0.3*traces_depth1_pos(1,3), 0.01, 0.5, 0.05];

 annotation('textbox', [traces_depth1_pos(1,1) traces_depth1_pos_top 0 0]+a_pos3,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Stimulation depth 2500 \mum', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [traces_depth2_pos(1,1) traces_depth2_pos_top 0 0]+a_pos3,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Stimulation depth 3500 \mum', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [traces_depth3_pos(1,1) traces_depth3_pos_top 0 0]+a_pos3,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Stimulation depth 4600 \mum', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 
annotation('textbox', [traces_depth1_pos(1,1) traces_depth1_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [traces_depth2_pos(1,1) traces_depth2_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [traces_depth3_pos(1,1) traces_depth3_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [PSD_depth1_pos(1,1) PSD_depth1_pos_top 0 0]+a_pos2,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [PSD_depth2_pos(1,1) PSD_depth2_pos_top 0 0]+a_pos2,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
  annotation('textbox', [PSD_depth3_pos(1,1) PSD_depth3_pos_top 0 0]+a_pos2,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 
%% 

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures'
filename='Fig 1 LFP traces+PSD_1-100Hz';
saveas(F,'Fig 1 LFP traces+PSD_1-100Hz.fig'); 
print(F,filename,'-dpng','-r600','-opengl') 
print(F, '-depsc2', filename);
