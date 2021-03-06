%% plotting a figure 

close all
clear all
ax_fontsize=12;
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\LFP+Blockers'
no_numbering_flag=1;
%opening saved figures:

%spontaneous:

% traces_21 = open('2015-07-29-002_LFP_traces_h1_t2  3  4  8  9.fig');    
traces_21 = open('2015-07-29-001_LFP_traces_h2_t2  3  4  5  7.fig');    
traces_21_ax = get(gcf, 'children');

traces_22 = open('2015-07-29-002_LFP_traces_h6_t2  3  4  5  7.fig');    
traces_22_ax = get(gcf, 'children');

%evoked

% PSD_21 = open('2015-07-29-002_LFP_PSD_1-100Hz_h1.fig');    %recording from end of experiment
PSD_21 = open('2015-07-29-001_LFP_PSD_1-100Hz_h2.fig');    %recording from beginning of experiment
PSD_21_ax = get(gcf, 'children');

PSD_22 = open('2015-07-29-002_LFP_PSD_1-100Hz_h6.fig');    
PSD_22_ax = get(gcf, 'children');


%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 21 15]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);


%% Positions:
traces_21_pos(1,:) = [0.05 , 0.51 , 0.58 , 0.4];
traces_22_pos(1,:) = [0.05 , 0.1 , traces_21_pos(1,3) , traces_21_pos(1,4)]; traces_22_pos(1,2)=traces_21_pos(1,2)-traces_21_pos(1,4)-0.04;

dist1=0.12;
h_dist1=traces_21_pos(1,1)+traces_21_pos(1,3)+dist1;
PSD_21_pos(1,:) = [h_dist1 , traces_21_pos(1,2)+0.08 ,  1-(2*traces_21_pos(1,1)+traces_21_pos(1,3)+dist1) ,  traces_21_pos(1,4)-0.1];
PSD_22_pos(1,:) = [h_dist1, traces_22_pos(1,2)+0.08 ,  PSD_21_pos(1,3) ,  traces_22_pos(1,4)-0.1];

traces_21_pos_top = traces_21_pos(1,2)+traces_21_pos(1,4);
traces_22_pos_top = traces_22_pos(1,2)+traces_22_pos(1,4);
PSD_21_pos_top = PSD_21_pos(1,2)+PSD_21_pos(1,4);
PSD_22_pos_top = PSD_22_pos(1,2)+PSD_22_pos(1,4);
%%
%Placing plots in the figure:
% LFP traces

traces_21_ax_copy = copyobj(traces_21_ax,F); % copy axes to new fig
set(traces_21_ax_copy,'position',traces_21_pos(1,:))

traces_22_ax_copy = copyobj(traces_22_ax,F); % copy axes to new fig
set(traces_22_ax_copy,'position',traces_22_pos(1,:))

     
% LFP PSD
PSD_21_ax_copy = copyobj(PSD_21_ax,F); % copy axes to new fig
set(PSD_21_ax_copy(2),'position',PSD_21_pos(1,:))
set(F, 'currentaxes', PSD_21_ax_copy(2)); 
xl=xlabel('');  
PSD_21_ax_copy(2).FontSize=ax_fontsize;
PSD_21_ax_copy(2).XTickLabel=[];
%position legend:
set(PSD_21_ax_copy(1),'position',[0.88 PSD_21_pos_top-0.025 0.08 0.05])
PSD_21_ax_copy(1).FontSize=ax_fontsize;
PSD_21_ax_copy(1).LineWidth=1.5;  
%   

PSD_22_ax_copy = copyobj(PSD_22_ax(2),F); % copy axes to new fig
set(PSD_22_ax_copy,'position',PSD_22_pos(1,:))
set(F, 'currentaxes', PSD_22_ax_copy); 
PSD_22_ax_copy.FontSize=ax_fontsize;
PSD_22_ax_copy.YLim=[0 0.05];
PSD_22_ax_copy.YMinorTick='on';
hold on
d=line([1 1],[0.0000001 0.0000001],'color',[1 1 1]);
hold off
delete(d)

%   annotation:
annotation('textbox', [traces_21_pos(1,1) traces_21_pos_top 0 0]+[0.17 0.03 0.8 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Barrel cortex LFP with Cholinergic blockers', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
%
 a_pos1=[0 -0.01 0.04 0.04];
 a_pos2=[-0.05 0 0.04 0.04];
 a_pos3=[0.5*traces_21_pos(1,3)-0.02, -0.01, 0.5, 0.05];
 a_pos4=[0.2*traces_21_pos(1,3)-0.06, -0.01, 0.5, 0.05];
 a_pos5=[0.5*traces_21_pos(1,3)-0.06, -0.01, 0.5, 0.05];


 
 annotation('textbox', [traces_21_pos(1,1) traces_21_pos_top 0 0]+a_pos3,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Control', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [traces_22_pos(1,1) traces_22_pos_top 0 0]+a_pos5,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Atropine+MEC', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 
 if no_numbering_flag==0;
annotation('textbox', [traces_21_pos(1,1) traces_21_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [traces_22_pos(1,1) traces_22_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')

 annotation('textbox', [PSD_21_pos(1,1) PSD_21_pos_top 0 0]+a_pos2,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [PSD_22_pos(1,1) PSD_22_pos_top 0 0]+a_pos2,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 end  
 
%% 
% 
if no_numbering_flag==1
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\No Numbering'
else
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures'
end
filename='Fig 1b LFP traces+PSD_1-100Hz+Blockers_v2';
saveas(F,filename,'fig'); 
print(F,filename,'-dpng','-r600','-opengl') 
print(F, '-depsc2', filename);
