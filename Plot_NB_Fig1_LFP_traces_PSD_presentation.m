%% plotting a figure 

close all
clear all
ax_fontsize=12;
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\LFP depth'
save_flag=1;
no_numbering_flag=1;
%opening saved figures:

%spontaneous:
% add header 34
% traces_depth1 = open('20150323-002_LFP_traces_h1_t3  6  7_2500um.fig');    
traces_depth1 = open('2016-10-20-001_LFP_traces_h33_t2  3  4  5  6_v3.fig');     
traces_depth1_ax = get(gcf, 'children');

% traces_depth2 = open('20150323-002_LFP_traces_h3_t2  4  8_3500um.fig');    
% traces_depth2 = open('2016-10-20-001_LFP_traces_h30_t2  3  4  5  6v2.fig');    
traces_depth2 = open('2016-10-20-001_LFP_traces_h31_t2  3  4  5  6_v3.fig');    
traces_depth2_ax = get(gcf, 'children');

% traces_depth3 = open('20150323-002_LFP_traces_h6_t3  6  7_4600um.fig');    
traces_depth3 = open('2016-10-20-001_LFP_traces_h28_t2  3  4  5  6_v3.fig');    
traces_depth3_ax = get(gcf, 'children');

%evoked
% PSD_depth1 = open('20150323-002_LFP_PSD_1-100Hz_h1_2500um.fig');   
PSD_depth1 = open('2016-10-20-001_LFP_PSD_1-100Hz_h33_v3.fig');
PSD_depth1_ax = get(gcf, 'children');

% PSD_depth2 = open('20150323-002_LFP_PSD_1-100Hz_h3_3500um.fig');    
% PSD_depth2 = open('2016-10-20-001_LFP_PSD_1-100Hz_h30.fig');    
PSD_depth2 = open('2016-10-20-001_LFP_PSD_1-100Hz_h31_v3.fig');    
PSD_depth2_ax = get(gcf, 'children');

% PSD_depth3 = open('20150323-002_LFP_PSD_1-100Hz_h6_4600um.fig');    
PSD_depth3 = open('2016-10-20-001_LFP_PSD_1-100Hz_h28_v3.fig');    
PSD_depth3_ax = get(gcf, 'children');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm Histograms and Powerspec\Powerspec LFP\Not Normalized';

%non-normalized power
total_power = open('Total_Power.fig');
total_power_ax = get(gcf, 'children');

delta = open('Lowfreq_Power.fig');
delta_ax = get(gcf, 'children');

% theta = open('Theta_Power.fig');
% theta_ax = get(gcf, 'children');
% 
% alpha = open('Alpha_Power.fig');
% alpha_ax = get(gcf, 'children');
% 
beta = open('Beta_Power.fig');
beta_ax = get(gcf, 'children');

gamma = open('Gamma_Power.fig');
gamma_ax = get(gcf, 'children');

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 27 29]); %[left, bottom, width, height] [1.2 1.2 18 20]
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);


%% Positions:
traces_depth1_pos(1,:) = [0.05 , 0.72 , 0.6 , 0.2];
traces_depth2_pos(1,:) = [0.05 , 0.6 , traces_depth1_pos(1,3) , traces_depth1_pos(1,4)]; traces_depth2_pos(1,2)=traces_depth1_pos(1,2)-traces_depth1_pos(1,4)-0.03;
traces_depth3_pos(1,:) = [0.05 , traces_depth2_pos(1,2)-traces_depth2_pos(1,4)-0.03 , traces_depth1_pos(1,3) , traces_depth1_pos(1,4)];

dist1=0.12;
dist2=0.08;
h_dist1=traces_depth1_pos(1,1)+traces_depth1_pos(1,3)+dist1;
PSD_depth1_pos(1,:) = [h_dist1 , traces_depth1_pos(1,2)+0.04 ,  1-(2*traces_depth1_pos(1,1)+traces_depth1_pos(1,3)+dist1) ,  traces_depth1_pos(1,4)-0.04];
PSD_depth2_pos(1,:) = [h_dist1, traces_depth2_pos(1,2)+0.04 ,  PSD_depth1_pos(1,3) ,  traces_depth2_pos(1,4)-0.04];
PSD_depth3_pos(1,:) = [h_dist1 , traces_depth3_pos(1,2)+0.04 ,  PSD_depth1_pos(1,3) ,  traces_depth3_pos(1,4)-0.04];

total_power_pos(1,:) = [0.14 , 0.5 , 0.12 , 0.14]; total_power_pos(1,2)=traces_depth3_pos(1,2)-total_power_pos(1,4)-0.08;
delta_pos(1,:) = [total_power_pos(1,1)+total_power_pos(1,3)+dist2 , total_power_pos(1,2) , total_power_pos(1,3) , total_power_pos(1,4)];
% theta_pos(1,:) =[delta_pos(1,1)+delta_pos(1,3)+dist2 , total_power_pos(1,2) , total_power_pos(1,3) , total_power_pos(1,4)];
% alpha_pos(1,:) =[theta_pos(1,1)+theta_pos(1,3)+dist2 , total_power_pos(1,2) , total_power_pos(1,3) , total_power_pos(1,4)];
% beta_pos(1,:) = [alpha_pos(1,1)+alpha_pos(1,3)+dist2 , total_power_pos(1,2) , total_power_pos(1,3) , total_power_pos(1,4)];
beta_pos(1,:) = [delta_pos(1,1)+delta_pos(1,3)+dist2 , total_power_pos(1,2) , total_power_pos(1,3) , total_power_pos(1,4)];
gamma_pos(1,:) = [beta_pos(1,1)+beta_pos(1,3)+dist2 , total_power_pos(1,2) , total_power_pos(1,3) , total_power_pos(1,4)];

%top positions
PSD_depth1_pos_top = PSD_depth1_pos(1,2)+PSD_depth1_pos(1,4);
PSD_depth2_pos_top = PSD_depth2_pos(1,2)+PSD_depth2_pos(1,4);
PSD_depth3_pos_top = PSD_depth3_pos(1,2)+PSD_depth3_pos(1,4);
traces_depth1_pos_top = traces_depth1_pos(1,2)+traces_depth1_pos(1,4);
traces_depth2_pos_top = traces_depth2_pos(1,2)+traces_depth2_pos(1,4);
traces_depth3_pos_top = traces_depth3_pos(1,2)+traces_depth3_pos(1,4);
total_power_pos_top = total_power_pos(1,2)+total_power_pos(1,4);
delta_pos_top = delta_pos(1,2)+delta_pos(1,4);
% theta_pos_top = theta_pos(1,2)+theta_pos(1,4);
% alpha_pos_top = alpha_pos(1,2)+alpha_pos(1,4);
beta_pos_top = beta_pos(1,2)+beta_pos(1,4);
gamma_pos_top = gamma_pos(1,2)+gamma_pos(1,4);

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

% theta_ax_copy = copyobj(theta_ax,F); % copy axes to new fig
% set(theta_ax_copy,'position',theta_pos(1,:))
% theta_ax_copy.FontSize=ax_fontsize;
% set(F, 'currentaxes', theta_ax_copy); 
% yl=ylabel(''); tl=title({'Theta Power'; ''},'fontweight','normal','fontsize',12);

% alpha_ax_copy = copyobj(alpha_ax,F); % copy axes to new fig
% set(alpha_ax_copy,'position',alpha_pos(1,:))
% alpha_ax_copy.FontSize=ax_fontsize;
% set(F, 'currentaxes', alpha_ax_copy); 
% yl=ylabel('');  tl=title({'Alpha Power'; ''},'fontweight','normal','fontsize',12);

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
annotation('textbox', [traces_depth1_pos(1,1) traces_depth1_pos_top 0 0]+[0.17 0.01 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Barrel cortex LFP traces', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
%  
 a_pos1=[0 -0.02 0.04 0.04];
 a_pos2=[-0.06 0 0.04 0.04];
 a_pos3=[0.4*traces_depth1_pos(1,3), -0.025, 0.5, 0.05];
 a_pos4=[-0.03, 0.02, 0.04, 0.04];
 

 annotation('textbox', [traces_depth1_pos(1,1) traces_depth1_pos_top 0 0]+a_pos3,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Stimulation depth 3000 \mum', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [traces_depth2_pos(1,1) traces_depth2_pos_top 0 0]+a_pos3,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Stimulation depth 3700 \mum', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [traces_depth3_pos(1,1) traces_depth3_pos_top 0 0]+a_pos3,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Stimulation depth 5000 \mum', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
%  
if no_numbering_flag==0
annotation('textbox', [traces_depth1_pos(1,1) traces_depth1_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [traces_depth2_pos(1,1) traces_depth2_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [traces_depth3_pos(1,1) traces_depth3_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [PSD_depth1_pos(1,1) PSD_depth1_pos_top 0 0]+a_pos2,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [PSD_depth2_pos(1,1) PSD_depth2_pos_top 0 0]+a_pos2,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
  annotation('textbox', [PSD_depth3_pos(1,1) PSD_depth3_pos_top 0 0]+a_pos2,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [total_power_pos(1,1) total_power_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold') 
annotation('textbox', [delta_pos(1,1) total_power_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')  
 annotation('textbox', [beta_pos(1,1) total_power_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold') 
  annotation('textbox', [gamma_pos(1,1) total_power_pos_top 0 0]+a_pos4,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold') 
end
 
%% 
if no_numbering_flag==1
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\No Numbering'
else 
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures'
end
if save_flag==1;
filename='Fig 1 LFP traces+PSD_1-100Hz_v4_not_normalized';
saveas(F,filename,'fig'); 
print(F,filename,'-dpng','-r600','-opengl') 
print(F, '-depsc2', filename);
end
