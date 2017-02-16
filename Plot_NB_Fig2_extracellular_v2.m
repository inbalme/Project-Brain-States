%% plotting a figure 

close all
clear all
save_flag=1;
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH'

%opening saved figures:

%spontaneous:
spont_trace_Off_f23 = open('file23 spontaneous activity NB-.fig');    
spont_trace_Off_f23_ax = get(gcf, 'children');
spont_trace_Off_f23_pos(1,:) = [0.1 , 0.83 , 0.2 , 0.07];
spont_trace_Off_f23_pos(2,:) = [0.1 , 0.5 , 0.2 , 0.07]; spont_trace_Off_f23_pos(2,2)=spont_trace_Off_f23_pos(1,2)-spont_trace_Off_f23_pos(2,4)-0.01;
spont_trace_Off_f23_pos_top = spont_trace_Off_f23_pos(1,2)+spont_trace_Off_f23_pos(1,4);

spont_trace_On_f23 = open('file23 spontaneous activity NB+.fig');    
spont_trace_On_f23_ax = get(gcf, 'children');
spont_trace_On_f23_pos(1,:) = [0.1 , 0.5 , 0.2 , 0.07];  spont_trace_On_f23_pos(1,2)=spont_trace_Off_f23_pos(2,2)-spont_trace_On_f23_pos(1,4)-0.06;
spont_trace_On_f23_pos(2,:) = [0.1 , 0.5 , 0.2 , 0.07]; spont_trace_On_f23_pos(2,2)=spont_trace_On_f23_pos(1,2)-spont_trace_On_f23_pos(2,4)-0.01;
spont_trace_On_f23_pos_top = spont_trace_On_f23_pos(1,2)+spont_trace_On_f23_pos(1,4);

%evoked
evoked_trace_Off_f23 = open('file23 whisker evoked activity NB-.fig');    
evoked_trace_Off_f23_ax = get(gcf, 'children');
h_dist=spont_trace_Off_f23_pos(1,1)+spont_trace_Off_f23_pos(1,3)+0.1;
evoked_trace_Off_f23_pos(1,:) = [h_dist , spont_trace_Off_f23_pos(1,2) ,  spont_trace_Off_f23_pos(1,3)+0.1 ,  spont_trace_Off_f23_pos(1,4)];
evoked_trace_Off_f23_pos(2,:) = [h_dist , spont_trace_Off_f23_pos(2,2) ,  spont_trace_Off_f23_pos(2,3)+0.1 ,  spont_trace_Off_f23_pos(2,4)];
evoked_trace_Off_f23_pos_top = evoked_trace_Off_f23_pos(1,2)+evoked_trace_Off_f23_pos(1,4);

evoked_trace_On_f23 = open('file23 whisker evoked activity NB+.fig');    
evoked_trace_On_f23_ax = get(gcf, 'children');
evoked_trace_On_f23_pos(1,:) = [h_dist , spont_trace_On_f23_pos(1,2) ,  evoked_trace_Off_f23_pos(1,3) ,  spont_trace_On_f23_pos(1,4)];
evoked_trace_On_f23_pos(2,:) = [h_dist , spont_trace_On_f23_pos(2,2) ,  evoked_trace_Off_f23_pos(2,3) ,  spont_trace_On_f23_pos(2,4)];
evoked_trace_On_f23_pos_top = evoked_trace_On_f23_pos(1,2)+evoked_trace_On_f23_pos(1,4);

% cell 33
%spontaneous:
spont_trace_Off_f33 = open('file33 spontaneous activity NB-.fig');    
spont_trace_Off_f33_ax = get(gcf, 'children');
spont_trace_Off_f33_pos(1,:) = [0.1 , spont_trace_On_f23_pos(2,2)-spont_trace_On_f23_pos(2,4)-0.08 , 0.2 , 0.07];
spont_trace_Off_f33_pos(2,:) = [0.1 , 0.5 , 0.2 , 0.07]; spont_trace_Off_f33_pos(2,2)=spont_trace_Off_f33_pos(1,2)-spont_trace_Off_f33_pos(2,4)-0.01;
spont_trace_Off_f33_pos_top = spont_trace_Off_f33_pos(1,2)+spont_trace_Off_f33_pos(1,4);

spont_trace_On_f33 = open('file33 spontaneous activity NB+.fig');    
spont_trace_On_f33_ax = get(gcf, 'children');
spont_trace_On_f33_pos(1,:) = [0.1 , 0.5 , 0.2 , 0.07];  spont_trace_On_f33_pos(1,2)=spont_trace_Off_f33_pos(2,2)-spont_trace_On_f33_pos(1,4)-0.07;
spont_trace_On_f33_pos(2,:) = [0.1 , 0.5 , 0.2 , 0.07]; spont_trace_On_f33_pos(2,2)=spont_trace_On_f33_pos(1,2)-spont_trace_On_f33_pos(2,4)-0.01;
spont_trace_On_f33_pos_top = spont_trace_On_f33_pos(1,2)+spont_trace_On_f33_pos(1,4);

%evoked
evoked_trace_Off_f33 = open('file33 whisker evoked activity NB-.fig');    
evoked_trace_Off_f33_ax = get(gcf, 'children');
h_dist=spont_trace_Off_f33_pos(1,1)+spont_trace_Off_f33_pos(1,3)+0.1;
evoked_trace_Off_f33_pos(1,:) = [h_dist , spont_trace_Off_f33_pos(1,2) ,  spont_trace_Off_f33_pos(1,3)+0.1 ,  spont_trace_Off_f33_pos(1,4)];
evoked_trace_Off_f33_pos(2,:) = [h_dist , spont_trace_Off_f33_pos(2,2) ,  spont_trace_Off_f33_pos(2,3)+0.1 ,  spont_trace_Off_f33_pos(2,4)];
evoked_trace_Off_f33_pos_top = evoked_trace_Off_f33_pos(1,2)+evoked_trace_Off_f33_pos(1,4);

evoked_trace_On_f33 = open('file33 whisker evoked activity NB+.fig');    
evoked_trace_On_f33_ax = get(gcf, 'children');
evoked_trace_On_f33_pos(1,:) = [h_dist , spont_trace_On_f33_pos(1,2) ,  evoked_trace_Off_f33_pos(1,3) ,  spont_trace_On_f33_pos(1,4)];
evoked_trace_On_f33_pos(2,:) = [h_dist , spont_trace_On_f33_pos(2,2) ,  evoked_trace_Off_f33_pos(2,3) ,  spont_trace_On_f33_pos(2,4)];
evoked_trace_On_f33_pos_top = evoked_trace_On_f33_pos(1,2)+evoked_trace_On_f33_pos(1,4);

%population parameters:
h_dist=evoked_trace_Off_f23_pos(1,1)+evoked_trace_Off_f23_pos(1,3)+ 0.14;
vert_dim=0.18;
%SNR Index
SNR = open('SNR.fig');    
SNR_ax = get(gcf, 'children');
SNR_pos(1,:) = [h_dist , evoked_trace_On_f33_pos(2,2)+0.02, 0.12 , vert_dim];
SNR_pos_top = SNR_pos(1,2)+SNR_pos(1,4);

% Modulation Index
Modulation_index = open('Modulation index.fig');    
Modulation_index_ax = get(gcf, 'children');
Modulation_index_pos(1,:) = [h_dist , SNR_pos(1,2)+vert_dim+0.08 , 0.12 , vert_dim];
Modulation_index_pos_top = Modulation_index_pos(1,2)+Modulation_index_pos(1,4);

%Response Modulation
Response_modulation = open('Response_modulation.fig');    
Response_modulation_ax = get(gcf, 'children');
Response_modulation_pos(1,:) = [h_dist , Modulation_index_pos(1,2)+vert_dim+0.1 , 0.12 , vert_dim];
Response_modulation_pos_top = Response_modulation_pos(1,2)+Response_modulation_pos(1,4);

%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 21 29.7]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
%% xlimits, y limits, ticks etc.
 x1limits = [0 3];
x1ticks = [0, 1, 2, 3];
x1ticklab=[0, 1, 2, 3]; %[0, 0.5, 1, 1.5, 2];
y1limits(1,:) = [0 25]; y1limits(2,:) = [0 70];
y1ticks(1,:) = [0 10 20];  y1ticks(2,:) = [0 30 60];  
        
        x2ticklab=[]; %[0, 0.5, 1, 1.5, 2];
       y2limits(1,:) = [0 25]; y2limits(2,:) = [0 200];
       y2ticks(1,:) = [0 10 20];  y2ticks(2,:) = [0 100 200];  
        y2ticklab=[];
              x3limits = [0 3];
               x3ticks = [0, 1, 2, 3];
               x3ticklab=[0, 1, 2, 3]; %[0, 0.5, 1, 1.5, 2];
               y3limits(1,:) = [0 30]; y3limits(2,:) = [0 20];
               y3ticks(1,:) = [0 15 30];  y3ticks(2,:) = [0 10 20];  

x4limits = [0 3];
x4ticks = [0, 1, 2, 3];
x4ticklab=[0, 1, 2, 3]; %[0, 0.5, 1, 1.5, 2];
y4limits(1,:) = [0 30]; y4limits(2,:) = [0 40];
y4ticks(1,:) = [0 15 30];  y4ticks(2,:) = [0 20 40];  

        x5limits = [6 8];
        x5ticks = [6, 6.5, 7, 7.5, 8];
        x5ticklab=[0, 0.5, 1, 1.5, 2];
        y5limits = [-20 30];
        y5ticks = [-20 0 20];
        y5ticklab=[];
             x6limits = [-0.2 0.2];
            x6ticks = [-0.2, 0, 0.2];
            x6ticklab=[-0.2, 0, 0.2];
            y6limits = [-0.5 0.5];
            y6ticks = [-0.5 0 0.5];
            %%
%Placing plots in the figure:
ax_fontsize=14;
%Cell 23 - spont
spont_trace_Off_f23_ax_copy = copyobj(spont_trace_Off_f23_ax,F); % copy axes to new fig
set(spont_trace_Off_f23_ax_copy(2),'position',spont_trace_Off_f23_pos(1,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:));
   set(F, 'currentaxes', spont_trace_Off_f23_ax_copy(2));  t=title(''); ylabel('Trial', 'fontsize', 13)
   set(spont_trace_Off_f23_ax_copy(1),'position',spont_trace_Off_f23_pos(2,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y1limits(2,:),'ytick', y1ticks(2,:)); %,'xlim',x1limits,'xtick',x1ticks, 'xticklabel',[]
 set(F, 'currentaxes', spont_trace_Off_f23_ax_copy(1)); xl=xlabel('');  t=title(''); ylabel('Spike/S', 'fontsize', 13); %h = get(gca, 'ylabel'); set(h, 'fontsize', 11)

spont_trace_On_f23_ax_copy = copyobj(spont_trace_On_f23_ax,F); % copy axes to new fig
set(spont_trace_On_f23_ax_copy(2),'position',spont_trace_On_f23_pos(1,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:));
   set(F, 'currentaxes', spont_trace_On_f23_ax_copy(2));   t=title(''); ylabel('Trial', 'fontsize', 13)
    set(spont_trace_On_f23_ax_copy(1),'position',spont_trace_On_f23_pos(2,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y1limits(2,:),'ytick', y1ticks(2,:)); %,'xlim',x1limits,'xtick',x1ticks, 'xticklabel',x1ticklab
  set(F, 'currentaxes', spont_trace_On_f23_ax_copy(1)); xl=xlabel(''); t=title('');  ylabel('Spike/S', 'fontsize', 13)
 
%    l=legend('Vm','LFP'); set(l,'box','off'), set(l,'FontSize', 8) ; set(l,'linewidth',1.5);set(l,'position',[0.73 0.89 0.03 0.03]);

   %Cell 23 evoked
   evoked_trace_Off_f23_ax_copy = copyobj(evoked_trace_Off_f23_ax,F); % copy axes to new fig
set(evoked_trace_Off_f23_ax_copy(2),'position',evoked_trace_Off_f23_pos(1,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:));
   set(F, 'currentaxes', evoked_trace_Off_f23_ax_copy(2)); t=title(''); yl=ylabel(''); xl=xlabel('');  
   set(evoked_trace_Off_f23_ax_copy(1),'position',evoked_trace_Off_f23_pos(2,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y2limits(2,:),'ytick', y2ticks(2,:) ); %'xlim',x1limits,'xtick',x1ticks, 'xticklabel',[]
  set(F, 'currentaxes', evoked_trace_Off_f23_ax_copy(1)); t=title(''); yl=ylabel(''); xl=xlabel('');  
  
evoked_trace_On_f23_ax_copy = copyobj(evoked_trace_On_f23_ax,F); % copy axes to new fig
set(evoked_trace_On_f23_ax_copy(2),'position',evoked_trace_On_f23_pos(1,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y1limits(1,:),'ytick', y1ticks(1,:));
   set(F, 'currentaxes', evoked_trace_On_f23_ax_copy(2)); t=title(''); yl=ylabel(''); xl=xlabel('');  
    set(evoked_trace_On_f23_ax_copy(1),'position',evoked_trace_On_f23_pos(2,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y2limits(2,:),'ytick', y2ticks(2,:)); %,'xlim',x1limits,'xtick',x1ticks, 'xticklabel',x1ticklab
 set(F, 'currentaxes', evoked_trace_On_f23_ax_copy(1)); t=title(''); yl=ylabel('');  xl=xlabel('');  
 
 %Cell 33 - spont
spont_trace_Off_f33_ax_copy = copyobj(spont_trace_Off_f33_ax,F); % copy axes to new fig
set(spont_trace_Off_f33_ax_copy(2),'position',spont_trace_Off_f33_pos(1,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y3limits(1,:),'ytick', y3ticks(1,:));
    set(F, 'currentaxes', spont_trace_Off_f33_ax_copy(2));  t=title(''); ylabel('Trial', 'fontsize', 13)
   set(spont_trace_Off_f33_ax_copy(1),'position',spont_trace_Off_f33_pos(2,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y3limits(2,:),'ytick', y3ticks(2,:)); %,'xlim',x1limits,'xtick',x1ticks, 'xticklabel',[]
 set(F, 'currentaxes', spont_trace_Off_f33_ax_copy(1)); t=title('');  xl=xlabel(''); ylabel('Spike/S', 'fontsize', 13)
 
spont_trace_On_f33_ax_copy = copyobj(spont_trace_On_f33_ax,F); % copy axes to new fig
set(spont_trace_On_f33_ax_copy(2),'position',spont_trace_On_f33_pos(1,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y3limits(1,:),'ytick', y3ticks(1,:));
    set(F, 'currentaxes', spont_trace_On_f33_ax_copy(2));  t=title(''); ylabel('Trial', 'fontsize', 13)
    set(spont_trace_On_f33_ax_copy(1),'position',spont_trace_On_f33_pos(2,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y3limits(2,:),'ytick', y3ticks(2,:)); %,'xlim',x1limits,'xtick',x1ticks, 'xticklabel',x1ticklab
   set(F, 'currentaxes', spont_trace_On_f33_ax_copy(1)); t=title('');ylabel('Spike/S', 'fontsize', 13); xlabel('Time [S]', 'fontsize', 13);
%   set(F, 'currentaxes', spont_trace_On_f33_ax_copy); t=title('NB+'); yl=ylabel(''); xl=xlabel('');
%    l=legend('Vm','LFP'); set(l,'box','off'), set(l,'FontSize', 8) ; set(l,'linewidth',1.5);set(l,'position',[0.73 0.89 0.03 0.03]);

   %Cell 33 evoked
   evoked_trace_Off_f33_ax_copy = copyobj(evoked_trace_Off_f33_ax,F); % copy axes to new fig
set(evoked_trace_Off_f33_ax_copy(2),'position',evoked_trace_Off_f33_pos(1,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y4limits(1,:),'ytick', y4ticks(1,:));
   set(F, 'currentaxes', evoked_trace_Off_f33_ax_copy(2)); t=title(''); yl=ylabel(''); xl=xlabel('');  
   set(evoked_trace_Off_f33_ax_copy(1),'position',evoked_trace_Off_f33_pos(2,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y4limits(2,:),'ytick', y4ticks(2,:)); % 'xlim',x1limits,'xtick',x1ticks, 'xticklabel',[] 
  set(F, 'currentaxes', evoked_trace_Off_f33_ax_copy(1)); t=title(''); yl=ylabel(''); xl=xlabel('');  
  
evoked_trace_On_f33_ax_copy = copyobj(evoked_trace_On_f33_ax,F); % copy axes to new fig
set(evoked_trace_On_f33_ax_copy(2),'position',evoked_trace_On_f33_pos(1,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y4limits(1,:),'ytick', y4ticks(1,:));
   set(F, 'currentaxes', evoked_trace_On_f33_ax_copy(2)); t=title(''); yl=ylabel(''); xl=xlabel('');  
    set(evoked_trace_On_f33_ax_copy(1),'position',evoked_trace_On_f33_pos(2,:),...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5, 'ylim', y4limits(2,:),'ytick', y4ticks(2,:)); %,'xlim',x1limits,'xtick',x1ticks, 'xticklabel',x1ticklab
 set(F, 'currentaxes', evoked_trace_On_f33_ax_copy(1)); t=title(''); yl=ylabel(''); xlabel('Time [S]', 'fontsize', 13);
 
 %Response modulation
Response_modulation_ax_copy = copyobj(Response_modulation_ax,F); % copy axes to new fig
set(Response_modulation_ax_copy,'position',Response_modulation_pos(1,:),'ylim',[-5 30] ,'ytick', [0 10 20],...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5,'box','off'); %, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks, 'xticklabel',x6ticklab);
 set(F, 'currentaxes', Response_modulation_ax_copy); t=title(''); yl=ylabel('Spikes/Stim. train','fontsize',13);

 %Modulation index
 Modulation_index_ax_copy = copyobj(Modulation_index_ax,F); % copy axes to new fig
set(Modulation_index_ax_copy,'position',Modulation_index_pos(1,:),'ylim',[-1 1] ,'ytick', [-1 0 1],...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5,'box','off'); %, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks, 'xticklabel',x6ticklab);
 set(F, 'currentaxes', Modulation_index_ax_copy); t=title(''); yl=ylabel('Modulation Index','fontsize',13);
 %SNR
SNR_ax_copy = copyobj(SNR_ax,F); % copy axes to new fig
set(SNR_ax_copy,'position',SNR_pos(1,:),'ylim',[-0.1 1.1] ,'ytick', [0 0.5 1],...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5,'box','off'); %, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks, 'xticklabel',x6ticklab);
 set(F, 'currentaxes', SNR_ax_copy); t=title(''); yl=ylabel('SNR Index','fontsize',13); %set Ytick to be on the left side and not on the right


%annotations:
annotation('textbox', [spont_trace_Off_f23_pos(1,1) spont_trace_Off_f23_pos_top 0 0]+[-0.02 0.02 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous Activity', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
annotation('textbox', [evoked_trace_Off_f23_pos(1,1) evoked_trace_Off_f23_pos_top 0 0]+[-0.02 0.02 0.5 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Sensory-evoked Responses', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
 
 a_pos1=[-0.03 0 0.04 0.04];
 
annotation('textbox', [spont_trace_Off_f23_pos(1,1) spont_trace_Off_f23_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
annotation('textbox', [spont_trace_On_f23_pos(1,1) spont_trace_On_f23_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [evoked_trace_Off_f23_pos(1,1) evoked_trace_Off_f23_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [evoked_trace_On_f23_pos(1,1) evoked_trace_On_f23_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_trace_Off_f33_pos(1,1) spont_trace_Off_f33_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_trace_On_f33_pos(1,1) spont_trace_On_f33_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [evoked_trace_Off_f33_pos(1,1) evoked_trace_Off_f33_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [evoked_trace_On_f33_pos(1,1) evoked_trace_On_f33_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'H', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [Response_modulation_pos(1,1) Response_modulation_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'I', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [Modulation_index_pos(1,1) Modulation_index_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'J', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
  annotation('textbox', [SNR_pos(1,1) SNR_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'K', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 
%%
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures'
if save_flag==1  
    filename='Fig 2 Extracellular_v2';
    saveas(F,'Fig 2 Extracellular.fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end