%% plotting a figure 

%% new version

close all
clear all
save_flag=0;
no_numbering_flag=0;
abslen=0.05; %[in cm]
ax_fontsize=12;

%opening saved figures:

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH'

% Modulation Index
Modulation_index = open('Modulation index.fig');    
Modulation_index_ax = get(gcf, 'children');    

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 6 10]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
 
%% Positions:
%NB
Modulation_index_pos(1,:) = [0.2,0.08, 0.7, 0.8];

%top positions:
Modulation_index_pos_top = Modulation_index_pos(1,2)+Modulation_index_pos(1,4);

% NB_illustration_pos(1,:)=[spont_amp_hist_pos(1,1)+0.1,spont_amp_hist_pos_top+0.05, 0.15, 0.15];
% NB_illustration_pos_top=NB_illustration_pos(1,2)+NB_illustration_pos(1,4);
%%
%Placing plots in the figure:

%Modulation index
 Modulation_index_ax_copy = copyobj(Modulation_index_ax,F); % copy axes to new fig
set(Modulation_index_ax_copy,'position',Modulation_index_pos(1,:),'ylim',[-1 1] ,'ytick', [-1 0 1],...
       'fontname', 'arial','fontsize',ax_fontsize,'linewidth',1.5,'box','off','tickdir','out'); %, 'ylim', y6limits,'ytick', y6ticks,'xlim',x6limits,'xtick',x6ticks, 'xticklabel',x6ticklab);
 set(F, 'currentaxes', Modulation_index_ax_copy); t=title(''); yl=ylabel('Modulation Index','fontsize',ax_fontsize);
 
% NB_illustration_ax = axes('position',NB_illustration_pos);
% imshow(NB_illustration, 'parent', NB_illustration_ax) 


    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'

if save_flag==1   
    filename='Fig S2';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end