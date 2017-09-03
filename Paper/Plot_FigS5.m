%% plotting figure S5 - event detection

close all
clear all
save_flag=0;
no_numbering_flag=0;
abslen=0.05; %[in cm]
ax_fontsize=10;
dotsize=15; %markersize of dots for detected events
color1=[255, 102,102]./256; %pink
color2=[102, 172,255]./256; %light blue
%opening saved figures:

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

        
cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous\detection f80'
    color_onset=[255,178,102]./256; color_top=[255,153,16]./256;
        detection_example1_ChAT= open('detection example - ongoing NB-, trace 5.fig');    
        a=findall(gca,'type','line'); delete(a(2));
        b=findall(gca,'type','text'); delete(b(2)); set(b(1),'rotation',0)
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
        b=findall(gca,'type','text'); set(b(1),'rotation',0) %delete(b(2));
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

            cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
NB_illustration = imread('NBES_Schematic_Illustration_one-electrode_tight','TIF');

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Illustrations+Histology'
ChAT_illustration = imread('ChAT_Schematic_Illustration_one-electrode_tight','TIF');

        %%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 18 15]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

%%
detection_example1_pos(1,:) = [0.05 , 0.68 , 0.45 , 0.25]; 
hdist3=detection_example1_pos(1,1)+detection_example1_pos(1,3)+0.05;
detection_example2_pos(1,:) = [hdist3 , detection_example1_pos(1,2) , detection_example1_pos(1,3) , detection_example1_pos(1,4)]; %detection_example2_pos(1,2)=detection_example1_pos(1,2)-detection_example1_pos(1,4)+0.05;
detection_example3_pos(1,:) = [detection_example1_pos(1,1) , 0.5 , detection_example1_pos(1,3) , detection_example1_pos(1,4)]; detection_example3_pos(1,2)=detection_example1_pos(1,2)-detection_example3_pos(1,4)+0.07;
detection_example4_pos(1,:) = [hdist3 , detection_example3_pos(1,2) , detection_example1_pos(1,3) , detection_example1_pos(1,4)]; %detection_example4_pos(1,2)=detection_example3_pos(1,2)-detection_example3_pos(1,4)+0.05;

detection_example1_ChAT_pos(1,:) = [0.05 , 0.21 , 0.45 , 0.25]; 
hdist3=detection_example1_ChAT_pos(1,1)+detection_example1_ChAT_pos(1,3)+0.05;
detection_example2_ChAT_pos(1,:) = [hdist3 , detection_example1_ChAT_pos(1,2) , detection_example1_ChAT_pos(1,3) , detection_example1_ChAT_pos(1,4)]; %detection_example2_pos(1,2)=detection_example1_pos(1,2)-detection_example1_pos(1,4)+0.05;
detection_example3_ChAT_pos(1,:) = [detection_example1_ChAT_pos(1,1) , 0.5 , detection_example1_ChAT_pos(1,3) , detection_example1_ChAT_pos(1,4)]; detection_example3_ChAT_pos(1,2)=detection_example1_ChAT_pos(1,2)-detection_example3_ChAT_pos(1,4)+0.09;
detection_example4_ChAT_pos(1,:) = [hdist3 , detection_example3_ChAT_pos(1,2) , detection_example1_ChAT_pos(1,3) , detection_example1_ChAT_pos(1,4)]; %detection_example4_pos(1,2)=detection_example3_pos(1,2)-detection_example3_pos(1,4)+0.05;

detection_example1_pos_top = detection_example1_pos(1,2)+detection_example1_pos(1,4);
detection_example2_pos_top = detection_example2_pos(1,2)+detection_example2_pos(1,4);
detection_example3_pos_top = detection_example3_pos(1,2)+detection_example3_pos(1,4);
detection_example4_pos_top = detection_example4_pos(1,2)+detection_example4_pos(1,4);
detection_example1_ChAT_pos_top = detection_example1_ChAT_pos(1,2)+detection_example1_ChAT_pos(1,4);
detection_example2_ChAT_pos_top = detection_example2_ChAT_pos(1,2)+detection_example2_ChAT_pos(1,4);
detection_example3_ChAT_pos_top = detection_example3_ChAT_pos(1,2)+detection_example3_ChAT_pos(1,4);
detection_example4_ChAT_pos_top = detection_example4_ChAT_pos(1,2)+detection_example4_ChAT_pos(1,4);

NB_illustration_pos(1,:)=[0.04,detection_example1_pos_top-0.04, 0.08, 0.09];
ChAT_illustration_pos(1,:)=[NB_illustration_pos(1,1),detection_example1_ChAT_pos_top-0.06, NB_illustration_pos(1,3), NB_illustration_pos(1,4)];

NB_illustration_pos_top=NB_illustration_pos(1,2)+NB_illustration_pos(1,4);
ChAT_illustration_pos_top=ChAT_illustration_pos(1,2)+ChAT_illustration_pos(1,4);

%% placing plots in the figure

ax_fontsize=10;
color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];   

detection_example1_ax_copy = copyobj(detection_example1_ax,F); % copy axes to new fig
set(detection_example1_ax_copy,'position',detection_example1_pos(1,:))

detection_example2_ax_copy = copyobj(detection_example2_ax,F); % copy axes to new fig
set(detection_example2_ax_copy,'position',detection_example2_pos(1,:))

detection_example3_ax_copy = copyobj(detection_example3_ax,F); % copy axes to new fig
set(detection_example3_ax_copy,'position',detection_example3_pos(1,:))

detection_example4_ax_copy = copyobj(detection_example4_ax,F); % copy axes to new fig
set(detection_example4_ax_copy,'position',detection_example4_pos(1,:))

%ChAT:
color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [102, 172,255]./256; [30,75,14]/256; [112,172,90]/256];  

detection_example1_ChAT_ax_copy = copyobj(detection_example1_ChAT_ax,F); % copy axes to new fig
set(detection_example1_ChAT_ax_copy,'position',detection_example1_ChAT_pos(1,:))

detection_example2_ChAT_ax_copy = copyobj(detection_example2_ChAT_ax,F); % copy axes to new fig
set(detection_example2_ChAT_ax_copy,'position',detection_example2_ChAT_pos(1,:))

detection_example3_ChAT_ax_copy = copyobj(detection_example3_ChAT_ax,F); % copy axes to new fig
set(detection_example3_ChAT_ax_copy,'position',detection_example3_ChAT_pos(1,:))

detection_example4_ChAT_ax_copy = copyobj(detection_example4_ChAT_ax,F); % copy axes to new fig
set(detection_example4_ChAT_ax_copy,'position',detection_example4_ChAT_pos(1,:))

set(findall(gcf,'type','scatter'),'sizedata',dotsize);

NB_illustration_ax = axes('position',NB_illustration_pos);
imshow(NB_illustration, 'parent', NB_illustration_ax) 

ChAT_illustration_ax = axes('position',ChAT_illustration_pos);
imshow(ChAT_illustration, 'parent', ChAT_illustration_ax) 

% Anotations
color_text=[102,51,0]./256;
 a_pos1=[-0.06 0 0.04 0.04];

annotation('textbox', [detection_example1_pos(1,1)+0.09,detection_example1_pos_top-0.05, 0.1 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'NB-', 'FontName','arial', 'fontsize', 12, 'color','k');
 annotation('textbox', [detection_example2_pos(1,1)+0.09,detection_example1_pos_top-0.05, 0.1 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'NB+', 'FontName','arial', 'fontsize', 12, 'color','k');

 annotation('textbox', [detection_example1_pos(1,1)-0.05,detection_example1_pos_top-0.13, 0.1 0.05],...
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

 if no_numbering_flag==0;
 annotation('textbox', [detection_example1_pos(1,1)-0.05, detection_example1_pos_top-0.09 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 annotation('textbox', [detection_example2_pos(1,1)-0.03 detection_example1_pos_top-0.08 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 %   annotation('textbox', [ an1pos(1,1), detection_example3_pos_top-0.12 0.04 0.04],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
%  annotation('textbox', [detection_example4_pos(1,1)-0.03 detection_example3_pos_top-0.12 0.04 0.04],...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [detection_example1_ChAT_pos(1,1)-0.05, detection_example1_ChAT_pos_top-0.11 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 annotation('textbox', [detection_example2_ChAT_pos(1,1)-0.03 detection_example1_ChAT_pos_top-0.11 0.04 0.04],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 10, 'fontweight', 'bold')
 end

 cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures\Neuron'

 if save_flag==1   
    filename='Fig S5';
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end