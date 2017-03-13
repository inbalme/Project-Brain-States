%% plotting a figure 

close all
clear all
save_flag=0;
exp_type=1; %1- NBES, 2- ChAT
%opening saved figures:
switch exp_type
    case 1
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm Histograms and Powerspec\Spontaneous histograms'
        spont_Vm_hist = open('f46_Vm_histogram_ongoing.fig');
        spont_Vm_hist_ax = get(gcf, 'children');
    case 2
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm Histograms and Powerspec\Spontaneous histograms'
        spont_Vm_hist = open('f80_Vm_histogram_ongoing.fig');
        spont_Vm_hist_ax = get(gcf, 'children');
end

spont_Vm_5prcentile = open('Vm_5prcentile_ongoing.fig');    
spont_Vm_5prcentile_ax = get(gcf, 'children');
switch exp_type
    case 1
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous\detection f46'
        detection_example1= open('detection example - ongoing NB-, trace 5.fig');    
        a=findall(gca,'type','line'); delete(a(2));
        b=findall(gca,'type','text'); delete(b(2));
        detection_example1_ax = get(gcf, 'children');
        detection_example2= open('detection example - ongoing NB+, trace 3.fig'); 
        a=findall(gca,'type','line'); delete(a(1:2));
        b=findall(gca,'type','text'); delete(b(1:2));
        detection_example2_ax = get(gcf, 'children');
         cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous\detection f16'
        detection_example3= open('detection example - ongoing NB-, trace 6.fig');    
        a=findall(gca,'type','line'); %delete(a(2));
        b=findall(gca,'type','text'); %delete(b(2));
        detection_example3_ax = get(gcf, 'children');
        detection_example4= open('detection example - ongoing NB+, trace 3.fig');    
        a=findall(gca,'type','line'); delete(a(1:2));
        b=findall(gca,'type','text'); delete(b(1:2));
        detection_example4_ax = get(gcf, 'children');
    case 2
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous\detection f80'
end

switch exp_type
    case 1
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous'
    case 2
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous'
end
spont_event_freq = open('Spontaneous event frequency.fig');    
spont_event_freq_ax = get(gcf, 'children');

spont_event_amp = open('Spontaneous event amplitude.fig');    
spont_event_amp_ax = get(gcf, 'children');

% spont_event_halfwidth = open('Spontaneous event half-width.fig');    
% spont_event_halfwidth_ax = get(gcf, 'children');
% 
% spont_event_risetime = open('Spontaneous event rise-time.fig');    
% spont_event_risetime_ax = get(gcf, 'children');

switch exp_type
    case 1
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Rin'
    case 2
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Rin'
end

spont_Rin = open('Rin.fig');    
spont_Rin_ax = get(gcf, 'children');

%%
%open a new figure:
F = figure;
set(gcf,'color','w');
set(gcf,'DefaultAxesFontSize',18);
set(gcf,'DefaultAxesFontName','arial');
set(gcf, 'PaperType', 'A4');
set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 25 25]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
 
%% Positions:

spont_Vm_hist_pos(1,:)=[0.1, 0.75, 0.4, 0.2];
h_dist1=spont_Vm_hist_pos(1,1)+spont_Vm_hist_pos(1,3)+0.2;
spont_Vm_5prcentile_pos(1,:) = [h_dist1 , spont_Vm_hist_pos(1,2)-0.02 ,  0.2 ,  spont_Vm_hist_pos(1,4)];

detection_example1_pos(1,:) = [spont_Vm_hist_pos(1,1) , 0.5 , 0.4 , 0.1]; detection_example1_pos(1,2)=spont_Vm_hist_pos(1,2)-spont_Vm_hist_pos(1,4)-0.02;
detection_example2_pos(1,:) = [detection_example1_pos(1,1) , 0.5 , detection_example1_pos(1,3) , detection_example1_pos(1,4)]; detection_example2_pos(1,2)=detection_example1_pos(1,2)-detection_example1_pos(1,4)+0.05;
detection_example3_pos(1,:) = [detection_example1_pos(1,1) , 0.5 , detection_example1_pos(1,3) , detection_example1_pos(1,4)]; detection_example3_pos(1,2)=detection_example2_pos(1,2)-detection_example2_pos(1,4)+0.05;
detection_example4_pos(1,:) = [detection_example1_pos(1,1) , 0.5 , detection_example1_pos(1,3) , detection_example1_pos(1,4)]; detection_example4_pos(1,2)=detection_example3_pos(1,2)-detection_example3_pos(1,4)+0.05;

spont_Rin_pos(1,:) = [spont_Vm_hist_pos(1,1), detection_example4_pos(1,2)-spont_Vm_5prcentile_pos(1,4)-0.07 ,  spont_Vm_5prcentile_pos(1,3) ,  spont_Vm_5prcentile_pos(1,4)];
h_dist2=spont_Rin_pos(1,1)+spont_Rin_pos(1,3)+0.1;
spont_event_freq_pos(1,:) =[h_dist2,spont_Rin_pos(1,2),spont_Rin_pos(1,3),spont_Rin_pos(1,4)];
spont_event_amp_pos(1,:)=[h_dist2+spont_event_freq_pos(1,3)+0.1,spont_event_freq_pos(1,2),spont_event_freq_pos(1,3),spont_event_freq_pos(1,4)];
%top positions:
spont_Vm_hist_pos_top = spont_Vm_hist_pos(1,2)+spont_Vm_hist_pos(1,4);
spont_Vm_5prcentile_pos_top = spont_Vm_5prcentile_pos(1,2)+spont_Vm_5prcentile_pos(1,4);
spont_event_freq_pos_top = spont_event_freq_pos(1,2)+spont_event_freq_pos(1,4);
spont_event_amp_pos_top = spont_event_amp_pos(1,2)+spont_event_amp_pos(1,4);
spont_Rin_pos_top = spont_Rin_pos(1,2)+spont_Rin_pos(1,4);
% spont_event_halfwidth_pos_top = spont_event_halfwidth_pos(1,2)+spont_event_halfwidth_pos(1,4);
% spont_event_risetime_pos_top = spont_event_risetime_pos(1,2)+spont_event_risetime_pos(1,4);

%%
%Placing plots in the figure:
spont_Vm_hist_ax_copy = copyobj(spont_Vm_hist_ax,F); % copy axes to new fig
set(spont_Vm_hist_ax_copy(2),'position',spont_Vm_hist_pos(1,:))
set(F, 'currentaxes', spont_Vm_hist_ax_copy(2));  tl=title('');
spont_Vm_hist_ax_copy(2).FontSize=12;
% spont_Vm_hist_ax_copy.XTickLabel=[];
%placing legend:
set(spont_Vm_hist_ax_copy(1),'position',[spont_Vm_hist_pos(1,1)+spont_Vm_hist_pos(1,3)-0.1,spont_Vm_hist_pos_top-0.05, 0.2,0.05])

spont_Vm_5prcentile_ax_copy = copyobj(spont_Vm_5prcentile_ax,F); % copy axes to new fig
set(spont_Vm_5prcentile_ax_copy,'position',spont_Vm_5prcentile_pos(1,:))
set(F, 'currentaxes', spont_Vm_5prcentile_ax_copy);   tl=title('');
yl=ylabel({'Lower 5th'; 'percentile [mV]'});
spont_Vm_5prcentile_ax_copy.FontSize=12;

detection_example1_ax_copy = copyobj(detection_example1_ax,F); % copy axes to new fig
set(detection_example1_ax_copy,'position',detection_example1_pos(1,:))

detection_example2_ax_copy = copyobj(detection_example2_ax,F); % copy axes to new fig
set(detection_example2_ax_copy,'position',detection_example2_pos(1,:))

detection_example3_ax_copy = copyobj(detection_example3_ax,F); % copy axes to new fig
set(detection_example3_ax_copy,'position',detection_example3_pos(1,:))

detection_example4_ax_copy = copyobj(detection_example4_ax,F); % copy axes to new fig
set(detection_example4_ax_copy,'position',detection_example4_pos(1,:))

set(findall(gcf,'type','scatter'),'sizedata',20);

spont_event_freq_ax_copy = copyobj(spont_event_freq_ax,F); % copy axes to new fig
set(spont_event_freq_ax_copy,'position',spont_event_freq_pos(1,:))
set(F, 'currentaxes', spont_event_freq_ax_copy);  xl=xlabel('');  tl=title('');
yl=ylabel('Frequency [Hz]');
spont_event_freq_ax_copy.FontSize=12;
% spont_event_freq_ax_copy.XTickLabel=[];

  spont_event_amp_ax_copy = copyobj(spont_event_amp_ax,F); % copy axes to new fig
set(spont_event_amp_ax_copy,'position',spont_event_amp_pos(1,:))
set(F, 'currentaxes', spont_event_amp_ax_copy);  xl=xlabel('');  tl=title(''); 
yl=ylabel('Amplitude [mV]');
spont_event_amp_ax_copy.FontSize=12;
% spont_event_amp_ax_copy.XTickLabel=[];

spont_Rin_ax_copy = copyobj(spont_Rin_ax,F); % copy axes to new fig
set(spont_Rin_ax_copy,'position',spont_Rin_pos(1,:))
set(F, 'currentaxes', spont_Rin_ax_copy);   tl=title(''); 
% yl=ylabel('Rin [MOhm]');
spont_Rin_ax_copy.FontSize=12;

% spont_event_halfwidth_ax_copy = copyobj(spont_event_halfwidth_ax,F); % copy axes to new fig
% set(spont_event_halfwidth_ax_copy,'position',spont_event_halfwidth_pos(1,:))
% set(F, 'currentaxes', spont_event_halfwidth_ax_copy);  tl=title('');
% yl=ylabel('Half-Width [mS]');
% spont_event_halfwidth_ax_copy.FontSize=12;

% spont_event_risetime_ax_copy = copyobj(spont_event_risetime_ax,F); % copy axes to new fig
% set(spont_event_risetime_ax_copy,'position',spont_event_risetime_pos(1,:))
% set(F, 'currentaxes', spont_event_risetime_ax_copy);   tl=title(''); 
% yl=ylabel('Rise-time [mS]');
% spont_event_risetime_ax_copy.FontSize=12;

 a_pos1=[-0.06 0.02 0.04 0.04];
 annotation('textbox', [spont_Vm_hist_pos(1,1),spont_Vm_hist_pos_top+0.03 0.4 0.05],...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'Spontaneous activity', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')

 annotation('textbox', [spont_Vm_hist_pos(1,1),spont_Vm_hist_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'A', 'FontName','arial', 'fontsize', 14, 'fontweight', 'bold')
annotation('textbox', [spont_Vm_5prcentile_pos(1,1),spont_Vm_5prcentile_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'B', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_event_freq_pos(1,1) spont_event_freq_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'C', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_event_amp_pos(1,1) spont_event_amp_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'D', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
 annotation('textbox', [spont_Rin_pos(1,1) spont_Rin_pos_top 0 0]+a_pos1,...
     'FitHeightToText', 'on', 'edgecolor', 'none','string', 'E', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
%  annotation('textbox', [spont_event_halfwidth_pos(1,1) spont_event_halfwidth_pos_top 0 0]+a_pos1,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'F', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')
%   annotation('textbox', [spont_event_risetime_pos(1,1) spont_event_risetime_pos_top 0 0]+a_pos1,...
%      'FitHeightToText', 'on', 'edgecolor', 'none','string', 'G', 'FontName','arial', 'fontsize', 12, 'fontweight', 'bold')


cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Paper Figures'
if save_flag==1
    switch exp_type
        case 1
            filename='Fig 4 spont events_v3';
        case 2
            filename='Fig 4 ChAT spont events_v3';
    end
    saveas(F,filename,'fig'); 
    print(F,filename,'-dpng','-r600','-opengl') 
    print(F, '-depsc2', filename);
end