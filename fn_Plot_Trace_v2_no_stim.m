%% Plotting a single trace
%stim_X is a two-row matrix, row 1:locations of laser pulse begin, row 2:
%locations of laser pulse end.
%if raw_data is a matrix, each column will be plotted separately on the
%same axes.
%F is the handle of the figure
%obj is the handle of the plotted traces, without the plot of the stimulus

function [F,obj] = fn_Plot_Trace_v2_no_stim(raw_data,dt_data, dt_stim1, stim1_X, dt_stim2, stim2_X, lineProps);


% cd 'C:\Users\inbalme\Dropbox\Inbal M.Sc\MATLAB\Project Brain States'
% cd 'D:\Inbal M.Sc\MATLAB\Project Brain States'

% max_time_axis = ceil(size(raw_data,1).*dt_data);
max_time_axis = size(raw_data,1).*dt_data;
max_data = max(max(raw_data));
min_data = min(min(raw_data));
five_percent = (max_data-min_data).*0.02;
min_data_wide = min_data-2.*five_percent;
max_data_wide = max_data+2.*five_percent;

 x1limits = [0 max_time_axis];
 x1ticks = [0 0.5.*max_time_axis max_time_axis];
 if abs(min_data_wide) < 1
     min_data_wide = (floor(min_data_wide.*10))./10;
 else
     min_data_wide = floor(min_data_wide);
 end
 if abs(max_data_wide) < 1 
      max_data_wide = (ceil(max_data_wide.*10))./10;
 else
      max_data_wide = ceil(max_data_wide);
 end
 y1limits = [min_data_wide max_data_wide]; %[-0.3 0.5]; [0.5 1.5]
 y1ticks =  [min_data_wide max_data_wide]; %[-0.3 0 0.5]; %y1limits; [0.5 1 1.5]

     color_table = [0 0 0; 0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0];
     if size(raw_data,2) > size(color_table,1)
         color_table = rand(size(raw_data,2),3);
     end

F=figure;
        xlabel('Time [sec]' ,'FontSize', 16);
        ylabel('Vm [mV]', 'FontSize', 16);        

hold on
for i = 1:size(raw_data,2);
    %Set default options
defaultProps={'LineWidth',1.5,'color', color_table(i,:)};
if nargin<7, lineProps=defaultProps; end
if isempty(lineProps), lineProps=defaultProps; end
if ~iscell(lineProps), lineProps={lineProps}; end
    obj(i)=plot((1:size(raw_data,1)).*dt_data,raw_data(:,i), lineProps{:});
    
end
set(gca, 'xlim', x1limits, 'ylim', y1limits)%,'xtick', x1ticks

if isempty(stim1_X); else
%     stim1_Y = ones(size(stim1_X)).*(max_data+five_percent); 
%     line(stim1_X.*dt_stim1,stim1_Y,'LineWidth',6,'Color','b')
patch_xdata=[stim1_X; flipud(stim1_X)];
ylim_data=[get(gca,'ylim')]';
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt_stim1,patch_ydata,patch_cdata,'faceColor','k','edgecolor','none','faceAlpha', 0.3);
set(gca,'linewidth',1.2)
end

if isempty(stim2_X); else
%     stim2_Y = ones(size(stim2_X)).*(max_data+five_percent); 
%     line(stim2_X.*dt_stim2,stim2_Y,'LineWidth',6,'Color','r')
patch_xdata=[stim2_X; flipud(stim2_X)];
ylim_data=[get(gca,'ylim')]';
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt_stim2,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
set(gca,'linewidth',1.2)
set(gca, 'xlim', x1limits, 'ylim', y1limits)%,'xtick', x1ticks
end
hold off
end