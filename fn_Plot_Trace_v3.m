%% Plotting a single trace
%stim_X is a two-row matrix, row 1:locations of laser pulse begin, row 2:
%locations of laser pulse end.
% raw_data is a single trace
%F is the handle of the figure
%obj is the handle of the plotted traces, without the plot of the stimulus

% for file 36 header 2: 

function [f] = fn_Plot_Trace_v3(raw_data,dt_data, dt_stim1, stim1_X, dt_stim2, stim2_X);

global Exp

y_light_stim = [];
y_whisker_stim = [];

 x1limits = [0 5];
%  x1ticks = [0 5];
 y1limits = [-70 10]; 
%  y1ticks =  [-0.3 0 0.5]; 
     
f=figure(1);
clf
set(f,'units','normalized','position',[0.2 0.3 0.4 0.4],'color',[ 1 1 1]);
set( gca, 'xlim', x1limits, 'ylim', y1limits...
        'XMinorTick'; 'on'; 'ticklength'; [0.010 0.010]; ...
        'YMinorTick'; 'on'; 'ticklength';  [0.010 0.010]; ...
        'XTickLabel'; 'on';  'YTickLabel';  'on'...
        'fontname', 'helvetica';  'fontweight';  'bold';  'box';  'off' );
        xlabel('Time [sec]' ,'FontSize', 16);
        ylabel('Vm [mV]', 'FontSize', 16);
        xlabel('Time [sec]' ,'FontSize', 16);
        ylabel('Vm [mV]', 'FontSize', 16);
hold on
    plot((1:size(raw_data,1)).*dt_data,raw_data(:,i), 'LineWidth',2,'k');

if nargin>=4
    stim1_Y = ones(size(stim1_X)).*1; 
    line(stim1_X.*dt_stim1,stim1_Y,'LineWidth',6,'Color','b')
end

if nargin==6
    stim2_Y = ones(size(stim2_X)).*1; 
    line(stim2_X.*dt_stim2,stim2_Y,'LineWidth',6,'Color','r')
end
hold off