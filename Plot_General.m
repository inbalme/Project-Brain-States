%% Plot_General
% this program creates a custumized figure, on which different axes can be
% presented.

figure
set(gcf,'DefaultAxesFontSize',12);
set(gcf,'DefaultAxesFontName','helvetica');
set(gcf,'PaperUnits','centimeters','PaperPosition',[0.2 0.2 8.3 12]); %[left, bottom, width, height] 
set(gcf,'PaperOrientation','portrait');
set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[4 5 0 0]);

%This property specifies the size and location on the screen of the figure 
%window, cancels Matlab’s automatic resizing of AXES.

%Placing plots in the figure:
axes('position', [0.09 , 0.38 , 0.28 , 0.24]); %[left, bottom, width, height]

%Specifying the axis:
xlimits = [0.7 4.3];
xticks = 1 : 4 ;
ylimits = [-28 2];
yticks = [-28 0];

set( gca, 'xlim', xlimits, 'xtick', xticks, 'ylim', ylimits, 'ytick',...
[ylimits(1) 0 ylimits(2)], 'ticklength', [0.030 0.030], 'box', 'off' );
% Set the limits and ticks you defined earlier
line( xlimits, [0 0], 'color', 'k', 'linewidth', 0.5 ); 
% Place line at y = 0
text( xlimits(1)-diff(xlimits)/2.8, ylimits(1)+diff(ylimits)/2.0,...
{'\Delta Information', '(bits/spike)'}, 'fontname', 'helvetica',...
'fontsize', 7, 'rotation', 90, 'HorizontalAlignment', 'center' );
% Instead of using ylabel – use a relative placement technique 
text(xlimits(1)-diff(xlimits)/1.8, ylimits(1)+diff(ylimits)*1.08,...
'A', 'fontsize', 12, 'fontname', 'helvetica', 'fontweight', 'bold');