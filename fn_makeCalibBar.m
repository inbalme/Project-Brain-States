function  [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes)
%This function was created by Inbal, based on the function of Ilan
%"makeCalibBar". An example for using it to make both horizontal and vertical bar, where x-axis is time in [Sec] and y-axis is voltage in [mV]:
        %plotting scale bar
% horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12;
%         [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
%  horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh), ' mV'];     c=[0,0,0];  fonsizes=12;
%         [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
% comment: I tried to make it such that the distance between the text and
% the scale bar would be the same for both x-axis bar and y-axis. This was
% not possible because the function "text" plots the text in units
% according to the axes and not in absolute units as "centimeters".


gca;
% set(gca,'units','centimeters');
% set(gcf,'Paperunits','centimeters');
hold on;
 perc1 = 0.02;%percentage of axis length to set distance of bar from axes
 perc2 = 0.02;%percentage of axis length to set distance of text from bar

xl = get(gca,'xlim');
yl = get(gca,'ylim');
if horiz_vert;
    barXstart=max(xl)+perc1*diff(xl)-lengthh;
    barXend=max(xl)+perc1*diff(xl);
    barYstart=min(yl)-perc1*diff(yl);
    barYend=barYstart; 
    pos = get(gca,'position');  
    y_units=perc2*diff(yl);
    p1 = plot([barXstart, barXend],[barYstart, barYend],'linewidth',2,'color',c);
    p2 = text(barXstart+lengthh/2,barYstart-y_units,textit,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',fonsizes);
    set(gca,'ylim',[barYstart-y_units, yl(2)]) %setting new ylim for the axes to contain the scale bar
else
    barXstart=max(xl)+perc1*diff(xl);
    barXend=barXstart;
    barYstart=min(yl);
    barYend=min(yl)+lengthh;
%     y_units=perc2*diff(yl);
x_units=perc2*diff(xl);
    %using the position of the axes to express units of x-axis and y-axis
    %in terms of absolut length of x and y axis
    pos = get(gca,'position');
%     y_abs=perc2*pos(4); %changing from y-axis units to absolute units
%     x_units=y_abs/pos(3)*diff(xl); %changing from absolute units to x-axis units
   p1 = plot([barXstart, barXend],[barYstart, barYend],'linewidth',2,'color',c);
    p2 = text(barXstart+x_units,barYstart+lengthh/2,textit,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',fonsizes);
set(gca,'xlim',[xl(1),barXend+x_units]);  %setting new ylim for the axes to contain the scale bar
end;



    

