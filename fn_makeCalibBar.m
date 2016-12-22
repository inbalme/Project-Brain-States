function  [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes)
gca;
hold on;
 perc1 = 0.02;%percentage of axis length to set distance of bar from axes
 perc2 = 0.02;%percentage of axis length to set distance of text from bar

xl = get(gca,'xlim');
yl = get(gca,'ylim');
if horiz_vert;
    barXstart=max(xl)+perc1*diff(xl)-lengthh;
    barXend=max(xl)+perc1*diff(xl);
    barYstart=min(yl)-perc1*diff(yl);
    barYend=barYstart; min(yl)-perc1*diff(yl)+lengthh;
    pos = get(gca,'position');  
    y_units=perc2*diff(yl);
    p1 = plot([barXstart, barXend],[barYstart, barYend],'linewidth',2,'color',c);
    p2 = text(barXstart+lengthh/2,barYstart-y_units,textit,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',fonsizes);
    set(gca,'ylim',[barYstart-y_units, yl(2)])
else
    barXstart=max(xl)+perc1*diff(xl);
    barXend=barXstart;
    barYstart=min(yl);
    barYend=min(yl)+lengthh;
    %using the position of the axes to express units of x-axis and y-axis in terms of absolut length of x and y axis
    pos = get(gca,'position');
    y_abs=perc2*pos(4)+0.01; %changing from y-axis units to absolute units
    x_units=y_abs/pos(3)*diff(xl); %changing from absolute units to x-axis units
   p1 = plot([barXstart, barXend],[barYstart, barYend],'linewidth',2,'color',c);
    p2 = text(barXstart+x_units,barYstart+lengthh/2,textit,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',fonsizes);
   set(gca,'xlim',[xl(1),barXend+x_units]); 
end;



%need to set the new ylim and xlim of the axes to contain the scale bar.   
    

