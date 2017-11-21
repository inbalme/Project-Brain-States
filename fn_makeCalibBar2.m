function  [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes, perc1,perc2)
%This function was created by Inbal, based on the function of Ilan "makeCalibBar". 
%The inputs to the function are:
%horiz_vert =1 for horizontal bar and 0 for vertical bar
%lengthh is the bar length. The units are same as the axis units.
%textit=the text next to the bar
%c is the bar color
%fontsize 
%perc1 = percentage of axis length to set distance of bar from axes
%perc2 = percentage of axis length to set distance of text from the bar
%if perc1 or perc2 =[] then the function will use default values.
%An example for using it to make both horizontal and vertical bar, where x-axis is time in [Sec] and y-axis is voltage in [mV]:
        %plotting scale bar
% horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12;perc1=[]; perc2=[];
%         [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
%  horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh), ' mV'];     c=[0,0,0];  fonsizes=12;perc1=[]; perc2=[];
%         [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
% comment: I tried to make it such that the distance between the text and
% the scale bar would be the same for both x-axis bar and y-axis. This was
% not possible because the function "text" plots the text in units
% according to the axes and not in absolute units as "centimeters". 

%The default perc1 and perc2:
if isempty(perc1)
    perc1 = 0.02;%percentage of axis length to set distance of bar from axes
end
if isempty(perc2)
perc2 = 0.02;%percentage of axis length to set distance of text from bar
end

gca;
% set(gca,'units','centimeters');
% set(gcf,'Paperunits','centimeters');
hold on;
 xl = get(gca,'xlim');
yl = get(gca,'ylim');
if horiz_vert; %for horizontal bar
    barXstart=max(xl)+perc1*diff(xl)-lengthh;
    barXend=max(xl)+perc1*diff(xl);
    barYstart=min(yl)-2*perc1*diff(yl);
    barYend=barYstart; 
    pos = get(gca,'position');  
    x_units=2*perc2*diff(xl);  
    y_units=perc2*diff(yl);
    p1 = plot([barXstart, barXend],[barYstart, barYend],'linewidth',1.5,'color',c);
    p2 = text(barXstart+lengthh/2,barYstart-y_units,textit,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',fonsizes);
    set(gca,'xlim',[xl(1),barXend+x_units]);  %setting new xlim for the axes to contain the scale bar
    set(gca,'ylim',[barYstart-y_units, yl(2)]) %setting new ylim for the axes to contain the scale bar
else %for vertical bar
    barXstart=max(xl)+perc1*diff(xl);
    barXend=barXstart;
    barYstart=min(yl);
    barYend=min(yl)+lengthh;
    x_units=4*perc2*diff(xl);  %if there is 90 degree rotation then use x_units=2;
     y_units=perc2*diff(yl);
    pos = get(gca,'position');
    set(gca,'xlim',[xl(1),barXend+x_units]);  %setting new xlim for the axes to contain the scale bar
    set(gca,'ylim',[barYstart-y_units, yl(2)]) %setting new ylim for the axes to contain the scale bar
   p1 = plot([barXstart, barXend],[barYstart, barYend],'linewidth',1.5,'color',c);
    p2 = text(barXstart+x_units,barYstart+lengthh/2,textit,'HorizontalAlignment', 'center','VerticalAlignment', 'middle', 'fontsize',fonsizes); %'rotation',90,
end;



    

