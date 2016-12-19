%% Load file
close all
clear all
save_flag= 1;
short_flag=0; %1- short trace, 0- long trace


%% Plotting traces+mean+STD
plot_stim_1=[1,0,1];
plot_stim_2=[0,1,1];
trace_ind = [2,5,6]; %1:size(plot_data,2);  %trace_ind is the numbers of traces to be plotted
x_axis=[];
% x_axis_long=0.5*sf{1}:10*sf{1}-1;
x_axis_long=(stim1_X{1}(2,1)+10):10*sf{1}-1; %the x_axis for the long trace
x_axis_short=stim2_X{2}(1,1)-0.5*sf{1}:stim2_X{2}(1,1)+2*sf{1}; %the x_axis for the "zoom-in" trace
if short_flag
     x_axis= x_axis_short; 
else
        x_axis= x_axis_long; 
end
clear color_table
    color_table=[0 0 0; [216 22 22]/256; [136 137 138]/256; [255 153 153]/256];
    rectangle_color=[239 239 239]/256;
    
plot_data=data_no_spikes{channel}; %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
plot_data_mean = mean(plot_data,2);  plot_data_mean(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_std =  std(data_no_spike_no_DC{channel},0,2); plot_data_std(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_var=var(data_no_spike_no_DC{channel},0,2); % plot_data_var=var(plot_data,0,2);
plot_data_ff = (-1).*plot_data_var./plot_data_mean; %fano factor

%for plotting several traces with spaces between them, I add a matrix DC to the matrix of the
%traces, where the for the next trace a larger DC is added.
       l=size(plot_data(:,trace_ind, 2),1);
        l=l/2-1;
        DC_space=20; %the space between traces
        DC=DC_space.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);        
              
%         x2lim=[stim2_X{x_value(1)}(1,1).*dt-0.5,stim2_X{2}(1,1).*dt+2];
 
%%
    %plot evoked traces no ES:
    trace_to_plot = plot_data(:,trace_ind,2)+DC ; %DC is added just to make space between the traces 
        f1=figure;
        hold on
htrace=plot(x_axis*dt,trace_to_plot(x_axis,:,:), 'LineWidth',1.2,'color', color_table(1,:));
for i=1:length(trace_ind);
    text(x_axis(1)*dt,trace_to_plot(x_axis(1),i),[num2str(floor(plot_data(x_axis(1),i,2))), ' mV '],'HorizontalAlignment','right','fontsize',13,'fontname','arial')
end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
rectangle('Position',[x_axis_short(1)*dt,ylim_data(1),(x_axis_short(end)-x_axis_short(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none')

% plotting stim_2:
if plot_stim_2(2);
patch_xdata=[stim2_X{2}; flipud(stim2_X{2})]; %stim2_X is contains the sample number of the beggining and end of the sensory stimulus. e.g. for a train at 10,20,30 ms with duration of 1 ms each stim: [10 20 30;11 21 31]*sf
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3); %facealpha is the transparency
end
%plotting scale bar
yline_start=ylim_data(1)-6; yline_end=yline_start+10;
xline_start=x_axis(1)*dt-1.5; xline_end=xline_start+1;
if isequal(x_axis,x_axis_short)
    xline_start=x_axis(1)*dt-0.5; xline_end=xline_start+0.5;
end
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',13); %'color',[0 0 0]
htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',13);
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ; %eliminates the axes

%%
%plot evoked traces with ES:
trace_to_plot = plot_data(:,trace_ind,3)+DC ; %DC is added just to make space between the traces 
trace_to_plot(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
        f2=figure;
        hold on
htrace=plot(x_axis*dt,trace_to_plot(x_axis,:,:), 'LineWidth',1.2,'color', color_table(2,:));
for i=1:length(trace_ind);
    text(x_axis(1)*dt,trace_to_plot(x_axis(1),i),[num2str(floor(plot_data(x_axis(1),i,3))), ' mV '],'HorizontalAlignment','right','fontsize',13,'fontname','arial')
end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');

% plotting zoom-in area:
rectangle('Position',[x_axis_short(1)*dt,ylim_data(1),(x_axis_short(end)-x_axis_short(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none')

%plotting stim_2:
if plot_stim_2(3);
patch_xdata=[stim2_X{3}; flipud(stim2_X{3})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
%plotting scale bar
yline_start=ylim_data(1)-6; yline_end=yline_start+10;
xline_start=x_axis(1)*dt-1.5; xline_end=xline_start+1;
if isequal(x_axis,x_axis_short)
    xline_start=x_axis(1)*dt-0.5; xline_end=xline_start+0.5;
end
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',13); %'color',[0 0 0]
htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',13);
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ;

 %%
 %plot evoked mean trace
        f3=figure;
        hold on
htrace1=plot(x_axis*dt,plot_data_mean(x_axis,:,2), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis*dt,plot_data_mean(x_axis,:,3), 'LineWidth',1.2,'color', color_table(2,:));
text(x_axis(1)*dt,plot_data_mean(x_axis(1),1,2),[num2str(floor(plot_data_mean(x_axis(1),1,2))), ' mV '],'HorizontalAlignment','right','fontsize',13,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
rectangle('Position',[x_axis_short(1)*dt,ylim_data(1),(x_axis_short(end)-x_axis_short(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none')

%for plotting stim_2:
if plot_stim_2(3);
patch_xdata=[stim2_X{2}; flipud(stim2_X{2})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
%plotting scale bar
yline_start=ylim_data(1)-2; yline_end=yline_start+5;
xline_start=x_axis(1)*dt-0.8; xline_end=xline_start+1;
if isequal(x_axis,x_axis_short)
    xline_start=x_axis(1)*dt-0.5; xline_end=xline_start+0.5;
end
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',13); %'color',[0 0 0]
htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',13);
hold off
ylabel('Mean Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ;

%%
%plot evoked SD trace
        f4=figure;
        hold on
htrace1=plot(x_axis*dt,plot_data_std(x_axis,:,2), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis*dt,plot_data_std(x_axis,:,3), 'LineWidth',1.2,'color', color_table(2,:));
text(x_axis(1)*dt,plot_data_std(x_axis(1),1,2),[num2str(floor(plot_data_std(x_axis(1),1,2))), ' mV '],'HorizontalAlignment','right','fontsize',13,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
rectangle('Position',[x_axis_short(1)*dt,ylim_data(1),(x_axis_short(end)-x_axis_short(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none')

%for plotting stim_2:
if plot_stim_2(3);
patch_xdata=[stim2_X{2}; flipud(stim2_X{2})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
%plotting scale bar
yline_start=ylim_data(1)-0.2;  yline_end=yline_start+2;
xline_start=x_axis(1)*dt-0.8; xline_end=xline_start+1;
if isequal(x_axis,x_axis_short)
    xline_start=x_axis(1)*dt-0.5; xline_end=xline_start+0.5;
end
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',13); %'color',[0 0 0]
htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',13);
hold off
ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ;

%%
if save_flag==1;
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\Long Trace'  
    if isequal(x_axis,x_axis_short)
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\Zoom-in Trace'  
    end
        saveas(f1,['f' num2str(files_to_analyze(fileind)) '_traces_x2.fig']) 
        print(f1,['f' num2str(files_to_analyze(fileind)) '_traces_x2'],'-dpng','-r600','-opengl') 
        saveas(f2,['f' num2str(files_to_analyze(fileind)) '_traces_x3.fig']) 
        print(f2,['f' num2str(files_to_analyze(fileind)) '_traces_x3'],'-dpng','-r600','-opengl') 
        saveas(f3,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3.fig']) 
        print(f3,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3'],'-dpng','-r600','-opengl') 
        saveas(f4,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt.fig']) 
        print(f4,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt'],'-dpng','-r600','-opengl') 
end