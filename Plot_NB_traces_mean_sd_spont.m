%% Load file
close all
clear all
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param 
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
save_flag= 1;
short_flag=0; %1- short trace, 0- long trace
files_to_analyze =16; %[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,44,46,48,50,52,56,58,62,72,75]; 
% for fileind=1;
    for fileind=1:length(files_to_analyze) ;
    close all
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    clear data_no_spike_no_DC
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
    color_table=[0 0 0;color_table(1:6,:)];
 %   
                sf{1} = Param.sf_Vm;
                sf{2} = Param.sf_I1;
                sf{3} = Param.sf_V2;
                sf{4} = Param.sf_I2;
                dt=1/sf{channel};
                             
                sf_airpuff = Param.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                sf_galvano = Param.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;  
% low-pass filtering below 300Hz                
                for xx=2:3
    for trace= 1:size(data_no_spikes{channel},2)    
            jj=data_no_spikes{channel}(:,trace,xx);
            data_no_spikes{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jj,dt,-1,300,0,0); 
    end
                end
% remove outlyers:
if files_to_analyze(fileind)==16;
    data_no_spikes{channel}(:,4,:)=[];
end
                % Subtract mean trace from data with or without spikes
                meansubtract_start_time = 0; %[sec]
                meansubtract_duration = 3; %[sec]
for x_value = 1:size(data.x_value,2) 
     
                 meansubtract_start_sample = meansubtract_start_time.*sf{channel};
                if meansubtract_start_time==0
                    meansubtract_start_sample = 1;
                end
                meansubtract_end_sample = meansubtract_start_sample+meansubtract_duration.*sf{channel}-1;
                meansubtract_interval = round(meansubtract_start_sample:meansubtract_end_sample);
               
                data_no_spike_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(data_no_spikes{channel}(:,:,x_value),meansubtract_interval);
%                 data_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(raw_data{channel}(:,:,x_value),meansubtract_interval);    
end
    end

%% Plotting traces+mean+STD
trace_fontsize=12;
scalebar_fontsize=12;
plot_stim_1=[1,0,1];
plot_stim_2=[0,1,1];
trace_ind = [1,3,4]; %1:size(plot_data,2);  %trace_ind is the numbers of traces to be plotted
x_axis=[];
x_axis_long=0.5*sf{1}:10*sf{1}-1;
x_axis_short(1,:)=0.4*sf{1}:2.9*sf{1};
x_axis_short(2,:)=stim2_X{2}(1,1)-0.5*sf{1}:stim2_X{2}(1,1)+2*sf{1};
% x_axis_short(2,:)=5.6*sf{1}:8.1*sf{1};

clear color_table
    color_table=[0 0 0; [216 22 22]/256; [136 137 138]/256; [255 153 153]/256];
    rectangle_color=[239 239 239]/256;
    
plot_data=data_no_spikes{channel}; %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
plot_data_mean = mean(plot_data,2);  plot_data_mean(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_std =  std(data_no_spike_no_DC{channel},0,2); plot_data_std(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_var=var(data_no_spike_no_DC{channel},0,2); % plot_data_var=var(plot_data,0,2);
plot_data_ff = (-1).*plot_data_var./plot_data_mean; 
       l=size(plot_data(:,trace_ind, 2),1);
        l=l/2-1;
        DC=20.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);        
              
%         x2lim=[stim2_X{x_value(1)}(1,1).*dt-0.5,stim2_X{2}(1,1).*dt+2];
 
%%
    %plot spont traces long:
    trace_to_plot = plot_data(:,trace_ind,1)+DC ; %DC is added just to make space between the traces 
        f1=figure;
        hold on
        rec1=rectangle;
        rec2=rectangle;
htrace1=plot(x_axis_long(1:2.5*sf{1})*dt,trace_to_plot(x_axis_long(1:2.5*sf{1}),:,:), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_long(3*sf{1}:end)*dt,trace_to_plot(x_axis_long(3*sf{1}:end),:,:), 'LineWidth',1.2,'color', color_table(2,:));

% for i=1:length(trace_ind);
%     text(x_axis_long(1)*dt,trace_to_plot(x_axis_long(1),i),[num2str(floor(plot_data(x_axis_long(1),i,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_short(1,1)*dt,ylim_data(1),(x_axis_short(1,end)-x_axis_short(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
set(rec2,'Position',[x_axis_short(2,1)*dt,ylim_data(1),(x_axis_short(2,end)-x_axis_short(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
rec3=rectangle('Position',[stim1_X{1}(1)*dt,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');

%plotting scale bar
yline_start=ylim_data(1)-7; yline_end=yline_start+10;
xline_start=x_axis_short(1,1)*dt-0.5; xline_end=xline_start+1;
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
% htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',10);
htext_v=text(xline(1,1),yline(2,1),stringv,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
hold off
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
% ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca, 'visible', 'off') ;
%%
    %plot spont traces short - no ES:
    trace_to_plot = plot_data(:,trace_ind,1)+DC ; %DC is added just to make space between the traces 
        f4=figure;
        hold on
         rec1=rectangle('Position',[x_axis_short(2,1)*dt,trace_to_plot(x_axis_short(2,1)),0.1,0.1]);
htrace1=plot(x_axis_short(1,:)*dt,trace_to_plot(x_axis_short(1,:),:,:), 'LineWidth',1.2,'color', color_table(1,:));
% htrace2=plot(x_axis_long(3*sf{1}:end)*dt,trace_to_plot(x_axis_long(3*sf{1}:end),:,:), 'LineWidth',1.2,'color', color_table(2,:));

% for i=1:length(trace_ind);
%     text(x_axis_short(1,1)*dt,trace_to_plot(x_axis_short(1,1),i),[num2str(floor(plot_data(x_axis_short(1,1),i,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_short(1,1)*dt,ylim_data(1),(x_axis_short(1,end)-x_axis_short(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
% set(rec2,'Position',[x_axis_short(2,1)*dt,ylim_data(1),(x_axis_short(2,end)-x_axis_short(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plotting scale bar
yline_start=ylim_data(1)-6; yline_end=yline_start+10;
xline_start=x_axis_short(1,1)*dt-0.2; xline_end=xline_start+0.5;
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ;

%%
    %plot spont traces short - with ES:
    trace_to_plot = plot_data(:,trace_ind,1)+DC ; %DC is added just to make space between the traces 
        f5=figure;
        hold on
        rec1=rectangle('Position',[x_axis_short(2,1)*dt,trace_to_plot(x_axis_short(2,1)),0.1,0.1]);
htrace1=plot(x_axis_short(2,:)*dt,trace_to_plot(x_axis_short(2,:),:,:), 'LineWidth',1.2,'color', color_table(2,:));

% for i=1:length(trace_ind);
%     text(x_axis_short(2,1)*dt,trace_to_plot(x_axis_short(2,1),i),[num2str(floor(plot_data(x_axis_short(2,1),i,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
% ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_short(2,1)*dt,ylim_data(1),(x_axis_short(2,end)-x_axis_short(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plotting scale bar
yline_start=ylim_data(1)-6; yline_end=yline_start+10;
xline_start=x_axis_short(2,1)*dt-0.2; xline_end=xline_start+0.5;
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ;

%%
    %plot mean trace short:   
        f6=figure;
        hold on
        rec1=rectangle('Position',[x_axis_short(2,1)*dt,plot_data_mean(x_axis_short(1,1),:,1),0.1,0.1]);
htrace1=plot(x_axis_short(1,:)*dt,plot_data_mean(x_axis_short(1,:),:,1), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_short(1,:)*dt,plot_data_mean(x_axis_short(2,:),:,1), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis_short(1,1)*dt,plot_data_mean(x_axis_short(1,1),1),[num2str(floor(plot_data_mean(x_axis_short(1,1),1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_short(1,1)*dt,ylim_data(1),(x_axis_short(1,end)-x_axis_short(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plotting scale bar
yline_start=ylim_data(1)-2; yline_end=yline_start+2;
xline_start=x_axis_short(1,1)*dt-0.2; xline_end=xline_start+0.5;
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ;

%%
    %plot std trace short:   
        f7=figure;
        hold on
         rec1=rectangle('Position',[x_axis_short(2,1)*dt,plot_data_std(x_axis_short(1,1),:,1),0.1,0.1]);
htrace1=plot(x_axis_short(1,:)*dt,plot_data_std(x_axis_short(1,:),:,1), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_short(1,:)*dt,plot_data_std(x_axis_short(2,:),:,1), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis_short(1,1)*dt,plot_data_std(x_axis_short(1,1),1),[num2str(floor(plot_data_std(x_axis_short(1,1),1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_short(1,1)*dt,ylim_data(1),(x_axis_short(1,end)-x_axis_short(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plotting scale bar
yline_start=ylim_data(1)-1; yline_end=yline_start+2;
xline_start=x_axis_short(1,1)*dt-0.2; xline_end=xline_start+0.5;
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
hold off
ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ;

%%
%plot mean trace
        f2=figure;
        hold on
       rec1=rectangle('Position',[x_axis_long(1,1)*dt,plot_data_mean(x_axis_long(1,1),:,1),0.1,0.1]);
         rec2=rectangle('Position',[x_axis_long(1,1)*dt,plot_data_mean(x_axis_long(1,1),:,1),0.1,0.1]);
htrace1=plot(x_axis_long(1:2.5*sf{1})*dt,plot_data_mean(x_axis_long(1:2.5*sf{1}),:,1), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_long(3*sf{1}:end)*dt,plot_data_mean(x_axis_long(3*sf{1}:end),:,1), 'LineWidth',1.2,'color', color_table(2,:));
% htrace1=plot(x_axis_long*dt,plot_data_mean(x_axis_long,:,1), 'LineWidth',1.2,'color', color_table(1,:));
% text(x_axis_long(1)*dt,plot_data_mean(x_axis_long(1),1,1),[num2str(floor(plot_data_mean(x_axis_long(1),1,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_short(1,1)*dt,ylim_data(1),(x_axis_short(1,end)-x_axis_short(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
set(rec2,'Position',[x_axis_short(2,1)*dt,ylim_data(1),(x_axis_short(2,end)-x_axis_short(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
rec3=rectangle('Position',[stim1_X{1}(1)*dt,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
% rec1=rectangle('Position',[x_axis_short(1,1)*dt,ylim_data(1),(x_axis_short(1,end)-x_axis_short(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
% %plotting scale bar
yline_start=ylim_data(1)-1; yline_end=yline_start+2;
xline_start=x_axis_long(1)*dt-0.5; xline_end=xline_start+1;
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
% htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
htext_v=text(xline(1,1),yline(2,1),stringv,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ;

%%
%plot std trace
        f3=figure;
        hold on
        rec1=rectangle('Position',[x_axis_long(1,1)*dt,plot_data_std(x_axis_long(1,1),:,1),0.1,0.1]);
         rec2=rectangle('Position',[x_axis_long(1,1)*dt,plot_data_std(x_axis_long(1,1),:,1),0.1,0.1]);
htrace1=plot(x_axis_long(1:2.5*sf{1})*dt,plot_data_std(x_axis_long(1:2.5*sf{1}),:,1), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_long(3*sf{1}:end)*dt,plot_data_std(x_axis_long(3*sf{1}:end),:,1), 'LineWidth',1.2,'color', color_table(2,:));
% htrace1=plot(x_axis_long*dt,plot_data_std(x_axis_long,:,1), 'LineWidth',1.2,'color', color_table(1,:));
% text(x_axis_long(1)*dt,plot_data_std(x_axis_long(1),1,1),[num2str(floor(plot_data_std(x_axis_long(1),1,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_short(1,1)*dt,ylim_data(1),(x_axis_short(1,end)-x_axis_short(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
set(rec2,'Position',[x_axis_short(2,1)*dt,ylim_data(1),(x_axis_short(2,end)-x_axis_short(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
rec3=rectangle('Position',[stim1_X{1}(1)*dt,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
% rec1=rectangle('Position',[x_axis_short(1,1)*dt,ylim_data(1),(x_axis_short(1,end)-x_axis_short(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
% %plotting scale bar
yline_start=ylim_data(1)-1; yline_end=yline_start+2;
xline_start=x_axis_long(1)*dt-0.5; xline_end=xline_start+1;
xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
stringh=[num2str(xline(2,1)-xline(1,1)), ' Sec'];
stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
% htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',10);
htext_v=text(xline(1,1),yline(2,1),stringv,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
hold off
ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
set(gca, 'visible', 'off') ;
%%
if save_flag==1;
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\Long Trace Presentation'              
        saveas(f1,['f' num2str(files_to_analyze(fileind)) '_traces_x1_v2.fig']) 
        print(f1,['f' num2str(files_to_analyze(fileind)) '_traces_x1_v2'],'-dpng','-r600','-opengl') 
        saveas(f2,['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2.fig']) 
        print(f2,['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2'],'-dpng','-r600','-opengl') 
        saveas(f3,['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_v2.fig']) 
        print(f3,['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_v2'],'-dpng','-r600','-opengl') 
        
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\Zoom-in Trace Presentation'  
       saveas(f4,['f' num2str(files_to_analyze(fileind)) '_traces_x1_noES_v2.fig']) 
        print(f4,['f' num2str(files_to_analyze(fileind)) '_traces_x1_noES_v2'],'-dpng','-r600','-opengl') 
        saveas(f5,['f' num2str(files_to_analyze(fileind)) '_traces_x1_ES_v2.fig']) 
        print(f5,['f' num2str(files_to_analyze(fileind)) '_traces_x1_ES_v2'],'-dpng','-r600','-opengl') 
        saveas(f6,['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2.fig']) 
        print(f6,['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2'],'-dpng','-r600','-opengl') 
        saveas(f7,['f' num2str(files_to_analyze(fileind)) '_std_1_mean-subt_v2.fig']) 
        print(f7,['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_v2'],'-dpng','-r600','-opengl') 
end