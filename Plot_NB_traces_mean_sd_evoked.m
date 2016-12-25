%% Load file
close all
clear all

 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data data_current
 global exp_type
exp_type=2; %1-NBES, 2-ChAT
save_flag= 0;
print_flag=0;
short_flag=1; %1- short trace, 0- long trace
LPF_flag=0;
baseline_flag=0; %adding dashed line under traces
bl=-70; %baseline value;

switch exp_type
    case 1
        files_to_analyze =46; %[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,44,46,48,50,52,56,58,62,72,75]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};

    case 2
        files_to_analyze =[76];
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light On', 'Light Off'};
end

%% 
    for fileind=1:length(files_to_analyze) ;
    close all
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
%     clear data_no_spike_no_DC
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
    data_current=data_no_spikes; %raw_data
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
if LPF_flag==1;
    for xx=1:3
        for trace= 1:size(data_current{channel},2)    
                jj=data_current{channel}(:,trace,xx);
                data_current{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jj,dt,-1,300,0,0); 
        end
    end
end
% remove outlyers:
if exp_type==1 && files_to_analyze(fileind)==16;
    data_current{channel}(:,4,:)=[];
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
               
%                 data_no_spike_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(data_no_spikes{channel}(:,:,x_value),meansubtract_interval);
                data_current_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(data_current{channel}(:,:,x_value),meansubtract_interval);    
end
    end

%% Plotting traces+mean+STD
trace_fontsize=12;
scalebar_fontsize=12;
plot_stim_1=[1,0,1];
plot_stim_2=[0,1,1];
trace_ind = [2,5,6]; %1:size(plot_data,2);  %trace_ind is the numbers of traces to be plotted [2,5,6];
x_axis=[];

clear color_table
    
 switch exp_type
    case 1
%         x_axis_long=0.5*sf{1}:12*sf{1}-1;
%         x_axis_short(1,:)=0.4*sf{1}:2.9*sf{1};
%         x_axis_short(2,:)=stim2_X{2}(1,1)-0.5*sf{1}:stim2_X{2}(1,1)+2*sf{1};
        x_axis_long=(stim1_X{1}(2,1)+10):10*sf{1}-1;
        x_axis_short=stim2_X{2}(1,1)-0.5*sf{1}:stim2_X{2}(1,1)+2*sf{1};
        color_table=[0 0 0; [216 22 22]/256; [136 137 138]/256; [255 153 153]/256];
        rectangle_color=[239 239 239]/256;
    case 2
        x_axis_long=1:12*sf{1};
        x_axis_short=stim2_X{2}(1,1)-0.05*sf{1}:stim2_X{2}(1,1)+2.45*sf{1};
        color_table=[0 0 0; [216 22 22]/256; [204 229 255]/256]; %[0 0 0; [216 22 22]/256; [136 137 138]/256]; [255 153 153]/256]; 
        rectangle_color= [239 239 239]/256;
 end
 
if short_flag
     x_axis= x_axis_short; 
else
        x_axis= x_axis_long; 
end

% segment1=x_axis_long(1):stim1_X{1}(1,1);
% segment2=stim1_X{1}(1,1)+1:x_axis_long(end);

plot_data=data_current{channel}; %data_current %data_current_no_DC 
plot_data_mean = mean(plot_data,2);  
% plot_data_mean(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_std =  std(data_current_no_DC{channel},0,2); 
% plot_data_std(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_var=var(data_current_no_DC{channel},0,2); % plot_data_var=var(plot_data,0,2);
plot_data_ff = (-1).*plot_data_var./plot_data_mean; 

       l=size(plot_data(:,trace_ind, 2),1);
        l=l/2-1;
        DC=20.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);        
              
%         x2lim=[stim2_X{x_value(1)}(1,1).*dt-0.5,stim2_X{2}(1,1).*dt+2];
 
%%
    %plot evoked traces no ES:
    trace_to_plot = plot_data(:,trace_ind,2)+DC ; %DC is added just to make space between the traces 
    y_bl=DC(1:2,:)'+bl;
    x_bl=[x_axis(1).*ones(size(y_bl,1),1), x_axis(end).*ones(size(y_bl,1),1)].*dt;  
  
        f1=figure;
        hold on
        rec1=rectangle('position',[x_axis(1)*dt ,min(trace_to_plot(x_axis,2,:)),0.1,0.1]);
htrace=plot(x_axis*dt,trace_to_plot(x_axis,:,:), 'LineWidth',1.2,'color', color_table(1,:));
if baseline_flag==1;
    h_baseline=line(x_bl',y_bl','linestyle','--','color',[136 137 138]./256,'linewidth',1);
end
% for i=1:length(trace_ind);
%     text(x_axis(1)*dt,trace_to_plot(x_axis(1),i),[num2str(floor(plot_data(x_axis(1),i,2))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
set(rec1,'Position',[x_axis_short(1)*dt,ylim_data(1),(x_axis_short(end)-x_axis_short(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
%plot light stimulus
% switch exp_type
%     case 1
%     case 2
%         stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
%         line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
% end
% plotting stim_2:
if plot_stim_2(2);
patch_xdata=[stim2_X{2}; flipud(stim2_X{2})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize;
        [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh), ' mV'];     c=[0,0,0];  fonsizes=scalebar_fontsize;
        [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);        
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca, 'visible', 'off') ;

%%
%plot evoked traces with ES:
trace_to_plot = plot_data(:,trace_ind,3)+DC ; %DC is added just to make space between the traces 
% trace_to_plot(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
        f2=figure;
        hold on
        rec1=rectangle('position',[x_axis(1)*dt ,min(trace_to_plot(x_axis,2,:)),0.1,0.1]);
htrace=plot(x_axis*dt,trace_to_plot(x_axis,:,:), 'LineWidth',1.2,'color', color_table(2,:));
if baseline_flag==1;
    h_baseline=line(x_bl',y_bl','linestyle','--','color',[136 137 138]./256,'linewidth',1);
end

% for i=1:length(trace_ind);
%     text(x_axis(1)*dt,trace_to_plot(x_axis(1),i),[num2str(floor(plot_data(x_axis(1),i,3))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');

% plotting zoom-in area:
set(rec1,'Position',[x_axis_short(1)*dt,ylim_data(1),(x_axis_short(end)-x_axis_short(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plot light stimulus
switch exp_type
    case 1
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
end

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
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize;
        [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh), ' mV'];     c=[0,0,0];  fonsizes=scalebar_fontsize;
        [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
        
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca, 'visible', 'off') ;

 %% plot evoked mean trace
        f3=figure;
        hold on
        rec1=rectangle('position',[x_axis(1)*dt ,min(plot_data_mean(x_axis,:,2)),0.1,0.1]);
htrace1=plot(x_axis*dt,plot_data_mean(x_axis,:,2), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis*dt,plot_data_mean(x_axis,:,3), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis(1)*dt,plot_data_mean(x_axis(1),1,2),[num2str(floor(plot_data_mean(x_axis(1),1,2))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
set(rec1,'Position',[x_axis_short(1)*dt,ylim_data(1),(x_axis_short(end)-x_axis_short(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plot light stimulus
switch exp_type
    case 1
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
end

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
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize;
        [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh), ' mV'];     c=[0,0,0];  fonsizes=scalebar_fontsize;
        [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
        
hold off
ylabel('Mean Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
% l=legend({'NB-','NB+'},'box','off','location','northeast');
set(gca, 'visible', 'off') ;

%%
%plot evoked SD trace
        f4=figure;
        hold on
        rec1=rectangle('position',[x_axis(1)*dt ,min(plot_data_std(x_axis,:,2)),0.1,0.1]);
htrace1=plot(x_axis*dt,plot_data_std(x_axis,:,2), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis*dt,plot_data_std(x_axis,:,3), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis(1)*dt,plot_data_std(x_axis(1),1,2),[num2str(floor(plot_data_std(x_axis(1),1,2))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
set(rec1,'Position',[x_axis_short(1)*dt,ylim_data(1),(x_axis_short(end)-x_axis_short(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plot light stimulus
switch exp_type
    case 1
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
end

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
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize;
        [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh), ' mV'];     c=[0,0,0];  fonsizes=scalebar_fontsize;
        [p1,p2] = fn_makeCalibBar(horiz_vert,lengthh,textit,c,fonsizes);hold off
ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca, 'visible', 'off') ;

%%
if save_flag==1;
    switch exp_type
        case 1
            cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\Long Trace Presentation'
             if isequal(x_axis,x_axis_short)
                cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\Zoom-in Trace Presentation'  
            end
        case 2
            cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary\Long Trace Presentation'
            if isequal(x_axis,x_axis_short)
                cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary\Zoom-in Trace Presentation'  
            end
    end
    
        saveas(f1,['f' num2str(files_to_analyze(fileind)) '_traces_x2_v2.fig']) 
        print(f1,['f' num2str(files_to_analyze(fileind)) '_traces_x2_v2'],'-dpng','-r600','-opengl') 
        saveas(f2,['f' num2str(files_to_analyze(fileind)) '_traces_x3_v2.fig']) 
        print(f2,['f' num2str(files_to_analyze(fileind)) '_traces_x3_v2'],'-dpng','-r600','-opengl') 
        saveas(f3,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_v2.fig']) 
        print(f3,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_v2'],'-dpng','-r600','-opengl') 
        saveas(f4,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt_v2.fig']) 
        print(f4,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt_v2'],'-dpng','-r600','-opengl') 
end