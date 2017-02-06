%% Load file
close all
clear all

 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data
 global exp_type
exp_type=2; %1-NBES, 2-ChAT
save_flag= 1;
print_flag=0;
clamp_flag=[]; %[]; %3; %clamp_flag=1 for hyperpolarization traces, clamp_flag=2 for depolarization traces and clamp_flag=3 for no current traces (only clamp to resting Vm)
% short_flag=0; %1- short trace, 0- long trace
baseline_flag=0; %adding dashed line under traces
bl=-70; %baseline value;
norm_flag=0;
BP50HzLFP_flag=1; %removing 50Hz noise from LFP signal
BP50HzVm_flag=1; %removing 50Hz noise from Vm signal
BPLFP_flag=1; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=[1,200]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=1; %filtering LFP and Vm same as LFP was filtered in the multiclamp
bp_manual_Vm=[0,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)

 %% set the path for saving figures and variables
% if BP50HzLFP_flag==1 && BP50HzVm_flag==1 && BPVm_flag==1 && BPLFP_flag==1
%     path_output=['LFP_50Hz+BP Vm_ 50Hz+BP\BP',num2str(bp_manual_Vm(1)),'-',num2str(bp_manual_Vm(2))];
% else if BP50HzLFP_flag==1 && BP50HzVm_flag==1 && BPLFP_flag==1
%      path_output='LFP_50Hz+BP Vm_ 50Hz';  
%     else if BP50HzLFP_flag==1 && BP50HzVm_flag==1 && BPVm_flag==1
%             path_output='LFP_50Hz Vm_50Hz+BP';  
%          else if BP50HzLFP_flag==1 && BP50HzVm_flag==1
%                  path_output='LFP_50Hz Vm_50Hz';  
%              else if BP50HzLFP_flag==1 
%                  path_output='LFP_50Hz';  
%                  else path_output='No offline filter';
%                  end
%              end
%         end
%     end
% end
% switch exp_type
%     case 1 
%         path_output=['D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\',path_output];   
%         a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
%         if a==0;
%             mkdir(path_output);
%         end
%         cd(path_output)       
%     case 2
%         path_output=['D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary\',path_output];        
%         a = exist(path_output,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
%         if a==0;
%             mkdir(path_output);
%         end
%         cd(path_output)
% end
%%
switch exp_type
    case 1
        files_to_analyze =84; %[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75,82,84]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};   y_ax_label={'Vm'}; y_ax_units={'mV'};   
         path_output1='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Long Trace Presentation';  
        path_output2='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Zoom-in Trace Presentation';  
        a = exist(path_output1,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        b = exist(path_output2,'dir'); %b will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output1);
        end 
       if b==0;
            mkdir(path_output2);
      end 

    case 2
        files_to_analyze =84; %[74,76,77,80,82,84,87];
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light On', 'Light Off'};    y_ax_label={'Vm'}; y_ax_units={'mV'}; 
         path_output1='D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Long Trace Presentation';  
        path_output2='D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Zoom-in Trace Presentation';  
        a = exist(path_output1,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        b = exist(path_output2,'dir'); %b will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output1);
        end 
       if b==0;
            mkdir(path_output2);
      end 
   case 3 
        files_to_analyze =74; %[31,38,42,51,61,64,67,69,71,74,77];
        clamp_flag=3; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};    y_ax_label={'Im'}; y_ax_units={'pA'};    
        path_output1='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean_VC\Long Trace Presentation';  
        path_output2='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean_VC\Zoom-in Trace Presentation';  
        a = exist(path_output1,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        b = exist(path_output2,'dir'); %b will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output1);
        end 
       if b==0;
            mkdir(path_output2);
      end 
end
%% 
    for fileind=1:length(files_to_analyze);
    close all
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
    %%
    Ch2_data= raw_data{3}./20; %dividing by the LFP gain       
      if isempty(data_no_spikes)
            current_data=raw_data{channel};
             if clamp_flag==1
                current_data=(-1).*raw_data{channel};
            end
            data_used='raw_data';
        else
        current_data=data_no_spikes{channel}; %raw_data{channel}; 
         if clamp_flag==1
            current_data=(-1).*data_no_spikes{channel}; %raw_data{channel};
        end
        data_used='data_no_spikes';
      end
        
   data_preprocessing 
   
    if isempty(current_data_filt)
     current_data_filt=current_data;
    end
    
        clear color_table    
        whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;
        switch exp_type
            case 1 
                color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
            case 2
                color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
            case 3 
                color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];       
        end
        
% remove outlyers:
if exp_type==1 && files_to_analyze(fileind)==16;
    current_data_filt(:,4,:)=[];
end
if exp_type==3 && files_to_analyze(fileind)==31;
    current_data_filt(:,2,:)=[];
end
                %% Subtract mean trace from data with or without spikes
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
                current_data_no_DC(:,:,x_value) = fn_Subtract_Mean(current_data_filt(:,:,x_value),meansubtract_interval);    
end
    end

%% Plotting traces+mean+STD
trace_fontsize=12;
scalebar_fontsize=12;
plot_stim_1=[1,0,1];
plot_stim_2=[0,1,1];
trace_ind =[1:size(current_data_filt,2)]; % [1,2,3,4,5,6]; %1:size(plot_data,2);  %trace_ind is the numbers of traces to be plotted
switch exp_type
    case 1
        stim2{1}=stim2_X{2};
        stim2{2}=stim2_X{2};
        x_value=1:3;   
        galvano_nstim = Param.facade(6);
        galvano_freq = Param.facade(7);
         if size(current_data_filt,1)>12*sf{1}
            x_axis_Slong=0.5*sf{1}:12*sf{1}-1;
            x_axis_Elong=round((stim1_X{1}(2,1)+10):10*sf{1}-1);
        else
            x_axis_Slong=0.5*sf{1}:8*sf{1}-1;
            x_axis_Elong=round((stim1_X{1}(2,1)+10):10*sf{1}-1);
        end
        x_axis_Sshort(1,:)=0.4*sf{1}:2.9*sf{1};
        x_axis_Sshort(2,:)=round(stim2{1}(1,1)-0.5*sf{1}:stim2{1}(1,1)+2*sf{1});
        
        x_axis_Eshort(1,:)=round(stim2{1}(1,1)-0.5*sf{1}:stim2{1}(1,1)+2*sf{1});
        x_axis_Eshort(2,:)=round(stim2{1}(1,1)-0.5*sf{1}:stim2{1}(1,1)+2*sf{1});
        rectangle_color=[239 239 239]/256;
        segment1=x_axis_Slong(1):stim1_X{1}(1,1);
        segment2=stim1_X{1}(1,1)+1:x_axis_Slong(end);

    case 2
        stim2{1}=stim2_X{2};
        stim2{2}=stim2_X{2};
        x_value=1:3;   
        galvano_nstim = Param.facade(6);
        galvano_freq = Param.facade(7);
        curr_inj_1=Param.facade(22); %in sec
        curr_inj_dur=Param.facade(20); %in sec
        curr_inj_2=Param.facade(22)+curr_inj_dur+Param.facade(23); %in sec        
        curr_inj_1_end=curr_inj_1+curr_inj_dur+0.15;
        curr_inj_2_end=curr_inj_2+curr_inj_dur+0.15;
        x_axis_Slong=round(1*sf{1}:12*sf{1}-1);
        x_axis_Sshort(1,:)=round(0.4*sf{1}:2.9*sf{1});
        x_axis_Sshort(2,:)=round(stim2{1}(1,1)-0.05*sf{1}:stim2{1}(1,1)+2.45*sf{1});
        x_axis_Elong=round(1*sf{1}:12*sf{1}-1);
        x_axis_Eshort(1,:)=round(stim2{1}(1,1)-0.05*sf{1}:stim2{1}(1,1)+2.45*sf{1});
        x_axis_Eshort(2,:)=round(stim2{1}(1,1)-0.05*sf{1}:stim2{1}(1,1)+2.45*sf{1});
        rectangle_color= [239 239 239]/256;
        segment1=x_axis_Slong(1):stim1_X{1}(1,1);
        segment2=stim1_X{1}(1,1)+1:x_axis_Slong(end);
    case 3
        x_value=[clamp_flag,clamp_flag+3,clamp_flag+3]; %x_value(1) is spontaneous and x_value(2) is evoked
        stim2{1}=stim2_X{4}(1:2,:);
        stim2{2}=stim2_X{4}(3:4,:);      
        galvano_nstim = Param.facade(15);
        galvano_freq = Param.facade(16);
        curr_inj_del = [Param.facade(27) Param.facade(30)] ;
            if length(Param.facade)>32
                curr_inj_del = [Param.facade(27) Param.facade(30) Param.facade(33)];
            end
        curr_inj_depo = Param.facade(1);
        curr_inj_hyper = Param.facade(2);
        curr_inj_dur = Param.facade(7);
        if size(current_data_filt,1)>12*sf{1}
            x_axis_Slong=0.5*sf{1}:12*sf{1}-1;
            x_axis_Elong=1:12*sf{1}-1;
        else
            x_axis_Slong=0.5*sf{1}:8*sf{1}-1;
            x_axis_Elong=1:8*sf{1}-1;
        end
        x_axis_Sshort(1,:)=round((curr_inj_del(1)+0.1)*sf{1}:(curr_inj_del(1)+curr_inj_dur-0.1)*sf{1});
        x_axis_Sshort(2,:)=round((curr_inj_del(2)+0.1)*sf{1}:(curr_inj_del(2)+curr_inj_dur-0.1)*sf{1});
        
        x_axis_Eshort(1,:)=round((curr_inj_del(1)+0.1)*sf{1}:(curr_inj_del(1)+curr_inj_dur-0.1)*sf{1});
        x_axis_Eshort(2,:)=round((curr_inj_del(2)+0.1)*sf{1}:(curr_inj_del(2)+curr_inj_dur-0.1)*sf{1});
        rectangle_color=[239 239 239]/256;
        segment1=x_axis_Slong(1):stim1_X{1}(1,1);
        segment2=stim1_X{1}(1,1)+1:x_axis_Slong(end);
end
 

plot_data=current_data_filt(:,:,x_value); %current_data_filt %current_data_no_DC 
% if exp_type==2;
%     plot_data(curr_inj_1*sf{1}:curr_inj_1_end*sf{1},:,:)=nan;
%     plot_data(curr_inj_2*sf{1}:curr_inj_2_end*sf{1},:,:)=nan;
% end
plot_data_mean = mean(plot_data,2);  
% plot_data_mean(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_std =  std(current_data_no_DC(:,:,x_value),0,2); 
% plot_data_std(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_var=var(current_data_no_DC(:,:,x_value),0,2); % plot_data_var=var(plot_data,0,2);
% plot_data_ff = (-1).*plot_data_var./plot_data_mean; 
       l=size(plot_data(:,trace_ind, 2),1);
        l=l/2-1;
        DC=20.*(0:length(trace_ind)-1);
        if exp_type==3
            DC=40.*(0:length(trace_ind)-1);
        end
        DC=wextend('addrow', 'per', DC, l);        
              
%         x2lim=[stim2_X{x_value(1)}(1,1).*dt-0.5,stim2_X{2}(1,1).*dt+2];
  
   
%% spontaneous figures    
            %plot spont traces long:
    trace_to_plot = plot_data(:,trace_ind,1)+DC ; %DC is added just to make space between the traces 
        s1=figure;
        hold on
        rec1=rectangle('Position',[segment1(1,1)*dt,trace_to_plot(segment1(1,1)),0.1,0.1]);
        rec2=rectangle('Position',[segment1(1,1)*dt,trace_to_plot(segment1(1,1)),0.1,0.1]);
htrace1=plot(segment1*dt,trace_to_plot(segment1,:,:), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(segment2*dt,trace_to_plot(segment2,:,:), 'LineWidth',1.2,'color', color_table(2,:));
%     text(x_axis_Slong(1)*dt,trace_to_plot(x_axis_Slong(1),1),[num2str(floor(plot_data(x_axis_Slong(1),1,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% for i=1:length(trace_ind); %adding the voltage value of the first data point to the left of the trace
%     text(x_axis_Slong(1)*dt,trace_to_plot(x_axis_Slong(1),i),[num2str(floor(plot_data(x_axis_Slong(1),i,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_Sshort(1,1)*dt,ylim_data(1),(x_axis_Sshort(1,end)-x_axis_Sshort(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
set(rec2,'Position',[x_axis_Sshort(2,1)*dt,ylim_data(1),(x_axis_Sshort(2,end)-x_axis_Sshort(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
switch exp_type
    case 1
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-500*dt,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
    case 3
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-500*dt,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

        %plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);

set(gca, 'visible', 'off') ;
%%
%plot mean trace
        s2=figure;
        hold on
       rec1=rectangle('Position',[segment1(1,1)*dt,plot_data_mean(segment1(1,1)),0.1,0.1]);
         rec2=rectangle('Position',[segment1(1,1)*dt,plot_data_mean(segment1(1,1)),0.1,0.1]);
htrace1=plot(segment1*dt,plot_data_mean(segment1,:,1), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(segment2*dt,plot_data_mean(segment2,:,1), 'LineWidth',1.2,'color', color_table(2,:));
% htrace1=plot(x_axis_Slong*dt,plot_data_mean(x_axis_Slong,:,1), 'LineWidth',1.2,'color', color_table(1,:));
% text(x_axis_Slong(1)*dt,plot_data_mean(x_axis_Slong(1),1,1),[num2str(floor(plot_data_mean(x_axis_Slong(1),1,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_Sshort(1,1)*dt,ylim_data(1),(x_axis_Sshort(1,end)-x_axis_Sshort(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
set(rec2,'Position',[x_axis_Sshort(2,1)*dt,ylim_data(1),(x_axis_Sshort(2,end)-x_axis_Sshort(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

switch exp_type
    case 1
        rec3=rectangle('Position',[stim1_X{1}(1)*dt,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b')
    case 3
        rec3=rectangle('Position',[stim1_X{1}(1)*dt,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');    
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)
% %plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=2;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12;perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
%plot std trace
        s3=figure;
        hold on
        rec1=rectangle('Position',[segment1(1,1)*dt,plot_data_std(segment1(1,1),1,1),0.1,0.1]);
         rec2=rectangle('Position',[segment1(1,1)*dt,plot_data_std(segment1(1,1)),0.1,0.1]);
htrace1=plot(segment1*dt,plot_data_std(segment1,:,1), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(segment2*dt,plot_data_std(segment2,:,1), 'LineWidth',1.2,'color', color_table(2,:));
% htrace1=plot(x_axis_Slong*dt,plot_data_std(x_axis_Slong,:,1), 'LineWidth',1.2,'color', color_table(1,:));

axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_Sshort(1,1)*dt,ylim_data(1),(x_axis_Sshort(1,end)-x_axis_Sshort(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
set(rec2,'Position',[x_axis_Sshort(2,1)*dt,ylim_data(1),(x_axis_Sshort(2,end)-x_axis_Sshort(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

switch exp_type
    case 1
        rec3=rectangle('Position',[stim1_X{1}(1)*dt,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
     case 3
        rec3=rectangle('Position',[stim1_X{1}(1)*dt,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');    
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)
% %plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=1;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
        ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;
%%
    %plot spont traces short - no ES:
    trace_to_plot = plot_data(:,trace_ind,1)+DC ; %DC is added just to make space between the traces 
        s4=figure;
        hold on
         rec1=rectangle('Position',[x_axis_Sshort(1,1)*dt,trace_to_plot(x_axis_Sshort(1,1),1),0.1,0.1]);
htrace1=plot(x_axis_Sshort(1,:)*dt,trace_to_plot(x_axis_Sshort(1,:),:), 'LineWidth',1.2,'color', color_table(1,:)); 

% for i=1:length(trace_ind);
%     text(x_axis_Sshort(1,1)*dt,trace_to_plot(x_axis_Sshort(1,1),i),[num2str(floor(plot_data(x_axis_Sshort(1,1),i,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_Sshort(1,1)*dt,ylim_data(1),(x_axis_Sshort(1,end)-x_axis_Sshort(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
    %plot spont traces short - with ES:
    trace_to_plot = plot_data(:,trace_ind,1)+DC ; %DC is added just to make space between the traces 
        s5=figure;
        hold on
        rec1=rectangle('Position',[x_axis_Sshort(2,1)*dt,trace_to_plot(x_axis_Sshort(2,1),1),0.1,0.1]);
htrace1=plot(x_axis_Sshort(2,:)*dt,trace_to_plot(x_axis_Sshort(2,:),:), 'LineWidth',1.2,'color', color_table(2,:));

axis tight
% ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_Sshort(2,1)*dt,ylim_data(1),(x_axis_Sshort(2,end)-x_axis_Sshort(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
if exp_type==2;
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line((stim1_X{1}).*dt,stim1_Y,'LineWidth',6,'Color','b') 
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
    %plot mean trace short:   
        s6=figure;
        hold on
        rec1=rectangle('Position',[x_axis_Sshort(1,1)*dt,plot_data_mean(x_axis_Sshort(1,1),:,1),0.1,0.1]);
htrace1=plot(x_axis_Sshort(1,:)*dt,plot_data_mean(x_axis_Sshort(1,:),:,1), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_Sshort(1,:)*dt,plot_data_mean(x_axis_Sshort(2,:),:,1), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis_Sshort(1,1)*dt,plot_data_mean(x_axis_Sshort(1,1),1),[num2str(floor(plot_data_mean(x_axis_Sshort(1,1),1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_Sshort(1,1)*dt,ylim_data(1),(x_axis_Sshort(1,end)-x_axis_Sshort(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
if exp_type==2;
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line((stim1_X{1}-x_axis_Sshort(2,1)+x_axis_Sshort(1,1)).*dt,stim1_Y,'LineWidth',6,'Color','b') 
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)
%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=2;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
    %plot std trace short:   
        s7=figure;
        hold on
         rec1=rectangle('Position',[x_axis_Sshort(1,1)*dt,plot_data_std(x_axis_Sshort(1,1),:,1),0.1,0.1]);
htrace1=plot(x_axis_Sshort(1,:)*dt,plot_data_std(x_axis_Sshort(1,:),:,1), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_Sshort(1,:)*dt,plot_data_std(x_axis_Sshort(2,:),:,1), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis_Sshort(1,1)*dt,plot_data_std(x_axis_Sshort(1,1),1),[num2str(floor(plot_data_std(x_axis_Sshort(1,1),1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_Sshort(1,1)*dt,ylim_data(1),(x_axis_Sshort(1,end)-x_axis_Sshort(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
if exp_type==2;
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line((stim1_X{1}-x_axis_Sshort(2,1)+x_axis_Sshort(1,1)).*dt,stim1_Y,'LineWidth',6,'Color','b') 
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=1;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);

ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;


%%
if save_flag==1;    
            cd(path_output1)
                          
        saveas(s1,['f' num2str(files_to_analyze(fileind)) '_traces_x1_v2'],'fig') 
        print(s1,['f' num2str(files_to_analyze(fileind)) '_traces_x1_v2'],'-dpng','-r600','-opengl') 
        saveas(s2,['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2'],'fig') 
        print(s2,['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2'],'-dpng','-r600','-opengl') 
        saveas(s3,['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_v2'],'fig') 
        print(s3,['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_v2'],'-dpng','-r600','-opengl') 
    
         cd(path_output2)
    
       saveas(s4,['f' num2str(files_to_analyze(fileind)) '_traces_x1_noES_v2'],'fig') 
        print(s4,['f' num2str(files_to_analyze(fileind)) '_traces_x1_noES_v2'],'-dpng','-r600','-opengl') 
        saveas(s5,['f' num2str(files_to_analyze(fileind)) '_traces_x1_ES_v2'],'fig') 
        print(s5,['f' num2str(files_to_analyze(fileind)) '_traces_x1_ES_v2'],'-dpng','-r600','-opengl') 
        saveas(s6,['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2'],'fig') 
        print(s6,['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2'],'-dpng','-r600','-opengl') 
        saveas(s7,['f' num2str(files_to_analyze(fileind)) '_std_1_mean-subt_v2'],'fig') 
        print(s7,['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_v2'],'-dpng','-r600','-opengl') 
end

%% Evoked figures

    %% plot evoked traces no ES:
    trace_to_plot = plot_data(:,trace_ind,2)+DC ; %DC is added just to make space between the traces 
    y_bl=DC(1:2,:)'+bl;
    x_bl=[x_axis_Elong(1).*ones(size(y_bl,1),1), x_axis_Elong(end).*ones(size(y_bl,1),1)].*dt;  
  
        e1=figure;
        hold on
        rec1=rectangle('position',[x_axis_Elong(1)*dt ,min(trace_to_plot(x_axis_Elong,2,:)),0.1,0.1]);
htrace=plot(x_axis_Elong*dt,trace_to_plot(x_axis_Elong,:,:), 'LineWidth',1.2,'color', color_table(1,:));
if baseline_flag==1;
    h_baseline=line(x_bl',y_bl','linestyle','--','color',[136 137 138]./256,'linewidth',1);
end
% for i=1:length(trace_ind);
%     text(x_axis_Elong(1)*dt,trace_to_plot(x_axis_Elong(1),i),[num2str(floor(plot_data(x_axis_Elong(1),i,2))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
set(rec1,'Position',[x_axis_Elong(1)*dt,ylim_data(1),(x_axis_Elong(end)-x_axis_Elong(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
%plot light stimulus
% switch exp_type
%     case 1
%     case 2
%         stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
%         line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
% end
% plotting stim_2:
if plot_stim_2(2);
patch_xdata=[stim2{1},stim2{2}; flipud(stim2{1}), flipud(stim2{2})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);        
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
%plot evoked traces with ES:
trace_to_plot = plot_data(:,trace_ind,3)+DC ; %DC is added just to make space between the traces 
% trace_to_plot(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
        e2=figure;
        hold on
        rec1=rectangle('position',[x_axis_Elong(1)*dt ,min(trace_to_plot(x_axis_Elong,2,:)),0.1,0.1]);
htrace=plot(x_axis_Elong*dt,trace_to_plot(x_axis_Elong,:,:), 'LineWidth',1.2,'color', color_table(2,:));
if baseline_flag==1;
    h_baseline=line(x_bl',y_bl','linestyle','--','color',[136 137 138]./256,'linewidth',1);
end

% for i=1:length(trace_ind);
%     text(x_axis_Elong(1)*dt,trace_to_plot(x_axis_Elong(1),i),[num2str(floor(plot_data(x_axis_Elong(1),i,3))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');

% plotting zoom-in area:
set(rec1,'Position',[x_axis_Elong(1)*dt,ylim_data(1),(x_axis_Elong(end)-x_axis_Elong(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plot light stimulus
switch exp_type
    case 1
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
end

%plotting stim_2:
if plot_stim_2(3);
patch_xdata=[stim2{1},stim2{2}; flipud(stim2{1}), flipud(stim2{2})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

 %% plot evoked mean trace
        e3=figure;
        hold on
        rec1=rectangle('position',[x_axis_Elong(1)*dt ,min(plot_data_mean(x_axis_Elong,:,2)),0.1,0.1]);
htrace1=plot(x_axis_Elong*dt,plot_data_mean(x_axis_Elong,:,2), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_Elong*dt,plot_data_mean(x_axis_Elong,:,3), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis_Elong(1)*dt,plot_data_mean(x_axis_Elong(1),1,2),[num2str(floor(plot_data_mean(x_axis_Elong(1),1,2))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
set(rec1,'Position',[x_axis_Elong(1)*dt,ylim_data(1),(x_axis_Elong(end)-x_axis_Elong(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plot light stimulus
switch exp_type
    case 1
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
end

%for plotting stim_2:
if plot_stim_2(3);
patch_xdata=[stim2{1},stim2{2}; flipud(stim2{1}), flipud(stim2{2})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=5;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
hold off
ylabel('Mean Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
%plot evoked SD trace
        e4=figure;
        hold on
        rec1=rectangle('position',[x_axis_Elong(1)*dt ,min(plot_data_std(x_axis_Elong,:,2)),0.1,0.1]);
htrace1=plot(x_axis_Elong*dt,plot_data_std(x_axis_Elong,:,2), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_Elong*dt,plot_data_std(x_axis_Elong,:,3), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis_Elong(1)*dt,plot_data_std(x_axis_Elong(1),1,2),[num2str(floor(plot_data_std(x_axis_Elong(1),1,2))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
set(rec1,'Position',[x_axis_Elong(1)*dt,ylim_data(1),(x_axis_Elong(end)-x_axis_Elong(1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plot light stimulus
switch exp_type
    case 1
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
end

%for plotting stim_2:
if plot_stim_2(3);
patch_xdata=[stim2{1},stim2{2}; flipud(stim2{1}), flipud(stim2{2})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=1;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);hold off
ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
if save_flag==1;
    cd(path_output1)
        saveas(e1,['f' num2str(files_to_analyze(fileind)) '_traces_x2_v2.fig']) 
        print(e1,['f' num2str(files_to_analyze(fileind)) '_traces_x2_v2'],'-dpng','-r600','-opengl') 
        saveas(e2,['f' num2str(files_to_analyze(fileind)) '_traces_x3_v2.fig']) 
        print(e2,['f' num2str(files_to_analyze(fileind)) '_traces_x3_v2'],'-dpng','-r600','-opengl') 
        saveas(e3,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_v2.fig']) 
        print(e3,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_v2'],'-dpng','-r600','-opengl') 
        saveas(e4,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt_v2.fig']) 
        print(e4,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt_v2'],'-dpng','-r600','-opengl') 
end
%%
%plot evoked short traces no ES:
    trace_to_plot = plot_data(:,trace_ind,2)+DC ; %DC is added just to make space between the traces 
    y_bl=DC(1:2,:)'+bl;
    x_bl=[x_axis_Eshort(1,1).*ones(size(y_bl,1),1), x_axis_Eshort(end,1).*ones(size(y_bl,1),1)].*dt;  
  
        e5=figure;
        hold on
        rec1=rectangle('position',[x_axis_Eshort(1,1)*dt ,min(trace_to_plot(x_axis_Eshort(1,:),2,:)),0.1,0.1]);
htrace=plot(x_axis_Eshort(1,:)*dt,trace_to_plot(x_axis_Eshort(1,:),:), 'LineWidth',1.2,'color', color_table(1,:));
if baseline_flag==1;
    h_baseline=line(x_bl',y_bl','linestyle','--','color',[136 137 138]./256,'linewidth',1);
end
% for i=1:length(trace_ind);
%     text(x_axis_Eshort(1,1)*dt,trace_to_plot(x_axis_Eshort(1,:),i),[num2str(floor(plot_data(x_axis_Eshort(1,1),i,2))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
set(rec1,'Position',[x_axis_Eshort(1,1)*dt,ylim_data(1),(x_axis_Eshort(1,end)-x_axis_Eshort(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
%plot light stimulus
% switch exp_type
%     case 1
%     case 2
%         stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
%         line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
% end
% plotting stim_2:
if plot_stim_2(2);
patch_xdata=[stim2{1}; flipud(stim2{1})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);         
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
%plot evoked traces with ES:
trace_to_plot = plot_data(:,trace_ind,3)+DC ; %DC is added just to make space between the traces 
% trace_to_plot(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
        e6=figure;
        hold on
        rec1=rectangle('position',[x_axis_Eshort(2,1)*dt ,min(trace_to_plot(x_axis_Eshort(2,:),2,:)),0.1,0.1]);
htrace=plot(x_axis_Eshort(2,:)*dt,trace_to_plot(x_axis_Eshort(2,:),:), 'LineWidth',1.2,'color', color_table(2,:));
if baseline_flag==1;
    h_baseline=line(x_bl',y_bl','linestyle','--','color',[136 137 138]./256,'linewidth',1);
end

% for i=1:length(trace_ind);
%     text(x_axis_Eshort(1,2)*dt,trace_to_plot(x_axis_Eshort(1,2),i),[num2str(floor(plot_data(x_axis_Eshort(1,2),i,3))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');

% plotting zoom-in area:
set(rec1,'Position',[x_axis_Eshort(2,1)*dt,ylim_data(1),(x_axis_Eshort(2,end)-x_axis_Eshort(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plot light stimulus
switch exp_type
    case 1
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
end

%plotting stim_2:
if plot_stim_2(3);
patch_xdata=[stim2{2}; flipud(stim2{2})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
hold off
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

 %% plot evoked mean trace
        e7=figure;
        hold on
        rec1=rectangle('position',[x_axis_Eshort(1,1)*dt ,min(plot_data_mean(x_axis_Eshort(1,:),:,2)),0.1,0.1]);
htrace1=plot(x_axis_Eshort(1,:)*dt,plot_data_mean(x_axis_Eshort(1,:),:,2), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_Eshort(1,:)*dt,plot_data_mean(x_axis_Eshort(2,:),:,3), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis_Eshort(1,1)*dt,plot_data_mean(x_axis_Eshort(1,1),1,2),[num2str(floor(plot_data_mean(x_axis_Eshort(1,1),1,2))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
set(rec1,'Position',[x_axis_Eshort(1,1)*dt,ylim_data(1),(x_axis_Eshort(1,end)-x_axis_Eshort(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plot light stimulus
switch exp_type
    case 1
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
end

%for plotting stim_2:
if plot_stim_2(3);
patch_xdata=[stim2{1}; flipud(stim2{1})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=2;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
hold off
ylabel('Mean Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
%plot evoked SD trace
        e8=figure;
        hold on
        rec1=rectangle('position',[x_axis_Eshort(1,1)*dt ,min(plot_data_std(x_axis_Eshort(1,:),:,2)),0.1,0.1]);
htrace1=plot(x_axis_Eshort(1,:)*dt,plot_data_std(x_axis_Eshort(1,:),:,2), 'LineWidth',1.2,'color', color_table(1,:));
htrace2=plot(x_axis_Eshort(1,:)*dt,plot_data_std(x_axis_Eshort(2,:),:,3), 'LineWidth',1.2,'color', color_table(2,:));
% text(x_axis_Eshort(1)*dt,plot_data_std(x_axis_Eshort(1),1,2),[num2str(floor(plot_data_std(x_axis_Eshort(1),1,2))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
% plotting zoom-in area:
set(rec1,'Position',[x_axis_Eshort(1,1)*dt,ylim_data(1),(x_axis_Eshort(1,end)-x_axis_Eshort(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');

%plot light stimulus
switch exp_type
    case 1
    case 2
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',6,'Color','b') 
end

%for plotting stim_2:
if plot_stim_2(3);
patch_xdata=[stim2{1}; flipud(stim2{1})];
yex=wextend('1D','sym',ylim_data,1);
l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
temp_y=wextend('ac','sym',yex,l);
patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=1;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);hold off
ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
if save_flag==1;
   cd(path_output2)
    
        saveas(e5,['f' num2str(files_to_analyze(fileind)) '_traces_x2_v2'],'fig') 
        print(e5,['f' num2str(files_to_analyze(fileind)) '_traces_x2_v2'],'-dpng','-r600','-opengl') 
        saveas(e6,['f' num2str(files_to_analyze(fileind)) '_traces_x3_v2'],'fig') 
        print(e6,['f' num2str(files_to_analyze(fileind)) '_traces_x3_v2'],'-dpng','-r600','-opengl') 
        saveas(e7,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_v2'],'fig') 
        print(e7,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_v2'],'-dpng','-r600','-opengl') 
        saveas(e8,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt_v2'],'fig') 
        print(e8,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt_v2'],'-dpng','-r600','-opengl') 
end