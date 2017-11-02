%% Load file
close all
clear all
%I changed line 243
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data
 global exp_type
exp_type=2; %1-NBES, 2-ChAT
data_type='Vm'; %'LFP', 'Vm'
save_flag= 0;
print_flag=1;
choose_traces=[]; %if choose_traces is empty the script takes all traces for analysis. If it is not empty the script takes the specified traces.
clamp_flag=[]; %[]; %3; %clamp_flag=1 for hyperpolarization traces, clamp_flag=2 for depolarization traces and clamp_flag=3 for no current traces (only clamp to resting Vm)
% short_flag=0; %1- short trace, 0- long trace
baseline_flag=0; %adding dashed line under traces
bl=-70; %baseline value;
norm_flag=0;
BP50HzLFP_flag=0; %removing 50Hz noise from LFP signal
BP50HzVm_flag=1; %removing 50Hz noise from Vm signal
BPLFP_flag=1; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=0; %[1,200]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=1; %filtering LFP and Vm same as LFP was filtered in the multiclamp
bp_manual_Vm=[0,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
lengthh_vert=[2,5,10; 2,5,10]; %[0.1,0.1,0.5];%f46:[2,5,10; 2,5,10]; f80:[1,2,10; 1,2,10] %lengthh_vert(1) is for STD, lengthh_vert(2) is for mean and (3) is for traces. row 1 spont. row 2 evoked
trace_ind =[]; %[2,3,4,5,6];%[2,3,4,5,6]; %if trace_ind is empty, the default would be to take all traces
DC_factor =16; %12; %sets the spacing between plotted traces. set DC_factor=1 for LFP [f44 13; f46 16] [f80 9; f74 12] 
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
        files_to_analyze =46; %[8,10,12,14,15,16,22,36,37,40,1,44,46,48,52,56,58,62,72,75,82,84]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB-', 'NB+'};   y_ax_label={'Vm'}; y_ax_units={'mV'};   
         path_output1='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Long Trace Presentation\style2';  
        path_output2='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean\Zoom-in Trace Presentation\style2';  
        a = exist(path_output1,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        b = exist(path_output2,'dir'); %b will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output1);
        end 
       if b==0;
            mkdir(path_output2);
      end 

    case 2
        files_to_analyze =129; %[74,76,77,80,82,84,87,90,92,112,114,115]; [94,96,98,101] [118,119,120,121,122,126,127,128,129]
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light Off', 'Light On'};    y_ax_label={'Vm'}; y_ax_units={'mV'}; 
         path_output1='D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Long Trace Presentation\style2';  
        path_output2='D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean\Zoom-in Trace Presentation\style2';  
        a = exist(path_output1,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        b = exist(path_output2,'dir'); %b will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output1);
        end 
       if b==0;
            mkdir(path_output2);
      end 
   case 3 
        files_to_analyze =80; %[31,38,42,51,61,64,67,69,71,74,77];
        clamp_flag=3; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB-', 'NB+'};    y_ax_label={'Im'}; y_ax_units={'pA'};    
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
%     data_no_spikes=[];
Ch2_data=[];
if length(raw_data)>=3
    Ch2_data= raw_data{3}./20; %dividing by the LFP gain       
%     Ch2_data= raw_data{3};
end
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
      
    if strcmp(data_type,'LFP')
    current_data=Ch2_data;
    data_used='Ch2_data';
    end

   data_preprocessing 
   
    if isempty(current_data_filt)
     current_data_filt=current_data;
    end
    
        clear color_table  current_data_no_DC  
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

if ~isempty(choose_traces)
    temp=current_data_filt;
    current_data_filt=[];
    current_data_filt=temp(:,choose_traces,:);
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
plot_stim_1=[1];
plot_stim_2=[0];
if isempty(trace_ind)
    trace_ind =[1:size(current_data_filt,2)]; % [1,2,3,4,5,6]; %1:size(plot_data,2);  %trace_ind is the numbers of traces to be plotted
end
switch exp_type
    case 1
        stim2{1}=stim2_X{2};
        stim2{2}=stim2_X{2};
        x_value=1:3;   
        galvano_nstim = Param.facade(6);
        galvano_freq = Param.facade(7);
        if length(Param.facade)>22
            curr_inj_1=Param.facade(22); %in sec
            curr_inj_dur=Param.facade(20); %in sec
            curr_inj_2=Param.facade(22)+curr_inj_dur+Param.facade(23); %in sec        
            curr_inj_1_end=curr_inj_1+curr_inj_dur+0.1;
            curr_inj_2_end=curr_inj_2+curr_inj_dur+0.1;
        end
         if size(current_data_filt,1)>12*sf{1}
            x_axis_Slong=0.4*sf{1}:12*sf{1}-1;
            x_axis_Elong=round((stim1_X{1}(2,1)+10):10*sf{1}-1);
        else
            x_axis_Slong=0.4*sf{1}:8*sf{1}-1;
            x_axis_Elong=round((stim1_X{1}(2,1)+10):10*sf{1}-1);
         end
             x_axis_Sshort(1,:)=0.4*sf{1}:1.4*sf{1};
            x_axis_Sshort(2,:)=round(stim2{1}(1,1)-0.5*sf{1}:stim2{1}(1,1)+0.5*sf{1});
%         x_axis_Sshort(1,:)=0.4*sf{1}:2.9*sf{1};
%         x_axis_Sshort(2,:)=round(stim2{1}(1,1)-0.5*sf{1}:stim2{1}(1,1)+2*sf{1});
        
%         x_axis_Eshort(1,:)=round(stim2{1}(1,1)-0.5*sf{1}:stim2{1}(1,1)+2*sf{1});
%         x_axis_Eshort(2,:)=round(stim2{1}(1,1)-0.5*sf{1}:stim2{1}(1,1)+2*sf{1});
        x_axis_Eshort(1,:)=round(stim2{1}(1,1)-0.18*sf{1}:stim2{1}(1,1)+1.82*sf{1});
        x_axis_Eshort(2,:)=round(stim2{1}(1,1)-0.18*sf{1}:stim2{1}(1,1)+1.82*sf{1});
        rectangle_color=[239 239 239]/256;
        segment1=x_axis_Slong(1):stim1_X{1}(1,1);
        segment2=stim1_X{1}(1,1)+1:x_axis_Slong(end);

    case 2
%         stim2{1}=stim2_X{2};
%         stim2{2}=stim2_X{2};
        x_value=1;   
%         galvano_nstim = Param.facade(6);
%         galvano_freq = Param.facade(7);
%         curr_inj_1=Param.facade(22); %in sec
%         curr_inj_dur=Param.facade(20); %in sec
%         curr_inj_2=Param.facade(22)+curr_inj_dur+Param.facade(23); %in sec        
%         curr_inj_1_end=curr_inj_1+curr_inj_dur+0.1;
%         curr_inj_2_end=curr_inj_2+curr_inj_dur+0.1;
        x_axis_Slong=round(1:12*sf{1}-1);
        x_axis_Sshort(1,:)=round(0.4*sf{1}:1.4*sf{1});
        x_axis_Sshort(2,:)=stim1_X{1}(1,1)+0.4*sf{1}:stim1_X{1}(1,1)+1.4*sf{1};
%         x_axis_Sshort(1,:)=round(0.4*sf{1}:2.9*sf{1});
%         x_axis_Sshort(2,:)=stim1_X{1}(1,1)+0.4*sf{1}:stim1_X{1}(1,1)+2.9*sf{1};

        x_axis_Elong=round(1*sf{1}:12*sf{1}-1);
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
            x_axis_Slong=0.4*sf{1}:12*sf{1}-1;
            x_axis_Elong=1:12*sf{1}-1;
        else
            x_axis_Slong=0.4*sf{1}:8*sf{1}-1;
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
%need to put nans in the place of the step current and then make a double-line
%according to the nans and write "step current".

plot_data=current_data_filt(:,:,x_value); %current_data_filt %current_data_no_DC 
%     mean_trace=mean(mean(plot_data((curr_inj_2-0.5)*sf{1}:curr_inj_2*sf{1},:,1),2));
%     plot_data(plot_data(curr_inj_2*sf{1}:curr_inj_2_end*sf{1},:,:)<(mean_trace-2),:,:)=nan;
if exp_type==1||exp_type==3
    plot_data(stim1_X{1}(1,1)-50:stim1_X{1}(2,1)+50,:,:)=nan;
    current_data_no_DC(stim1_X{1}(1,1)-50:stim1_X{1}(2,1)+50,:,:)=nan;
end
plot_data_mean = nanmean(plot_data,2);  
% plot_data_mean(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_std =  nanstd(plot_data(:,:,x_value),0,2);  %current_data_no_DC
% plot_data_std(stim1_X{1}(1):stim1_X{1}(2),:)=nan;
plot_data_var=nanvar(current_data_no_DC(:,:,x_value),0,2); % plot_data_var=var(plot_data,0,2);
% plot_data_ff = (-1).*plot_data_var./plot_data_mean; 
       l=size(plot_data(:,trace_ind),1);
        l=l/2-1;
        DC=DC_factor.*(0:length(trace_ind)-1);
        if exp_type==3
            DC=40.*(0:length(trace_ind)-1);
        end
%         DC=wextend('addrow', 'per', DC, l);        
              
%         x2lim=[stim2_X{x_value(1)}(1,1).*dt-0.5,stim2_X{2}(1,1).*dt+2];
  
   
%% spontaneous figures    
            %plot spont traces long:
    trace_to_plot = bsxfun(@plus,plot_data(:,trace_ind,1),DC); %DC is added just to make space between the traces 
    
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
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-0.005,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt+0.005,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
    case 2
%         stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2)+0.1.*abs(ylim_data(2)-ylim_data(1)); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',3,'Color','b') 
    case 3
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-0.005,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt+0.005,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

        %plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];  c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=lengthh_vert(1,3);     textit=[num2str(lengthh),y_ax_units{1}];  c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
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
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-0.005,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt+0.005,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
    case 2
%         stim1_Y = ones(size(stim1_X{1})).*ylim_data(2);
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2)+0.1.*abs(ylim_data(2)-ylim_data(1)); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',3,'Color','b')
    case 3
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-0.005,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt+0.005,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)
% %plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=lengthh_vert(1,2);     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12;perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

l=legend([htrace1(1),htrace2(1)],legend_string, 'box', 'off','position',[0.8,0.9,0.2,0.1]); 
l.FontSize=11;
for ix=1:length(l.String)
  str = l.String{ix};
  h = findobj(gcf,'DisplayName',str);
  h.LineWidth =1.5;
end
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
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-0.005,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt+0.005,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
    case 2
%         stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2)+0.1.*abs(ylim_data(2)-ylim_data(1)); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',3,'Color','b') 
     case 3
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-0.005,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt+0.005,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)
% %plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=lengthh_vert(1,1);     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
        ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;
%%
    %plot spont traces short - no ES:
    trace_to_plot = bsxfun(@plus,plot_data(:,trace_ind,1),DC) ; %DC is added just to make space between the traces 
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
 horiz_vert=0;        lengthh=lengthh_vert(1,3);     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;

%%
    %plot spont traces short - with ES:
    trace_to_plot = bsxfun(@plus,plot_data(:,trace_ind,1),DC); %DC is added just to make space between the traces 
        s5=figure;
        hold on
        rec1=rectangle('Position',[x_axis_Sshort(2,1)*dt,trace_to_plot(x_axis_Sshort(2,1)+100,1),0.1,0.1]);
htrace1=plot(x_axis_Sshort(2,:)*dt,trace_to_plot(x_axis_Sshort(2,:),:), 'LineWidth',1.2,'color', color_table(2,:));

axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_Sshort(2,1)*dt,ylim_data(1),(x_axis_Sshort(2,end)-x_axis_Sshort(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
if exp_type==2;
%         stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2)+0.1.*abs(ylim_data(2)-ylim_data(1)); 
        line((stim1_X{1}).*dt,stim1_Y,'LineWidth',3,'Color','b') 
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=lengthh_vert(1,3);     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
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
%         stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2)+0.1.*abs(ylim_data(2)-ylim_data(1)); 
        line((stim1_X{1}-x_axis_Sshort(2,1)+x_axis_Sshort(1,1)).*dt,stim1_Y,'LineWidth',3,'Color','b') 
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)
%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=lengthh_vert(1,2);     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
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
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2)+0.1.*abs(ylim_data(2)-ylim_data(1)); 
        line((stim1_X{1}-x_axis_Sshort(2,1)+x_axis_Sshort(1,1)).*dt,stim1_Y,'LineWidth',3,'Color','b') 
end
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[],'xlim',xlim_data)

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' S'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
horiz_vert=0;        lengthh=lengthh_vert(1,1);     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);

ylabel('Vm SD [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca, 'visible', 'off') ;
%%
if save_flag==1;    
            cd(path_output1)
    if strcmp(data_type,'LFP')
        filename1= ['f' num2str(files_to_analyze(fileind)) '_traces_x1_LFP'];
         filename2= ['f' num2str(files_to_analyze(fileind)) '_mean_x1_LFP'];
         filename3= ['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_LFP'];
         filename4=['f' num2str(files_to_analyze(fileind)) '_traces_x1_noES_LFP'];
        filename5=['f' num2str(files_to_analyze(fileind)) '_traces_x1_ES_LFP'];
        filename6=['f' num2str(files_to_analyze(fileind)) '_mean_x1_LFP'];
        filename7=['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_LFP'];
        else
         filename1= ['f' num2str(files_to_analyze(fileind)) '_traces_x1_v2'];
         filename2= ['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2'];
         filename3= ['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_v2'];
         filename4=['f' num2str(files_to_analyze(fileind)) '_traces_x1_noES_v2'];
        filename5=['f' num2str(files_to_analyze(fileind)) '_traces_x1_ES_v2'];
        filename6=['f' num2str(files_to_analyze(fileind)) '_mean_x1_v2'];
        filename7=['f' num2str(files_to_analyze(fileind)) '_std_x1_mean-subt_v2'];
    end
         
        saveas(s1,filename1,'fig') 
        print(s1,filename1,'-dpng','-r600','-opengl') 
        saveas(s2,filename2,'fig') 
        print(s2,filename2,'-dpng','-r600','-opengl') 
        saveas(s3,filename3,'fig') 
        print(s3,filename2,'-dpng','-r600','-opengl') 
    
         cd(path_output2)
            
         saveas(s4,filename4,'fig') 
        print(s4,filename4,'-dpng','-r600','-opengl') 
        saveas(s5,filename5,'fig') 
        print(s5,filename5,'-dpng','-r600','-opengl') 
        saveas(s6,filename6,'fig') 
        print(s6,filename6,'-dpng','-r600','-opengl') 
        saveas(s7,filename7,'fig') 
        print(s7,filename7,'-dpng','-r600','-opengl') 
end
