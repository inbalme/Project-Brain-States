
clear all
close all
cc_stat=[]; cc_spont=[]; cc_evoked=[]; cc=[]; cc_shuffled_it=[]; cc_shuff_sub=[];
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data Ch2_data stim2_X stim1_X
 
 global exp_type
exp_type=2; %1-NBES, 2-ChAT
trace_type_input=[1]; %[3,2] for exp_type=1; %for exp_type=2 or 3 use [1,2]
analyze_time_before_train=0;
analyze_train_only_flag=1; %use analyze_train_only_flag=1
add_to_plot=0.08; %seconds from each side of the trace. use 0.1 for NBES and 0.080 for ChAT
 plot_trace=[2,4,5];%1:size(current_data,2); %5;for f46: [2,3,4], [3 4 6]; for f80: [2 4 5]
save_flag=0;
print_flag=1;
norm_flag=0;
clamp_flag=[]; %[]; %3; %clamp_flag=1 for hyperpolarization traces, clamp_flag=2 for depolarization traces and clamp_flag=3 for no current traces (only clamp to resting Vm)
BP50HzLFP_flag=1; %removing 50Hz noise from LFP signal
BP50HzVm_flag=1; %removing 50Hz noise from Vm signal
BPLFP_flag=1; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=[1,200];%[0.1,200]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=0; %filtering LFP and Vm same as LFP was filtered in the multiclamp
bp_manual_Vm=[0,300];%[0,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
    
switch exp_type
    case 1
        files_to_analyze =[44,46,48,52,56,58,62,72,75,82,84]; %[44,46,48,50,52,56,58,62,72,75]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB-', 'NB+'}; y_ax_label={'Vm'}; y_ax_units={'mV'};   
        legend_string_shuff={'NB- shuffled', 'NB+ shuffled'};       

    case 2
        files_to_analyze =80; %[76,77,80,82,84,87,90,92,112, 114, 115]; %[92,94,96,98,101]
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={ 'Light Off','Light On'}; y_ax_label={'Vm'}; y_ax_units={'mV'};   
        legend_string_shuff={'Light Off shuffled','Light On shuffled'};
        
    case 3 
        files_to_analyze =[31,38,42,51,69,71,74]; %[31,38,42,51,61,64,67,69,71,74,77]; [51,67];
        clamp_flag=3;
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};    y_ax_label={'Im'}; y_ax_units={'pA'};     
        legend_string_shuff={'NB+ shuffled', 'NB- shuffled'};      
end
    
    for fileind=1:length(files_to_analyze) ;    
        close all
    clearvars -except cc_stat cc_spont cc_evoked  files_to_analyze fileind files cc_spont_for_xls_mean...
        cc_evoked_for_xls_mean cc lags cc_shuffled_mean cc_shuffled_it cc_mean cc_shuff_sub_mean save_flag print_flag...
        cc_lag0_mat cc_lag0_shuff_mat cc_max_mat cc_max_time_mat cc_maxdiff_mat  cc_max_shuff_mat...
        norm_flag  BP50HzLFP_flag BP50HzVm_flag BPLFP_flag bp_manual_LFP BPVm_flag bp_manual_Vm exp_type trace_type_input...
        legend_string legend_string_shuff analyze_time_before_train analyze_train_only_flag clamp_flag y_ax_label y_ax_units add_to_plot...
        plot_trace
   
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
 
 path_output='VmLFP traces';  
 
 clear color_table    
    whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;

switch exp_type
    case 1
        color_table=[0 0 0; [30,75,14]/256; [136 137 138]/256; [112,172,90]/256; [216 22 22]/256; [255 153 153]/256];
        path_output=['D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\',path_output];   
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations'
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        cd(path_output)       
    case 2
%         color_table=[0 0 0; [30,75,14]/256; [136 137 138]/256; [112,172,90]/256; [0 0 204]/256; [255 153 153]/256];
        color_table=[0 0 0; [160,160,160]./256; [136 137 138]/256; [102, 172,255]./256; [0 0 204]/256; [255 153 153]/256];
        path_output=['D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm-LFP correlations\',path_output]; 
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm-LFP correlations'
        a = exist(path_output,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        cd(path_output)
    case 3 
        color_table=[0 0 0; [30,75,14]/256; [136 137 138]/256; [112,172,90]/256; [216 22 22]/256; [255 153 153]/256];
        path_output=['D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations_VC\',path_output]; 
%         cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations_VC'
        a = exist(path_output,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        cd(path_output)       
end
 
%%
for trace_type=trace_type_input; %1:2; 
        norm_LFP=[]; norm_Vm=[]; 
        %%
           interval=[];   interval_temp=[];     
%         intervals_to_analyze

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
        stim2{1}=stim2_X{2};
        stim2{2}=stim2_X{2};
        x_value=1:3;   
        galvano_nstim = Param.facade(6);
        galvano_freq = Param.facade(7);
        curr_inj_1=Param.facade(22); %in sec
        curr_inj_dur=Param.facade(20); %in sec
        curr_inj_2=Param.facade(22)+curr_inj_dur+Param.facade(23); %in sec        
        curr_inj_1_end=curr_inj_1+curr_inj_dur+0.1;
        curr_inj_2_end=curr_inj_2+curr_inj_dur+0.1;
        x_axis_Slong=round(1:12*sf{1}-1);
        x_axis_Sshort(1,:)=round(0.4*sf{1}:1.4*sf{1});
        x_axis_Sshort(2,:)=stim1_X{1}(1,1)+0.4*sf{1}:stim1_X{1}(1,1)+1.4*sf{1};
%         x_axis_Sshort(1,:)=round(0.4*sf{1}:2.9*sf{1});
%         x_axis_Sshort(2,:)=stim1_X{1}(1,1)+0.4*sf{1}:stim1_X{1}(1,1)+2.9*sf{1};

        x_axis_Elong=round(1*sf{1}:12*sf{1}-1);
%         x_axis_Eshort(1,:)=round(stim2{1}(1,1)-0.05*sf{1}:stim2{1}(1,1)+2.45*sf{1});
%         x_axis_Eshort(2,:)=round(stim2{1}(1,1)-0.05*sf{1}:stim2{1}(1,1)+2.45*sf{1});
        x_axis_Eshort(1,:)=round(stim2{1}(1,1)-0.18*sf{1}:stim2{1}(1,1)+1.82*sf{1});
        x_axis_Eshort(2,:)=round(stim2{1}(1,1)-0.18*sf{1}:stim2{1}(1,1)+1.82*sf{1});
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
%%
%% Subtract mean from the interval (fragment of trace)
% data_LFP{t}=Ch2_data_filt(interval(:,t),:,x_value(t));
% data_LFP_noDC{t} = fn_Subtract_Mean(data_LFP{t});
% data_Vm{t}=current_data_filt(interval(:,t),:,x_value(t));   
% data_Vm_noDC{t} = fn_Subtract_Mean(data_Vm{t});
% data_Vm_filt{t}=current_data_filt(interval(:,t),:,x_value(t));
% data_Vm_filt_noDC{t} = fn_Subtract_Mean(data_Vm_filt{t});
% for plots:
%      interval_plot(:,t)=[interval(1,t)-round(add_to_plot*sf{1}):interval(end,t)+round(add_to_plot*sf{1})]; 
%     data_Vm_plot{t}=current_data_filt(interval_plot(:,t),:,x_value(t));
%     data_LFP_plot{t}=Ch2_data_filt(interval_plot(:,t),:,x_value(t)).*20;
    
    plot_data_Vm=current_data_filt(:,:,x_value); %current_data_filt %current_data_no_DC 
    plot_data_LFP=Ch2_data_filt(:,:,x_value).*20; %current_data_filt %current_data_no_DC 
if exp_type==2 || any(files_to_analyze(fileind)==[62,72,75,82,84]);
    for x=x_value
        plot_data_Vm_no_step=plot_data_Vm(:,:,x); 
        plot_data_Vm_no_step(curr_inj_1*sf{1}-100:curr_inj_1_end*sf{1}+100,:)=nan;
        plot_data_Vm_no_step(curr_inj_2*sf{1}-100:curr_inj_2_end*sf{1}+100,:)=nan;
        min_trace=min(plot_data_Vm_no_step);
        diffmat=bsxfun(@minus,plot_data_Vm(:,:,x),min_trace);
        plot_data_Vm_copy=plot_data_Vm(:,:,x); 
        plot_data_Vm_copy(diffmat<0)=nan;
        plot_data_Vm(:,:,x)=plot_data_Vm_copy;
        
        plot_data_LFP_no_step=plot_data_LFP(:,:,x); 
        plot_data_LFP_no_step(curr_inj_1*sf{1}-100:curr_inj_1_end*sf{1}+100,:)=nan;
        plot_data_LFP_no_step(curr_inj_2*sf{1}-100:curr_inj_2_end*sf{1}+100,:)=nan;
        min_trace=min(plot_data_LFP_no_step);
        diffmat=bsxfun(@minus,plot_data_LFP(:,:,x),min_trace);
        plot_data_LFP_copy=plot_data_LFP(:,:,x); 
        plot_data_LFP_copy(diffmat<0)=nan;
        plot_data_LFP(:,:,x)=plot_data_LFP_copy;
    end
%     mean_trace=mean(mean(plot_data((curr_inj_2-0.5)*sf{1}:curr_inj_2*sf{1},:,1),2));
%     plot_data(plot_data(curr_inj_2*sf{1}:curr_inj_2_end*sf{1},:,:)<(mean_trace-2),:,:)=nan;
end
if exp_type==1||exp_type==3
    plot_data_Vm(stim1_X{1}(1,1)-50:stim1_X{1}(2,1)+50,:,:)=nan;
    current_data_Vm_no_DC(stim1_X{1}(1,1)-50:stim1_X{1}(2,1)+50,:,:)=nan;
end

%% Plots
trace_fontsize=12;
scalebar_fontsize=11;

    if print_flag==1;
% plotting one trace of data and LFP against each other -
lengthh_vert=[5,10,5]; %[1,2,10]; %lengthh_vert(1) is for spontaneous, lengthh_vert(2) is for evoked

%          interval_plot(:,1)=[interval(1,1)-add_to_plot*sf{1}:interval(end,1)+add_to_plot*sf{1}];  
%          data_Vm_plot=data_Vm_filt{1}(:,plot_trace);
%          data_LFP_plot=data_LFP{1}(:,plot_trace).*20;

Fig{fileind}(trace_type)=figure;
ymax=max([plot_data_Vm(:,plot_trace),plot_data_Vm(:,plot_trace)]);
ymin=min([plot_data_Vm(:,plot_trace),plot_data_Vm(:,plot_trace)]);
ydiff=ymax-ymin;
yrange=ceil(1.2*max(ydiff));
ymaxLFP=max([plot_data_LFP(:,plot_trace),plot_data_LFP(:,plot_trace)]);
yminLFP=min([plot_data_LFP(:,plot_trace),plot_data_LFP(:,plot_trace)]);
ydiffLFP=ymaxLFP-yminLFP;
yrangeLFP=ceil(max(ydiffLFP));

for tr_ind=1:length(plot_trace)
subplot(2*length(plot_trace),1,2*tr_ind)
    hold on
%         p1=plot(interval_plot(:,1).*dt*1000,plot_data_Vm{1}(:,plot_trace(tr_ind)), 'color',color_table(1,:),'LineWidth',1.5);
       htrace1=plot(segment1*dt,plot_data_Vm(segment1,plot_trace(tr_ind)), 'LineWidth',1.5,'color', color_table(1,:));
       htrace2=plot(segment2*dt,plot_data_Vm(segment2,plot_trace(tr_ind)), 'LineWidth',1.5,'color', color_table(5,:));

        axis tight
     ylim_data=[get(gca,'ylim')]';
     ylim_data(2)=ylim_data(1)+(1.1).*yrange;
     xlim_data=[get(gca,'xlim')]';
  set(gca,'ylim',ylim_data, 'xlim',xlim_data)  
    
          if trace_type==2;   
%             fn_plot_sensory_stim(dt, stim2{1}(:,1:size(stim2{1},2)-1),whisker_stim_color)
        %dashed lines:
                    yline=ones(size(stim2{2})); 
                    yline(1,:)=yline(1,:)*ylim_data(1);
                    yline(2,:)=yline(2,:)*ylim_data(2);
                    xline=[stim2{2}(1,:).*dt.*1000;stim2{2}(1,:).*dt.*1000];
                    line(xline,yline,'linestyle','--','LineWidth',1,'Color',[0 0 0 0.6])
       end
            
        %plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' s'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=0.03; perc2=[];
        [b1,b2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=lengthh_vert(trace_type);     textit=[num2str(lengthh), ' mV'];     c=[0,0,0];  fonsizes=scalebar_fontsize; %perc1=0.03; perc2=[];
        [b3,b4] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
         if tr_ind~=length(plot_trace)-1;
             delete(b1); delete(b2); delete(b3); delete(b4);
        end
    hold off
%     uistack(htrace1,'top');
    set(gca, 'visible', 'off') ;
    
  subplot(2*length(plot_trace),1,2*tr_ind-1)     
 	hold on
%         p2=plot(interval_plot(:,1).*dt*1000,data_LFP_plot{1}(:,plot_trace(tr_ind)),'color',color_table(2,:),'LineWidth',1.5);
        htrace3=plot(segment1*dt,plot_data_LFP(segment1,plot_trace(tr_ind)), 'LineWidth',1.5,'color', color_table(3,:));
        htrace4=plot(segment2*dt,plot_data_LFP(segment2,plot_trace(tr_ind)), 'LineWidth',1.5,'color', color_table(4,:));
        
        axis tight            
     ylim_data=[get(gca,'ylim')]';
     ylim_data(1)=ylim_data(2)-yrangeLFP;   
     xlim_data=[get(gca,'xlim')]';
     ylim_data(2)=ylim_data(2)+(0.1).*yrangeLFP; 
     set(gca,'ylim',ylim_data, 'xlim',xlim_data)  
          
if tr_ind==1;
    switch exp_type
    case 1
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-0.005,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt+0.005,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
    case 2
%         stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        stim1_Y = ones(size(stim1_X{1})).*ylim_data(2); 
        line(stim1_X{1}.*dt,stim1_Y,'LineWidth',3,'Color','b') 
    case 3
        rec3=rectangle('Position',[stim1_X{1}(1)*dt-0.005,ylim_data(1),(stim1_X{1}(2)-stim1_X{1}(1))*dt+0.005,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');
    end
end


if trace_type==2;
%              fn_plot_sensory_stim(dt, stim2{1}(:,1:size(stim2{1},2)-1),whisker_stim_color)
        %dashed lines:
                    yline=ones(size(stim2{2})); 
                    yline(1,:)=yline(1,:)*ylim_data(1);
                    yline(2,:)=yline(2,:)*ylim_data(2);
                    xline=[stim2{2}(1,:).*dt.*1000;stim2{2}(1,:).*dt.*1000];
                    line(xline,yline,'linestyle','--','LineWidth',1,'Color',[0 0 0 0.6])
       end
        
        %plotting scale bar
             horiz_vert=1;        lengthh=200;     textit=[num2str(lengthh), ' mS'];     c=[0,0,0];  fonsizes=scalebar_fontsize; %perc1=[]; perc2=[];
                [b1,b2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
             horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh./20), ' mV'];     c=[0,0,0];  fonsizes=scalebar_fontsize; %perc1=[]; perc2=[];
                [b3,b4] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
           delete(b1); delete(b2); 
         if tr_ind~=length(plot_trace);
             delete(b3); delete(b4); 
        end
    hold off  
%     uistack(htrace3,'top');
    set(gca, 'visible', 'off') ;

end
        line1=findall(gcf,'type','line');
        linexdata=get(line1(8),'xdata');
        set(line1(8),'xdata',linexdata-0.3);
        set(line1(1),'xdata',linexdata-0.3);
    title('Vm-LFP single trace ES Off','FontSize', 16); 
    ylabel('Potential [mV]' ,'FontSize', 14); xlabel('Time (s)','FontSize', 14)
l=legend([htrace3,htrace4,htrace1,htrace2],{'LFP light off','LFP light on','Vm light off','Vm light on'}, 'box', 'off','position',[0.8,0.85,0.2,0.1]); 
    l.FontSize=11;
        end
    end
    end
    
cd(path_output)
   if save_flag==1
        saveas(Fig{fileind}(trace_type),['Vm-LFP_spont_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'fig') 
        print(Fig{fileind}(trace_type),['Vm-LFP_spont_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'-dpng','-r600','-opengl') 
   end