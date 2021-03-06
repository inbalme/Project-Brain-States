%Analyze_VmLFPcc_galvano_v3
% This file was created on 27/12/2016, based on Analyze NBES protocol
% ES+galvano v3.
%This file is used for the analysis of files created with
%extract_NBES_Data_v2 or extract_ChAT_Data_v3

%       protocol 5 (NBES+galvnao (5 x-values))
%       protocol 6 (ChAT intracellular current injections light+galvnao (15 x-values))
%       protocol 8 - ChAT intracellular: light only (3 x-values): 1) train
%       20 Hz, 10ms ondur, 50 pulses. 2) single pulse 5 sec ondur. 3)single pulse 10 ms
%       protocol 9 - IC: I-V curve
%       protocol 10 (NBES+galvnao train+test (3 x-values)) or 
%                 11 - NBES: 3 currents+galvano train+test (6 x-values)

clear all
close all
cc_stat=[]; cc_spont=[]; cc_evoked=[]; cc=[]; cc_shuffled_it=[]; cc_shuff_sub=[];
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data Ch2_data stim2_X stim1_X
 
 global exp_type
exp_type=2; %1-NBES, 2-ChAT
trace_type_input=[1,2]; %[3,2] for exp_type=1; %for exp_type=2 or 3 use [1,2]
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
        files_to_analyze =[76,77,80,82,84,87,90,92,112, 114, 115]; %[92,94,96,98,101]
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
%% set the path for saving figures and variables
if BP50HzLFP_flag==1 && BP50HzVm_flag==1 && BPVm_flag==1 && BPLFP_flag==1
    path_output=['LFP_50Hz+BP Vm_50Hz+BP',num2str(bp_manual_Vm(1)),'-',num2str(bp_manual_Vm(2))];
else if BP50HzLFP_flag==1 && BP50HzVm_flag==1 && BPLFP_flag==1
     path_output='LFP_50Hz+BP Vm_50Hz';  
    else if BP50HzLFP_flag==1 && BP50HzVm_flag==1 && BPVm_flag==1
            path_output='LFP_50Hz Vm_50Hz+BP';  
         else if BP50HzLFP_flag==1 && BP50HzVm_flag==1
                 path_output='LFP_50Hz Vm_50Hz';  
             else if BP50HzLFP_flag==1 
                 path_output='LFP_50Hz';  
                 else path_output='No offline filter';
                 end
             end
        end
    end
end

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
        color_table=[0 0 0; [30,75,14]/256; [136 137 138]/256; [112,172,90]/256; [0 0 204]/256; [255 153 153]/256];
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
        intervals_to_analyze
%%
for t=1:2
%% Subtract mean from the interval (fragment of trace)
data_LFP{t}=Ch2_data_filt(interval(:,t),:,x_value(t));
data_LFP_noDC{t} = fn_Subtract_Mean(data_LFP{t});
data_Vm{t}=current_data_filt(interval(:,t),:,x_value(t));   
data_Vm_noDC{t} = fn_Subtract_Mean(data_Vm{t});
data_Vm_filt{t}=current_data_filt(interval(:,t),:,x_value(t));
data_Vm_filt_noDC{t} = fn_Subtract_Mean(data_Vm_filt{t});
% for plots:
     interval_plot(:,t)=[interval(1,t)-round(add_to_plot*sf{1}):interval(end,t)+round(add_to_plot*sf{1})]; 
    data_Vm_plot{t}=current_data_filt(interval_plot(:,t),:,x_value(t));
    data_LFP_plot{t}=Ch2_data_filt(interval_plot(:,t),:,x_value(t)).*20;

%normalize data (so that a general decrease in the power of the trace will not affect the correlations:
 if norm_flag==1;
      trials=size(current_data,2);    
        for trace = 1:trials
            norm_LFP{t}(:,trace)=data_LFP_noDC{t}(:,trace)./sum((data_LFP_noDC{t}(:,trace)).^2);
            norm_Vm{t}(:,trace)=data_Vm_filt_noDC{t}(:,trace)./sum((data_Vm_filt_noDC{t}(:,trace)).^2);
        end
       data_LFP_noDC{t}=norm_LFP{t};
       data_Vm_filt_noDC{t}=norm_Vm{t};
 end 
       
%% bootstrap
 % bootstrap to calculate shuffeled cc
    trials=size(current_data,2);
    iterations=60;  
    DC_sub=0; %put 1 if the mean trace was not subtracted from the data.
        
    for i=1:iterations
        cc_shuffled_it{fileind,t}(:,i)=fn_bootstrap_cc(data_LFP_noDC{t}(:,:),data_Vm_filt_noDC{t}, 1, DC_sub);
    end
cc_shuffled_mean{fileind}(:,t)=mean(cc_shuffled_it{fileind,t}(:,:),2);
        % Calculate Correlations trace by trace
  for trace = 1:trials
        [cc{fileind,t}(:,trace),lags{fileind,t}(:,1)] = xcorr(data_LFP_noDC{t}(:,trace),data_Vm_filt_noDC{t}(:,trace),'coeff') ; 
%         [cc_bound{fileind,t}(:,trace), lags_bound{fileind,t}(:,1), bounds{fileind,t}(:,trace)] = crosscorr(data_Vm_filt{t}(:,trace),data_LFP{t}(:,trace));
        cc_shuff_sub{fileind,t}(:,trace)=cc{fileind,t}(:,trace)-cc_shuffled_mean{fileind}(:,t);     
  end
  
% %using shuffled-subtracted cc:
cc_shuff_sub_mean{fileind}(:,t) = mean(cc_shuff_sub{fileind,t}(:,:),2);
l_cc=length(cc_shuff_sub_mean{fileind}(:,t));
[c_mean_max_abs_val(t), c_mean_max_r(t)]=max(abs(cc_shuff_sub_mean{fileind}(ceil(l_cc./2)-50.*sf{1}./1000:ceil(l_cc./2)+50.*sf{1}./1000,t))); %look for the maximum in an interval of 50ms around the lag-zero point
% c_mean_max_r(t) = c_mean_max_r(t)+ceil(l_cc./2)-50.*sf{1}./1000;
c_mean_max_val(t)=cc_shuff_sub_mean{fileind}(c_mean_max_r(1)+ceil(l_cc./2)-50.*sf{1}./1000,t);
c_mean_max_val_shuff(t)=cc_shuffled_mean{fileind}(c_mean_max_r(1)+ceil(l_cc./2)-50.*sf{1}./1000,t); %not sure if this is correct to take it
c_mean_max_lag(t)=c_mean_max_r(t)-50.*sf{1}./1000;
cc_mean{fileind}(:,t) = mean(cc{fileind,t}(:,:),2);
%  [c_mean_max_abs_val(t), c_mean_max_r(t)]=max(abs(cc_mean{fileind}(:,t)));
% c_mean_max_val(t)=cc_mean{fileind}(c_mean_max_r(1),t);
c_max(:,t)=cc_shuff_sub{fileind,t}(c_mean_max_r(1),:);
c_lag0(:,t)=cc_shuff_sub{fileind,t}(lags{fileind,1}==0,:);
end
[max_diff_val, max_diff_loc]=max(abs(cc_shuff_sub_mean{fileind}(:,1)-cc_shuff_sub_mean{fileind}(:,2)));
% [max_diff_val, max_diff_loc]=max(abs(cc_mean{fileind}(:,1)-cc_mean{fileind}(:,2)));
c_maxdiff(:,1)=cc_shuff_sub{fileind,1}(max_diff_loc,:);
c_maxdiff(:,2)=cc_shuff_sub{fileind,2}(max_diff_loc,:);
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
ymax=max([data_Vm_plot{1}(:,plot_trace),data_Vm_plot{2}(:,plot_trace)]);
ymin=min([data_Vm_plot{1}(:,plot_trace),data_Vm_plot{2}(:,plot_trace)]);
ydiff=ymax-ymin;
yrange=ceil(1.2*max(ydiff));
ymaxLFP=max([data_LFP_plot{1}(:,plot_trace),data_LFP_plot{2}(:,plot_trace)]);
yminLFP=min([data_LFP_plot{1}(:,plot_trace),data_LFP_plot{2}(:,plot_trace)]);
ydiffLFP=ymaxLFP-yminLFP;
yrangeLFP=ceil(max(ydiffLFP));

for tr_ind=1:length(plot_trace)
subplot(2*length(plot_trace),1,2*tr_ind)
    hold on
        p1=plot(interval_plot(:,1).*dt*1000,data_Vm_plot{1}(:,plot_trace(tr_ind)), 'color',color_table(1,:),'LineWidth',1.5);
                           
        axis tight
     ylim_data=[get(gca,'ylim')]';
     ylim_data(2)=ylim_data(1)+yrange;
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
horiz_vert=1;        lengthh=200;     textit=[num2str(lengthh), ' mS'];     c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=0.03; perc2=[];
        [b1,b2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=lengthh_vert(trace_type);     textit=[num2str(lengthh), ' mV'];     c=[0,0,0];  fonsizes=scalebar_fontsize; %perc1=0.03; perc2=[];
        [b3,b4] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
         if tr_ind~=length(plot_trace)-1;
             delete(b1); delete(b2); delete(b3); delete(b4);
        end
    hold off
    uistack(p1,'top');
    set(gca, 'visible', 'off') ;
    
  subplot(2*length(plot_trace),1,2*tr_ind-1)     
 	hold on
        p2=plot(interval_plot(:,1).*dt*1000,data_LFP_plot{1}(:,plot_trace(tr_ind)),'color',color_table(2,:),'LineWidth',1.5);
%         text(interval_plot(1,1).*dt*1000,data_Vm_plot{1}(1),[num2str(floor(data_Vm_plot{1}(1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')  

        axis tight            
     ylim_data=[get(gca,'ylim')]';
     ylim_data(1)=ylim_data(2)-yrangeLFP;
     xlim_data=[get(gca,'xlim')]';
     set(gca,'ylim',ylim_data, 'xlim',xlim_data)  
     
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
    uistack(p2,'top');
    set(gca, 'visible', 'off') ;

end
    title('Vm-LFP single trace ES Off','FontSize', 16); 
    ylabel('Potential [mV]' ,'FontSize', 14); xlabel('Time [S]','FontSize', 14)
l=legend([p2,p1],{'LFP','Vm'}, 'box', 'off','position',[0.8,0.9,0.2,0.1]); 
    l.FontSize=11;
%     for ix=1:length(l.String)
%       str = l.String{ix};
%       h = findobj(gcf,'DisplayName',str);
%       h.LineWidth =2;
%     end

% pause
% close(gcf);
% end %end of trace for loop

% after ES
       
% ymax=max(data_Vm_plot{2}(:,plot_trace));
% ymin=min(data_Vm_plot{2}(:,plot_trace));
% ydiff=ymax-ymin;
% yrange=ceil(1.4*max(ydiff));
% ymaxLFP=max(data_LFP_plot{2}(:,plot_trace));
% yminLFP=min(data_LFP_plot{2}(:,plot_trace));
% ydiffLFP=ymaxLFP-yminLFP;
% yrangeLFP=ceil(max(ydiffLFP));
          
Fig{fileind}(trace_type+3)=figure;
for tr_ind=1:length(plot_trace)
subplot(2*length(plot_trace),1,2*tr_ind)
    hold on
        p1=plot(interval_plot(:,1).*dt*1000,data_Vm_plot{2}(:,plot_trace(tr_ind)),'color',color_table(5,:),'LineWidth',1.5);
                           
        axis tight
    ylim_data=[get(gca,'ylim')]';
     ylim_data(2)=ylim_data(1)+yrange;
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
            horiz_vert=1;        lengthh=200;     textit=[num2str(lengthh), ' mS'];     c=[0,0,0];  fonsizes=scalebar_fontsize; %perc1=0.03; perc2=[];
                    [b1,b2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
             horiz_vert=0;        lengthh=lengthh_vert(trace_type);     textit=[num2str(lengthh), ' mV'];     c=[0,0,0];  fonsizes=scalebar_fontsize; %perc1=0.03; perc2=[];
                    [b3,b4] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
       if tr_ind~=length(plot_trace)-1;
             delete(b1); delete(b2); delete(b3); delete(b4);
        end
    hold off
    uistack(p1,'top');
    set(gca, 'visible', 'off') ;
    
  subplot(2*length(plot_trace),1,2*tr_ind-1)
 	hold on
       p2=plot(interval_plot(:,1).*dt*1000,data_LFP_plot{2}(:,plot_trace(tr_ind)),'color',color_table(4,:),'LineWidth',1.5);
%         text(interval_plot(1,1).*dt*1000,data_Vm_plot{2}(1),[num2str(floor(data_Vm_plot{2}(1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')  

        axis tight            
     ylim_data=[get(gca,'ylim')]';
     ylim_data(1)=ylim_data(2)-yrangeLFP;
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
            horiz_vert=1;        lengthh=200;     textit=[num2str(lengthh), ' mS'];     c=[0,0,0];  fonsizes=scalebar_fontsize; %perc1=[]; perc2=[];
                [b1,b2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
             horiz_vert=0;        lengthh=10;     textit=[num2str(lengthh./20), ' mV'];     c=[0,0,0];  fonsizes=scalebar_fontsize; %perc1=[]; perc2=[];
                [b3,b4] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
                delete(b1); delete(b2);
                
         if tr_ind~=length(plot_trace);
             delete(b3); delete(b4);
        end
    hold off  
   uistack(p2,'top');
    set(gca, 'visible', 'off') ;
end
title('Vm-LFP single trace ES On','FontSize', 16); 
ylabel('Potential [mV]' ,'FontSize', 14); xlabel('Time [S]','FontSize', 14)
l=legend([p2,p1],{'LFP','Vm'}, 'box', 'off','position',[0.8,0.9,0.2,0.1]); 
l.FontSize=11;
% for ix=1:length(l.String)
%   str = l.String{ix};
%   h = findobj(gcf,'DisplayName',str);
%   h.LineWidth =1.5;
% end

%% clear prcntile1_off prcntile2_off ci_off patch_xdata patch_ydata
%plotting the crosscorrelation for a single trace+the mean
%calculating the critical t.
alpha=0.05;
df=size(cc_shuff_sub{fileind,1}(:,:),2)-1;
alphaup = 1-alpha/2;
alphalow = alpha/2;
tupp = tinv(alphaup,df);
tlow = tinv(alphalow,df); 
%  [prcntile1_off, prcntile2_off]=fn_get_CI_w_bootstrap(cc_shuff_sub{fileind,1}(:,:),0,5000);
%  [prcntile1_on, prcntile2_on]=fn_get_CI_w_bootstrap(cc_shuff_sub{fileind,2}(:,:),0,5000);

% %for CI:
%  cc_shuff_sub_sem1=tupp.*std(cc_shuff_sub{fileind,1}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,1}(:,:),2)); %this is ci and not SEM
%  cc_shuff_sub_sem2=tupp.*std(cc_shuff_sub{fileind,2}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,2}(:,:),2));
%  cc_shuffled_sem1=tupp.*std(cc_shuffled_it{fileind,1}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,1}(:,:),2));
%  cc_shuffled_sem2=tupp.*std(cc_shuffled_it{fileind,2}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,2}(:,:),2));

 %for SEM
 cc_shuff_sub_sem1=std(cc_shuff_sub{fileind,1}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,1}(:,:),2));
 cc_shuff_sub_sem2=std(cc_shuff_sub{fileind,2}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,2}(:,:),2));
 cc_shuffled_sem1=std(cc_shuffled_it{fileind,1}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,1}(:,:),2));
 cc_shuffled_sem2=std(cc_shuffled_it{fileind,2}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,2}(:,:),2));

 Fig{fileind}(trace_type+6)=figure;
 hold on
% sem_xdata(:,1)=[lags{1}(:); fliplr(lags{1}(:))];
% sem_ydata(:,1)=[prcntile1_off; fliplr(prcntile2_off)];
% sem_cdata=ones(size(sem_xdata,1),1);
% sem_xdata(:,2)=[lags{2}(:); fliplr(lags{2}(:))];
% sem_ydata(:,2)=[prcntile1_on; fliplr(prcntile2_on)];
%plotting mean+bootstrap CI:
% plot( lags{fileind,1}.*dt,cc_shuff_sub_mean{fileind}(:,1), 'color',color_table(1,:), 'LineWidth',1.5)
% patch(sem_xdata(:,1).*dt,sem_ydata(:,1),sem_cdata,'faceColor',color_table(3,:),'edgecolor','none','faceAlpha', 0.3)
% plot( lags{fileind,2}.*dt,cc_shuff_sub_mean{fileind}(:,2), 'color',color_table(5,:), 'LineWidth',1.5)
% patch(sem_xdata(:,2).*dt,sem_ydata(:,2),sem_cdata,'faceColor',color_table(5,:),'edgecolor','none','faceAlpha', 0.3)
%alternative - SEM:
h1=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuff_sub_mean{fileind}(:,1),cc_shuff_sub_sem1,{'LineWidth',1.5,'color', color_table(1,:)},1);
h2=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuff_sub_mean{fileind}(:,2),cc_shuff_sub_sem2,{'LineWidth',1.5,'color', color_table(5,:)},1);
h3=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuffled_mean{fileind}(:,1),cc_shuffled_sem1,{'LineWidth',1.5,'linestyle','--','color', color_table(1,:)},1);
h4=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuffled_mean{fileind}(:,2),cc_shuffled_sem2,{'LineWidth',1.5,'linestyle','--','color', color_table(5,:)},1);

ylim=get(gca,'ylim');
ylim(2)=0.8;
yticks=[-0.8, -0.4, 0, 0.4, 0.8];
hline_zero=line([0 0],[ylim(1) ylim(2)],'linestyle','-.','color',[136 137 138]./256,'linewidth',1);
%     plot( lags{fileind,1}.*dt,cc_shuffled_mean{fileind}(:,1),'k-.', 'LineWidth',1.5) %shuffled correlations that were subtracted
%     plot( lags{fileind,1}.*dt,cc_shuffled_mean{fileind}(:,2), 'b-.','LineWidth',1.5)

% title('Vm-LFP crosscorrelation','FontSize', 16); 
l=legend([h1.mainLine h2.mainLine h3.mainLine h4.mainLine],{'NB-','NB+','NB- shuffled', 'NB+ shuffled'},'fontsize',12,'Location','northeast', 'box', 'off');
l.FontSize=11;
% legend('ES off mean cc','ES on mean cc'); legend('boxoff');%legend('boxoff','fontsize',6,'linewidth',1.5,'Position',[0.39,0.45,0.25,0.1]);
xlabel('Lags (s)' ,'FontSize', 14); ylabel('CC' ,'FontSize', 14);
set(gca,'xlim',[-0.5 0.5],'xtick',[-0.5 0 0.5],'ytick',yticks,'fontsize',14); %set(gca,'xlim',[-0.2 0.2]); set(gca,'xtick',[-0.2 0 0.2]);
hold off
% axis tight 

% plot( lags{1,1}.*dt,cc{1}(:,plot_trace),'k-.', 'LineWidth',1)
% plot( lags{1,1}.*dt,cc{2}(:,plot_trace), 'b-.','LineWidth',1)
% plot( lags{fileind,1}.*dt,cc_mean{fileind}(:,1), 'k-', 'LineWidth',1.5)
% plot( lags{fileind,1}.*dt,cc_mean{fileind}(:,2), 'color',[13 49 133]./256, 'LineWidth',1.5)
% plot( lags{fileind,1}.*dt,cc_shuff_sub_mean{fileind}(:,2), 'color',color_table(2,:), 'LineWidth',1.5)
%  h1=scatter(c_mean_max_lag(1).*dt,c_mean_max_val(1),'y','fill');
% h2=scatter(c_mean_max_lag(1).*dt,c_mean_max_val(2),'y','fill');
%% Plot actual and shuffled data on separate figures

Fig{fileind}(trace_type+9)=figure;
 hold on
 h1=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuff_sub_mean{fileind}(:,1),cc_shuff_sub_sem1,{'LineWidth',1.5,'color', color_table(1,:)},1);
h2=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuff_sub_mean{fileind}(:,2),cc_shuff_sub_sem2,{'LineWidth',1.5,'color', color_table(5,:)},1);
ylim=get(gca,'ylim');
% ylim=[-0.4 0.4];
% yticks=[-0.4, -0.2 0, 0.2];
%  if trace_type==2
% ylim=[-0.8 0.8];
% yticks=[-0.8, -0.4, 0, 0.4, 0.8];
%  end
% hline_zero=line([0 0],[ylim(1) ylim(2)],'linestyle','-.','color',[136 137 138]./256,'linewidth',1);
% title('Vm-LFP crosscorrelation','FontSize', 16); 
l=legend([h1.mainLine h2.mainLine ],legend_string,'Location','northeast', 'box', 'off');
l.FontSize=11;
for ix=1:length(l.String)
  str = l.String{ix};
  h = findobj(gcf,'DisplayName',str);
  h.LineWidth =1.5;
end
% legend('ES off mean cc','ES on mean cc'); legend('boxoff');%legend('boxoff','fontsize',6,'linewidth',1.5,'Position',[0.39,0.45,0.25,0.1]);
xlabel('Lags (s)' ,'FontSize', 14); ylabel('CC' ,'FontSize', 14);
set(gca,'xlim',[-0.5 0.5],'xtick',[-0.5 0 0.5],'fontsize',14); 
if trace_type~=2
    set(gca,'ylim',ylim);
end
hold off

Fig{fileind}(trace_type+12)=figure;
 hold on
h3=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuffled_mean{fileind}(:,1),cc_shuffled_sem1,{'LineWidth',1.5,'linestyle','--','color', color_table(1,:)},1);
h4=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuffled_mean{fileind}(:,2),cc_shuffled_sem2,{'LineWidth',1.5,'linestyle','--','color', color_table(5,:)},1);
%  ylim=get(gca,'ylim');
% ylim=[-0.2, 0.2];
%  yticks=[-0.2,0,0.2];
%  if trace_type==2
% ylim=[-0.8 0.8];
% yticks=[-0.8, -0.4, 0, 0.4, 0.8];
%  end
% hline_zero=line([0 0],[ylim(1) ylim(2)],'linestyle','-.','color',[136 137 138]./256,'linewidth',1);
% title('Vm-LFP crosscorrelation','FontSize', 16); 
l=legend([h3.mainLine h4.mainLine ],legend_string_shuff,'Location','northeast', 'box', 'off');
l.FontSize=11;
for ix=1:length(l.String)
  str = l.String{ix};
  h = findobj(gcf,'DisplayName',str);
  h.LineWidth =1.5;
end
xlabel('Lags (s)' ,'FontSize', 14); ylabel('CC' ,'FontSize', 14);
set(gca,'xlim',[-0.5 0.5],'xtick',[-0.5 0 0.5],'fontsize',14);
if trace_type~=2
    set(gca,'ylim',ylim);
end
hold off
% pause
        %% Save figures    
        cd(path_output)
   if save_flag==1
       if trace_type==1 || trace_type==3         
                saveas(Fig{fileind}(trace_type),['Vm-LFP_spont_stim1_Off_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'fig') 
                print(Fig{fileind}(trace_type),['Vm-LFP_spont_stim1_Off_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'-dpng','-r600','-opengl') 
                saveas(Fig{fileind}(trace_type+3),['Vm-LFP_spont_stim1_On_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'fig') 
                print(Fig{fileind}(trace_type+3),['Vm-LFP_spont_stim1_On_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'-dpng','-r600','-opengl') 
                saveas(Fig{fileind}(trace_type+6),['Vm-LFPcc_spont_f' num2str(files_to_analyze(fileind))],'fig') 
                print(Fig{fileind}(trace_type+6),['Vm-LFPcc_spont_f' num2str(files_to_analyze(fileind))],'-dpng','-r600','-opengl') 
                saveas(Fig{fileind}(trace_type+9),['Vm-LFPcc_spont_f' num2str(files_to_analyze(fileind)), 'actual_data'],'fig') 
                print(Fig{fileind}(trace_type+9),['Vm-LFPcc_spont_f' num2str(files_to_analyze(fileind)), 'actual_data'],'-dpng','-r600','-opengl') 
                saveas(Fig{fileind}(trace_type+12),['Vm-LFPcc_spont_f' num2str(files_to_analyze(fileind)), 'shuffled_data'],'fig') 
                print(Fig{fileind}(trace_type+12),['Vm-LFPcc_spont_f' num2str(files_to_analyze(fileind)), 'shuffled_data'],'-dpng','-r600','-opengl') 
       end

        if trace_type==2
                saveas(Fig{fileind}(trace_type),['Vm-LFP_evoked_stim1_Off_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'fig') 
                print(Fig{fileind}(trace_type),['Vm-LFP_evoked_stim1_Off_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'-dpng','-r600','-opengl') 
                saveas(Fig{fileind}(trace_type+3),['Vm-LFP_evoked_stim1_On_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'fig') 
                print(Fig{fileind}(trace_type+3),['Vm-LFP_evoked_stim1_On_f' num2str(files_to_analyze(fileind)),'_t', num2str(plot_trace)],'-dpng','-r600','-opengl') 
                saveas(Fig{fileind}(trace_type+6),['Vm-LFPcc_evoked_f' num2str(files_to_analyze(fileind))],'fig') 
                print(Fig{fileind}(trace_type+6),['Vm-LFPcc_evoked_f' num2str(files_to_analyze(fileind))],'-dpng','-r600','-opengl') 
                saveas(Fig{fileind}(trace_type+9),['Vm-LFPcc_evoked_f' num2str(files_to_analyze(fileind)), 'actual_data'],'fig') 
                print(Fig{fileind}(trace_type+9),['Vm-LFPcc_evoked_f' num2str(files_to_analyze(fileind)), 'actual_data'],'-dpng','-r600','-opengl') 
                saveas(Fig{fileind}(trace_type+12),['Vm-LFPcc_evoked_f' num2str(files_to_analyze(fileind)), 'shuffled_data'],'fig') 
                print(Fig{fileind}(trace_type+12),['Vm-LFPcc_evoked_f' num2str(files_to_analyze(fileind)), 'shuffled_data'],'-dpng','-r600','-opengl') 
                % saveas(Fig10,['Vm-LFPcc_ongoing+evoked_file_' num2str(files_to_analyze(fileind)) '.fig']) 
                % print(Fig10,['Vm-LFPcc_ongoing+evoked_file_' num2str(files_to_analyze(fileind))],'-dpng','-r600','-opengl') 
        end
    end
    end
%%

    if trace_type==1 || trace_type==3
        
        cc_spont(fileind).fname = fname;
        cc_spont(fileind).win_duration = duration;
        cc_spont(fileind).trace_type_input = trace_type_input;
        cc_spont(fileind).analyze_time_before_train = analyze_time_before_train;
        cc_spont(fileind).V = analyze_train_only_flag;
        cc_spont(fileind).BP50HzLFP=BP50HzLFP_flag; %removing 50Hz noise from LFP signal
        cc_spont(fileind).BP50HzVm=BP50HzVm_flag; %removing 50Hz noise from Vm signal
        cc_spont(fileind).BPLFP_flag=BPLFP_flag; %filtering LFP 
        cc_spont(fileind).BPLFP=bp_filt_LFP; %the BP frequency filter for  LFP, used if BPLFP_flag=1
        cc_spont(fileind).BPVm_flag=BPVm_flag; %filtering Vm 
        cc_spont(fileind).BPVm=bp_filt_Vm; %the BP frequency filter for Vm, used if BPVm_flag=1
        cc_spont(fileind).cc_max=c_mean_max_val; %value of maximal correlation 
        cc_spont(fileind).cc_max_lag=lags{fileind,1}(c_mean_max_r)';         
        cc_spont(fileind).cc_max_time=lags{fileind,1}(c_mean_max_r).*dt;
        cc_spont(fileind).cc_max_shuff=c_mean_max_val_shuff; %value of cc_shuffled at the maximal correlation point at cc_shuff_sub
        cc_spont(fileind).cc_lag0_shuff=cc_shuffled_mean{fileind}(lags{fileind,1}==0,:); %mean value of correlation at lag zero
        cc_spont(fileind).cc_lag0=cc_shuff_sub_mean{fileind}(lags{fileind,1}==0,:); %mean value of correlation at lag zero
%         cc_spont(fileind).cc_maxdiff=cc_mean{fileind}(max_diff_loc,:); %location (row) of maximal difference in correlation (the lag)
        cc_spont(fileind).cc_maxdiff=cc_shuff_sub_mean{fileind}(max_diff_loc,:); %location (row) of maximal difference in correlation (the lag)
        cc_spont(fileind).cc_max_trace = c_max; %location (row) of maximal correlation (the lag)
        cc_spont(fileind).cc_lag0_trace = c_lag0; %mean values of correlation at lag zero (for all traces)
        cc_spont(fileind).cc_maxdiff_trace = c_maxdiff; %mean values of correlation at largest difference point (for all traces)

         cc_lag0_mat{trace_type}(fileind,:)=cc_spont(fileind).cc_lag0;
         cc_lag0_shuff_mat{trace_type}(fileind,:)=cc_spont(fileind).cc_lag0_shuff;
         cc_max_mat{trace_type}(fileind,:)=cc_spont(fileind).cc_max;
         cc_max_shuff_mat{trace_type}(fileind,:)=cc_spont(fileind).cc_max_shuff;
         cc_max_time_mat{trace_type}(fileind,:)=cc_spont(fileind).cc_max_time;
         cc_maxdiff_mat{trace_type}(fileind,:)=cc_spont(fileind).cc_maxdiff;
    end
    if trace_type==2
        cc_evoked(fileind).fname = fname;
        cc_evoked(fileind).win_duration = duration;
        cc_evoked(fileind).trace_type_input = trace_type_input;
        cc_evoked(fileind).analyze_time_before_train = analyze_time_before_train;
        cc_evoked(fileind).V = analyze_train_only_flag;
        cc_evoked(fileind).BP50HzLFP=BP50HzLFP_flag; %removing 50Hz noise from LFP signal
        cc_evoked(fileind).BP50HzVm=BP50HzVm_flag; %removing 50Hz noise from Vm signal
        cc_evoked(fileind).BPLFP_flag=BPLFP_flag; %filtering LFP 
        cc_evoked(fileind).BPLFP=bp_filt_LFP; %the BP frequency filter for  LFP, used if BPLFP_flag=1
        cc_evoked(fileind).BPVm_flag=BPVm_flag; %filtering Vm 
        cc_evoked(fileind).BPVm=bp_filt_Vm; %the BP frequency filter for Vm, used if BPVm_flag=1
        cc_evoked(fileind).cc_max=c_mean_max_val; %value of maximal correlation
        cc_evoked(fileind).cc_max_shuff=c_mean_max_val_shuff; %value of cc_shuffled at the maximal correlation point at cc_shuff_sub
        cc_evoked(fileind).cc_max_lag=lags{fileind,1}(c_mean_max_r);
        cc_evoked(fileind).cc_max_time=lags{fileind,1}(c_mean_max_r).*dt;
        cc_evoked(fileind).cc_lag0_shuff=cc_shuffled_mean{fileind}(lags{fileind,1}==0,:); %mean value of correlation at lag zero
        cc_evoked(fileind).cc_lag0=cc_shuff_sub_mean{fileind}(lags{fileind,1}==0,:); %mean value of correlation at lag zero
%         cc_evoked(fileind).cc_maxdiff=cc_mean{fileind}(max_diff_loc,:); %location (row) of maximal difference in correlation (the lag)
        cc_evoked(fileind).cc_maxdiff=cc_shuff_sub_mean{fileind}(max_diff_loc,:); %location (row) of maximal difference in correlation (the lag)
        cc_evoked(fileind).cc_max_trace = c_max; %location (row) of maximal correlation (the lag)
        cc_evoked(fileind).cc_lag0_trace = c_lag0; %mean values of correlation at lag zero (for all traces)
        cc_evoked(fileind).cc_maxdiff_trace = c_maxdiff; %mean values of correlation at largest difference point (for all traces)

         cc_lag0_mat{trace_type}(fileind,:)=cc_evoked(fileind).cc_lag0;
         cc_lag0_shuff_mat{trace_type}(fileind,:)=cc_evoked(fileind).cc_lag0_shuff;
         cc_max_mat{trace_type}(fileind,:)=cc_evoked(fileind).cc_max;
         cc_max_shuff_mat{trace_type}(fileind,:)=cc_evoked(fileind).cc_max_shuff;
         cc_max_time_mat{trace_type}(fileind,:)=cc_evoked(fileind).cc_max_time;
         cc_maxdiff_mat{trace_type}(fileind,:)=cc_evoked(fileind).cc_maxdiff;         
    end
clear c_max c_lag0 temp_mean c_mean_max_r c_mean_max_c...% c_mean_max_val...
    max_diff_val max_diff_loc c_maxdiff cc_shuff_sub cc_shuffled_mean  cc_shuffled_it cc lags cc_mean cc_shuff_sub_mean

end  
% close all
    end
%% Statistics for cross-correlation
%t-tests between cells  
%Spontaneous
        %normalize cc max values and perform ttest   
        cc_stat.spont.max_val=cc_max_mat{trace_type_input(1)}(:,:);  
        cc_stat.spont.max_val_sign=sign(cc_stat.spont.max_val);
        cc_stat.spont.max_val_sign_change=cc_stat.spont.max_val_sign(:,1).*cc_stat.spont.max_val_sign(:,2)<0;
        cc_stat.spont.max_val=abs(cc_stat.spont.max_val);
        cc_stat.spont.max_val_m=mean(cc_stat.spont.max_val,1);
        cc_stat.spont.max_val_std=std(cc_stat.spont.max_val,0,1);
        cc_stat.spont.change_max_val=[(cc_max_mat{trace_type_input(1)}(:,2)-cc_max_mat{trace_type_input(1)}(:,1))./abs(cc_max_mat{trace_type_input(1)}(:,1))].*100; %percent change
        cc_stat.spont.change_max_val_m=mean(cc_stat.spont.change_max_val,1);
        cc_stat.spont.change_max_val_std=std(cc_stat.spont.change_max_val,0,1);
        change_max_val_mat= cc_stat.spont.change_max_val;
    if length(files_to_analyze)>=5
        %testing for normal distribution       
        [cc_stat.spont.lillietest_h_max_val, cc_stat.spont.lillietest_p_max_val] = lillietest(cc_stat.spont.max_val(:,2)- cc_stat.spont.max_val(:,1));
        [cc_stat.spont.lillietest_h_change_max_val, cc_stat.spont.lillietest_p_change_max_val] = lillietest( cc_stat.spont.change_max_val);
    end
        %paired ttest 
        [cc_stat.spont.ttest_h_max_val, cc_stat.spont.ttest_p_max_val]= ttest(cc_stat.spont.max_val(:,1),cc_stat.spont.max_val(:,2));   
        [cc_stat.spont.wilcoxon_p_max_val, cc_stat.spont.wilcoxon_h_max_val]= signrank(cc_stat.spont.max_val(:,1),cc_stat.spont.max_val(:,2));
        [cc_stat.spont.ttest_h_change_max_val, cc_stat.spont.ttest_p_change_max_val]= ttest(cc_stat.spont.change_max_val);
        [cc_stat.spont.wilcoxon_p_change_max_val, cc_stat.spont.wilcoxon_h_change_max_val]= signrank(cc_stat.spont.change_max_val);
  
        %normalize cc max values and perform ttest   
        cc_stat.spont.max_val_shuff=cc_max_shuff_mat{trace_type_input(1)}(:,:);  
        cc_stat.spont.max_val_shuff_sign=sign(cc_stat.spont.max_val_shuff);
        cc_stat.spont.max_val_shuff_sign_change=cc_stat.spont.max_val_shuff_sign(:,1).*cc_stat.spont.max_val_shuff_sign(:,2)<0;
        cc_stat.spont.max_val_shuff=abs(cc_stat.spont.max_val_shuff);
        cc_stat.spont.max_val_shuff_m=mean(cc_stat.spont.max_val_shuff,1);
        cc_stat.spont.max_val_shuff_std=std(cc_stat.spont.max_val_shuff,0,1);
        cc_stat.spont.change_max_val_shuff=[(cc_max_shuff_mat{trace_type_input(1)}(:,2)-cc_max_shuff_mat{trace_type_input(1)}(:,1))./abs(cc_max_shuff_mat{trace_type_input(1)}(:,1))].*100; %percent change
        cc_stat.spont.change_max_val_shuff_m=mean(cc_stat.spont.change_max_val_shuff,1);
        cc_stat.spont.change_max_val_shuff_std=std(cc_stat.spont.change_max_val_shuff,0,1);
        change_max_val_shuff_mat= cc_stat.spont.change_max_val_shuff;
        %testing for normal distribution     
    if length(files_to_analyze)>=5
        [cc_stat.spont.lillietest_h_max_val_shuff, cc_stat.spont.lillietest_p_max_val_shuff] = lillietest(cc_stat.spont.max_val_shuff(:,2)- cc_stat.spont.max_val_shuff(:,1));
        [cc_stat.spont.lillietest_h_change_max_val_shuff, cc_stat.spont.lillietest_p_change_max_val_shuff] = lillietest( cc_stat.spont.change_max_val_shuff);
    end
        %paired ttest 
        [cc_stat.spont.ttest_h_max_val_shuff, cc_stat.spont.ttest_p_max_val_shuff]= ttest(cc_stat.spont.max_val_shuff(:,1),cc_stat.spont.max_val_shuff(:,2));   
        [cc_stat.spont.wilcoxon_p_max_val_shuff, cc_stat.spont.wilcoxon_h_max_val_shuff]= signrank(cc_stat.spont.max_val_shuff(:,1),cc_stat.spont.max_val_shuff(:,2));
        [cc_stat.spont.ttest_h_change_max_val_shuff, cc_stat.spont.ttest_p_change_max_val_shuff]= ttest(cc_stat.spont.change_max_val_shuff);
        [cc_stat.spont.wilcoxon_p_change_max_val_shuff, cc_stat.spont.wilcoxon_h_change_max_val_shuff]= signrank(cc_stat.spont.change_max_val_shuff);

         %normalize cc zero-lag values and perform ttest
        cc_stat.spont.lag0=cc_lag0_mat{trace_type_input(1)}(:,:);
        cc_stat.spont.lag0_sign=sign(cc_stat.spont.lag0);
        cc_stat.spont.lag0_sign_change=cc_stat.spont.lag0_sign(:,1).*cc_stat.spont.lag0_sign(:,2)<0;
        cc_stat.spont.lag0=abs(cc_stat.spont.lag0);
        cc_stat.spont.lag0_m=mean(cc_stat.spont.lag0,1);
        cc_stat.spont.lag0_std=std(cc_stat.spont.lag0,0,1);
        cc_stat.spont.change_lag0=[(cc_lag0_mat{trace_type_input(1)}(:,2)-cc_lag0_mat{trace_type_input(1)}(:,1))./abs(cc_lag0_mat{trace_type_input(1)}(:,1))].*100; %percent change
        cc_stat.spont.change_lag0_m=mean(cc_stat.spont.change_lag0,1);
        cc_stat.spont.change_lag0_std=std(cc_stat.spont.change_lag0,0,1);
        change_lag0_mat= cc_stat.spont.change_lag0;
        %testing for normal distribution       
    if length(files_to_analyze)>=5
        [cc_stat.spont.lillietest_h_lag0, cc_stat.spont.lillietest_p_lag0] = lillietest(cc_stat.spont.lag0(:,2)- cc_stat.spont.lag0(:,1));
        [cc_stat.spont.lillietest_h_change_lag0, cc_stat.spont.lillietest_p_change_lag0] = lillietest( cc_stat.spont.change_lag0);
    end
        %paired ttest 
        [cc_stat.spont.ttest_h_lag0, cc_stat.spont.ttest_p_lag0]= ttest(cc_stat.spont.lag0(:,1),cc_stat.spont.lag0(:,2));   
        [cc_stat.spont.wilcoxon_p_lag0, cc_stat.spont.wilcoxon_h_lag0]= signrank(cc_stat.spont.lag0(:,1),cc_stat.spont.lag0(:,2));
        [cc_stat.spont.ttest_h_change_lag0, cc_stat.spont.ttest_p_change_lag0]= ttest(cc_stat.spont.change_lag0);
        [cc_stat.spont.wilcoxon_p_change_lag0, cc_stat.spont.wilcoxon_h_change_lag0]= signrank(cc_stat.spont.change_lag0);
 
          %normalize cc zero-lag shuffled values and perform ttest
        cc_stat.spont.lag0_shuff=cc_lag0_shuff_mat{trace_type_input(1)}(:,:);
        cc_stat.spont.lag0_shuff_sign=sign(cc_stat.spont.lag0_shuff);
        cc_stat.spont.lag0_shuff_sign_change=cc_stat.spont.lag0_shuff_sign(:,1).*cc_stat.spont.lag0_shuff_sign(:,2)<0;
        cc_stat.spont.lag0_shuff=abs(cc_stat.spont.lag0_shuff);
        cc_stat.spont.lag0_shuff_m=mean(cc_stat.spont.lag0_shuff,1);
        cc_stat.spont.lag0_shuff_std=std(cc_stat.spont.lag0_shuff,0,1);
        cc_stat.spont.change_lag0_shuff=[(cc_lag0_shuff_mat{trace_type_input(1)}(:,2)-cc_lag0_shuff_mat{trace_type_input(1)}(:,1))./abs(cc_lag0_shuff_mat{trace_type_input(1)}(:,1))].*100; %percent change
        cc_stat.spont.change_lag0_shuff_m=mean(cc_stat.spont.change_lag0_shuff,1);
        cc_stat.spont.change_lag0_shuff_std=std(cc_stat.spont.change_lag0_shuff,0,1);
        change_lag0_shuff_mat= cc_stat.spont.change_lag0_shuff;
        %testing for normal distribution       
     if length(files_to_analyze)>=5
        [cc_stat.spont.lillietest_h_lag0_shuff, cc_stat.spont.lillietest_p_lag0_shuff] = lillietest(cc_stat.spont.lag0_shuff(:,2)- cc_stat.spont.lag0_shuff(:,1));
        [cc_stat.spont.lillietest_h_change_lag0_shuff, cc_stat.spont.lillietest_p_change_lag0_shuff] = lillietest( cc_stat.spont.change_lag0_shuff);
     end
        %paired ttest 
        [cc_stat.spont.ttest_h_lag0_shuff, cc_stat.spont.ttest_p_lag0_shuff]= ttest(cc_stat.spont.lag0_shuff(:,1),cc_stat.spont.lag0_shuff(:,2));   
        [cc_stat.spont.wilcoxon_p_lag0_shuff, cc_stat.spont.wilcoxon_h_lag0_shuff]= signrank(cc_stat.spont.lag0_shuff(:,1),cc_stat.spont.lag0_shuff(:,2));
        [cc_stat.spont.ttest_h_change_lag0_shuff, cc_stat.spont.ttest_p_change_lag0_shuff]= ttest(cc_stat.spont.change_lag0_shuff);
        [cc_stat.spont.wilcoxon_p_change_lag0_shuff, cc_stat.spont.wilcoxon_h_change_lag0_shuff]= signrank(cc_stat.spont.change_lag0_shuff);
        
         %normalize cc maxdiff values and perform ttest
        cc_stat.spont.maxdiff=cc_maxdiff_mat{trace_type_input(1)}(:,:);
        cc_stat.spont.maxdiff_sign=sign(cc_stat.spont.maxdiff);
        cc_stat.spont.maxdiff_sign_change=cc_stat.spont.maxdiff_sign(:,1).*cc_stat.spont.maxdiff_sign(:,2)<0;
        cc_stat.spont.maxdiff=abs(cc_stat.spont.maxdiff);
        cc_stat.spont.maxdiff_m=mean(cc_stat.spont.maxdiff,1);
        cc_stat.spont.maxdiff_std=std(cc_stat.spont.maxdiff,0,1);
        cc_stat.spont.change_maxdiff=[(cc_maxdiff_mat{trace_type_input(1)}(:,2)-cc_maxdiff_mat{trace_type_input(1)}(:,1))./abs(cc_maxdiff_mat{trace_type_input(1)}(:,1))].*100; %percent change
        cc_stat.spont.change_maxdiff_m=mean(cc_stat.spont.change_maxdiff,1);
        cc_stat.spont.change_maxdiff_std=std(cc_stat.spont.change_maxdiff,0,1);
        change_maxdiff_mat= cc_stat.spont.change_maxdiff;
        %testing for normal distribution       
     if length(files_to_analyze)>=5
        [cc_stat.spont.lillietest_h_maxdiff, cc_stat.spont.lillietest_p_maxdiff] = lillietest(cc_stat.spont.maxdiff(:,2)- cc_stat.spont.maxdiff(:,1));
        [cc_stat.spont.lillietest_h_change_maxdiff, cc_stat.spont.lillietest_p_change_maxdiff] = lillietest( cc_stat.spont.change_maxdiff);
     end
        %paired ttest 
        [cc_stat.spont.ttest_h_maxdiff, cc_stat.spont.ttest_p_maxdiff]= ttest(cc_stat.spont.maxdiff(:,1),cc_stat.spont.maxdiff(:,2));   
        [cc_stat.spont.wilcoxon_p_maxdiff, cc_stat.spont.wilcoxon_h_maxdiff]= signrank(cc_stat.spont.maxdiff(:,1),cc_stat.spont.maxdiff(:,2));
        [cc_stat.spont.ttest_h_change_maxdiff, cc_stat.spont.ttest_p_change_maxdiff]= ttest(cc_stat.spont.change_maxdiff);
        [cc_stat.spont.wilcoxon_p_change_maxdiff, cc_stat.spont.wilcoxon_h_change_maxdiff]= signrank(cc_stat.spont.change_maxdiff);
        
%   Evoked
        %normalize cc max values and perform ttest   
        cc_stat.evoked.max_val=cc_max_mat{2}(:,:);
        cc_stat.evoked.max_val_sign=sign(cc_stat.evoked.max_val);
        cc_stat.evoked.max_val_sign_change=cc_stat.evoked.max_val_sign(:,1).*cc_stat.evoked.max_val_sign(:,2)<0;
        cc_stat.evoked.max_val=abs(cc_stat.evoked.max_val);
        cc_stat.evoked.max_val_m=mean(cc_stat.evoked.max_val,1);
        cc_stat.evoked.max_val_std=std(cc_stat.evoked.max_val,0,1);
        cc_stat.evoked.change_max_val=[(cc_max_mat{2}(:,2)-cc_max_mat{2}(:,1))./abs(cc_max_mat{2}(:,1))].*100; %percent change
        cc_stat.evoked.change_max_val_m=mean(cc_stat.evoked.change_max_val,1);
        cc_stat.evoked.change_max_val_std=std(cc_stat.evoked.change_max_val,0,1);
        change_max_val_mat= cc_stat.evoked.change_max_val;
        %testing for normal distribution       
     if length(files_to_analyze)>=5
        [cc_stat.evoked.lillietest_h_max_val, cc_stat.evoked.lillietest_p_max_val] = lillietest(cc_stat.evoked.max_val(:,2)- cc_stat.evoked.max_val(:,1));
        [cc_stat.evoked.lillietest_h_change_max_val, cc_stat.evoked.lillietest_p_change_max_val] = lillietest( cc_stat.evoked.change_max_val);
     end
        %paired ttest 
        [cc_stat.evoked.ttest_h_max_val, cc_stat.evoked.ttest_p_max_val]= ttest(cc_stat.evoked.max_val(:,1),cc_stat.evoked.max_val(:,2));   
        [cc_stat.evoked.wilcoxon_p_max_val, cc_stat.evoked.wilcoxon_h_max_val]= signrank(cc_stat.evoked.max_val(:,1),cc_stat.evoked.max_val(:,2));
        [cc_stat.evoked.ttest_h_change_max_val, cc_stat.evoked.ttest_p_change_max_val]= ttest(cc_stat.evoked.change_max_val);
        [cc_stat.evoked.wilcoxon_p_change_max_val, cc_stat.evoked.wilcoxon_h_change_max_val]= signrank(cc_stat.evoked.change_max_val);
  
        %normalize cc max values and perform ttest   
        cc_stat.evoked.max_val_shuff=cc_max_shuff_mat{2}(:,:);  
        cc_stat.evoked.max_val_shuff_sign=sign(cc_stat.evoked.max_val_shuff);
        cc_stat.evoked.max_val_shuff_sign_change=cc_stat.evoked.max_val_shuff_sign(:,1).*cc_stat.evoked.max_val_shuff_sign(:,2)<0;
        cc_stat.evoked.max_val_shuff=abs(cc_stat.evoked.max_val_shuff);
        cc_stat.evoked.max_val_shuff_m=mean(cc_stat.evoked.max_val_shuff,1);
        cc_stat.evoked.max_val_shuff_std=std(cc_stat.evoked.max_val_shuff,0,1);
        cc_stat.evoked.change_max_val_shuff=[(cc_max_shuff_mat{2}(:,2)-cc_max_shuff_mat{2}(:,1))./abs(cc_max_shuff_mat{2}(:,1))].*100; %percent change
        cc_stat.evoked.change_max_val_shuff_m=mean(cc_stat.evoked.change_max_val_shuff,1);
        cc_stat.evoked.change_max_val_shuff_std=std(cc_stat.evoked.change_max_val_shuff,0,1);
        change_max_val_shuff_mat= cc_stat.evoked.change_max_val_shuff;
        %testing for normal distribution       
     if length(files_to_analyze)>=5
        [cc_stat.evoked.lillietest_h_max_val_shuff, cc_stat.evoked.lillietest_p_max_val_shuff] = lillietest(cc_stat.evoked.max_val_shuff(:,2)- cc_stat.evoked.max_val_shuff(:,1));
        [cc_stat.evoked.lillietest_h_change_max_val_shuff, cc_stat.evoked.lillietest_p_change_max_val_shuff] = lillietest( cc_stat.evoked.change_max_val_shuff);
     end
        %paired ttest 
        [cc_stat.evoked.ttest_h_max_val_shuff, cc_stat.evoked.ttest_p_max_val_shuff]= ttest(cc_stat.evoked.max_val_shuff(:,1),cc_stat.evoked.max_val_shuff(:,2));   
        [cc_stat.evoked.wilcoxon_p_max_val_shuff, cc_stat.evoked.wilcoxon_h_max_val_shuff]= signrank(cc_stat.evoked.max_val_shuff(:,1),cc_stat.evoked.max_val_shuff(:,2));
        [cc_stat.evoked.ttest_h_change_max_val_shuff, cc_stat.evoked.ttest_p_change_max_val_shuff]= ttest(cc_stat.evoked.change_max_val_shuff);
        [cc_stat.evoked.wilcoxon_p_change_max_val_shuff, cc_stat.evoked.wilcoxon_h_change_max_val_shuff]= signrank(cc_stat.evoked.change_max_val_shuff);

         %normalize cc zero-lag values and perform ttest
        cc_stat.evoked.lag0=cc_lag0_mat{2}(:,:);
        cc_stat.evoked.lag0_sign=sign(cc_stat.evoked.lag0);
        cc_stat.evoked.lag0_sign_change=cc_stat.evoked.lag0_sign(:,1).*cc_stat.evoked.lag0_sign(:,2)<0;
        cc_stat.evoked.lag0=abs(cc_stat.evoked.lag0);
        cc_stat.evoked.lag0_m=mean(cc_stat.evoked.lag0,1);
        cc_stat.evoked.lag0_std=std(cc_stat.evoked.lag0,0,1);
        cc_stat.evoked.change_lag0=[(cc_lag0_mat{2}(:,2)-cc_lag0_mat{2}(:,1))./abs(cc_lag0_mat{2}(:,1))].*100; %percent change
        cc_stat.evoked.change_lag0_m=mean(cc_stat.evoked.change_lag0,1);
        cc_stat.evoked.change_lag0_std=std(cc_stat.evoked.change_lag0,0,1);
        change_lag0_mat= cc_stat.evoked.change_lag0;
        %testing for normal distribution       
     if length(files_to_analyze)>=5
        [cc_stat.evoked.lillietest_h_lag0, cc_stat.evoked.lillietest_p_lag0] = lillietest(cc_stat.evoked.lag0(:,2)- cc_stat.evoked.lag0(:,1));
        [cc_stat.evoked.lillietest_h_change_lag0, cc_stat.evoked.lillietest_p_change_lag0] = lillietest( cc_stat.evoked.change_lag0);
     end
        %paired ttest 
        [cc_stat.evoked.ttest_h_lag0, cc_stat.evoked.ttest_p_lag0]= ttest(cc_stat.evoked.lag0(:,1),cc_stat.evoked.lag0(:,2));   
        [cc_stat.evoked.wilcoxon_p_lag0, cc_stat.evoked.wilcoxon_h_lag0]= signrank(cc_stat.evoked.lag0(:,1),cc_stat.evoked.lag0(:,2));
        [cc_stat.evoked.ttest_h_change_lag0, cc_stat.evoked.ttest_p_change_lag0]= ttest(cc_stat.evoked.change_lag0);
        [cc_stat.evoked.wilcoxon_p_change_lag0, cc_stat.evoked.wilcoxon_h_change_lag0]= signrank(cc_stat.evoked.change_lag0);
 
          %normalize cc zero-lag shuffled values and perform ttest
        cc_stat.evoked.lag0_shuff=cc_lag0_shuff_mat{2}(:,:);
        cc_stat.evoked.lag0_shuff_sign=sign(cc_stat.evoked.lag0_shuff);
        cc_stat.evoked.lag0_shuff_sign_change=cc_stat.evoked.lag0_shuff_sign(:,1).*cc_stat.evoked.lag0_shuff_sign(:,2)<0;
        cc_stat.evoked.lag0_shuff=abs(cc_stat.evoked.lag0_shuff);
        cc_stat.evoked.lag0_shuff_m=mean(cc_stat.evoked.lag0_shuff,1);
        cc_stat.evoked.lag0_shuff_std=std(cc_stat.evoked.lag0_shuff,0,1);
        cc_stat.evoked.change_lag0_shuff=[(cc_lag0_shuff_mat{2}(:,2)-cc_lag0_shuff_mat{2}(:,1))./abs(cc_lag0_shuff_mat{2}(:,1))].*100; %percent change
        cc_stat.evoked.change_lag0_shuff_m=mean(cc_stat.evoked.change_lag0_shuff,1);
        cc_stat.evoked.change_lag0_shuff_std=std(cc_stat.evoked.change_lag0_shuff,0,1);
        change_lag0_shuff_mat= cc_stat.evoked.change_lag0_shuff;
        %testing for normal distribution       
     if length(files_to_analyze)>=5
        [cc_stat.evoked.lillietest_h_lag0_shuff, cc_stat.evoked.lillietest_p_lag0_shuff] = lillietest(cc_stat.evoked.lag0_shuff(:,2)- cc_stat.evoked.lag0_shuff(:,1));
        [cc_stat.evoked.lillietest_h_change_lag0_shuff, cc_stat.evoked.lillietest_p_change_lag0_shuff] = lillietest( cc_stat.evoked.change_lag0_shuff);
     end
        %paired ttest 
        [cc_stat.evoked.ttest_h_lag0_shuff, cc_stat.evoked.ttest_p_lag0_shuff]= ttest(cc_stat.evoked.lag0_shuff(:,1),cc_stat.evoked.lag0_shuff(:,2));   
        [cc_stat.evoked.wilcoxon_p_lag0_shuff, cc_stat.evoked.wilcoxon_h_lag0_shuff]= signrank(cc_stat.evoked.lag0_shuff(:,1),cc_stat.evoked.lag0_shuff(:,2));
        [cc_stat.evoked.ttest_h_change_lag0_shuff, cc_stat.evoked.ttest_p_change_lag0_shuff]= ttest(cc_stat.evoked.change_lag0_shuff);
        [cc_stat.evoked.wilcoxon_p_change_lag0_shuff, cc_stat.evoked.wilcoxon_h_change_lag0_shuff]= signrank(cc_stat.evoked.change_lag0_shuff);
        
         %normalize cc maxdiff values and perform ttest
        cc_stat.evoked.maxdiff=cc_maxdiff_mat{2}(:,:);
        cc_stat.evoked.maxdiff_sign=sign(cc_stat.evoked.maxdiff);
        cc_stat.evoked.maxdiff_sign_change=cc_stat.evoked.maxdiff_sign(:,1).*cc_stat.evoked.maxdiff_sign(:,2)<0;
        cc_stat.evoked.maxdiff=abs(cc_stat.evoked.maxdiff);
        cc_stat.evoked.maxdiff_m=mean(cc_stat.evoked.maxdiff,1);
        cc_stat.evoked.maxdiff_std=std(cc_stat.evoked.maxdiff,0,1);
        cc_stat.evoked.change_maxdiff=[(cc_maxdiff_mat{2}(:,2)-cc_maxdiff_mat{2}(:,1))./abs(cc_maxdiff_mat{2}(:,1))].*100; %percent change
        cc_stat.evoked.change_maxdiff_m=mean(cc_stat.evoked.change_maxdiff,1);
        cc_stat.evoked.change_maxdiff_std=std(cc_stat.evoked.change_maxdiff,0,1);
        change_maxdiff_mat= cc_stat.evoked.change_maxdiff;
        %testing for normal distribution       
     if length(files_to_analyze)>=5
        [cc_stat.evoked.lillietest_h_maxdiff, cc_stat.evoked.lillietest_p_maxdiff] = lillietest(cc_stat.evoked.maxdiff(:,2)- cc_stat.evoked.maxdiff(:,1));
        [cc_stat.evoked.lillietest_h_change_maxdiff, cc_stat.evoked.lillietest_p_change_maxdiff] = lillietest( cc_stat.evoked.change_maxdiff);
     end
        %paired ttest 
        [cc_stat.evoked.ttest_h_maxdiff, cc_stat.evoked.ttest_p_maxdiff]= ttest(cc_stat.evoked.maxdiff(:,1),cc_stat.evoked.maxdiff(:,2));   
        [cc_stat.evoked.wilcoxon_p_maxdiff, cc_stat.evoked.wilcoxon_h_maxdiff]= signrank(cc_stat.evoked.maxdiff(:,1),cc_stat.evoked.maxdiff(:,2));
        [cc_stat.evoked.ttest_h_change_maxdiff, cc_stat.evoked.ttest_p_change_maxdiff]= ttest(cc_stat.evoked.change_maxdiff);
        [cc_stat.evoked.wilcoxon_p_change_maxdiff, cc_stat.evoked.wilcoxon_h_change_maxdiff]= signrank(cc_stat.evoked.change_maxdiff);
 
%% rmANOVA with built-in matlab function - lag0
% clear all
% load cc_analysis
lag0_spont_noNB(:,1)=cc_stat.spont.lag0(:,1);
lag0_spont_NB(:,1)=cc_stat.spont.lag0(:,2);
lag0_evoked_noNB(:,1)=cc_stat.evoked.lag0(:,1);
lag0_evoked_NB(:,1)=cc_stat.evoked.lag0(:,2);
lag0_spont(:,1)=[lag0_spont_noNB;lag0_spont_NB];
lag0_evoked(:,1)=[lag0_evoked_noNB;lag0_evoked_NB];
 S_1(:,1)=1:length(files_to_analyze);

ta_vector_names={'lag0_spont_noNB','lag0_spont_NB','lag0_evoked_noNB','lag0_evoked_NB'};
ta=table(S_1,lag0_spont_noNB,lag0_spont_NB,lag0_evoked_noNB,lag0_evoked_NB,...
    'variablenames', {'cells','Y1','Y2','Y3','Y4'});
factorNames = {'Sensory_stim','NB'};
within = table({'Sp';'Sp';'Ev';'Ev'},{'N';'Y';'N';'Y'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% fit the repeated measures model
rm_lag0 = fitrm(ta,'Y1-Y4~1','WithinDesign',within);
%test for sphericity
mauchly_lag0=mauchly(rm_lag0);
% run my repeated measures anova here
[ranovatbl_lag0] = ranova(rm_lag0, 'WithinModel','Sensory_stim*NB');
%multiple comparisons:
NBbySS_lag0 = multcompare(rm_lag0,'NB','By','Sensory_stim');
SSbyNB_lag0 = multcompare(rm_lag0,'Sensory_stim','By','NB');
lag0_rmANOVA_p= table2cell(ranovatbl_lag0(5,5));

eps_lag0 = epsilon(rm_lag0);
eps=table2cell(eps_lag0(1,2));
%in case where the sphericity assumption is violated (mauchly test p<0.05):
%if eps<0.75 - use the Greenhouse-Geisser correction.
%if eps>0.75 - use the Huynh-Feldt correction
if eps{1}<0.75
    lag0_p_evoked=table2cell(NBbySS_lag0(1,6));
    lag0_p_spont=table2cell(NBbySS_lag0(3,6));
else
     lag0_p_evoked=table2cell(NBbySS_lag0(1,6));
     lag0_p_spont=table2cell(NBbySS_lag0(3,6));
end

lag0_stat.table=ta;
lag0_stat.table_data_vecs=ta_vector_names;
lag0_stat.within_design=within;
lag0_stat.rm=rm_lag0;
lag0_stat.mauchly=mauchly_lag0;
lag0_stat.eps=eps_lag0;
lag0_stat.ANOVA=ranovatbl_lag0;
lag0_stat.multcomp_NBbySS=NBbySS_lag0;
lag0_stat.multcomp_SSbyNB=SSbyNB_lag0;
lag0_stat.rmANOVA_p=lag0_rmANOVA_p;
lag0_stat.evoked_p=lag0_p_evoked;
lag0_stat.spont_p=lag0_p_spont;
%% rmANOVA with matlab function - lag0 shuffled
lag0_shuff_spont_noNB(:,1)=cc_stat.spont.lag0_shuff(:,1);
lag0_shuff_spont_NB(:,1)=cc_stat.spont.lag0_shuff(:,2);
lag0_shuff_evoked_noNB(:,1)=cc_stat.evoked.lag0_shuff(:,1);
lag0_shuff_evoked_NB(:,1)=cc_stat.evoked.lag0_shuff(:,2);
lag0_shuff_spont(:,1)=[lag0_shuff_spont_noNB;lag0_shuff_spont_NB];
lag0_shuff_evoked(:,1)=[lag0_shuff_evoked_noNB;lag0_shuff_evoked_NB];
 S_1(:,1)=1:length(files_to_analyze);
 
clear rm ranovatbl ta
ta_vector_names={'lag0_shuff_spont_noNB','lag0_shuff_spont_NB','lag0_shuff_evoked_noNB','lag0_shuff_evoked_NB'};
ta=table(S_1,lag0_shuff_spont_noNB,lag0_shuff_spont_NB,lag0_shuff_evoked_noNB,lag0_shuff_evoked_NB,...
    'variablenames', {'cells','Y1','Y2','Y3','Y4'});
factorNames = {'Sensory_stim','NB'};
within = table({'Sp';'Sp';'Ev';'Ev'},{'N';'Y';'N';'Y'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% fit the repeated measures model
rm_lag0_shuff = fitrm(ta,'Y1-Y4~1','WithinDesign',within);
%test for sphericity
mauchly_lag0_shuff=mauchly(rm_lag0_shuff);
% run my repeated measures anova here
[ranovatbl_lag0_shuff] = ranova(rm_lag0_shuff, 'WithinModel','Sensory_stim*NB');
%multiple comparisons:
NBbySS_lag0_shuff = multcompare(rm_lag0_shuff,'NB','By','Sensory_stim');
SSbyNB_lag0_shuff = multcompare(rm_lag0_shuff,'Sensory_stim','By','NB');
lag0_shuff_rmANOVA_p= table2cell(ranovatbl_lag0_shuff(5,5));

eps_lag0_shuff = epsilon(rm_lag0_shuff);
eps=table2cell(eps_lag0_shuff(1,2));
%in case where the sphericity assumption is violated (mauchly test p<0.05):
%if eps<0.75 - use the Greenhouse-Geisser correction.
%if eps>0.75 - use the Huynh-Feldt correction
if eps{1}<0.75
    lag0_shuff_p_evoked=table2cell(NBbySS_lag0_shuff(1,6));
    lag0_shuff_p_spont=table2cell(NBbySS_lag0_shuff(3,6));
else
     lag0_shuff_p_evoked=table2cell(NBbySS_lag0_shuff(1,7));
     lag0_shuff_p_spont=table2cell(NBbySS_lag0_shuff(3,7));
end

lag0_shuff_stat.table=ta;
lag0_shuff_stat.table_data_vecs=ta_vector_names;
lag0_shuff_stat.within_design=within;
lag0_shuff_stat.rm=rm_lag0_shuff;
lag0_shuff_stat.mauchly=mauchly_lag0_shuff;
lag0_shuff_stat.eps=eps_lag0_shuff;
lag0_shuff_stat.ANOVA=ranovatbl_lag0_shuff;
lag0_shuff_stat.multcomp_NBbySS=NBbySS_lag0_shuff;
lag0_shuff_stat.multcomp_SSbyNB=SSbyNB_lag0_shuff;
lag0_shuff_stat.rmANOVA_p=lag0_shuff_rmANOVA_p;
lag0_shuff_stat.evoked_p=lag0_shuff_p_evoked;
lag0_shuff_stat.spont_p=lag0_shuff_p_spont;
%% rmANOVA with matlab function - max val
max_val_spont_noNB(:,1)=cc_stat.spont.max_val(:,1);
max_val_spont_NB(:,1)=cc_stat.spont.max_val(:,2);
max_val_evoked_noNB(:,1)=cc_stat.evoked.max_val(:,1);
max_val_evoked_NB(:,1)=cc_stat.evoked.max_val(:,2);
max_val_spont(:,1)=[max_val_spont_noNB;max_val_spont_NB];
max_val_evoked(:,1)=[max_val_evoked_noNB;max_val_evoked_NB];
 S_1(:,1)=1:length(files_to_analyze);

 ta_vector_names={'max_val_spont_noNB','max_val_spont_NB','max_val_evoked_noNB','max_val_evoked_NB'};
clear rm ranovatbl ta
ta=table(S_1,max_val_spont_noNB,max_val_spont_NB,max_val_evoked_noNB,max_val_evoked_NB,...
    'variablenames', {'cells','Y1','Y2','Y3','Y4'});
factorNames = {'Sensory_stim','NB'};
within = table({'Sp';'Sp';'Ev';'Ev'},{'N';'Y';'N';'Y'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% fit the repeated measures model
rm_mv = fitrm(ta,'Y1-Y4~1','WithinDesign',within);
%test for sphericity
mauchly_mv=mauchly(rm_mv);
% run my repeated measures anova here
[ranovatbl_mv] = ranova(rm_mv, 'WithinModel','Sensory_stim*NB');
%multiple comparisons:
NBbySS_mv = multcompare(rm_mv,'NB','By','Sensory_stim');
SSbyNB_mv = multcompare(rm_mv,'Sensory_stim','By','NB');
max_val_rmANOVA_p= table2cell(ranovatbl_mv(5,5));

eps_mv = epsilon(rm_mv);
eps=table2cell(eps_mv(1,2));
%in case where the sphericity assumption is violated (mauchly test p<0.05):
%if eps<0.75 - use the Greenhouse-Geisser correction.
%if eps>0.75 - use the Huynh-Feldt correction
if eps{1}<0.75
    max_val_p_evoked=table2cell(NBbySS_mv(1,6));
    max_val_p_spont=table2cell(NBbySS_mv(3,6));
else
     max_val_p_evoked=table2cell(NBbySS_mv(1,7));
    max_val_p_spont=table2cell(NBbySS_mv(3,7));
end
mv_stat.table=ta;
mv_stat.table_data_vecs=ta_vector_names;
mv_stat.within_design=within;
mv_stat.rm=rm_mv;
mv_stat.mauchly=mauchly_mv;
mv_stat.eps=eps_mv;
mv_stat.ANOVA=ranovatbl_mv;
mv_stat.multcomp_NBbySS=NBbySS_mv;
mv_stat.multcomp_SSbyNB=SSbyNB_mv;
mv_stat.rmANOVA_p=max_val_rmANOVA_p;
mv_stat.evoked_p=max_val_p_evoked;
mv_stat.spont_p=max_val_p_spont;

%%
max_val_shuff_spont_noNB(:,1)= cc_stat.spont.max_val_shuff(:,1);
max_val_shuff_spont_NB(:,1)=cc_stat.spont.max_val_shuff(:,2);
max_val_shuff_evoked_noNB(:,1)=cc_stat.evoked.max_val_shuff(:,1);
max_val_shuff_evoked_NB(:,1)=cc_stat.evoked.max_val_shuff(:,2);
 S_1(:,1)=1:length(files_to_analyze);
 
 ta_vector_names={'max_val_shuff_spont_noNB','max_val_shuff_spont_NB','max_val_shuff_evoked_noNB','max_val_shuff_evoked_NB'};

clear rm ranovatbl ta
ta=table(S_1,max_val_shuff_spont_noNB,max_val_shuff_spont_NB,max_val_shuff_evoked_noNB,max_val_shuff_evoked_NB,...
    'variablenames', {'cells','Y1','Y2','Y3','Y4'});
factorNames = {'Sensory_stim','NB'};
within = table({'Sp';'Sp';'Ev';'Ev'},{'N';'Y';'N';'Y'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% fit the repeated measures model
rm_mv_shuff = fitrm(ta,'Y1-Y4~1','WithinDesign',within);
%test for sphericity
mauchly_mv_shuff=mauchly(rm_mv_shuff);
% run my repeated measures anova here
[ranovatbl_mv_shuff] = ranova(rm_mv_shuff, 'WithinModel','Sensory_stim*NB');
%multiple comparisons:
NBbySS_mv_shuff = multcompare(rm_mv_shuff,'NB','By','Sensory_stim');
SSbyNB_mv_shuff = multcompare(rm_mv_shuff,'Sensory_stim','By','NB');
max_val_shuff_rmANOVA_p= table2cell(ranovatbl_mv_shuff(5,5));

eps_mv_shuff = epsilon(rm_mv_shuff);
eps=table2cell(eps_mv_shuff(1,2));
%in case where the sphericity assumption is violated (mauchly test p<0.05):
%if eps<0.75 - use the Greenhouse-Geisser correction.
%if eps>0.75 - use the Huynh-Feldt correction
if eps{1}<0.75
    max_val_shuff_p_evoked=table2cell(NBbySS_mv_shuff(1,6));
    max_val_shuff_p_spont=table2cell(NBbySS_mv_shuff(3,6));
else
     max_val_shuff_p_evoked=table2cell(NBbySS_mv_shuff(1,7));
    max_val_shuff_p_spont=table2cell(NBbySS_mv_shuff(3,7));
end

mv_shuff_stat.table=ta;
mv_shuff_stat.table_data_vecs=ta_vector_names;
mv_shuff_stat.within_design=within;
mv_shuff_stat.rm=rm_mv_shuff;
mv_shuff_stat.mauchly=mauchly_mv_shuff;
mv_shuff_stat.eps=eps_mv_shuff;
mv_shuff_stat.ANOVA=ranovatbl_mv_shuff;
mv_shuff_stat.multcomp_NBbySS=NBbySS_mv_shuff;
mv_shuff_stat.multcomp_SSbyNB=SSbyNB_mv;
mv_shuff_stat.rmANOVA_p=max_val_rmANOVA_p;
mv_shuff_stat.evoked_p=max_val_shuff_p_evoked;
mv_shuff_stat.spont_p=max_val_shuff_p_spont;
         %% Paired plot of spont+evoked non-normalized lag0 absolute values
% CC_Y_max_val_sp=cc_stat.spont.max_val';
CC_Y_max_val_sp=cc_stat.spont.lag0';
CC_X_sp(1,:)=ones(1,size(CC_Y_max_val_sp,2));
CC_X_sp(2,:)=2*ones(1,size(CC_Y_max_val_sp,2));
E_sp = std(CC_Y_max_val_sp,0,2);
% CC_Y_max_val_ev=cc_stat.evoked.max_val';
CC_Y_max_val_ev=cc_stat.evoked.lag0';
CC_X_ev(1,:)=3*ones(1,size(CC_Y_max_val_ev,2));
CC_X_ev(2,:)=4*ones(1,size(CC_Y_max_val_ev,2));
E_ev = std(CC_Y_max_val_ev,0,2);
linex=[1 3;2 4];
% my=max(max([cc_stat.spont.max_val,cc_stat.evoked.max_val]))*1.1; 
my=0.8;
liney=[my my;my my];
x1limits = [0.75 4.25];
x1ticks = [1,2,3,4];
 y1limits = [-0.1 1];
y1ticks = [0,0.5,1]; 
% p_spont=mv_stat.spont_p;
% p_evoked=mv_stat.evoked_p;
p_spont=lag0_stat.spont_p;
p_evoked=lag0_stat.evoked_p;
if p_spont{1,1} >0.05 
    asterisk_sp='n.s.';
else if p_spont{1,1}<0.05 && p_spont{1,1}>0.01
    asterisk_sp='*';
    else if p_spont{1,1}<0.01 && p_spont{1,1}>0.001
            asterisk_sp='**';
    else if p_spont{1,1}<0.001
             asterisk_sp='***';
        end
        end
    end
end

if p_evoked{1,1} >0.05 
    asterisk_ev='n.s.';
else if p_evoked{1,1}<0.05 && p_evoked{1,1}>0.01
    asterisk_ev='*';
    else if p_evoked{1,1}<0.01 && p_evoked{1,1}>0.001
            asterisk_ev='**';
    else if p_evoked{1,1}<0.001
             asterisk_ev='***';
        end
        end
    end
end

figure
hold on
set(gca,'ylim', y1limits)
line(CC_X_sp,CC_Y_max_val_sp,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_sp(:,1), mean(CC_Y_max_val_sp,2),E_sp,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(CC_X_ev,CC_Y_max_val_ev,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_ev(:,1), mean(CC_Y_max_val_ev,2),E_ev,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(3.5,my,asterisk_ev,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(1.5,my+0.15,'Spontaneous','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',10)
text(3.5,my+0.15,'Sensory evoked','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',10)
hold off

set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',22,'linewidth',1, 'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',[legend_string, legend_string] ,'box', 'off'); %'fontweight', 'bold',  
         ylabel('CC', 'FontSize', 20,'fontname', 'arial');
%            title(['Spontaneous and Sensory-evoked Max-Peak Cross-Correlation, n=' num2str(length(files_to_analyze))] ,'fontname', 'arial','FontSize', 20);   
 if save_flag==1;
 cd(path_output)
%  filename1='Vm-LFPcc_spont+evoked_max-peak_paired_population';
filename1='Vm-LFPcc_spont+evoked_lag0_paired_population';
saveas(gcf,filename1,'fig') 
print(gcf,filename1,'-r600','-opengl') 
 end
         %% Paired plot of shuffled spont+evoked non-normalized max peak absolute values
% CC_Y_max_val_sp_shuff=cc_stat.spont.max_val_shuff';
CC_Y_max_val_sp_shuff=cc_stat.spont.lag0_shuff';
CC_X_sp_shuff(1,:)=ones(1,size(CC_Y_max_val_sp_shuff,2));
CC_X_sp_shuff(2,:)=2*ones(1,size(CC_Y_max_val_sp_shuff,2));
E_sp_shuff = std(CC_Y_max_val_sp_shuff,0,2);
% CC_Y_max_val_ev_shuff=cc_stat.evoked.max_val_shuff';
CC_Y_max_val_ev_shuff=cc_stat.evoked.lag0_shuff';
CC_X_ev_shuff(1,:)=3*ones(1,size(CC_Y_max_val_ev_shuff,2));
CC_X_ev_shuff(2,:)=4*ones(1,size(CC_Y_max_val_ev_shuff,2));
E_ev_shuff = std(CC_Y_max_val_ev_shuff,0,2);
linex=[1 3;2 4];
% my=max(max([cc_stat.spont.max_val_shuff,cc_stat.evoked.max_val_shuff]))*1.1; 
liney=[my my;my my];
x1limits = [0.75 4.25];
x1ticks = [1,2,3,4];
 y1limits = [-0.1 1];
y1ticks = [0,0.5,1];
% shuff_p_spont=mv_shuff_stat.spont_p;
% shuff_p_evoked=mv_shuff_stat.evoked_p;
shuff_p_spont=lag0_shuff_stat.spont_p;
shuff_p_evoked=lag0_shuff_stat.evoked_p;
if shuff_p_spont{1,1} >0.05 
    asterisk_sp='n.s.';
else if shuff_p_spont{1,1}<0.05 && shuff_p_spont{1,1}>0.01
    asterisk_sp='*';
    else if shuff_p_spont{1,1}<0.01 && shuff_p_spont{1,1}>0.001
            asterisk_sp='**';
    else if shuff_p_spont{1,1}<0.001
             asterisk_sp='***';
        end
        end
    end
end

if shuff_p_evoked{1,1} >0.05 
    asterisk_ev='n.s.';
else if shuff_p_evoked{1,1}<0.05 && shuff_p_evoked{1,1}>0.01
    asterisk_ev='*';
    else if shuff_p_evoked{1,1}<0.01 && shuff_p_evoked{1,1}>0.001
            asterisk_ev='**';
    else if shuff_p_evoked{1,1}<0.001
             asterisk_ev='***';
        end
        end
    end
end

figure
set(gca,'ylim', y1limits)
hold on
line(CC_X_sp_shuff,CC_Y_max_val_sp_shuff,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_sp_shuff(:,1), mean(CC_Y_max_val_sp_shuff,2),E_sp_shuff,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(CC_X_ev_shuff,CC_Y_max_val_ev_shuff,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_ev_shuff(:,1), mean(CC_Y_max_val_ev_shuff,2),E_ev_shuff,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(3.5,my,asterisk_ev,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(1.5,my+0.15,'Spontaneous','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',10)
text(3.5,my+0.15,'Sensory evoked','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',10)
hold off

set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',22,'linewidth',1, 'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',[legend_string, legend_string] ,'box', 'off'); %'fontweight', 'bold',  
         ylabel('CC', 'FontSize', 20,'fontname', 'arial');
%            title(['Spontaneous and Sensory-evoked Max-Peak Cross-Correlation, n=' num2str(length(files_to_analyze))] ,'fontname', 'arial','FontSize', 20);   
 if save_flag==1;
 cd(path_output)
%  filename1='Vm-LFPcc_spont+evoked_max-peak_paired_population_shuff';
filename2='Vm-LFPcc_spont+evoked_lag0_paired_population_shuff';
saveas(gcf,filename2,'fig') 
print(gcf,filename2,'-dpng','-r600','-opengl') 
 end
  %% Paired plot of spont non-normalized lag0 absolute values: real data+shuffled data
CC_Y_sp_shuff=cc_stat.spont.lag0_shuff';
CC_X_sp_shuff(1,:)=ones(1,size(CC_Y_sp_shuff,2));
CC_X_sp_shuff(2,:)=2*ones(1,size(CC_Y_sp_shuff,2));
E_sp_shuff = std(CC_Y_sp_shuff,0,2);
CC_Y_sp=cc_stat.spont.lag0';
CC_X_sp(1,:)=3*ones(1,size(CC_Y_sp,2));
CC_X_sp(2,:)=4*ones(1,size(CC_Y_sp,2));
E_sp = std(CC_Y_sp,0,2);

linex=[1 3;2 4];
% my=max(max([cc_stat.spont.max_val,cc_stat.evoked.max_val]))*1.1; 
my=0.8;
liney=[my my;my my];
x1limits = [0.75 4.25];
x1ticks = [1,2,3,4];
 y1limits = [-0.1 1];
y1ticks = [0,0.5,1]; 

p_spont_shuff=lag0_shuff_stat.spont_p; 
p_spont=lag0_stat.spont_p;
if p_spont_shuff{1,1} >0.05 
    asterisk_sp_shuff='n.s.';
else if p_spont_shuff{1,1}<0.05 && p_spont_shuff{1,1}>0.01
    asterisk_sp_shuff='*';
    else if p_spont_shuff{1,1}<0.01 && p_spont_shuff{1,1}>0.001
            asterisk_sp_shuff='**';
    else if p_spont_shuff{1,1}<0.001
             asterisk_sp_shuff='***';
        end
        end
    end
end

if p_spont{1,1} >0.05 
    asterisk_sp='n.s.';
else if p_spont{1,1}<0.05 && p_spont{1,1}>0.01
    asterisk_sp='*';
    else if p_spont{1,1}<0.01 && p_spont{1,1}>0.001
            asterisk_sp='**';
    else if p_spont{1,1}<0.001
             asterisk_sp='***';
        end
        end
    end
end

figure
hold on
set(gca,'ylim', y1limits)
line(CC_X_sp_shuff,CC_Y_sp_shuff,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_sp_shuff(:,1), mean(CC_Y_sp_shuff,2),E_sp_shuff,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(CC_X_sp,CC_Y_sp,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_sp(:,1), mean(CC_Y_sp,2),E_sp,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp_shuff,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(3.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(1.5,my+0.15,'Shuffled','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',11)
text(3.5,my+0.15,'Real','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',11)
hold off

set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',22,'linewidth',1, 'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',[legend_string, legend_string] ,'box', 'off'); %'fontweight', 'bold',  
         ylabel('CC', 'FontSize', 20,'fontname', 'arial');
%            title(['Spontaneous and Sensory-evoked Max-Peak Cross-Correlation, n=' num2str(length(files_to_analyze))] ,'fontname', 'arial','FontSize', 20);   
 if save_flag==1;
 cd(path_output)
%  filename1='Vm-LFPcc_spont+evoked_max-peak_paired_population';
filename1='Vm-LFPcc_spont_lag0_paired_population_real-shuff';
saveas(gcf,filename1,'fig') 
print(gcf,filename1,'-r600','-opengl') 
 end
 %% Paired plot of evoked non-normalized max peak absolute values: real data+shuffled data
CC_Y_ev_shuff=cc_stat.evoked.lag0_shuff';
CC_X_ev_shuff(1,:)=ones(1,size(CC_Y_ev_shuff,2));
CC_X_ev_shuff(2,:)=2*ones(1,size(CC_Y_ev_shuff,2));
E_ev_shuff = std(CC_Y_ev_shuff,0,2);
CC_Y_ev=cc_stat.evoked.lag0';
CC_X_ev(1,:)=3*ones(1,size(CC_Y_ev,2));
CC_X_ev(2,:)=4*ones(1,size(CC_Y_ev,2));
E_ev = std(CC_Y_ev,0,2);

linex=[1 3;2 4];
% my=max(max([cc_stat.spont.max_val_shuff,cc_stat.evoked.max_val_shuff]))*1.1; 
liney=[my my;my my];
x1limits = [0.75 4.25];
x1ticks = [1,2,3,4];
 y1limits = [-0.1 1];
y1ticks = [0,0.5,1];

p_evoked_shuff=lag0_shuff_stat.evoked_p;
p_evoked=lag0_stat.evoked_p;
if p_evoked_shuff{1,1} >0.05 
    asterisk_ev_shuff='n.s.';
else if p_evoked_shuff{1,1}<0.05 && p_evoked_shuff{1,1}>0.01
    asterisk_ev_shuff='*';
    else if p_evoked_shuff{1,1}<0.01 && p_evoked_shuff{1,1}>0.001
            asterisk_ev_shuff='**';
    else if p_evoked_shuff{1,1}<0.001
             asterisk_ev_shuff='***';
        end
        end
    end
end

if p_evoked{1,1} >0.05 
    asterisk_ev='n.s.';
else if p_evoked{1,1}<0.05 && p_evoked{1,1}>0.01
    asterisk_ev='*';
    else if p_evoked{1,1}<0.01 && p_evoked{1,1}>0.001
            asterisk_ev='**';
    else if p_evoked{1,1}<0.001
             asterisk_ev='***';
        end
        end
    end
end

figure
set(gca,'ylim', y1limits)
hold on
line(CC_X_ev_shuff,CC_Y_ev_shuff,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_ev_shuff(:,1), mean(CC_Y_ev_shuff,2),E_ev_shuff,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(CC_X_ev,CC_Y_ev,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_ev(:,1), mean(CC_Y_ev,2),E_ev,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_ev_shuff,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(3.5,my,asterisk_ev,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(1.5,my+0.15,'Shuffled','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',11)
text(3.5,my+0.15,'Real','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',11)
hold off

set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',22,'linewidth',1, 'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',[legend_string, legend_string] ,'box', 'off'); %'fontweight', 'bold',  
         ylabel('CC', 'FontSize', 20,'fontname', 'arial');
%            title(['Spontaneous and Sensory-evoked Max-Peak Cross-Correlation, n=' num2str(length(files_to_analyze))] ,'fontname', 'arial','FontSize', 20);   
 if save_flag==1;
 cd(path_output)
%  filename1='Vm-LFPcc_spont+evoked_max-peak_paired_population_shuff';
filename2='Vm-LFPcc_evoked_lag0_paired_population_real-shuff';
saveas(gcf,filename2,'fig') 
print(gcf,filename2,'-dpng','-r600','-opengl') 
 end
          %% Paired plot of spont+evoked non-normalized max peak real values
CC_Y_max_val_sp=(cc_stat.spont.max_val.*cc_stat.spont.max_val_sign)';
CC_X_sp(1,:)=ones(1,size(CC_Y_max_val_sp,2));
CC_X_sp(2,:)=2*ones(1,size(CC_Y_max_val_sp,2));
E_sp = std(CC_Y_max_val_sp,0,2);
CC_Y_max_val_ev=(cc_stat.evoked.max_val.*cc_stat.evoked.max_val_sign)';
CC_X_ev(1,:)=3*ones(1,size(CC_Y_max_val_ev,2));
CC_X_ev(2,:)=4*ones(1,size(CC_Y_max_val_ev,2));
E_ev = std(CC_Y_max_val_ev,0,2);
linex=[1 3;2 4];
% my=max(max([cc_stat.spont.max_val.*cc_stat.spont.max_val_sign,cc_stat.evoked.max_val.*cc_stat.evoked.max_val_sign]))*1.1; 
liney=[my my;my my];
x1limits = [0.75 4.25];
x1ticks = [1,2,3,4];
 y1limits = [-1 1];
y1ticks = [-1,-0.5,0,0.5,1];

figure
hold on
set(gca,'ylim', y1limits)
line(CC_X_sp,CC_Y_max_val_sp,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_sp(:,1), mean(CC_Y_max_val_sp,2),E_sp,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(CC_X_ev,CC_Y_max_val_ev,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_ev(:,1), mean(CC_Y_max_val_ev,2),E_ev,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
% text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
% text(3.5,my,asterisk_ev,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(1.5,my+0.15,'Spontaneous','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',10)
text(3.5,my+0.15,'Sensory evoked','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',10)
hold off

set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',22,'linewidth',1,'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',[legend_string, legend_string] ,'box', 'off'); %'fontweight', 'bold',  
         ylabel('Peak CC', 'FontSize', 20,'fontname', 'arial');
%            title(['Spontaneous and Sensory-evoked Max-Peak Cross-Correlation, n=' num2str(length(files_to_analyze))] ,'fontname', 'arial','FontSize', 20);   
 if save_flag==1;
 cd(path_output)
saveas(gcf,'Vm-LFPcc_spont+evoked_max-peak_paired_population_real_value','fig') 
print(gcf,'Vm-LFPcc_spont+evoked_max-peak_paired_population_real_value','-dpng','-r600','-opengl') 
 end
         %% Paired plot of shuffled spont+evoked non-normalized max peak real values
CC_Y_max_val_sp_shuff=(cc_stat.spont.max_val_shuff.*cc_stat.spont.max_val_shuff_sign)';
CC_X_sp_shuff(1,:)=ones(1,size(CC_Y_max_val_sp_shuff,2));
CC_X_sp_shuff(2,:)=2*ones(1,size(CC_Y_max_val_sp_shuff,2));
E_sp_shuff = std(CC_Y_max_val_sp_shuff,0,2);
CC_Y_max_val_ev_shuff=(cc_stat.evoked.max_val_shuff.*cc_stat.evoked.max_val_shuff_sign)';
CC_X_ev_shuff(1,:)=3*ones(1,size(CC_Y_max_val_ev_shuff,2));
CC_X_ev_shuff(2,:)=4*ones(1,size(CC_Y_max_val_ev_shuff,2));
E_ev_shuff = std(CC_Y_max_val_ev_shuff,0,2);
linex=[1 3;2 4];
% my=max(max([cc_stat.spont.max_val_shuff.*cc_stat.spont.max_val_shuff_sign,cc_stat.evoked.max_val_shuff.*cc_stat.evoked.max_val_shuff_sign]))*1.1; 
liney=[my my;my my];
x1limits = [0.75 4.25];
x1ticks = [1,2,3,4];
 y1limits = [-1 1];
y1ticks = [-1,-0.5,0,0.5];

figure
hold on
set(gca,'ylim', y1limits)
line(CC_X_sp_shuff,CC_Y_max_val_sp_shuff,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_sp_shuff(:,1), mean(CC_Y_max_val_sp_shuff,2),E_sp_shuff,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(CC_X_ev_shuff,CC_Y_max_val_ev_shuff,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(CC_X_ev_shuff(:,1), mean(CC_Y_max_val_ev_shuff,2),E_ev_shuff,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
% text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
% text(3.5,my,asterisk_ev,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(1.5,my+0.15,'Spontaneous','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',10)
text(3.5,my+0.15,'Sensory evoked','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',10)
hold off

set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',22,'linewidth',1, 'ylim', y1limits,'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',[legend_string, legend_string] ,'box', 'off'); %'fontweight', 'bold',  
         ylabel('Peak CC', 'FontSize', 20,'fontname', 'arial');
%            title(['Spontaneous and Sensory-evoked Max-Peak Cross-Correlation, n=' num2str(length(files_to_analyze))] ,'fontname', 'arial','FontSize', 20);   
 if save_flag==1;
 cd(path_output)
saveas(gcf,'Vm-LFPcc_spont+evoked_max-peak_paired_population_shuff_real_value','fig') 
print(gcf,'Vm-LFPcc_spont+evoked_max-peak_paired_population_shuff_real_value','-dpng','-r600','-opengl') 
 end
%% Save analysis
 if save_flag==1;
    cd(path_output)
     filename='cc_analysis'; 
    save(filename, 'files_to_analyze', 'cc_spont', 'cc_evoked', 'cc_stat','mv_stat','mv_shuff_stat','lag0_stat','lag0_shuff_stat')
 end
