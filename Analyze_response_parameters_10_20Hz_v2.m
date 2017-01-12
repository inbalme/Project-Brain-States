%% Analyze_response_parameters_10_20Hz_v2
% This file was created on 19/12/2016 based on Analyze_NBES_response_parameters_10_20Hz_v2
%This file is used for the analysis of files created with
%extract_NBES_Data_v3 or extract_ChAT_Data_v3

%       protocol 5 (NBES+galvnao (5 x-values))
%       protocol 6 (ChAT intracellular current injections light+galvnao (15 x-values))
%       protocol 8 - ChAT intracellular: light only (3 x-values): 1) train
%       20 Hz, 10ms ondur, 50 pulses. 2) single pulse 5 sec ondur. 3)single pulse 10 ms
%       protocol 9 - IC: I-V curve
%       protocol 10 (NBES:galvnao train+test (3 x-values) or ChAT:galvnao train+test )
%                 11 - NBES: 3 currents+galvano train+test (6 x-values), or ChAT: 3 currents+galvano train+test (6 x-values)
%% for opening workspace saved 
clear all
global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data Ch2_data stim2_X stim1_X 
 global exp_type
 channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
exp_type=1; %1-NBES, 2-ChAT
% trace_type_input=[3,2]; %1:3
save_flag= 0;
print_flag=1;
norm_flag=0;
BP50HzLFP_flag=1; %removing 50Hz noise from LFP signal
BP50HzVm_flag=1; %removing 50Hz noise from Vm signal
BPLFP_flag=1; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=[0.1,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=1; %filtering LFP and Vm same as LFP was filtered in the multiclamp
bp_manual_Vm=[0,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)

 %% set the path for saving figures and variables
if BP50HzLFP_flag==1 && BP50HzVm_flag==1 && BPVm_flag==1 && BPLFP_flag==1
    path_output='LFP_50Hz+BP Vm_ 50Hz+BP';
else if BP50HzLFP_flag==1 && BP50HzVm_flag==1 && BPLFP_flag==1
     path_output='LFP_50Hz+BP Vm_ 50Hz';  
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
switch exp_type
    case 1 
        path_output=['D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\',path_output];   
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        cd(path_output)       
    case 2
        path_output=['D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary\',path_output];        
        a = exist(path_output,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        cd(path_output)
end

switch exp_type
    case 1
        files_to_analyze =[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,44,46,48,50,52,56,58,62,72,75]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};       

    case 2
        files_to_analyze =[74,76,77,80,82,84,87];
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light On', 'Light Off'};
end
        
% for fileind=1;
    for fileind=1:length(files_to_analyze) ;
    close all
    clear data_no_spike_no_DC
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
    %%
        Ch2_data= raw_data{3}./20; %dividing by the LFP gain           
        current_data=data_no_spikes{channel};
        galvano_nstim = Param.facade(6);
        galvano_freq = Param.facade(7);
   data_preprocessing 
        clear color_table    
        whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;
        switch exp_type
            case 1 
                color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
            case 2
                color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
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
               
                data_no_spike_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(data_no_spikes{channel}(:,:,x_value),meansubtract_interval);
%                 data_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(raw_data{channel}(:,:,x_value),meansubtract_interval);    
end
%% visualizing raw data traces for ongoing and evoked NB- and NB+ together  
%  figure
%  hold on
%  plot(interval(:,1).*dt,raw_data{1}(interval(:,1),1:5,1),'k');
%   plot(interval(:,1).*dt,raw_data{1}(interval(:,2),1:5,1),'b');
%   hold off
%   pause
% 
%  figure
%  hold on
%  plot((stim2_X{2}(1,1):(stim2_X{2}(1,1)+sf{1}-1)).*dt,raw_data{1}(stim2_X{2}(1,1):stim2_X{2}(1,1)+sf{1}-1,1:3,2),'k');
%   plot((stim2_X{2}(1,1):(stim2_X{2}(1,1)+sf{1}-1)).*dt,raw_data{1}(stim2_X{2}(1,1):stim2_X{2}(1,1)+sf{1}-1,1:3,3),'b');
%   hold off
%   pause
%     end %temp
%     for fileind=1; %temp
%% Variance plots - two x_values on each plot, for flexible data
x_value = 2:3;

plot_data=current_data; %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
plot_data_mean = mean(plot_data,2);
plot_data_std =  std(data_no_spike_no_DC{channel},0,2);
% plot_data_var=var(plot_data,0,2);
plot_data_var=var(data_no_spike_no_DC{channel},0,2);
plot_data_ff = (-1).*plot_data_var./plot_data_mean; 
    if exp_type==1;
        for i=1:length(x_value)
            if ES_flag(x_value(i))==1
                plot_data(29910:35010,:,x_value(i)) =nan; %Ignoring the ES artifact
                plot_data_mean(29910:35010,:,x_value(i)) =nan;
                plot_data_std(29910:35010,:,x_value(i)) =nan;
                plot_data_CV(29910:35010,:,x_value(i)) =nan;
            end
        end
    end 

%        trace_ind = [2,4,5]; %1:size(plot_data,2);  %trace_ind is the numbers of traces to be plotted
%        l=size(plot_data(:,trace_ind, x_value(1)),1);
%         l=l/2-1;
%         DC=20.*(0:length(trace_ind)-1);
%         DC=wextend('addrow', 'per', DC, l);
%         trace_to_plot = plot_data(:,trace_ind, x_value(1))+DC ; %DC is added just to make space between the traces
%         tmp1=[]; %spacer instead of stim1_X, so that the light stim. will not be plotted
%         lineprops1={'LineWidth',1.2,'color', [0 0 0]};
%         [Fig1,h1]= fn_Plot_Trace_v2(trace_to_plot, dt, dt,tmp1 , dt, stim2_X{x_value(1)},lineprops1);       
%         hline = findobj(gcf, 'type', 'line'); set(hline(1:end),'color',[0 0 0]);
%         y1lim=[get(gca,'ylim')]';
%         trace_to_plot = plot_data(:,trace_ind, x_value(2))+DC ;
%         lineprops2={'LineWidth',1.2,'color', [0 0 1]};
%         [Fig2,h2]= fn_Plot_Trace_v2(trace_to_plot, dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)},lineprops2);
%          hline = findobj(gcf, 'type', 'line'); set(hline(1:end),'color',[0 0 1]);
%         set(gca,'ylim',y1lim);
% %         [Fig3,h3]= fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
% %         ylabel('mean Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
% % %         hline = findobj(gcf, 'type', 'line'); set(hline(13),'color',[0 0 0]);  set(hline(14),'color',[0 0 1]);
% %         [Fig4,h4] = fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
% %         ylabel('std Vm [mV]', 'FontSize', 16);          
%         [Fig5,h5]= fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});        
%         ylabel('mean Vm [mV]', 'FontSize', 16);  legend(legend_string); legend('boxoff','Location','northeast'); 
% %         hline = findobj(gcf, 'type', 'line'); set(hline(13),'color',[0 0 0]);  set(hline(14),'color',[0 0 1]);
%         [Fig6,h6] = fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
%         ylabel('std Vm [mV]', 'FontSize', 16);  legend(legend_string); legend('boxoff','Location','northeast');

     
%% saving figures
% if save_flag==1;
%     cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10Hz'  
%         saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x2.fig']) 
%         print(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x2'],'-dpng','-r600','-opengl') 
%         saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x3.fig']) 
%         print(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x3'],'-dpng','-r600','-opengl') 
%         saveas(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt.fig']) 
%         print(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt'],'-dpng','-r600','-opengl') 
%         saveas(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3.fig']) 
%         print(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3'],'-dpng','-r600','-opengl') 

%         set(ax1,'xlim',x2lim); set(ax2,'xlim',x2lim); set(ax5,'xlim',x2lim); set(ax6,'xlim',x2lim); set(ax7,'xlim',x2lim); set(ax9,'xlim',x2lim); 
%         
%         saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x2_zoom-in.fig']) 
%         print(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x2_zoom-in'],'-dpng','-r600','-opengl') 
%         saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x3_zoom-in.fig']) 
%         print(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x3_zoom-in'],'-dpng','-r600','-opengl') 
%         saveas(Fig5,['f' num2str(files_to_analyze(fileind)) '_Vm std_x2+3_mean-subt_zoom-in.fig']) 
%         print(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt_zoom-in'],'-dpng','-r600','-opengl') 
%         saveas(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_zoom-in.fig']) 
%         print(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_zoom-in'],'-dpng','-r600','-opengl') 
        
% end        
% close all

%% Response parameters - amplitude, latency and more...
%find peak values and locations in interval following each whisker stim.

    clear data_mean data_std data_mean_smooth mean_peak_val_all mean_peak_loc_all mean_peak_val mean_peak_loc baseline_interval baseline_local_M baseline_local_average...
        mean_peak_half_peak_rise mean_peak_half_peak_fall mean_peak_half_peak_width temp temp_mean half_peak_val response_rise_interval...
        response_fall_interval baseline_global_M mean_10per_val mean_10per_loc mean_90per_val mean_90per_loc...
        baseline_local_std baseline_global baseline_global_resid baseline_global_STD baseline_global_average...
        data_pre_response data_pre_response_resid data_train_response data_train_response_resid data_post_train_response data_post_train_response_resid...
        data_response data_response_resid train_interval signal signal_noDC signal_noDC_squared noise1 noise1_noDC noise1_noDC_squared noise2 noise2_noDC noise2_noDC_squared
t=0; 
    int_time=50; %[ms]
    min_res_del=4; %[ms]
    baseline_time=3; %[ms] %local baseline before each peak
    baseline_time_total=100; %[ms] %global baseline for the whole trace
Y_abs = []; f = []; Y_abs_exp=[];  f_exp = [];
for x_value = 2:3;
    t=t+1; %index for the power spectrum
    data_mean = mean(current_data(:,:,x_value),2); %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC 
    data_mean_smooth=sgolayfilt(data_mean, 1, 29);
%     data_deriv1=diff(data_mean_smooth);
%     data_deriv2=diff(data_mean_smooth,2);
  
    
    %taking the 50 ms prior to the first stim:
                data_pre_response(:,:,x_value) = current_data(stim2_X{2}(1,1)-int_time.*sf{1}./1000:stim2_X{2}(1,1)-1,:,x_value);
                data_pre_response_resid(:,:,x_value) = fn_Subtract_Mean(data_pre_response(:,:,x_value));
                pre_response_M(:,1) =mean(mean(data_pre_response(:,:,x_value),2));
                pre_response_STD(:,1) = mean(std(data_pre_response(:,:,x_value),0,2)); %mean std across traces (trial-to-trial)
                pre_response_VAR(:,1) = mean(var(data_pre_response(:,:,x_value),0,2)); %mean var across traces (trial-to-trial)
                pre_response_VAR_time(:,1) = mean(var(data_pre_response(:,:,x_value),0,1)); %mean var across time (along trace)
                pre_response_CV(:,1) = mean(std(data_pre_response(:,:,x_value),0,2)./mean(data_pre_response(:,:,x_value),2)); %mean CV across traces (trial-to-trial)

    %taking the response as the whole train:  
                start_sample =stim2_X{x_value(1)}(1,1);  %start with sensory stim
                galvano_nstim = Param.facade(6);
                galvano_freq = Param.facade(7);
                duration = galvano_nstim./galvano_freq+0.05;
                end_sample = start_sample+duration.*sf{1}-1;
                train_interval(:,1) = round(start_sample:end_sample);           
                data_train_response(:,:,x_value) = current_data(train_interval,:,x_value);
                data_train_response_resid(:,:,x_value) = fn_Subtract_Mean(data_train_response(:,:,x_value));  
                Vm_res_train_M(:,1) = mean(mean(data_train_response(:,:,x_value),2)); %mean along the trace and across traces (trial-to-trial)
                Vm_res_train_STD(:,1) = mean(std(data_train_response(:,:,x_value),0,2)); %mean std across traces (trial-to-trial)
                Vm_res_train_VAR(:,1) = mean(var(data_train_response(:,:,x_value),0,2)); %mean var across traces (trial-to-trial)
                Vm_res_train_VAR_time(:,1) = mean(var(data_train_response(:,:,x_value),0,1)); %mean var across time
                Vm_res_train_CV(:,1) = mean(std(data_train_response(:,:,x_value),0,2)./mean(data_train_response(:,:,x_value),2)); %mean CV across traces (trial-to-trial)
                Vm_res_train_M_resid(:,1) = mean(mean(data_train_response_resid(:,:,x_value),2)); %mean across traces (trial-to-trial)
        
        %looking at the 300ms after the train:
                post_train_interval(:,1) = round(end_sample+0.1*sf{1}:end_sample+0.3*sf{1}-1);           
                data_post_train_response(:,:,x_value) = current_data(post_train_interval,:,x_value);
                data_post_train_response_resid(:,:,x_value) = fn_Subtract_Mean(data_post_train_response(:,:,x_value));  
                post_train_M(:,1) = mean(mean(data_post_train_response(:,:,x_value),2)); %mean along the trace and across traces (trial-to-trial)
                post_train_STD(:,1) = mean(std(data_post_train_response(:,:,x_value),0,2)); %mean std across traces (trial-to-trial)
                post_train_VAR(:,1) = mean(var(data_post_train_response(:,:,x_value),0,2)); %mean var across traces (trial-to-trial)
                post_train_VAR_time(:,1) = mean(var(data_post_train_response(:,:,x_value),0,1)); %mean var across time (along trace)
                post_train_CV(:,1) = mean(std(data_post_train_response(:,:,x_value),0,2)./mean(data_post_train_response(:,:,x_value),2)); %mean CV across traces (trial-to-trial)
                post_train_M_resid(:,1) = mean(mean(data_post_train_response_resid(:,:,x_value),2)); %mean along the trace and across traces (trial-to-trial)
%computing SNR:
signal(:,1)=mean(data_train_response(:,:,x_value),2); %mean across trials
signal_noDC=signal(:,1)-mean(signal(:,1)); %subtracting the mean over time
signal_noDC_squared(:,1)=signal_noDC.^2;
signal_noDC_squared_m=mean(signal_noDC_squared); %mean over time
signal_noDC_squared_m_sqrt=sqrt(signal_noDC_squared_m);
Amplitude_signal=signal_noDC_squared_m_sqrt;
for trace=1:size(data_train_response,2);
    noise1(:,trace)=data_train_response(:,trace,x_value)-signal;
    noise1_noDC(:,trace)=noise1(:,trace)-mean(noise1(:,trace));
    noise1_noDC_squared(:,trace)=noise1_noDC(:,trace).^2;
    noise1_noDC_squared_m(:,trace)=mean(noise1_noDC_squared(:,trace));
    noise1_noDC_squared_m_sqrt(:,trace)=sqrt(noise1_noDC_squared_m(:,trace));
end
Amplitude_noise1=mean(noise1_noDC_squared_m_sqrt);
%alternative calculation of noise: taking the ongoing activity before the sensory stim as the noise
for trace=1:size(data_train_response,2);
    noise2(:,trace)=current_data(start_sample-duration.*sf{1}:start_sample-1,trace,x_value);
    noise2_noDC(:,trace)=noise2(:,trace)-mean(noise2(:,trace));
    noise2_noDC_squared(:,trace)=noise2_noDC(:,trace).^2;
    noise2_noDC_squared_m(:,trace)=mean(noise2_noDC_squared(:,trace));
    noise2_noDC_squared_m_sqrt(:,trace)=sqrt(noise2_noDC_squared_m(:,trace));
end
Amplitude_noise2=mean(noise2_noDC_squared_m_sqrt);
SNR1=(Amplitude_signal/Amplitude_noise1)^2;
SNR2=(Amplitude_signal/Amplitude_noise2)^2;


                
        % global baseline (100ms before the onset of stim. train):
                baseline_global_interval(:,1) = stim2_X{x_value}(1,1)-(baseline_time_total./1000)*sf{1}:stim2_X{x_value}(1,1)-1;
                baseline_global(:,:,x_value) = current_data(baseline_global_interval,:,x_value);
                baseline_global_resid(:,:,x_value) = fn_Subtract_Mean(baseline_global(:,:,x_value));  
                baseline_global_average(1,:) = mean(baseline_global(:,:,x_value),1); %mean along time, within the trace, over 100ms prior to the first stim in the train
                baseline_global_M(:,1) = mean(baseline_global_average); %mean across traces
                baseline_global_STD(:,1) = mean(std(baseline_global(:,:,x_value),0,2)); %%mean across time of trial-to-trial variability (std)
                baseline_global_VAR(:,1) = mean(var(baseline_global(:,:,x_value),0,2)); %%mean across time of trial-to-trial variability (var)
                baseline_global_VAR_time(:,1) = mean(var(baseline_global(:,:,x_value),0,1)); %%mean of var across time (alonng the trace)
                
           for stim_num=1:size(stim2_X{x_value},2);
               
                response_begin=stim2_X{x_value}(1,stim_num)+(min_res_del./1000).*sf{1}; %time onset of whisker stim.+ minimum response latency.
                response_end=response_begin+(int_time-min_res_del).*sf{1}./1000-1;  
                %taking 50 ms of response according to stim. number in the train:
                data_response(:,:,x_value) = current_data(response_begin:response_end,:,x_value);
                data_response_resid(:,:,x_value) = fn_Subtract_Mean(data_response(:,:,x_value));  
                 Vm_res_M(:,stim_num) = mean(mean(data_response(:,:,x_value),2));  %mean along the trace and across traces (trial-to-trial)
                Vm_res_STD(:,stim_num) = mean(std(data_response(:,:,x_value),0,2)); %mean std across traces (trial-to-trial)
                Vm_res_VAR(:,stim_num) = mean(var(data_response(:,:,x_value),0,2)); %mean var across traces (trial-to-trial)
                Vm_res_VAR_time(:,stim_num) = mean(var(data_response(:,:,x_value),0,1)); %mean var across time (along trace)
                Vm_res_CV(:,stim_num) = mean(std(data_response(:,:,x_value),0,2)./mean(data_response(:,:,x_value),2)); %mean CV across traces (trial-to-trial)
                Vm_res_M_resid(:,stim_num) = mean(mean(data_response_resid(:,:,x_value),2)); %mean across traces (trial-to-trial)
                
                % local baseline (3ms before each peak):
                baseline_interval(:,1)=[stim2_X{x_value}(1,stim_num)-(baseline_time./1000).*sf{1}:stim2_X{x_value}(1,stim_num)-1];
                baseline_local_std(stim_num,:) = std(current_data(baseline_interval,:,x_value),0,1); %std along time, within the trace
                baseline_local_average(1,:) = mean(current_data(baseline_interval,:,x_value),1); %mean along time, within the trace
                baseline_local_M(stim_num,1)=mean(baseline_local_average); %mean across traces
                
                
%                 minpeak=mean(data_mean_smooth(baseline_interval))+2.*mean(baseline_local_std(stim_num,1)); %threshold criterion for minimal peak height
                [mean_peak_val_all{x_value,stim_num}, mean_peak_loc_all{x_value,stim_num}] =  findpeaks(data_mean_smooth(response_begin:response_end,1),'sortstr','descend');
%                 [mean_peak_val_all{x_value,stim_num}, mean_peak_loc_all{x_value,stim_num}] =  findpeaks(data_mean_smooth(response_begin:response_end,1),'minpeakheight',minpeak);               
%                 mean_peak_val{x_value}(stim_num,1)=mean_peak_val_all{x_value,stim_num}(1);
                mean_peak_loc{x_value}(stim_num,1)=mean_peak_loc_all{x_value,stim_num}(1);
                mean_peak_loc{x_value}(stim_num,1)=mean_peak_loc{x_value}(stim_num,1)+response_begin;
                mean_peak_val{x_value}(stim_num,1)=data_mean(mean_peak_loc{x_value}(stim_num,1));
                 mean_peak_amp_abs{x_value}(stim_num,1)=mean_peak_val{x_value}(stim_num,1)-baseline_global_M;
                mean_peak_amp{x_value}(stim_num,1)=mean_peak_val{x_value}(stim_num,1)-baseline_local_M(stim_num,1);
                
%                 max_deriv2_loc(:,stim_num)=find(data_deriv2(response_begin:response_end,1)==0);
                if mean_peak_amp{x_value}(stim_num,1)<0
%                      mean_peak_amp{x_value}(stim_num,1)=0;
                    mean_peak_amp{x_value}(stim_num,1)=nan;
                    mean_peak_rise_loc{x_value}(stim_num,1)=nan;
                    half_peak_val(stim_num,1)=nan;
                    mean_10per_loc{x_value}(stim_num,1) = nan;  
                    mean_90per_loc{x_value}(stim_num,1) = nan; 
                    mean_10per_val{x_value}(stim_num,1) = nan;
                    mean_90per_val{x_value}(stim_num,1) = nan;
                    mean_peak_half_peak_rise{x_value}(stim_num,1)=nan;
                    mean_peak_half_peak_fall{x_value}(stim_num,1)=nan;
                    mean_peak_half_peak_width{x_value}(stim_num)=nan;
                else
                                               
                %finding half-peak width and 10%peak value and latency from the mean trace:
                half_peak_val(stim_num,1)=baseline_local_M(stim_num,1)+mean_peak_amp{1,x_value}(stim_num)./2;
                mean_10per_val{x_value}(stim_num,1)=baseline_local_M(stim_num,1)+mean_peak_amp{1,x_value}(stim_num).*0.1;
                mean_90per_val{x_value}(stim_num,1)=baseline_local_M(stim_num,1)+mean_peak_amp{1,x_value}(stim_num).*0.9;
%                 mean_peak_rise_loc{x_value}(stim_num,1)=find(data_deriv1(response_begin:response_end,1)==max(data_deriv1(response_begin:response_end,1)));
%                 mean_peak_rise_loc{x_value}(stim_num,1)=mean_peak_rise_loc{x_value}(stim_num,1)+response_begin;
%                 mean_peak_rise_val{x_value}(stim_num,1)=data_mean(mean_peak_rise_loc{x_value}(stim_num,1),1);
                response_rise_interval=[response_begin:mean_peak_loc{x_value}(stim_num,1)];
                response_fall_interval=[mean_peak_loc{x_value}(stim_num,1):stim2_X{x_value}(1,stim_num)+(100./1000).*sf{1}];        
%                  if isempty(find(data_mean(response_rise_interval)>= mean_10per_val{x_value}(stim_num,1),1,'first'))
%                     mean_10per_loc{x_value}(stim_num,1) = nan;                   
%                  else                    
                    mean_10per_loc{x_value}(stim_num,1)=response_begin+find(data_mean(response_rise_interval)>= mean_10per_val{x_value}(stim_num,1),1,'first');
%                  end
%                 if isempty(find(data_mean(response_rise_interval)>= mean_10per_val{x_value}(stim_num,1),1,'first'))
%                     mean_10per_loc{x_value}(stim_num,1) = nan;                   
%                 else                     
                    mean_90per_loc{x_value}(stim_num,1)=response_begin+find(data_mean(response_rise_interval)>= mean_90per_val{x_value}(stim_num,1),1,'first');
%                 end
                mean_peak_half_peak_rise{x_value}(stim_num,1)=response_begin+find(data_mean(response_rise_interval)>=half_peak_val(stim_num,1),1,'first');
                    if isempty(find(data_mean(response_fall_interval)<=half_peak_val(stim_num,1),1,'first'))
                        mean_peak_half_peak_fall{x_value}(stim_num)=mean_peak_loc{x_value}(stim_num)+find(data_mean(response_fall_interval)==min(data_mean(response_fall_interval)),1,'first');
                    else                    
                        mean_peak_half_peak_fall{x_value}(stim_num,1)=mean_peak_loc{x_value}(stim_num,1)+find(data_mean(response_fall_interval)<=half_peak_val(stim_num,1),1,'first');
                    end
                mean_peak_half_peak_width{x_value}(stim_num)= (mean_peak_half_peak_fall{x_value}(stim_num)-mean_peak_half_peak_rise{x_value}(stim_num)).*dt; 
                end
                %putting nans in the place of the 11th stim if this cell doesn't have this data:
                if stim_num>size(stim2_X{x_value},2)
                    half_peak_val(stim_num,1)=nan;
                    mean_10per_val{x_value}(stim_num,1)=nan;
                    mean_90per_val{x_value}(stim_num,1)=nan;
                    mean_10per_loc{x_value}(stim_num,1)=nan;
                    mean_90per_loc{x_value}(stim_num,1)=nan;
                    mean_peak_half_peak_rise{x_value}(stim_num,1)=nan;
                    mean_peak_half_peak_fall{x_value}(stim_num,1)=nan;
                    mean_peak_half_peak_width{x_value}(stim_num)=nan;
                end
           end
 % parameters from single traces:               
                %try keti's way of finding max response in single traces.
%                loop on all traces, findpeaks on smoothed data without spikes
 %               [mean_peak_val_all{x_value,stim_num}, mean_peak_loc_all{x_value,stim_num}] =  findpeaks(data_mean_smooth(response_begin:response_end,1),'sortstr','descend');                        
%%    % verification of detected peaks in the raw data:
%      [Fig6,h6] = fn_Plot_Trace_v2([data_mean_smooth], dt, dt, stim1_X{channel}, dt, stim2_X{x_value}); %add to plot: raw_data{channel}(:,:,x_value)
%         ylabel('mean Vm [mV]', 'FontSize', 16);  legend('NB-', 'NB+'); legend('boxoff','Location','northeast');
%                     hold on 
% %                         for trace = 1:size(data_no_spikes{1},2)         
% %                              h1=plot([1:size(current_data,1)].*dt, current_data(:,trace,x_value),'b');
% %                              h2=scatter(mean_peak_loc{x_value}.*dt,current_data(mean_peak_loc{x_value},trace,x_value),'r','fill'); %mark the mean peaks on the single traces
%                              h2=scatter(mean_peak_loc{x_value}.*dt,data_mean(mean_peak_loc{x_value}),'r','fill'); %mark the mean peaks on the mean traces
% %                              h2=scatter(mean_peak_loc{x_value}.*dt,data_mean_smooth(mean_peak_loc{x_value}),'r', 'fill'); %mark the peaks on the smoothed mean trace
% %                              h3=scatter(mean_peak_half_peak_rise{x_value}.*dt,current_data(mean_peak_half_peak_rise{x_value},trace,x_value),'g','fill');
% %                              h4=scatter(mean_peak_half_peak_fall{x_value}.*dt,current_data(mean_peak_half_peak_fall{x_value},trace,x_value),'g','fill');                       
%                              h3=scatter(mean_peak_half_peak_rise{x_value}.*dt,data_mean(mean_peak_half_peak_rise{x_value}),'g','fill');
%                              h4=scatter(mean_peak_half_peak_fall{x_value}.*dt,data_mean(mean_peak_half_peak_fall{x_value}),'g','fill');                       
%                                 h5=scatter(stim2_X{x_value}(1,:).*dt,baseline_local_M,'c','fill');
%                                  h6=scatter(stim2_X{x_value}(1,:).*dt,baseline_global_M.*ones(size(stim2_X{x_value}(1,:))),'y','fill');
% %                                  h7=scatter(mean_peak_rise_loc{x_value}.*dt,data_mean_smooth(mean_peak_rise_loc{x_value}),'r', 'fill'); %mark the peaks on the mean trace
%                                 h7=scatter(mean_10per_loc{x_value}.*dt,data_mean(mean_10per_loc{x_value}),'r','fill'); %mark the mean peaks on the mean traces
%                                  pause
% %                              delete(h1)
%                              delete(h2)
%                              delete(h3)
%                              delete(h4)
%                              delete(h5)
%                              delete(h6)
%                              delete(h7)
%                              figure(Fig6) 
% %                          end
%                       hold off
%                         close (gcf)
%% power spectrum of sensory response for adaptation F1
fstim=Param.facade(7); %train stim frequency
nstim=Param.facade(6); %number of stim in the train
stim1_duration=1000./fstim./1000; %1000ms/galvano frequency converted to [sec]
% duration = stim1_duration.*nstim; %stim duration multiplied by the number of stim
DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = []; interval = []; expected_mat=[];

start_sample = stim2_X{x_value}(1,1); %[sec]
end_sample = start_sample+duration.*sf{channel}-1;
interval = start_sample:end_sample; %taking the interval of the sensory stim train
spec_mat = current_data(interval,:,x_value);
expected_mat=repmat(spec_mat(1:stim1_duration.*sf{channel},:),nstim,1); %replicating the response for the 1st stim in the train
 [spec_mat_noDC] = fn_Subtract_Mean(spec_mat);
 [expected_mat_noDC] = fn_Subtract_Mean(expected_mat);
% DC= mean(spec_mat,1);
% l=size(spec_mat,1);
% l=l/2-1;
% DC=wextend('addrow', 'per', DC, l);
% spec_mat_noDC = spec_mat-DC;
[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf{channel},0,0);
[Y_abs_exp(:,:,t),f_exp(:,t)] = fn_Amp_Spectrum(expected_mat_noDC,sf{channel},0,0);
max_res(1,t)=max(mean(Y_abs(f(:,t)>(fstim-1) &f(:,t)<(fstim+1),:,t),2));
max_res_exp(1,t)=max(mean(Y_abs_exp(f_exp(:,t)>(fstim-1) &f_exp(:,t)<(fstim+1) ,:,t),2));
F1(1,t)=max_res(1,t)./max_res_exp(1,t);

% figure(10) %(fileind)
%       hold on
%         plot(f(:,t),mean(Y_abs(:,:,t),2),'color', color_table(t,:)) 
%         plot(f_exp(:,t),mean(Y_abs_exp(:,:,t),2),'color', color_table(t+2,:)) 
%       hold off
%         xlim([0 50]); ylim([-0.01 25])
% 
%         set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2)
%         xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
%         ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
%         title(['Response Power file ', num2str(files_to_analyze(fileind))],'FontSize', 12); 
%         leg=legend('NB-','expected NB-','NB+','expected NB+');
%         set(leg,'box','off', 'FontSize', 12);
%%
%Parameters from mean trace                    peak_10per_lat    
peaks(fileind).fname = fname; 
peaks(fileind).cells=files_to_analyze;
peaks(fileind).analysis_mfile='Analyze_response_parameters_10_20Hz_v2.m'; 
peaks(fileind).BP50HzLFP=BP50HzLFP_flag; %removing 50Hz noise from LFP signal
peaks(fileind).BP50HzVm=BP50HzVm_flag; %removing 50Hz noise from Vm signal
peaks(fileind).BPLFP_flag=BPLFP_flag; %filtering LFP 
peaks(fileind).BPLFP=bp_filt_LFP; %the BP frequency filter for  LFP, used if BPLFP_flag=1
peaks(fileind).BPVm_flag=BPVm_flag; %filtering Vm 
peaks(fileind).BPVm=bp_filt_Vm; %the BP frequency filter for Vm, used if BPVm_flag=1
peaks(fileind).sf=sf{1};
peaks(fileind).min_res_del=min_res_del; %start looking for response after this min delay from stim onset.
peaks(fileind).baseline_time_local=baseline_time; %baseline before each whisker stim.
peaks(fileind).baseline_time_global=baseline_time_total; %global baseline for the whole mean trace
peaks(fileind).baseline_local(:,x_value)=baseline_local_M; 
peaks(fileind).baseline_global(:,x_value)=baseline_global_M; 
peaks(fileind).baseline_global_std(:,x_value)=baseline_global_STD;
peaks(fileind).baseline_global_VAR_time(:,x_value)=baseline_global_VAR_time;
peaks(fileind).pre_response_M(:,x_value) = pre_response_M(:,1); 
peaks(fileind).pre_response_STD(:,x_value) = pre_response_STD(:,1); %mean std across traces (trial-to-trial) of 50ms prior to train
peaks(fileind).pre_response_VAR(:,x_value) = pre_response_VAR(:,1); %mean var across traces (trial-to-trial) of 50ms prior to train
peaks(fileind).pre_response_VAR_time(:,x_value) = pre_response_VAR_time(:,1);
peaks(fileind).pre_response_CV(:,x_value) = pre_response_CV(:,1); 
peaks(fileind).Vm_res_M(:,x_value)=Vm_res_M;
peaks(fileind).Vm_res_M_meanstim(:,x_value)=mean(peaks(fileind).Vm_res_M(:,x_value)); %trial-to-trial variability averaged over stim.
peaks(fileind).Vm_res_STD(:,x_value)=Vm_res_STD;
peaks(fileind).Vm_res_VAR(:,x_value)=Vm_res_VAR;
peaks(fileind).Vm_res_VAR_meanstim(:,x_value)=mean(peaks(fileind).Vm_res_VAR(:,x_value)); %trial-to-trial variability averaged over stim.
peaks(fileind).Vm_res_VAR_time(:,x_value)=Vm_res_VAR_time; %mean across time (trace)
peaks(fileind).Vm_res_STD_meanstim(:,x_value)=mean(peaks(fileind).Vm_res_STD(:,x_value));
peaks(fileind).Vm_res_CV(:,x_value)=Vm_res_CV;
peaks(fileind).Vm_res_M_resid(:,x_value)=Vm_res_M_resid;
peaks(fileind).Vm_train_interval(:,x_value)=train_interval;
peaks(fileind).Vm_res_train_M(:,x_value)=Vm_res_train_M;
peaks(fileind).Vm_res_train_STD(:,x_value)=Vm_res_train_STD;
peaks(fileind).Vm_res_train_VAR(:,x_value)=Vm_res_train_VAR;
peaks(fileind).Vm_res_train_VAR_time(:,x_value)=Vm_res_train_VAR_time;
peaks(fileind).Vm_res_train_CV(:,x_value)=Vm_res_train_CV;
peaks(fileind).Vm_res_train_M_resid(:,x_value)=Vm_res_train_M_resid;
peaks(fileind).post_train_interval(:,x_value)=post_train_interval;
peaks(fileind).post_train_M(:,x_value)=post_train_M;
peaks(fileind).post_train_STD(:,x_value)=post_train_STD;
peaks(fileind).post_train_VAR(:,x_value)=post_train_VAR;
peaks(fileind).post_train_VAR_time(:,x_value)=post_train_VAR_time;
peaks(fileind).post_train_CV(:,x_value)=post_train_CV;
peaks(fileind).post_train_M_resid(:,x_value)=post_train_M_resid;
peaks(fileind).val(:,x_value) = mean_peak_val{1,x_value};
peaks(fileind).amp(:,x_value) = mean_peak_amp{1,x_value};
peaks(fileind).amp_abs(:,x_value) = mean_peak_amp_abs{1,x_value};
peaks(fileind).loc(:,x_value) = mean_peak_loc{1,x_value};
peaks(fileind).latency(:,x_value)=(mean_peak_loc{1,x_value}-(stim2_X{x_value}(1,:))').*dt.*1000; %[msec]
peaks(fileind).per10_val(:,x_value) = mean_10per_val{x_value}(1,:);
peaks(fileind).per10_lat(:,x_value)= (mean_10per_loc{x_value}-(stim2_X{x_value}(1,:))').*dt.*1000; %[msec]
peaks(fileind).per90_val(:,x_value) = mean_90per_val{x_value}(1,:);
peaks(fileind).per90_lat(:,x_value)= (mean_90per_loc{x_value}-(stim2_X{x_value}(1,:))').*dt.*1000; %[msec]
peaks(fileind).half_width(:,x_value) = mean_peak_half_peak_width{x_value}.*1000;
peaks(fileind).adapt_amp(:,x_value)=mean(mean_peak_amp{1,x_value}(8:10))./mean_peak_amp{1,x_value}(1);
peaks(fileind).SNR_var(:,x_value)=mean(peaks(fileind).Vm_res_VAR_time(:,x_value))/peaks(fileind).pre_response_VAR_time(:,x_value); %SNR=SIGMAsignal/SIGMAnoise. SIGMAsignal is the mean of the variances in each response in the train
peaks(fileind).SNR1(:,x_value)=SNR1;
peaks(fileind).SNR2(:,x_value)=SNR2;
peaks(fileind).Amplitude_signal(:,x_value)=Amplitude_signal;
peaks(fileind).Amplitude_noise1(:,x_value)=Amplitude_noise1;
peaks(fileind).Amplitude_noise2(:,x_value)=Amplitude_noise2;

peaks(fileind).fft_max_res(:,x_value)=max_res(1,t);
peaks(fileind).fft_max_res_exp(:,x_value)=max_res_exp(1,t);
peaks(fileind).F1(:,x_value)=F1(1,t);
                
         for stim_num=1:size(stim2_X{x_value},2);
            peak_val{x_value}(fileind,stim_num)=peaks(fileind).val(stim_num,x_value);
            peak_amp{x_value}(fileind,stim_num)=peaks(fileind).amp(stim_num,x_value);
            peak_amp_abs{x_value}(fileind,stim_num)=peaks(fileind).amp_abs(stim_num,x_value);
            peak_lat{x_value}(fileind,stim_num)=peaks(fileind).latency(stim_num,x_value);
            peak_half_width{x_value}(fileind,stim_num)=peaks(fileind).half_width(stim_num,x_value);
            peak_per10_lat{x_value}(fileind,stim_num)=peaks(fileind).per10_lat(stim_num,x_value);
            peak_per90_lat{x_value}(fileind,stim_num)=peaks(fileind).per90_lat(stim_num,x_value);
            peak_baseline_local{x_value}(fileind,stim_num)=peaks(fileind).baseline_local(stim_num,x_value);
            peak_pre_response_STD{x_value}(fileind,stim_num)=peaks(fileind).pre_response_STD(:,x_value); %mean std across traces (trial-to-trial) of 50ms prior to train
            peak_pre_response_VAR{x_value}(fileind,stim_num)=peaks(fileind).pre_response_VAR(:,x_value); % mean variance across traces
            peak_pre_response_VAR_time{x_value}(fileind,stim_num)=peaks(fileind).pre_response_VAR_time(:,x_value); 
            peak_pre_response_CV{x_value}(fileind,stim_num)=peaks(fileind).pre_response_CV(:,x_value); 
            peak_pre_response_M{x_value}(fileind,stim_num)=peaks(fileind).pre_response_M(:,x_value); 
            peak_Vm_res_M{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_M(stim_num,x_value);
             peak_Vm_res_M_meanstim{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_M_meanstim(:,x_value); %mean across time and stim.
            peak_Vm_res_STD{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_STD(stim_num,x_value);
            peak_Vm_res_STD_meanstim{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_STD_meanstim(:,x_value); %mean across trials and stim
            peak_Vm_res_VAR{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_VAR(stim_num,x_value); %mean across time and stim.
            peak_Vm_res_VAR_meanstim{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_VAR_meanstim(:,x_value); %mean across time and stim.
             peak_Vm_res_VAR_time{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_VAR_time(stim_num,x_value);
            peak_Vm_res_CV{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_CV(stim_num,x_value);
            peak_Vm_res_M_resid{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_M_resid(stim_num,x_value);
            peak_Vm_res_train_M{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_train_M(:,x_value);
            peak_Vm_res_train_STD{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_train_STD(:,x_value);
            peak_Vm_res_train_VAR{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_train_VAR(:,x_value);
            peak_Vm_res_train_VAR_time{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_train_VAR_time(:,x_value);
            peak_Vm_res_train_CV{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_train_CV(:,x_value);
            peak_Vm_res_train_M_resid{x_value}(fileind,stim_num)=peaks(fileind).Vm_res_train_M_resid(:,x_value);
            peak_post_train_M{x_value}(fileind,stim_num)=peaks(fileind).post_train_M(:,x_value);
            peak_post_train_STD{x_value}(fileind,stim_num)=peaks(fileind).post_train_STD(:,x_value);
            peak_post_train_VAR{x_value}(fileind,stim_num)=peaks(fileind).post_train_VAR(:,x_value);
            peak_post_train_VAR_time{x_value}(fileind,stim_num)=peaks(fileind).post_train_VAR_time(:,x_value);
            peak_post_train_CV{x_value}(fileind,stim_num)=peaks(fileind).post_train_CV(:,x_value);
            peak_post_train_M_resid{x_value}(fileind,stim_num)=peaks(fileind).post_train_M_resid(:,x_value);
            peak_SNR_var{x_value}(fileind,stim_num)=peaks(fileind).SNR_var(:,x_value);
            peak_SNR1{x_value}(fileind,stim_num)=peaks(fileind).SNR1(:,x_value);
            peak_SNR2{x_value}(fileind,stim_num)=peaks(fileind).SNR2(:,x_value);
            peak_Amplitude_signal{x_value}(fileind,stim_num)=peaks(fileind).Amplitude_signal(:,x_value);
            peak_Amplitude_noise1{x_value}(fileind,stim_num)=peaks(fileind).Amplitude_noise1(:,x_value);
            peak_Amplitude_noise2{x_value}(fileind,stim_num)=peaks(fileind).Amplitude_noise2(:,x_value);
            
            peak_adapt_amp{x_value}(fileind,stim_num)=peaks(fileind).adapt_amp(1,x_value);
            peak_F1{x_value}(fileind,stim_num)=peaks(fileind).F1(1,x_value);
            peak_baseline_global{x_value}(fileind,stim_num)=peaks(fileind).baseline_global(x_value);
            peak_baseline_global_std{x_value}(fileind,stim_num)=peaks(fileind).baseline_global_std(x_value);
            peak_baseline_global_VAR_time{x_value}(fileind,stim_num)=peaks(fileind).baseline_global_VAR_time(x_value);
           
         end
end
% close all
clear data_no_spike_no_DC
    end %end of files loop
 
    %% Statistics for peaks
%t-tests between cells      
%relative change=(X-Xreference)/|Xreference|
    for stim_num=1:11; %1:10;
        %peak values - absolute values + relative change
        peaks_stat(stim_num).val=[peak_val{2}(:,stim_num),  peak_val{3}(:,stim_num)];
        peaks_stat(stim_num).val_m=nanmean(peaks_stat(stim_num).val,1);
        peaks_stat(stim_num).val_std=nanstd(peaks_stat(stim_num).val,0,1);
        peaks_stat(stim_num).change_val=[(peak_val{3}(:,stim_num)-peak_val{2}(:,stim_num))./abs(peak_val{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_val_m=nanmean(peaks_stat(stim_num).change_val,1);
        peaks_stat(stim_num).change_val_std=nanstd(peaks_stat(stim_num).change_val,0,1);
        change_val_mat(:,stim_num)= peaks_stat(stim_num).change_val;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_val, peaks_stat(stim_num).lillietest_p_val] = lillietest(peaks_stat(stim_num).val(:,2)- peaks_stat(stim_num).val(:,1));
        [peaks_stat(stim_num).lillietest_h_change_val, peaks_stat(stim_num).lillietest_p_change_val] = lillietest( peaks_stat(stim_num).change_val);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_val, peaks_stat(stim_num).ttest_p_val]= ttest(peaks_stat(stim_num).val(:,1),peaks_stat(stim_num).val(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_val, peaks_stat(stim_num).wilcoxon_h_val]= signrank(peaks_stat(stim_num).val(:,1),peaks_stat(stim_num).val(:,2));
        [peaks_stat(stim_num).ttest_h_change_val, peaks_stat(stim_num).ttest_p_change_val]= ttest(peaks_stat(stim_num).change_val);
        [peaks_stat(stim_num).wilcoxon_p_change_val, peaks_stat(stim_num).wilcoxon_h_change_val]= signrank(peaks_stat(stim_num).change_val);
        
%peak amplitudes - absolute values + relative change
        
%         if stim_num==9  %outlyer
%             peak_amp{2}(1,stim_num)=nan;
%             peak_amp{3}(1,stim_num)=nan;
%         end
%           if stim_num==8  %outlyer
%             peak_amp{2}(8,stim_num)=nan;
%             peak_amp{3}(8,stim_num)=nan;
%         end
        peaks_stat(stim_num).amp=[peak_amp{2}(:,stim_num),  peak_amp{3}(:,stim_num)];
        peaks_stat(stim_num).amp_m=nanmean(peaks_stat(stim_num).amp,1);
        peaks_stat(stim_num).amp_std=nanstd(peaks_stat(stim_num).amp,0,1);
        peaks_stat(stim_num).change_amp=[(peak_amp{3}(:,stim_num)-peak_amp{2}(:,stim_num))./abs(peak_amp{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_amp_m=nanmean(peaks_stat(stim_num).change_amp,1);
        peaks_stat(stim_num).change_amp_std=nanstd(peaks_stat(stim_num).change_amp,0,1);
        change_amp_mat(:,stim_num)= peaks_stat(stim_num).change_amp;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_amp, peaks_stat(stim_num).lillietest_p_amp] = lillietest(peaks_stat(stim_num).amp(:,2)- peaks_stat(stim_num).amp(:,1));
        [peaks_stat(stim_num).lillietest_h_change_amp, peaks_stat(stim_num).lillietest_p_change_amp] = lillietest( peaks_stat(stim_num).change_amp);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_amp, peaks_stat(stim_num).ttest_p_amp]= ttest(peaks_stat(stim_num).amp(:,1),peaks_stat(stim_num).amp(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_amp, peaks_stat(stim_num).wilcoxon_h_amp]= signrank(peaks_stat(stim_num).amp(:,1),peaks_stat(stim_num).amp(:,2));
        [peaks_stat(stim_num).ttest_h_change_amp, peaks_stat(stim_num).ttest_p_change_amp]= ttest(peaks_stat(stim_num).change_amp);
        [peaks_stat(stim_num).wilcoxon_p_change_amp, peaks_stat(stim_num).wilcoxon_h_change_amp]= signrank(peaks_stat(stim_num).change_amp);
       
        %peak absolute amplitude - absolute values + relative change
        peaks_stat(stim_num).amp_abs=[peak_amp_abs{2}(:,stim_num),  peak_amp_abs{3}(:,stim_num)];
        peaks_stat(stim_num).amp_abs_m=nanmean(peaks_stat(stim_num).amp_abs,1);
        peaks_stat(stim_num).amp_abs_std=nanstd(peaks_stat(stim_num).amp_abs,0,1);
        peaks_stat(stim_num).change_amp_abs=[(peak_amp_abs{3}(:,stim_num)-peak_amp_abs{2}(:,stim_num))./abs(peak_amp_abs{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_amp_abs_m=nanmean(peaks_stat(stim_num).change_amp_abs,1);
        peaks_stat(stim_num).change_amp_abs_std=nanstd(peaks_stat(stim_num).change_amp_abs,0,1);
         change_amp_abs_mat(:,stim_num)= peaks_stat(stim_num).change_amp_abs;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_amp_abs, peaks_stat(stim_num).lillietest_p_amp_abs] = lillietest(peaks_stat(stim_num).amp_abs(:,2)- peaks_stat(stim_num).amp_abs(:,1));
        [peaks_stat(stim_num).lillietest_h_change_amp_abs, peaks_stat(stim_num).lillietest_p_change_amp_abs] = lillietest( peaks_stat(stim_num).change_amp_abs);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_amp_abs, peaks_stat(stim_num).ttest_p_amp_abs]= ttest(peaks_stat(stim_num).amp_abs(:,1),peaks_stat(stim_num).amp_abs(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_amp_abs, peaks_stat(stim_num).wilcoxon_h_amp_abs]= signrank(peaks_stat(stim_num).amp_abs(:,1),peaks_stat(stim_num).amp_abs(:,2));
        [peaks_stat(stim_num).ttest_h_change_amp_abs, peaks_stat(stim_num).ttest_p_change_amp_abs]= ttest(peaks_stat(stim_num).change_amp_abs);
        [peaks_stat(stim_num).wilcoxon_p_change_amp_abs, peaks_stat(stim_num).wilcoxon_h_change_amp_abs]= signrank(peaks_stat(stim_num).change_amp_abs);
        
        
        %peak latency - absolute values + relative change
        peaks_stat(stim_num).latency=[peak_lat{2}(:,stim_num),  peak_lat{3}(:,stim_num)];
        peaks_stat(stim_num).latency_m=nanmean(peaks_stat(stim_num).latency,1);
        peaks_stat(stim_num).latency_std=nanstd(peaks_stat(stim_num).latency,0,1);
        peaks_stat(stim_num).change_latency=[(peak_lat{3}(:,stim_num)-peak_lat{2}(:,stim_num))./abs(peak_lat{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_latency_m=nanmean(peaks_stat(stim_num).change_latency,1);
        peaks_stat(stim_num).change_latency_std=nanstd(peaks_stat(stim_num).change_latency,0,1);
         change_latency_mat(:,stim_num)= peaks_stat(stim_num).change_latency;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_latency, peaks_stat(stim_num).lillietest_p_latency] = lillietest(peaks_stat(stim_num).latency(:,2)- peaks_stat(stim_num).latency(:,1));
        [peaks_stat(stim_num).lillietest_h_change_latency, peaks_stat(stim_num).lillietest_p_change_latency] = lillietest( peaks_stat(stim_num).change_latency);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_latency, peaks_stat(stim_num).ttest_p_latency]= ttest(peaks_stat(stim_num).latency(:,1),peaks_stat(stim_num).latency(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_latency, peaks_stat(stim_num).wilcoxon_h_latency]= signrank(peaks_stat(stim_num).latency(:,1),peaks_stat(stim_num).latency(:,2));
        [peaks_stat(stim_num).ttest_h_change_latency, peaks_stat(stim_num).ttest_p_change_latency]= ttest(peaks_stat(stim_num).change_latency);
        [peaks_stat(stim_num).wilcoxon_p_change_latency, peaks_stat(stim_num).wilcoxon_h_change_latency]= signrank(peaks_stat(stim_num).change_latency);
               
       %peak half-width - absolute values + relative change
        peaks_stat(stim_num).half_width=[peak_half_width{2}(:,stim_num),  peak_half_width{3}(:,stim_num)];
        peaks_stat(stim_num).half_width_m=nanmean(peaks_stat(stim_num).half_width,1);
        peaks_stat(stim_num).half_width_std=nanstd(peaks_stat(stim_num).half_width,0,1);
        peaks_stat(stim_num).change_half_width=[(peak_half_width{3}(:,stim_num)-peak_half_width{2}(:,stim_num))./abs(peak_half_width{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_half_width_m=nanmean(peaks_stat(stim_num).change_half_width,1);
        peaks_stat(stim_num).change_half_width_std=nanstd(peaks_stat(stim_num).change_half_width,0,1);
         change_half_width_mat(:,stim_num)= peaks_stat(stim_num).change_half_width;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_half_width, peaks_stat(stim_num).lillietest_p_half_width] = lillietest(peaks_stat(stim_num).half_width(:,2)- peaks_stat(stim_num).half_width(:,1));
        [peaks_stat(stim_num).lillietest_h_change_half_width, peaks_stat(stim_num).lillietest_p_change_half_width] = lillietest( peaks_stat(stim_num).change_half_width);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_half_width, peaks_stat(stim_num).ttest_p_half_width]= ttest(peaks_stat(stim_num).half_width(:,1),peaks_stat(stim_num).half_width(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_half_width, peaks_stat(stim_num).wilcoxon_h_half_width]= signrank(peaks_stat(stim_num).half_width(:,1),peaks_stat(stim_num).half_width(:,2));
        [peaks_stat(stim_num).ttest_h_change_half_width, peaks_stat(stim_num).ttest_p_change_half_width]= ttest(peaks_stat(stim_num).change_half_width);
        [peaks_stat(stim_num).wilcoxon_p_change_half_width, peaks_stat(stim_num).wilcoxon_h_change_half_width]= signrank(peaks_stat(stim_num).change_half_width);
        
%peak 10% latency (onset latency) - absolute values + relative change
        peaks_stat(stim_num).per10_lat=[peak_per10_lat{2}(:,stim_num),  peak_per10_lat{3}(:,stim_num)];
        peaks_stat(stim_num).per10_lat_m=nanmean(peaks_stat(stim_num).per10_lat,1);
        peaks_stat(stim_num).per10_lat_std=nanstd(peaks_stat(stim_num).per10_lat,0,1);
        peaks_stat(stim_num).change_per10_lat=[(peak_per10_lat{3}(:,stim_num)-peak_per10_lat{2}(:,stim_num))./abs(peak_per10_lat{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_per10_lat_m=nanmean(peaks_stat(stim_num).change_per10_lat,1);
        peaks_stat(stim_num).change_per10_lat_std=nanstd(peaks_stat(stim_num).change_per10_lat,0,1);
         change_per10_lat_mat(:,stim_num)= peaks_stat(stim_num).change_per10_lat;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_per10_lat, peaks_stat(stim_num).lillietest_p_per10_lat] = lillietest(peaks_stat(stim_num).per10_lat(:,2)- peaks_stat(stim_num).per10_lat(:,1));
        [peaks_stat(stim_num).lillietest_h_change_per10_lat, peaks_stat(stim_num).lillietest_p_change_per10_lat] = lillietest( peaks_stat(stim_num).change_per10_lat);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_per10_lat, peaks_stat(stim_num).ttest_p_per10_lat]= ttest(peaks_stat(stim_num).per10_lat(:,1),peaks_stat(stim_num).per10_lat(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_per10_lat, peaks_stat(stim_num).wilcoxon_h_per10_lat]= signrank(peaks_stat(stim_num).per10_lat(:,1),peaks_stat(stim_num).per10_lat(:,2));
        [peaks_stat(stim_num).ttest_h_change_per10_lat, peaks_stat(stim_num).ttest_p_change_per10_lat]= ttest(peaks_stat(stim_num).change_per10_lat);
        [peaks_stat(stim_num).wilcoxon_p_change_per10_lat, peaks_stat(stim_num).wilcoxon_h_change_per10_lat]= signrank(peaks_stat(stim_num).change_per10_lat);
                
        %peak local baseline - absolute values + relative change
        peaks_stat(stim_num).baseline_local=[peak_baseline_local{2}(:,stim_num),  peak_baseline_local{3}(:,stim_num)];
        peaks_stat(stim_num).baseline_local_m=nanmean(peaks_stat(stim_num).baseline_local,1);
        peaks_stat(stim_num).baseline_local_std=nanstd(peaks_stat(stim_num).baseline_local,0,1);
        peaks_stat(stim_num).change_baseline_local=[(peak_baseline_local{3}(:,stim_num)-peak_baseline_local{2}(:,stim_num))./abs(peak_baseline_local{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_baseline_local_m=nanmean(peaks_stat(stim_num).change_baseline_local,1);
        peaks_stat(stim_num).change_baseline_local_std=nanstd(peaks_stat(stim_num).change_baseline_local,0,1);
         change_baseline_local_mat(:,stim_num)= peaks_stat(stim_num).change_baseline_local;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_baseline_local, peaks_stat(stim_num).lillietest_p_baseline_local] = lillietest(peaks_stat(stim_num).baseline_local(:,2)- peaks_stat(stim_num).baseline_local(:,1));
        [peaks_stat(stim_num).lillietest_h_change_baseline_local, peaks_stat(stim_num).lillietest_p_change_baseline_local] = lillietest( peaks_stat(stim_num).change_baseline_local);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_baseline_local, peaks_stat(stim_num).ttest_p_baseline_local]= ttest(peaks_stat(stim_num).baseline_local(:,1),peaks_stat(stim_num).baseline_local(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_baseline_local, peaks_stat(stim_num).wilcoxon_h_baseline_local]= signrank(peaks_stat(stim_num).baseline_local(:,1),peaks_stat(stim_num).baseline_local(:,2));
        [peaks_stat(stim_num).ttest_h_change_baseline_local, peaks_stat(stim_num).ttest_p_change_baseline_local]= ttest(peaks_stat(stim_num).change_baseline_local);
        [peaks_stat(stim_num).wilcoxon_p_change_baseline_local, peaks_stat(stim_num).wilcoxon_h_change_baseline_local]= signrank(peaks_stat(stim_num).change_baseline_local);
        
        %global baseline - absolute values + relative change
        peaks_stat(stim_num).baseline_global=[peak_baseline_global{2}(:,stim_num),  peak_baseline_global{3}(:,stim_num)];
        peaks_stat(stim_num).baseline_global_m=nanmean(peaks_stat(stim_num).baseline_global,1);
        peaks_stat(stim_num).baseline_global_std=nanstd(peaks_stat(stim_num).baseline_global,0,1);
        peaks_stat(stim_num).change_baseline_global=[(peak_baseline_global{3}(:,stim_num)-peak_baseline_global{2}(:,stim_num))./abs(peak_baseline_global{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_baseline_global_m=nanmean(peaks_stat(stim_num).change_baseline_global,1);
        peaks_stat(stim_num).change_baseline_global_std=nanstd(peaks_stat(stim_num).change_baseline_global,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_baseline_global, peaks_stat(stim_num).lillietest_p_baseline_global] = lillietest(peaks_stat(stim_num).baseline_global(:,2)- peaks_stat(stim_num).baseline_global(:,1));
        [peaks_stat(stim_num).lillietest_h_change_baseline_global, peaks_stat(stim_num).lillietest_p_change_baseline_global] = lillietest( peaks_stat(stim_num).change_baseline_global);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_baseline_global, peaks_stat(stim_num).ttest_p_baseline_global]= ttest(peaks_stat(stim_num).baseline_global(:,1),peaks_stat(stim_num).baseline_global(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_baseline_global, peaks_stat(stim_num).wilcoxon_h_baseline_global]= signrank(peaks_stat(stim_num).baseline_global(:,1),peaks_stat(stim_num).baseline_global(:,2));
        [peaks_stat(stim_num).ttest_h_change_baseline_global, peaks_stat(stim_num).ttest_p_change_baseline_global]= ttest(peaks_stat(stim_num).change_baseline_global);
        [peaks_stat(stim_num).wilcoxon_p_change_baseline_global, peaks_stat(stim_num).wilcoxon_h_change_baseline_global]= signrank(peaks_stat(stim_num).change_baseline_global);
        
        %F1 - absolute values + relative change
        peaks_stat(stim_num).F1=[peak_F1{2}(:,stim_num),  peak_F1{3}(:,stim_num)];
        peaks_stat(stim_num).F1_m=nanmean(peaks_stat(stim_num).F1,1);
        peaks_stat(stim_num).F1_std=nanstd(peaks_stat(stim_num).F1,0,1);
        peaks_stat(stim_num).change_F1=[(peak_F1{3}(:,stim_num)-peak_F1{2}(:,stim_num))./abs(peak_F1{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_F1_m=nanmean(peaks_stat(stim_num).change_F1,1);
        peaks_stat(stim_num).change_F1_std=nanstd(peaks_stat(stim_num).change_F1,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_F1, peaks_stat(stim_num).lillietest_p_F1] = lillietest(peaks_stat(stim_num).F1(:,2)- peaks_stat(stim_num).F1(:,1));
        [peaks_stat(stim_num).lillietest_h_change_F1, peaks_stat(stim_num).lillietest_p_change_F1] = lillietest( peaks_stat(stim_num).change_F1);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_F1, peaks_stat(stim_num).ttest_p_F1]= ttest(peaks_stat(stim_num).F1(:,1),peaks_stat(stim_num).F1(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_F1, peaks_stat(stim_num).wilcoxon_h_F1]= signrank(peaks_stat(stim_num).F1(:,1),peaks_stat(stim_num).F1(:,2));
        [peaks_stat(stim_num).ttest_h_change_F1, peaks_stat(stim_num).ttest_p_change_F1]= ttest(peaks_stat(stim_num).change_F1);
        [peaks_stat(stim_num).wilcoxon_p_change_F1, peaks_stat(stim_num).wilcoxon_h_change_F1]= signrank(peaks_stat(stim_num).change_F1);
                        
        %peak adaptation amplitude ratio - absolute values + relative change
        peaks_stat(stim_num).adapt_amp=[peak_adapt_amp{2}(:,stim_num),  peak_adapt_amp{3}(:,stim_num)];
        if exp_type==1
            peaks_stat(stim_num).adapt_amp([1,3,5,6,8,9],:)=[];
        end
        peaks_stat(stim_num).adapt_amp_m=nanmean(peaks_stat(stim_num).adapt_amp,1);
        peaks_stat(stim_num).adapt_amp_std=nanstd(peaks_stat(stim_num).adapt_amp,0,1);
        peaks_stat(stim_num).change_adapt_amp=[(peak_adapt_amp{3}(:,stim_num)-peak_adapt_amp{2}(:,stim_num))./abs(peak_adapt_amp{2}(:,stim_num))].*100; %percent change
        if exp_type==1
            peaks_stat(stim_num).change_adapt_amp([1,3,5,6,8,9],:)=[];
        end
        peaks_stat(stim_num).change_adapt_amp_m=nanmean(peaks_stat(stim_num).change_adapt_amp,1);
        peaks_stat(stim_num).change_adapt_amp_std=nanstd(peaks_stat(stim_num).change_adapt_amp,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_adapt_amp, peaks_stat(stim_num).lillietest_p_adapt_amp] = lillietest(peaks_stat(stim_num).adapt_amp(:,2)- peaks_stat(stim_num).adapt_amp(:,1));
        [peaks_stat(stim_num).lillietest_h_change_adapt_amp, peaks_stat(stim_num).lillietest_p_change_adapt_amp] = lillietest( peaks_stat(stim_num).change_adapt_amp);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_adapt_amp, peaks_stat(stim_num).ttest_p_adapt_amp]= ttest(peaks_stat(stim_num).adapt_amp(:,1),peaks_stat(stim_num).adapt_amp(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_adapt_amp, peaks_stat(stim_num).wilcoxon_h_adapt_amp]= signrank(peaks_stat(stim_num).adapt_amp(:,1),peaks_stat(stim_num).adapt_amp(:,2));
        [peaks_stat(stim_num).ttest_h_change_adapt_amp, peaks_stat(stim_num).ttest_p_change_adapt_amp]= ttest(peaks_stat(stim_num).change_adapt_amp);
        [peaks_stat(stim_num).wilcoxon_p_change_adapt_amp, peaks_stat(stim_num).wilcoxon_h_change_adapt_amp]= signrank(peaks_stat(stim_num).change_adapt_amp);
        
 %global baseline std - absolute values + relative change
        peaks_stat(stim_num).baseline_global_std=[peak_baseline_global_std{2}(:,stim_num),  peak_baseline_global_std{3}(:,stim_num)];
        peaks_stat(stim_num).baseline_global_std_m=nanmean(peaks_stat(stim_num).baseline_global_std,1);
        peaks_stat(stim_num).baseline_global_std_std=nanstd(peaks_stat(stim_num).baseline_global_std,0,1);
        peaks_stat(stim_num).change_baseline_global_std=[(peak_baseline_global_std{3}(:,stim_num)-peak_baseline_global_std{2}(:,stim_num))./abs(peak_baseline_global_std{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_baseline_global_std_m=nanmean(peaks_stat(stim_num).change_baseline_global_std,1);
        peaks_stat(stim_num).change_baseline_global_std_std=nanstd(peaks_stat(stim_num).change_baseline_global_std,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_baseline_global_std, peaks_stat(stim_num).lillietest_p_baseline_global_std] = lillietest(peaks_stat(stim_num).baseline_global_std(:,2)- peaks_stat(stim_num).baseline_global_std(:,1));
        [peaks_stat(stim_num).lillietest_h_change_baseline_global_std, peaks_stat(stim_num).lillietest_p_change_baseline_global_std] = lillietest( peaks_stat(stim_num).change_baseline_global_std);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_baseline_global_std, peaks_stat(stim_num).ttest_p_baseline_global_std]= ttest(peaks_stat(stim_num).baseline_global_std(:,1),peaks_stat(stim_num).baseline_global_std(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_baseline_global_std, peaks_stat(stim_num).wilcoxon_h_baseline_global_std]= signrank(peaks_stat(stim_num).baseline_global_std(:,1),peaks_stat(stim_num).baseline_global_std(:,2));
        [peaks_stat(stim_num).ttest_h_change_baseline_global_std, peaks_stat(stim_num).ttest_p_change_baseline_global_std]= ttest(peaks_stat(stim_num).change_baseline_global_std);
        [peaks_stat(stim_num).wilcoxon_p_change_baseline_global_std, peaks_stat(stim_num).wilcoxon_h_change_baseline_global_std]= signrank(peaks_stat(stim_num).change_baseline_global_std);

        %pre-train response M - absolute values + relative change
        peaks_stat(stim_num).pre_response_M=[peak_pre_response_M{2}(:,stim_num),  peak_pre_response_M{3}(:,stim_num)];
        peaks_stat(stim_num).pre_response_M_m=nanmean(peaks_stat(stim_num).pre_response_M,1);
        peaks_stat(stim_num).pre_response_M_std=nanstd(peaks_stat(stim_num).pre_response_M,0,1);

%pre-train response STD - absolute values + relative change
        peaks_stat(stim_num).pre_response_STD=[peak_pre_response_STD{2}(:,stim_num),  peak_pre_response_STD{3}(:,stim_num)];
        peaks_stat(stim_num).pre_response_STD_m=nanmean(peaks_stat(stim_num).pre_response_STD,1);
        peaks_stat(stim_num).pre_response_STD_std=nanstd(peaks_stat(stim_num).pre_response_STD,0,1);
        peaks_stat(stim_num).change_pre_response_STD=[(peak_pre_response_STD{3}(:,stim_num)-peak_pre_response_STD{2}(:,stim_num))./abs(peak_pre_response_STD{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_pre_response_STD_m=nanmean(peaks_stat(stim_num).change_pre_response_STD,1);
        peaks_stat(stim_num).change_pre_response_STD_std=nanstd(peaks_stat(stim_num).change_pre_response_STD,0,1);        
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_pre_response_STD, peaks_stat(stim_num).lillietest_p_pre_response_STD] = lillietest(peaks_stat(stim_num).pre_response_STD(:,2)- peaks_stat(stim_num).pre_response_STD(:,1));
        [peaks_stat(stim_num).lillietest_h_change_pre_response_STD, peaks_stat(stim_num).lillietest_p_change_pre_response_STD] = lillietest( peaks_stat(stim_num).change_pre_response_STD);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_pre_response_STD, peaks_stat(stim_num).ttest_p_pre_response_STD]= ttest(peaks_stat(stim_num).pre_response_STD(:,1),peaks_stat(stim_num).pre_response_STD(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_pre_response_STD, peaks_stat(stim_num).wilcoxon_h_pre_response_STD]= signrank(peaks_stat(stim_num).pre_response_STD(:,1),peaks_stat(stim_num).pre_response_STD(:,2));
        [peaks_stat(stim_num).ttest_h_change_pre_response_STD, peaks_stat(stim_num).ttest_p_change_pre_response_STD]= ttest(peaks_stat(stim_num).change_pre_response_STD);
        [peaks_stat(stim_num).wilcoxon_p_change_pre_response_STD, peaks_stat(stim_num).wilcoxon_h_change_pre_response_STD]= signrank(peaks_stat(stim_num).change_pre_response_STD);

                %pre-train response variance across trials - absolute values + relative change
        peaks_stat(stim_num).pre_response_VAR=[peak_pre_response_VAR{2}(:,stim_num),  peak_pre_response_VAR{3}(:,stim_num)];
        peaks_stat(stim_num).pre_response_VAR_m=nanmean(peaks_stat(stim_num).pre_response_VAR,1);
        peaks_stat(stim_num).pre_response_VAR_std=nanstd(peaks_stat(stim_num).pre_response_VAR,0,1);
        peaks_stat(stim_num).change_pre_response_VAR=[(peak_pre_response_VAR{3}(:,stim_num)-peak_pre_response_VAR{2}(:,stim_num))./abs(peak_pre_response_VAR{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_pre_response_VAR_m=nanmean(peaks_stat(stim_num).change_pre_response_VAR,1);
        peaks_stat(stim_num).change_pre_response_VAR_std=nanstd(peaks_stat(stim_num).change_pre_response_VAR,0,1);        
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_pre_response_VAR, peaks_stat(stim_num).lillietest_p_pre_response_VAR] = lillietest(peaks_stat(stim_num).pre_response_VAR(:,2)- peaks_stat(stim_num).pre_response_VAR(:,1));
        [peaks_stat(stim_num).lillietest_h_change_pre_response_VAR, peaks_stat(stim_num).lillietest_p_change_pre_response_VAR] = lillietest( peaks_stat(stim_num).change_pre_response_VAR);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_pre_response_VAR, peaks_stat(stim_num).ttest_p_pre_response_VAR]= ttest(peaks_stat(stim_num).pre_response_VAR(:,1),peaks_stat(stim_num).pre_response_VAR(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_pre_response_VAR, peaks_stat(stim_num).wilcoxon_h_pre_response_VAR]= signrank(peaks_stat(stim_num).pre_response_VAR(:,1),peaks_stat(stim_num).pre_response_VAR(:,2));
        [peaks_stat(stim_num).ttest_h_change_pre_response_VAR, peaks_stat(stim_num).ttest_p_change_pre_response_VAR]= ttest(peaks_stat(stim_num).change_pre_response_VAR);
        [peaks_stat(stim_num).wilcoxon_p_change_pre_response_VAR, peaks_stat(stim_num).wilcoxon_h_change_pre_response_VAR]= signrank(peaks_stat(stim_num).change_pre_response_VAR);

        %pre-train response variance in time - absolute values + relative change
        peaks_stat(stim_num).pre_response_VAR_time=[peak_pre_response_VAR_time{2}(:,stim_num),  peak_pre_response_VAR_time{3}(:,stim_num)];
        peaks_stat(stim_num).pre_response_VAR_time_m=nanmean(peaks_stat(stim_num).pre_response_VAR_time,1);
        peaks_stat(stim_num).pre_response_VAR_time_std=nanstd(peaks_stat(stim_num).pre_response_VAR_time,0,1);
        peaks_stat(stim_num).change_pre_response_VAR_time=[(peak_pre_response_VAR_time{3}(:,stim_num)-peak_pre_response_VAR_time{2}(:,stim_num))./abs(peak_pre_response_VAR_time{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_pre_response_VAR_time_m=nanmean(peaks_stat(stim_num).change_pre_response_VAR_time,1);
        peaks_stat(stim_num).change_pre_response_VAR_time_std=nanstd(peaks_stat(stim_num).change_pre_response_VAR_time,0,1);        
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_pre_response_VAR_time, peaks_stat(stim_num).lillietest_p_pre_response_VAR_time] = lillietest(peaks_stat(stim_num).pre_response_VAR_time(:,2)- peaks_stat(stim_num).pre_response_VAR_time(:,1));
        [peaks_stat(stim_num).lillietest_h_change_pre_response_VAR_time, peaks_stat(stim_num).lillietest_p_change_pre_response_VAR_time] = lillietest( peaks_stat(stim_num).change_pre_response_VAR_time);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_pre_response_VAR_time, peaks_stat(stim_num).ttest_p_pre_response_VAR_time]= ttest(peaks_stat(stim_num).pre_response_VAR_time(:,1),peaks_stat(stim_num).pre_response_VAR_time(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_pre_response_VAR_time, peaks_stat(stim_num).wilcoxon_h_pre_response_VAR_time]= signrank(peaks_stat(stim_num).pre_response_VAR_time(:,1),peaks_stat(stim_num).pre_response_VAR_time(:,2));
        [peaks_stat(stim_num).ttest_h_change_pre_response_VAR_time, peaks_stat(stim_num).ttest_p_change_pre_response_VAR_time]= ttest(peaks_stat(stim_num).change_pre_response_VAR_time);
        [peaks_stat(stim_num).wilcoxon_p_change_pre_response_VAR_time, peaks_stat(stim_num).wilcoxon_h_change_pre_response_VAR_time]= signrank(peaks_stat(stim_num).change_pre_response_VAR_time);

        %pre-train response CV - absolute values + relative change
        peaks_stat(stim_num).pre_response_CV=[peak_pre_response_CV{2}(:,stim_num),  peak_pre_response_CV{3}(:,stim_num)];
        peaks_stat(stim_num).pre_response_CV_m=nanmean(peaks_stat(stim_num).pre_response_CV,1);
        peaks_stat(stim_num).pre_response_CV_std=nanstd(peaks_stat(stim_num).pre_response_CV,0,1);
        peaks_stat(stim_num).change_pre_response_CV=[(peak_pre_response_CV{3}(:,stim_num)-peak_pre_response_CV{2}(:,stim_num))./abs(peak_pre_response_CV{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_pre_response_CV_m=nanmean(peaks_stat(stim_num).change_pre_response_CV,1);
        peaks_stat(stim_num).change_pre_response_CV_std=nanstd(peaks_stat(stim_num).change_pre_response_CV,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_pre_response_CV, peaks_stat(stim_num).lillietest_p_pre_response_CV] = lillietest(peaks_stat(stim_num).pre_response_CV(:,2)- peaks_stat(stim_num).pre_response_CV(:,1));
        [peaks_stat(stim_num).lillietest_h_change_pre_response_CV, peaks_stat(stim_num).lillietest_p_change_pre_response_CV] = lillietest( peaks_stat(stim_num).change_pre_response_CV);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_pre_response_CV, peaks_stat(stim_num).ttest_p_pre_response_CV]= ttest(peaks_stat(stim_num).pre_response_CV(:,1),peaks_stat(stim_num).pre_response_CV(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_pre_response_CV, peaks_stat(stim_num).wilcoxon_h_pre_response_CV]= signrank(peaks_stat(stim_num).pre_response_CV(:,1),peaks_stat(stim_num).pre_response_CV(:,2));
        [peaks_stat(stim_num).ttest_h_change_pre_response_CV, peaks_stat(stim_num).ttest_p_change_pre_response_CV]= ttest(peaks_stat(stim_num).change_pre_response_CV);
        [peaks_stat(stim_num).wilcoxon_p_change_pre_response_CV, peaks_stat(stim_num).wilcoxon_h_change_pre_response_CV]= signrank(peaks_stat(stim_num).change_pre_response_CV);
        
        %response (Vm over 50ms) STD - absolute values + relative change
        peaks_stat(stim_num).Vm_res_STD=[peak_Vm_res_STD{2}(:,stim_num),  peak_Vm_res_STD{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_STD_m=nanmean(peaks_stat(stim_num).Vm_res_STD,1);
        peaks_stat(stim_num).Vm_res_STD_std=nanstd(peaks_stat(stim_num).Vm_res_STD,0,1);
        peaks_stat(stim_num).change_Vm_res_STD=[(peak_Vm_res_STD{3}(:,stim_num)-peak_Vm_res_STD{2}(:,stim_num))./abs(peak_Vm_res_STD{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_STD_m=nanmean(peaks_stat(stim_num).change_Vm_res_STD,1);
        peaks_stat(stim_num).change_Vm_res_STD_std=nanstd(peaks_stat(stim_num).change_Vm_res_STD,0,1);
        change_Vm_res_STD_mat(:,stim_num)= peaks_stat(stim_num).change_Vm_res_STD;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_STD, peaks_stat(stim_num).lillietest_p_Vm_res_STD] = lillietest(peaks_stat(stim_num).Vm_res_STD(:,2)- peaks_stat(stim_num).Vm_res_STD(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_STD, peaks_stat(stim_num).lillietest_p_change_Vm_res_STD] = lillietest( peaks_stat(stim_num).change_Vm_res_STD);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_STD, peaks_stat(stim_num).ttest_p_Vm_res_STD]= ttest(peaks_stat(stim_num).Vm_res_STD(:,1),peaks_stat(stim_num).Vm_res_STD(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_STD, peaks_stat(stim_num).wilcoxon_h_Vm_res_STD]= signrank(peaks_stat(stim_num).Vm_res_STD(:,1),peaks_stat(stim_num).Vm_res_STD(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_STD, peaks_stat(stim_num).ttest_p_change_Vm_res_STD]= ttest(peaks_stat(stim_num).change_Vm_res_STD);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_STD, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_STD]= signrank(peaks_stat(stim_num).change_Vm_res_STD);

         %response (Vm over 50ms) variance over time - absolute values + relative change
        peaks_stat(stim_num).Vm_res_VAR_time=[peak_Vm_res_VAR_time{2}(:,stim_num),  peak_Vm_res_VAR_time{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_VAR_time_m=nanmean(peaks_stat(stim_num).Vm_res_VAR_time,1);
        peaks_stat(stim_num).Vm_res_VAR_time_std=nanstd(peaks_stat(stim_num).Vm_res_VAR_time,0,1);
        peaks_stat(stim_num).change_Vm_res_VAR_time=[(peak_Vm_res_VAR_time{3}(:,stim_num)-peak_Vm_res_VAR_time{2}(:,stim_num))./abs(peak_Vm_res_VAR_time{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_VAR_time_m=nanmean(peaks_stat(stim_num).change_Vm_res_VAR_time,1);
        peaks_stat(stim_num).change_Vm_res_VAR_time_std=nanstd(peaks_stat(stim_num).change_Vm_res_VAR_time,0,1);
        change_Vm_res_VAR_time_mat(:,stim_num)= peaks_stat(stim_num).change_Vm_res_VAR_time;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_VAR_time, peaks_stat(stim_num).lillietest_p_Vm_res_VAR_time] = lillietest(peaks_stat(stim_num).Vm_res_VAR_time(:,2)- peaks_stat(stim_num).Vm_res_VAR_time(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_VAR_time, peaks_stat(stim_num).lillietest_p_change_Vm_res_VAR_time] = lillietest( peaks_stat(stim_num).change_Vm_res_VAR_time);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_VAR_time, peaks_stat(stim_num).ttest_p_Vm_res_VAR_time]= ttest(peaks_stat(stim_num).Vm_res_VAR_time(:,1),peaks_stat(stim_num).Vm_res_VAR_time(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_VAR_time, peaks_stat(stim_num).wilcoxon_h_Vm_res_VAR_time]= signrank(peaks_stat(stim_num).Vm_res_VAR_time(:,1),peaks_stat(stim_num).Vm_res_VAR_time(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_VAR_time, peaks_stat(stim_num).ttest_p_change_Vm_res_VAR_time]= ttest(peaks_stat(stim_num).change_Vm_res_VAR_time);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_VAR_time, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_VAR_time]= signrank(peaks_stat(stim_num).change_Vm_res_VAR_time);

        %response (Vm over 50ms) STD averaged on all stim. - absolute values + relative change
        peaks_stat(stim_num).Vm_res_STD_meanstim=[peak_Vm_res_STD_meanstim{2}(:,stim_num),  peak_Vm_res_STD_meanstim{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_STD_meanstim_m=nanmean(peaks_stat(stim_num).Vm_res_STD_meanstim,1);
        peaks_stat(stim_num).Vm_res_STD_meanstim_std=nanstd(peaks_stat(stim_num).Vm_res_STD_meanstim,0,1);
        peaks_stat(stim_num).change_Vm_res_STD_meanstim=[(peak_Vm_res_STD_meanstim{3}(:,stim_num)-peak_Vm_res_STD_meanstim{2}(:,stim_num))./abs(peak_Vm_res_STD_meanstim{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_STD_meanstim_m=nanmean(peaks_stat(stim_num).change_Vm_res_STD_meanstim,1);
        peaks_stat(stim_num).change_Vm_res_STD_meanstim_std=nanstd(peaks_stat(stim_num).change_Vm_res_STD_meanstim,0,1);
        change_Vm_res_STD_meanstim_mat(:,stim_num)= peaks_stat(stim_num).change_Vm_res_STD_meanstim;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_STD_meanstim, peaks_stat(stim_num).lillietest_p_Vm_res_STD_meanstim] = lillietest(peaks_stat(stim_num).Vm_res_STD_meanstim(:,2)- peaks_stat(stim_num).Vm_res_STD_meanstim(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_STD_meanstim, peaks_stat(stim_num).lillietest_p_change_Vm_res_STD_meanstim] = lillietest( peaks_stat(stim_num).change_Vm_res_STD_meanstim);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_STD_meanstim, peaks_stat(stim_num).ttest_p_Vm_res_STD_meanstim]= ttest(peaks_stat(stim_num).Vm_res_STD_meanstim(:,1),peaks_stat(stim_num).Vm_res_STD_meanstim(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_STD_meanstim, peaks_stat(stim_num).wilcoxon_h_Vm_res_STD_meanstim]= signrank(peaks_stat(stim_num).Vm_res_STD_meanstim(:,1),peaks_stat(stim_num).Vm_res_STD_meanstim(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_STD_meanstim, peaks_stat(stim_num).ttest_p_change_Vm_res_STD_meanstim]= ttest(peaks_stat(stim_num).change_Vm_res_STD_meanstim);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_STD_meanstim, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_STD_meanstim]= signrank(peaks_stat(stim_num).change_Vm_res_STD_meanstim);
 
        %response (Vm over 50ms) VAR averaged on all stim. - absolute values + relative change
        peaks_stat(stim_num).Vm_res_VAR_meanstim=[peak_Vm_res_VAR_meanstim{2}(:,stim_num),  peak_Vm_res_VAR_meanstim{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_VAR_meanstim_m=nanmean(peaks_stat(stim_num).Vm_res_VAR_meanstim,1);
        peaks_stat(stim_num).Vm_res_VAR_meanstim_std=nanstd(peaks_stat(stim_num).Vm_res_VAR_meanstim,0,1);
        peaks_stat(stim_num).change_Vm_res_VAR_meanstim=[(peak_Vm_res_VAR_meanstim{3}(:,stim_num)-peak_Vm_res_VAR_meanstim{2}(:,stim_num))./abs(peak_Vm_res_VAR_meanstim{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_VAR_meanstim_m=nanmean(peaks_stat(stim_num).change_Vm_res_VAR_meanstim,1);
        peaks_stat(stim_num).change_Vm_res_VAR_meanstim_std=nanstd(peaks_stat(stim_num).change_Vm_res_VAR_meanstim,0,1);
        change_Vm_res_VAR_meanstim_mat(:,stim_num)= peaks_stat(stim_num).change_Vm_res_VAR_meanstim;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_VAR_meanstim, peaks_stat(stim_num).lillietest_p_Vm_res_VAR_meanstim] = lillietest(peaks_stat(stim_num).Vm_res_VAR_meanstim(:,2)- peaks_stat(stim_num).Vm_res_VAR_meanstim(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_VAR_meanstim, peaks_stat(stim_num).lillietest_p_change_Vm_res_VAR_meanstim] = lillietest( peaks_stat(stim_num).change_Vm_res_VAR_meanstim);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_VAR_meanstim, peaks_stat(stim_num).ttest_p_Vm_res_VAR_meanstim]= ttest(peaks_stat(stim_num).Vm_res_VAR_meanstim(:,1),peaks_stat(stim_num).Vm_res_VAR_meanstim(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_VAR_meanstim, peaks_stat(stim_num).wilcoxon_h_Vm_res_VAR_meanstim]= signrank(peaks_stat(stim_num).Vm_res_VAR_meanstim(:,1),peaks_stat(stim_num).Vm_res_VAR_meanstim(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_VAR_meanstim, peaks_stat(stim_num).ttest_p_change_Vm_res_VAR_meanstim]= ttest(peaks_stat(stim_num).change_Vm_res_VAR_meanstim);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_VAR_meanstim, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_VAR_meanstim]= signrank(peaks_stat(stim_num).change_Vm_res_VAR_meanstim);
       
   %response CV - absolute values + relative change
        peaks_stat(stim_num).Vm_res_CV=[peak_Vm_res_CV{2}(:,stim_num),  peak_Vm_res_CV{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_CV_m=nanmean(peaks_stat(stim_num).Vm_res_CV,1);
        peaks_stat(stim_num).Vm_res_CV_std=nanstd(peaks_stat(stim_num).Vm_res_CV,0,1);
        peaks_stat(stim_num).change_Vm_res_CV=[(peak_Vm_res_CV{3}(:,stim_num)-peak_Vm_res_CV{2}(:,stim_num))./abs(peak_Vm_res_CV{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_CV_m=nanmean(peaks_stat(stim_num).change_Vm_res_CV,1);
        peaks_stat(stim_num).change_Vm_res_CV_std=nanstd(peaks_stat(stim_num).change_Vm_res_CV,0,1);
        change_Vm_res_CV_mat(:,stim_num)= peaks_stat(stim_num).change_Vm_res_CV;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_CV, peaks_stat(stim_num).lillietest_p_Vm_res_CV] = lillietest(peaks_stat(stim_num).Vm_res_CV(:,2)- peaks_stat(stim_num).Vm_res_CV(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_CV, peaks_stat(stim_num).lillietest_p_change_Vm_res_CV] = lillietest( peaks_stat(stim_num).change_Vm_res_CV);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_CV, peaks_stat(stim_num).ttest_p_Vm_res_CV]= ttest(peaks_stat(stim_num).Vm_res_CV(:,1),peaks_stat(stim_num).Vm_res_CV(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_CV, peaks_stat(stim_num).wilcoxon_h_Vm_res_CV]= signrank(peaks_stat(stim_num).Vm_res_CV(:,1),peaks_stat(stim_num).Vm_res_CV(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_CV, peaks_stat(stim_num).ttest_p_change_Vm_res_CV]= ttest(peaks_stat(stim_num).change_Vm_res_CV);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_CV, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_CV]= signrank(peaks_stat(stim_num).change_Vm_res_CV);
        
  %response M - absolute values + relative change
        peaks_stat(stim_num).Vm_res_M=[peak_Vm_res_M{2}(:,stim_num),  peak_Vm_res_M{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_M_m=nanmean(peaks_stat(stim_num).Vm_res_M,1);
        peaks_stat(stim_num).Vm_res_M_std=nanstd(peaks_stat(stim_num).Vm_res_M,0,1);
        peaks_stat(stim_num).change_Vm_res_M=[(peak_Vm_res_M{3}(:,stim_num)-peak_Vm_res_M{2}(:,stim_num))./abs(peak_Vm_res_M{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_M_m=nanmean(peaks_stat(stim_num).change_Vm_res_M,1);
        peaks_stat(stim_num).change_Vm_res_M_std=nanstd(peaks_stat(stim_num).change_Vm_res_M,0,1);
        change_Vm_res_M_mat(:,stim_num)= peaks_stat(stim_num).change_Vm_res_M;
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_M, peaks_stat(stim_num).lillietest_p_Vm_res_M] = lillietest(peaks_stat(stim_num).Vm_res_M(:,2)- peaks_stat(stim_num).Vm_res_M(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_M, peaks_stat(stim_num).lillietest_p_change_Vm_res_M] = lillietest( peaks_stat(stim_num).change_Vm_res_M);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_M, peaks_stat(stim_num).ttest_p_Vm_res_M]= ttest(peaks_stat(stim_num).Vm_res_M(:,1),peaks_stat(stim_num).Vm_res_M(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_M, peaks_stat(stim_num).wilcoxon_h_Vm_res_M]= signrank(peaks_stat(stim_num).Vm_res_M(:,1),peaks_stat(stim_num).Vm_res_M(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_M, peaks_stat(stim_num).ttest_p_change_Vm_res_M]= ttest(peaks_stat(stim_num).change_Vm_res_M);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_M, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_M]= signrank(peaks_stat(stim_num).change_Vm_res_M);
         
          %response (Vm over 50ms) M averaged on all stim. - absolute values + relative change
        peaks_stat(stim_num).Vm_res_M_meanstim=[peak_Vm_res_M_meanstim{2}(:,stim_num),  peak_Vm_res_M_meanstim{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_M_meanstim_m=nanmean(peaks_stat(stim_num).Vm_res_M_meanstim,1);
        peaks_stat(stim_num).Vm_res_M_meanstim_std=nanstd(peaks_stat(stim_num).Vm_res_M_meanstim,0,1);
        
      %train response STD - absolute values + relative change
        peaks_stat(stim_num).Vm_res_train_STD=[peak_Vm_res_train_STD{2}(:,stim_num),  peak_Vm_res_train_STD{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_train_STD_m=nanmean(peaks_stat(stim_num).Vm_res_train_STD,1);
        peaks_stat(stim_num).Vm_res_train_STD_std=nanstd(peaks_stat(stim_num).Vm_res_train_STD,0,1);
        peaks_stat(stim_num).change_Vm_res_train_STD=[(peak_Vm_res_train_STD{3}(:,stim_num)-peak_Vm_res_train_STD{2}(:,stim_num))./abs(peak_Vm_res_train_STD{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_train_STD_m=nanmean(peaks_stat(stim_num).change_Vm_res_train_STD,1);
        peaks_stat(stim_num).change_Vm_res_train_STD_std=nanstd(peaks_stat(stim_num).change_Vm_res_train_STD,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_train_STD, peaks_stat(stim_num).lillietest_p_Vm_res_train_STD] = lillietest(peaks_stat(stim_num).Vm_res_train_STD(:,2)- peaks_stat(stim_num).Vm_res_train_STD(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_train_STD, peaks_stat(stim_num).lillietest_p_change_Vm_res_train_STD] = lillietest( peaks_stat(stim_num).change_Vm_res_train_STD);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_train_STD, peaks_stat(stim_num).ttest_p_Vm_res_train_STD]= ttest(peaks_stat(stim_num).Vm_res_train_STD(:,1),peaks_stat(stim_num).Vm_res_train_STD(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_train_STD, peaks_stat(stim_num).wilcoxon_h_Vm_res_train_STD]= signrank(peaks_stat(stim_num).Vm_res_train_STD(:,1),peaks_stat(stim_num).Vm_res_train_STD(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_train_STD, peaks_stat(stim_num).ttest_p_change_Vm_res_train_STD]= ttest(peaks_stat(stim_num).change_Vm_res_train_STD);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_train_STD, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_train_STD]= signrank(peaks_stat(stim_num).change_Vm_res_train_STD);
 
        %train response variance (trial to trial) - absolute values + relative change
        peaks_stat(stim_num).Vm_res_train_VAR=[peak_Vm_res_train_VAR{2}(:,stim_num),  peak_Vm_res_train_VAR{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_train_VAR_m=nanmean(peaks_stat(stim_num).Vm_res_train_VAR,1);
        peaks_stat(stim_num).Vm_res_train_VAR_std=nanstd(peaks_stat(stim_num).Vm_res_train_VAR,0,1);
        peaks_stat(stim_num).change_Vm_res_train_VAR=[(peak_Vm_res_train_VAR{3}(:,stim_num)-peak_Vm_res_train_VAR{2}(:,stim_num))./abs(peak_Vm_res_train_VAR{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_train_VAR_m=nanmean(peaks_stat(stim_num).change_Vm_res_train_VAR,1);
        peaks_stat(stim_num).change_Vm_res_train_VAR_std=nanstd(peaks_stat(stim_num).change_Vm_res_train_VAR,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_train_VAR, peaks_stat(stim_num).lillietest_p_Vm_res_train_VAR] = lillietest(peaks_stat(stim_num).Vm_res_train_VAR(:,2)- peaks_stat(stim_num).Vm_res_train_VAR(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_train_VAR, peaks_stat(stim_num).lillietest_p_change_Vm_res_train_VAR] = lillietest( peaks_stat(stim_num).change_Vm_res_train_VAR);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_train_VAR, peaks_stat(stim_num).ttest_p_Vm_res_train_VAR]= ttest(peaks_stat(stim_num).Vm_res_train_VAR(:,1),peaks_stat(stim_num).Vm_res_train_VAR(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_train_VAR, peaks_stat(stim_num).wilcoxon_h_Vm_res_train_VAR]= signrank(peaks_stat(stim_num).Vm_res_train_VAR(:,1),peaks_stat(stim_num).Vm_res_train_VAR(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_train_VAR, peaks_stat(stim_num).ttest_p_change_Vm_res_train_VAR]= ttest(peaks_stat(stim_num).change_Vm_res_train_VAR);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_train_VAR, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_train_VAR]= signrank(peaks_stat(stim_num).change_Vm_res_train_VAR);

         %train response variance along time - absolute values + relative change
        peaks_stat(stim_num).Vm_res_train_VAR_time=[peak_Vm_res_train_VAR_time{2}(:,stim_num),  peak_Vm_res_train_VAR_time{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_train_VAR_time_m=nanmean(peaks_stat(stim_num).Vm_res_train_VAR_time,1);
        peaks_stat(stim_num).Vm_res_train_VAR_time_std=nanstd(peaks_stat(stim_num).Vm_res_train_VAR_time,0,1);
        peaks_stat(stim_num).change_Vm_res_train_VAR_time=[(peak_Vm_res_train_VAR_time{3}(:,stim_num)-peak_Vm_res_train_VAR_time{2}(:,stim_num))./abs(peak_Vm_res_train_VAR_time{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_train_VAR_time_m=nanmean(peaks_stat(stim_num).change_Vm_res_train_VAR_time,1);
        peaks_stat(stim_num).change_Vm_res_train_VAR_time_std=nanstd(peaks_stat(stim_num).change_Vm_res_train_VAR_time,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_train_VAR_time, peaks_stat(stim_num).lillietest_p_Vm_res_train_VAR_time] = lillietest(peaks_stat(stim_num).Vm_res_train_VAR_time(:,2)- peaks_stat(stim_num).Vm_res_train_VAR_time(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_train_VAR_time, peaks_stat(stim_num).lillietest_p_change_Vm_res_train_VAR_time] = lillietest( peaks_stat(stim_num).change_Vm_res_train_VAR_time);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_train_VAR_time, peaks_stat(stim_num).ttest_p_Vm_res_train_VAR_time]= ttest(peaks_stat(stim_num).Vm_res_train_VAR_time(:,1),peaks_stat(stim_num).Vm_res_train_VAR_time(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_train_VAR_time, peaks_stat(stim_num).wilcoxon_h_Vm_res_train_VAR_time]= signrank(peaks_stat(stim_num).Vm_res_train_VAR_time(:,1),peaks_stat(stim_num).Vm_res_train_VAR_time(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_train_VAR_time, peaks_stat(stim_num).ttest_p_change_Vm_res_train_VAR_time]= ttest(peaks_stat(stim_num).change_Vm_res_train_VAR_time);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_train_VAR_time, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_train_VAR_time]= signrank(peaks_stat(stim_num).change_Vm_res_train_VAR_time);
      
  %train response M - absolute values + relative change
        peaks_stat(stim_num).Vm_res_train_M=[peak_Vm_res_train_M{2}(:,stim_num),  peak_Vm_res_train_M{3}(:,stim_num)];
        peaks_stat(stim_num).Vm_res_train_M_m=nanmean(peaks_stat(stim_num).Vm_res_train_M,1);
        peaks_stat(stim_num).Vm_res_train_M_std=nanstd(peaks_stat(stim_num).Vm_res_train_M,0,1);
        peaks_stat(stim_num).change_Vm_res_train_M=[(peak_Vm_res_train_M{3}(:,stim_num)-peak_Vm_res_train_M{2}(:,stim_num))./abs(peak_Vm_res_train_M{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Vm_res_train_M_m=nanmean(peaks_stat(stim_num).change_Vm_res_train_M,1);
        peaks_stat(stim_num).change_Vm_res_train_M_std=nanstd(peaks_stat(stim_num).change_Vm_res_train_M,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Vm_res_train_M, peaks_stat(stim_num).lillietest_p_Vm_res_train_M] = lillietest(peaks_stat(stim_num).Vm_res_train_M(:,2)- peaks_stat(stim_num).Vm_res_train_M(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Vm_res_train_M, peaks_stat(stim_num).lillietest_p_change_Vm_res_train_M] = lillietest( peaks_stat(stim_num).change_Vm_res_train_M);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Vm_res_train_M, peaks_stat(stim_num).ttest_p_Vm_res_train_M]= ttest(peaks_stat(stim_num).Vm_res_train_M(:,1),peaks_stat(stim_num).Vm_res_train_M(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Vm_res_train_M, peaks_stat(stim_num).wilcoxon_h_Vm_res_train_M]= signrank(peaks_stat(stim_num).Vm_res_train_M(:,1),peaks_stat(stim_num).Vm_res_train_M(:,2));
        [peaks_stat(stim_num).ttest_h_change_Vm_res_train_M, peaks_stat(stim_num).ttest_p_change_Vm_res_train_M]= ttest(peaks_stat(stim_num).change_Vm_res_train_M);
        [peaks_stat(stim_num).wilcoxon_p_change_Vm_res_train_M, peaks_stat(stim_num).wilcoxon_h_change_Vm_res_train_M]= signrank(peaks_stat(stim_num).change_Vm_res_train_M);
        
 %post train response STD - absolute values + relative change
        peaks_stat(stim_num).post_train_STD=[peak_post_train_STD{2}(:,stim_num),  peak_post_train_STD{3}(:,stim_num)];
        peaks_stat(stim_num).post_train_STD_m=nanmean(peaks_stat(stim_num).post_train_STD,1);
        peaks_stat(stim_num).post_train_STD_std=nanstd(peaks_stat(stim_num).post_train_STD,0,1);
        peaks_stat(stim_num).change_post_train_STD=[(peak_post_train_STD{3}(:,stim_num)-peak_post_train_STD{2}(:,stim_num))./abs(peak_post_train_STD{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_post_train_STD_m=nanmean(peaks_stat(stim_num).change_post_train_STD,1);
        peaks_stat(stim_num).change_post_train_STD_std=nanstd(peaks_stat(stim_num).change_post_train_STD,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_post_train_STD, peaks_stat(stim_num).lillietest_p_post_train_STD] = lillietest(peaks_stat(stim_num).post_train_STD(:,2)- peaks_stat(stim_num).post_train_STD(:,1));
        [peaks_stat(stim_num).lillietest_h_change_post_train_STD, peaks_stat(stim_num).lillietest_p_change_post_train_STD] = lillietest( peaks_stat(stim_num).change_post_train_STD);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_post_train_STD, peaks_stat(stim_num).ttest_p_post_train_STD]= ttest(peaks_stat(stim_num).post_train_STD(:,1),peaks_stat(stim_num).post_train_STD(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_post_train_STD, peaks_stat(stim_num).wilcoxon_h_post_train_STD]= signrank(peaks_stat(stim_num).post_train_STD(:,1),peaks_stat(stim_num).post_train_STD(:,2));
        [peaks_stat(stim_num).ttest_h_change_post_train_STD, peaks_stat(stim_num).ttest_p_change_post_train_STD]= ttest(peaks_stat(stim_num).change_post_train_STD);
        [peaks_stat(stim_num).wilcoxon_p_change_post_train_STD, peaks_stat(stim_num).wilcoxon_h_change_post_train_STD]= signrank(peaks_stat(stim_num).change_post_train_STD);

        %post train response variance (trial to trial) - absolute values + relative change
        peaks_stat(stim_num).post_train_VAR=[peak_post_train_VAR{2}(:,stim_num),  peak_post_train_VAR{3}(:,stim_num)];
        peaks_stat(stim_num).post_train_VAR_m=nanmean(peaks_stat(stim_num).post_train_VAR,1);
        peaks_stat(stim_num).post_train_VAR_std=nanstd(peaks_stat(stim_num).post_train_VAR,0,1);
        peaks_stat(stim_num).change_post_train_VAR=[(peak_post_train_VAR{3}(:,stim_num)-peak_post_train_VAR{2}(:,stim_num))./abs(peak_post_train_VAR{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_post_train_VAR_m=nanmean(peaks_stat(stim_num).change_post_train_VAR,1);
        peaks_stat(stim_num).change_post_train_VAR_std=nanstd(peaks_stat(stim_num).change_post_train_VAR,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_post_train_VAR, peaks_stat(stim_num).lillietest_p_post_train_VAR] = lillietest(peaks_stat(stim_num).post_train_VAR(:,2)- peaks_stat(stim_num).post_train_VAR(:,1));
        [peaks_stat(stim_num).lillietest_h_change_post_train_VAR, peaks_stat(stim_num).lillietest_p_change_post_train_VAR] = lillietest( peaks_stat(stim_num).change_post_train_VAR);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_post_train_VAR, peaks_stat(stim_num).ttest_p_post_train_VAR]= ttest(peaks_stat(stim_num).post_train_VAR(:,1),peaks_stat(stim_num).post_train_VAR(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_post_train_VAR, peaks_stat(stim_num).wilcoxon_h_post_train_VAR]= signrank(peaks_stat(stim_num).post_train_VAR(:,1),peaks_stat(stim_num).post_train_VAR(:,2));
        [peaks_stat(stim_num).ttest_h_change_post_train_VAR, peaks_stat(stim_num).ttest_p_change_post_train_VAR]= ttest(peaks_stat(stim_num).change_post_train_VAR);
        [peaks_stat(stim_num).wilcoxon_p_change_post_train_VAR, peaks_stat(stim_num).wilcoxon_h_change_post_train_VAR]= signrank(peaks_stat(stim_num).change_post_train_VAR);
       
        %post train response variance across time - absolute values + relative change
        peaks_stat(stim_num).post_train_VAR_time=[peak_post_train_VAR_time{2}(:,stim_num),  peak_post_train_VAR_time{3}(:,stim_num)];
        peaks_stat(stim_num).post_train_VAR_time_m=nanmean(peaks_stat(stim_num).post_train_VAR_time,1);
        peaks_stat(stim_num).post_train_VAR_time_std=nanstd(peaks_stat(stim_num).post_train_VAR_time,0,1);
        peaks_stat(stim_num).change_post_train_VAR_time=[(peak_post_train_VAR_time{3}(:,stim_num)-peak_post_train_VAR_time{2}(:,stim_num))./abs(peak_post_train_VAR_time{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_post_train_VAR_time_m=nanmean(peaks_stat(stim_num).change_post_train_VAR_time,1);
        peaks_stat(stim_num).change_post_train_VAR_time_std=nanstd(peaks_stat(stim_num).change_post_train_VAR_time,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_post_train_VAR_time, peaks_stat(stim_num).lillietest_p_post_train_VAR_time] = lillietest(peaks_stat(stim_num).post_train_VAR_time(:,2)- peaks_stat(stim_num).post_train_VAR_time(:,1));
        [peaks_stat(stim_num).lillietest_h_change_post_train_VAR_time, peaks_stat(stim_num).lillietest_p_change_post_train_VAR_time] = lillietest( peaks_stat(stim_num).change_post_train_VAR_time);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_post_train_VAR_time, peaks_stat(stim_num).ttest_p_post_train_VAR_time]= ttest(peaks_stat(stim_num).post_train_VAR_time(:,1),peaks_stat(stim_num).post_train_VAR_time(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_post_train_VAR_time, peaks_stat(stim_num).wilcoxon_h_post_train_VAR_time]= signrank(peaks_stat(stim_num).post_train_VAR_time(:,1),peaks_stat(stim_num).post_train_VAR_time(:,2));
        [peaks_stat(stim_num).ttest_h_change_post_train_VAR_time, peaks_stat(stim_num).ttest_p_change_post_train_VAR_time]= ttest(peaks_stat(stim_num).change_post_train_VAR_time);
        [peaks_stat(stim_num).wilcoxon_p_change_post_train_VAR_time, peaks_stat(stim_num).wilcoxon_h_change_post_train_VAR_time]= signrank(peaks_stat(stim_num).change_post_train_VAR_time);
       
   %post train response CV - absolute values + relative change
        peaks_stat(stim_num).post_train_CV=[peak_post_train_CV{2}(:,stim_num),  peak_post_train_CV{3}(:,stim_num)];
        peaks_stat(stim_num).post_train_CV_m=nanmean(peaks_stat(stim_num).post_train_CV,1);
        peaks_stat(stim_num).post_train_CV_std=nanstd(peaks_stat(stim_num).post_train_CV,0,1);
        peaks_stat(stim_num).change_post_train_CV=[(peak_post_train_CV{3}(:,stim_num)-peak_post_train_CV{2}(:,stim_num))./abs(peak_post_train_CV{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_post_train_CV_m=nanmean(peaks_stat(stim_num).change_post_train_CV,1);
        peaks_stat(stim_num).change_post_train_CV_std=nanstd(peaks_stat(stim_num).change_post_train_CV,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_post_train_CV, peaks_stat(stim_num).lillietest_p_post_train_CV] = lillietest(peaks_stat(stim_num).post_train_CV(:,2)- peaks_stat(stim_num).post_train_CV(:,1));
        [peaks_stat(stim_num).lillietest_h_change_post_train_CV, peaks_stat(stim_num).lillietest_p_change_post_train_CV] = lillietest( peaks_stat(stim_num).change_post_train_CV);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_post_train_CV, peaks_stat(stim_num).ttest_p_post_train_CV]= ttest(peaks_stat(stim_num).post_train_CV(:,1),peaks_stat(stim_num).post_train_CV(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_post_train_CV, peaks_stat(stim_num).wilcoxon_h_post_train_CV]= signrank(peaks_stat(stim_num).post_train_CV(:,1),peaks_stat(stim_num).post_train_CV(:,2));
        [peaks_stat(stim_num).ttest_h_change_post_train_CV, peaks_stat(stim_num).ttest_p_change_post_train_CV]= ttest(peaks_stat(stim_num).change_post_train_CV);
        [peaks_stat(stim_num).wilcoxon_p_change_post_train_CV, peaks_stat(stim_num).wilcoxon_h_change_post_train_CV]= signrank(peaks_stat(stim_num).change_post_train_CV);
        
  %post train response M - absolute values + relative change
        peaks_stat(stim_num).post_train_M=[peak_post_train_M{2}(:,stim_num),  peak_post_train_M{3}(:,stim_num)];
        peaks_stat(stim_num).post_train_M_m=nanmean(peaks_stat(stim_num).post_train_M,1);
        peaks_stat(stim_num).post_train_M_std=nanstd(peaks_stat(stim_num).post_train_M,0,1);
        peaks_stat(stim_num).change_post_train_M=[(peak_post_train_M{3}(:,stim_num)-peak_post_train_M{2}(:,stim_num))./abs(peak_post_train_M{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_post_train_M_m=nanmean(peaks_stat(stim_num).change_post_train_M,1);
        peaks_stat(stim_num).change_post_train_M_std=nanstd(peaks_stat(stim_num).change_post_train_M,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_post_train_M, peaks_stat(stim_num).lillietest_p_post_train_M] = lillietest(peaks_stat(stim_num).post_train_M(:,2)- peaks_stat(stim_num).post_train_M(:,1));
        [peaks_stat(stim_num).lillietest_h_change_post_train_M, peaks_stat(stim_num).lillietest_p_change_post_train_M] = lillietest( peaks_stat(stim_num).change_post_train_M);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_post_train_M, peaks_stat(stim_num).ttest_p_post_train_M]= ttest(peaks_stat(stim_num).post_train_M(:,1),peaks_stat(stim_num).post_train_M(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_post_train_M, peaks_stat(stim_num).wilcoxon_h_post_train_M]= signrank(peaks_stat(stim_num).post_train_M(:,1),peaks_stat(stim_num).post_train_M(:,2));
        [peaks_stat(stim_num).ttest_h_change_post_train_M, peaks_stat(stim_num).ttest_p_change_post_train_M]= ttest(peaks_stat(stim_num).change_post_train_M);
        [peaks_stat(stim_num).wilcoxon_p_change_post_train_M, peaks_stat(stim_num).wilcoxon_h_change_post_train_M]= signrank(peaks_stat(stim_num).change_post_train_M);
        
%SNR (variance) - absolute values + relative change
        peaks_stat(stim_num).SNR_var=[peak_SNR_var{2}(:,stim_num),  peak_SNR_var{3}(:,stim_num)];
        peaks_stat(stim_num).SNR_var_m=nanmean(peaks_stat(stim_num).SNR_var,1);
        peaks_stat(stim_num).SNR_var_std=nanstd(peaks_stat(stim_num).SNR_var,0,1);
        peaks_stat(stim_num).change_SNR_var=[(peak_SNR_var{3}(:,stim_num)-peak_SNR_var{2}(:,stim_num))./abs(peak_SNR_var{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_SNR_var_m=nanmean(peaks_stat(stim_num).change_SNR_var,1);
        peaks_stat(stim_num).change_SNR_var_std=nanstd(peaks_stat(stim_num).change_SNR_var,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_SNR_var, peaks_stat(stim_num).lillietest_p_SNR_var] = lillietest(peaks_stat(stim_num).SNR_var(:,2)- peaks_stat(stim_num).SNR_var(:,1));
        [peaks_stat(stim_num).lillietest_h_change_SNR_var, peaks_stat(stim_num).lillietest_p_change_SNR_var] = lillietest( peaks_stat(stim_num).change_SNR_var);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_SNR_var, peaks_stat(stim_num).ttest_p_SNR_var]= ttest(peaks_stat(stim_num).SNR_var(:,1),peaks_stat(stim_num).SNR_var(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_SNR_var, peaks_stat(stim_num).wilcoxon_h_SNR_var]= signrank(peaks_stat(stim_num).SNR_var(:,1),peaks_stat(stim_num).SNR_var(:,2));
        [peaks_stat(stim_num).ttest_h_change_SNR_var, peaks_stat(stim_num).ttest_p_change_SNR_var]= ttest(peaks_stat(stim_num).change_SNR_var);
        [peaks_stat(stim_num).wilcoxon_p_change_SNR_var, peaks_stat(stim_num).wilcoxon_h_change_SNR_var]= signrank(peaks_stat(stim_num).change_SNR_var);

        %SNR1 (noise is during sensory train) - absolute values + relative change
        peaks_stat(stim_num).SNR1=[peak_SNR1{2}(:,stim_num),  peak_SNR1{3}(:,stim_num)];
        peaks_stat(stim_num).SNR1_m=nanmean(peaks_stat(stim_num).SNR1,1);
        peaks_stat(stim_num).SNR1_std=nanstd(peaks_stat(stim_num).SNR1,0,1);
        peaks_stat(stim_num).change_SNR1=[(peak_SNR1{3}(:,stim_num)-peak_SNR1{2}(:,stim_num))./abs(peak_SNR1{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_SNR1_m=nanmean(peaks_stat(stim_num).change_SNR1,1);
        peaks_stat(stim_num).change_SNR1_std=nanstd(peaks_stat(stim_num).change_SNR1,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_SNR1, peaks_stat(stim_num).lillietest_p_SNR1] = lillietest(peaks_stat(stim_num).SNR1(:,2)- peaks_stat(stim_num).SNR1(:,1));
        [peaks_stat(stim_num).lillietest_h_change_SNR1, peaks_stat(stim_num).lillietest_p_change_SNR1] = lillietest( peaks_stat(stim_num).change_SNR1);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_SNR1, peaks_stat(stim_num).ttest_p_SNR1]= ttest(peaks_stat(stim_num).SNR1(:,1),peaks_stat(stim_num).SNR1(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_SNR1, peaks_stat(stim_num).wilcoxon_h_SNR1]= signrank(peaks_stat(stim_num).SNR1(:,1),peaks_stat(stim_num).SNR1(:,2));
        [peaks_stat(stim_num).ttest_h_change_SNR1, peaks_stat(stim_num).ttest_p_change_SNR1]= ttest(peaks_stat(stim_num).change_SNR1);
        [peaks_stat(stim_num).wilcoxon_p_change_SNR1, peaks_stat(stim_num).wilcoxon_h_change_SNR1]= signrank(peaks_stat(stim_num).change_SNR1);

        %SNR2 (noise is the ongoing Vm before the sensory train) - absolute values + relative change
        peaks_stat(stim_num).SNR2=[peak_SNR2{2}(:,stim_num),  peak_SNR2{3}(:,stim_num)];
        peaks_stat(stim_num).SNR2_m=nanmean(peaks_stat(stim_num).SNR2,1);
        peaks_stat(stim_num).SNR2_std=nanstd(peaks_stat(stim_num).SNR2,0,1);
        peaks_stat(stim_num).change_SNR2=[(peak_SNR2{3}(:,stim_num)-peak_SNR2{2}(:,stim_num))./abs(peak_SNR2{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_SNR2_m=nanmean(peaks_stat(stim_num).change_SNR2,1);
        peaks_stat(stim_num).change_SNR2_std=nanstd(peaks_stat(stim_num).change_SNR2,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_SNR2, peaks_stat(stim_num).lillietest_p_SNR2] = lillietest(peaks_stat(stim_num).SNR2(:,2)- peaks_stat(stim_num).SNR2(:,1));
        [peaks_stat(stim_num).lillietest_h_change_SNR2, peaks_stat(stim_num).lillietest_p_change_SNR2] = lillietest( peaks_stat(stim_num).change_SNR2);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_SNR2, peaks_stat(stim_num).ttest_p_SNR2]= ttest(peaks_stat(stim_num).SNR2(:,1),peaks_stat(stim_num).SNR2(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_SNR2, peaks_stat(stim_num).wilcoxon_h_SNR2]= signrank(peaks_stat(stim_num).SNR2(:,1),peaks_stat(stim_num).SNR2(:,2));
        [peaks_stat(stim_num).ttest_h_change_SNR2, peaks_stat(stim_num).ttest_p_change_SNR2]= ttest(peaks_stat(stim_num).change_SNR2);
        [peaks_stat(stim_num).wilcoxon_p_change_SNR2, peaks_stat(stim_num).wilcoxon_h_change_SNR2]= signrank(peaks_stat(stim_num).change_SNR2);
  
        %Amplitude signal - absolute values + relative change
        peaks_stat(stim_num).Amplitude_signal=[peak_Amplitude_signal{2}(:,stim_num),  peak_Amplitude_signal{3}(:,stim_num)];
        peaks_stat(stim_num).Amplitude_signal_m=nanmean(peaks_stat(stim_num).Amplitude_signal,1);
        peaks_stat(stim_num).Amplitude_signal_std=nanstd(peaks_stat(stim_num).Amplitude_signal,0,1);
        peaks_stat(stim_num).change_Amplitude_signal=[(peak_Amplitude_signal{3}(:,stim_num)-peak_Amplitude_signal{2}(:,stim_num))./abs(peak_Amplitude_signal{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Amplitude_signal_m=nanmean(peaks_stat(stim_num).change_Amplitude_signal,1);
        peaks_stat(stim_num).change_Amplitude_signal_std=nanstd(peaks_stat(stim_num).change_Amplitude_signal,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Amplitude_signal, peaks_stat(stim_num).lillietest_p_Amplitude_signal] = lillietest(peaks_stat(stim_num).Amplitude_signal(:,2)- peaks_stat(stim_num).Amplitude_signal(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Amplitude_signal, peaks_stat(stim_num).lillietest_p_change_Amplitude_signal] = lillietest( peaks_stat(stim_num).change_Amplitude_signal);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Amplitude_signal, peaks_stat(stim_num).ttest_p_Amplitude_signal]= ttest(peaks_stat(stim_num).Amplitude_signal(:,1),peaks_stat(stim_num).Amplitude_signal(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Amplitude_signal, peaks_stat(stim_num).wilcoxon_h_Amplitude_signal]= signrank(peaks_stat(stim_num).Amplitude_signal(:,1),peaks_stat(stim_num).Amplitude_signal(:,2));
        [peaks_stat(stim_num).ttest_h_change_Amplitude_signal, peaks_stat(stim_num).ttest_p_change_Amplitude_signal]= ttest(peaks_stat(stim_num).change_Amplitude_signal);
        [peaks_stat(stim_num).wilcoxon_p_change_Amplitude_signal, peaks_stat(stim_num).wilcoxon_h_change_Amplitude_signal]= signrank(peaks_stat(stim_num).change_Amplitude_signal);

        %Amplitude noise1 (noise is during sensory train) - absolute values + relative change
        peaks_stat(stim_num).Amplitude_noise1=[peak_Amplitude_noise1{2}(:,stim_num),  peak_Amplitude_noise1{3}(:,stim_num)];
        peaks_stat(stim_num).Amplitude_noise1_m=nanmean(peaks_stat(stim_num).Amplitude_noise1,1);
        peaks_stat(stim_num).Amplitude_noise1_std=nanstd(peaks_stat(stim_num).Amplitude_noise1,0,1);
        peaks_stat(stim_num).change_Amplitude_noise1=[(peak_Amplitude_noise1{3}(:,stim_num)-peak_Amplitude_noise1{2}(:,stim_num))./abs(peak_Amplitude_noise1{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Amplitude_noise1_m=nanmean(peaks_stat(stim_num).change_Amplitude_noise1,1);
        peaks_stat(stim_num).change_Amplitude_noise1_std=nanstd(peaks_stat(stim_num).change_Amplitude_noise1,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Amplitude_noise1, peaks_stat(stim_num).lillietest_p_Amplitude_noise1] = lillietest(peaks_stat(stim_num).Amplitude_noise1(:,2)- peaks_stat(stim_num).Amplitude_noise1(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Amplitude_noise1, peaks_stat(stim_num).lillietest_p_change_Amplitude_noise1] = lillietest( peaks_stat(stim_num).change_Amplitude_noise1);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Amplitude_noise1, peaks_stat(stim_num).ttest_p_Amplitude_noise1]= ttest(peaks_stat(stim_num).Amplitude_noise1(:,1),peaks_stat(stim_num).Amplitude_noise1(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Amplitude_noise1, peaks_stat(stim_num).wilcoxon_h_Amplitude_noise1]= signrank(peaks_stat(stim_num).Amplitude_noise1(:,1),peaks_stat(stim_num).Amplitude_noise1(:,2));
        [peaks_stat(stim_num).ttest_h_change_Amplitude_noise1, peaks_stat(stim_num).ttest_p_change_Amplitude_noise1]= ttest(peaks_stat(stim_num).change_Amplitude_noise1);
        [peaks_stat(stim_num).wilcoxon_p_change_Amplitude_noise1, peaks_stat(stim_num).wilcoxon_h_change_Amplitude_noise1]= signrank(peaks_stat(stim_num).change_Amplitude_noise1);

        %Amplitude noise2 (noise is the ongoing Vm before the sensory train) - absolute values + relative change
        peaks_stat(stim_num).Amplitude_noise2=[peak_Amplitude_noise2{2}(:,stim_num),  peak_Amplitude_noise2{3}(:,stim_num)];
        peaks_stat(stim_num).Amplitude_noise2_m=nanmean(peaks_stat(stim_num).Amplitude_noise2,1);
        peaks_stat(stim_num).Amplitude_noise2_std=nanstd(peaks_stat(stim_num).Amplitude_noise2,0,1);
        peaks_stat(stim_num).change_Amplitude_noise2=[(peak_Amplitude_noise2{3}(:,stim_num)-peak_Amplitude_noise2{2}(:,stim_num))./abs(peak_Amplitude_noise2{2}(:,stim_num))].*100; %percent change
        peaks_stat(stim_num).change_Amplitude_noise2_m=nanmean(peaks_stat(stim_num).change_Amplitude_noise2,1);
        peaks_stat(stim_num).change_Amplitude_noise2_std=nanstd(peaks_stat(stim_num).change_Amplitude_noise2,0,1);
        %testing for normal distribution       
        [peaks_stat(stim_num).lillietest_h_Amplitude_noise2, peaks_stat(stim_num).lillietest_p_Amplitude_noise2] = lillietest(peaks_stat(stim_num).Amplitude_noise2(:,2)- peaks_stat(stim_num).Amplitude_noise2(:,1));
        [peaks_stat(stim_num).lillietest_h_change_Amplitude_noise2, peaks_stat(stim_num).lillietest_p_change_Amplitude_noise2] = lillietest( peaks_stat(stim_num).change_Amplitude_noise2);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_Amplitude_noise2, peaks_stat(stim_num).ttest_p_Amplitude_noise2]= ttest(peaks_stat(stim_num).Amplitude_noise2(:,1),peaks_stat(stim_num).Amplitude_noise2(:,2));   
        [peaks_stat(stim_num).wilcoxon_p_Amplitude_noise2, peaks_stat(stim_num).wilcoxon_h_Amplitude_noise2]= signrank(peaks_stat(stim_num).Amplitude_noise2(:,1),peaks_stat(stim_num).Amplitude_noise2(:,2));
        [peaks_stat(stim_num).ttest_h_change_Amplitude_noise2, peaks_stat(stim_num).ttest_p_change_Amplitude_noise2]= ttest(peaks_stat(stim_num).change_Amplitude_noise2);
        [peaks_stat(stim_num).wilcoxon_p_change_Amplitude_noise2, peaks_stat(stim_num).wilcoxon_h_change_Amplitude_noise2]= signrank(peaks_stat(stim_num).change_Amplitude_noise2);
%%        
  % creating matrices with the value of h to plot asterisks when significant according to the normalized values
  if peaks_stat(stim_num).lillietest_h_change_per10_lat==0
        h_per10_lat(1,stim_num)=peaks_stat(stim_num).ttest_h_change_per10_lat; 
  else  h_per10_lat(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_per10_lat; 
  end
  if peaks_stat(stim_num).lillietest_h_change_baseline_local==0
        h_baseline_local(1,stim_num)=peaks_stat(stim_num).ttest_h_change_baseline_local; 
  else h_baseline_local(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_baseline_local; 
  end
  if peaks_stat(stim_num).lillietest_h_change_Vm_res_STD==0      
        h_Vm_res_STD(1,stim_num)=peaks_stat(stim_num).ttest_h_change_Vm_res_STD; 
  else h_Vm_res_STD(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_Vm_res_STD; 
  end
  if peaks_stat(stim_num).lillietest_h_change_Vm_res_M==0      
        h_Vm_res_M(1,stim_num)=peaks_stat(stim_num).ttest_h_change_Vm_res_M; 
  else h_Vm_res_M(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_Vm_res_M; 
  end
  if peaks_stat(stim_num).lillietest_h_change_val==0      
        h_val(1,stim_num)=peaks_stat(stim_num).ttest_h_change_val; 
  else h_val(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_val; 
  end
  if peaks_stat(stim_num).lillietest_h_change_amp==0      
        h_amp(1,stim_num)=peaks_stat(stim_num).ttest_h_change_amp; 
  else h_amp(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_amp; 
  end
  if peaks_stat(stim_num).lillietest_h_change_amp_abs==0     
        h_amp_abs(1,stim_num)=peaks_stat(stim_num).ttest_h_change_amp_abs; 
  else h_amp_abs(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_amp_abs; 
  end
  if peaks_stat(stim_num).lillietest_h_change_latency==0      
        h_latency(1,stim_num)=peaks_stat(stim_num).ttest_h_change_latency; 
  else  h_latency(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_latency; 
  end
  if peaks_stat(stim_num).lillietest_h_change_half_width==0      
        h_half_width(1,stim_num)=peaks_stat(stim_num).ttest_h_change_half_width; 
  else h_half_width(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_half_width; 
  end
  
    %% creating matrices with the value of h to plot asterisks when significant  
%   if peaks_stat(stim_num).lillietest_h_per10_lat==0
%         h_per10_lat(1,stim_num)=peaks_stat(stim_num).ttest_h_per10_lat; 
%   else  h_per10_lat(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_per10_lat; 
%   end
%   if peaks_stat(stim_num).lillietest_h_baseline_local==0
%         h_baseline_local(1,stim_num)=peaks_stat(stim_num).ttest_h_baseline_local; 
%   else h_baseline_local(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_baseline_local; 
%   end
%   if peaks_stat(stim_num).lillietest_h_Vm_res_STD==0      
%         h_Vm_res_STD(1,stim_num)=peaks_stat(stim_num).ttest_h_Vm_res_STD; 
%   else h_Vm_res_STD(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_Vm_res_STD; 
%   end
%   if peaks_stat(stim_num).lillietest_h_Vm_res_M==0      
%         h_Vm_res_M(1,stim_num)=peaks_stat(stim_num).ttest_h_Vm_res_M; 
%   else h_Vm_res_M(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_Vm_res_M; 
%   end
%   if peaks_stat(stim_num).lillietest_h_val==0      
%         h_val(1,stim_num)=peaks_stat(stim_num).ttest_h_val; 
%   else h_val(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_val; 
%   end
%   if peaks_stat(stim_num).lillietest_h_amp==0      
%         h_amp(1,stim_num)=peaks_stat(stim_num).ttest_h_amp; 
%   else h_amp(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_amp; 
%   end
%   if peaks_stat(stim_num).lillietest_h_amp_abs==0     
%         h_amp_abs(1,stim_num)=peaks_stat(stim_num).ttest_h_amp_abs; 
%   else h_amp_abs(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_amp_abs; 
%   end
%   if peaks_stat(stim_num).lillietest_h_latency==0      
%         h_latency(1,stim_num)=peaks_stat(stim_num).ttest_h_latency; 
%   else  h_latency(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_latency; 
%   end
%   if peaks_stat(stim_num).lillietest_h_half_width==0      
%         h_half_width(1,stim_num)=peaks_stat(stim_num).ttest_h_half_width; 
%   else h_half_width(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_half_width; 
%   end
    end
    %getting the x-axis (stimulus number) and y-axis for the asterisks of significance
     s_per10_lat=find(h_per10_lat==1);
     s_baseline_local=find(h_baseline_local==1);
     s_Vm_res_STD=find(h_Vm_res_STD==1);
     s_Vm_res_M=find(h_Vm_res_M==1);
     s_val=find(h_val==1);
     s_amp=find(h_amp==1);
     s_amp_abs=find(h_amp_abs==1);
     s_latency=find(h_latency==1);
     s_half_width=find(h_half_width==1);
%% Plot response amplitude Vs. local baseline
%     h1=figure;
%     cell_num=1;
%     for cell_num=1:length(files_to_analyze)
%         linfit_ES_Off=polyfit(peaks(cell_num).baseline_local(:,2),peaks(cell_num).amp(:,2),1);
%         linval_ES_Off=polyval(linfit_ES_Off,peaks(cell_num).baseline_local(:,2));
%         linfit_ES_On=polyfit(peaks(cell_num).baseline_local(:,3),peaks(cell_num).amp(:,3),1);
%         linval_ES_On=polyval(linfit_ES_On,peaks(cell_num).baseline_local(:,3));
%        
% 
%         h(cell_num)=figure;
%         color(cell_num,:)=rand(1,3);
%  hold on
% scatter(peaks(cell_num).baseline_local(:,2),peaks(cell_num).amp(:,2),30,'markeredgecolor',color(cell_num,:))
% scatter(peaks(cell_num).baseline_local(:,3),peaks(cell_num).amp(:,3),30,'markerfacecolor',color(cell_num,:),'markeredgecolor',color(cell_num,:))
% plot(peaks(cell_num).baseline_local(:,2),linval_ES_Off,'k-')
% plot(peaks(cell_num).baseline_local(:,3),linval_ES_On,'b-')
% 
% hold off
%         xlabel('Baseline' ,'FontSize', 16);  ylabel('Amplitude','FontSize', 16); 
%         title(['Response amplitude Vs. Baseline, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
% %         set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line]);
%     end   
%     
    %% Plot response amplitude Vs. local baseline - all cells on same plot
%      h=figure;
%     for cell_num=1:length(files_to_analyze)    
%         color(cell_num,:)=rand(1,3);
%  hold on
% scatter(peaks(cell_num).baseline_local(:,2),peaks(cell_num).amp(:,2),30,'markeredgecolor',color(cell_num,:))
% scatter(peaks(cell_num).baseline_local(:,3),peaks(cell_num).amp(:,3),30,'markerfacecolor',color(cell_num,:),'markeredgecolor',color(cell_num,:))
% hold off
%         xlabel('Baseline' ,'FontSize', 16);  ylabel('Amplitude','FontSize', 16); 
%         title(['Response amplitude Vs. Baseline, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
% %         set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line]);
%     end   
      %% Plot response amplitude to first stim Vs. local baseline - all cells on same plot
%     baseline_vec=[]; amp_vec=[];
%     for cell_num=1:length(files_to_analyze)   
%         baseline_vec = [baseline_vec; peaks(cell_num).baseline_local(:,2:3)];
%         amp_vec = [amp_vec; peaks(cell_num).amp(:,2:3)];
%     end
%     if print_flag==1;
%          h=figure;
%  hold on
% scatter(baseline_vec(:,1),amp_vec(:,1),50,'markeredgecolor',color(cell_num,:))
% scatter(baseline_vec(:,2),amp_vec(:,2),50,'markerfacecolor',color(cell_num,:),'markeredgecolor',color(cell_num,:))
% hold off
%         xlabel('Baseline' ,'FontSize', 16);  ylabel('Amplitude','FontSize', 16); 
%         title(['Response amplitude Vs. Baseline, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
% %         set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line]);
%     end     
%% Plot parameters along the train stim - version 1: line+errorbars of percent change
  %creatig matrices with the value of h to plot asterisks when significant  
%   for stim_num=1:11;
%   if peaks_stat(stim_num).lillietest_h_change_amp==0      
%         h_amp(1,stim_num)=peaks_stat(stim_num).ttest_h_change_amp; 
%   else h_amp(1,stim_num)=peaks_stat(stim_num).wilcoxon_h_change_amp; 
%   end
%     end
%     %getting the x-axis (stimulus number) and y-axis for the asterisks of significance
%      s_amp=find(h_amp==1);

% close all
%     if print_flag==1;
% h1=figure;        
%         hold on
%         errorbar([1:size(change_amp_mat,2)], nanmean(change_amp_mat,1),nanstd(change_amp_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_amp=max(ylim_data).*ones(size(s_amp));
%         plot(s_amp,p_amp,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b') %change line to zero
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Mean Response Amplitude [%change]' ,'FontSize', 16);    
%         title(['Mean Response Amplitude - local baseline,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);     
%     end
    %% Plot parameters along the train stim - version 3:bars+error bars of the mean values   
    h2=figure;     
    barwidth1=0.3;
    nstim=1:11;
        hold on
        errbar_h=errorbar([1:size(peak_amp{2}(:,nstim),2)]-0.3,nanmean(peak_amp{2}(:,nstim),1),zeros(1,size(peak_amp{2}(:,nstim),2)), nanstd(peak_amp{2}(:,nstim),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak_amp{3}(:,nstim),2)],nanmean(peak_amp{3}(:,nstim),1),zeros(1,size(peak_amp{3}(:,nstim),2)), nanstd(peak_amp{3}(:,nstim),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak_amp{2}(:,nstim),2)]-0.3, nanmean(peak_amp{2}(:,nstim),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak_amp{3}(:,nstim),2)], nanmean(peak_amp{3}(:,nstim),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
  set(gca,'xlim',[0, length(nstim)+1],'xtick',[nstim], 'fontsize',16);
        ylim_data=[get(gca,'ylim')]';
        p_amp=max(ylim_data).*ones(size(s_amp));
        plot(s_amp-0.15,p_amp,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Amplitude [mV]' ,'FontSize', 16);    
        title(['Mean Response Amplitude - local baseline,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);    
%%        
         h3=figure;  
    barwidth1=0.3;
        hold on
        errbar_h=errorbar([1:size(peak_lat{2}(:,:),2)]-0.3,nanmean(peak_lat{2}(:,:),1),zeros(1,size(peak_lat{2}(:,:),2)), nanstd(peak_lat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak_lat{3}(:,:),2)],nanmean(peak_lat{3}(:,:),1),zeros(1,size(peak_lat{3}(:,:),2)), nanstd(peak_lat{3}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak_lat{2}(:,:),2)]-0.3, nanmean(peak_lat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak_lat{3}(:,:),2)], nanmean(peak_lat{3}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_latency=max(ylim_data).*ones(size(s_latency));
        plot(s_latency-0.15,p_latency,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b') %change line to zero
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Latency [mS]' ,'FontSize', 16);    
        title(['Mean Peak Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);    
        
                 h4=figure;
    barwidth1=0.3;
        hold on
        errbar_h=errorbar([1:size(peak_per10_lat{2}(:,:),2)]-0.3,nanmean(peak_per10_lat{2}(:,:),1),zeros(1,size(peak_per10_lat{2}(:,:),2)), nanstd(peak_per10_lat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak_per10_lat{3}(:,:),2)],nanmean(peak_per10_lat{3}(:,:),1),zeros(1,size(peak_per10_lat{3}(:,:),2)), nanstd(peak_per10_lat{3}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak_per10_lat{2}(:,:),2)]-0.3, nanmean(peak_per10_lat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak_per10_lat{3}(:,:),2)], nanmean(peak_per10_lat{3}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_per10_lat=max(ylim_data).*ones(size(s_per10_lat));
        plot(s_per10_lat-0.15,p_per10_lat,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b') %change line to zero
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Latency [mS]' ,'FontSize', 16);    
        title(['Mean Onset Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);    
        
                       h9=figure;  
    barwidth1=0.3;
        hold on
        errbar_h=errorbar([1:size(peak_Vm_res_STD{2}(:,:),2)]-0.3,nanmean(peak_Vm_res_STD{2}(:,:),1),zeros(1,size(peak_Vm_res_STD{2}(:,:),2)), nanstd(peak_Vm_res_STD{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak_Vm_res_STD{3}(:,:),2)],nanmean(peak_Vm_res_STD{3}(:,:),1),zeros(1,size(peak_Vm_res_STD{3}(:,:),2)), nanstd(peak_Vm_res_STD{3}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak_Vm_res_STD{2}(:,:),2)]-0.3, nanmean(peak_Vm_res_STD{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak_Vm_res_STD{3}(:,:),2)], nanmean(peak_Vm_res_STD{3}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_Vm_res_STD=max(ylim_data).*ones(size(s_Vm_res_STD));
        plot(s_Vm_res_STD-0.15,p_Vm_res_STD,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b') %change line to zero
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean Vm STD [mV]' ,'FontSize', 16);    
        title(['Mean Vm STD, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);    
    %%  Paired plots for global mean trace parameters
stim_num=1;
        % Plot global baseline
baseline_global_Y=peaks_stat(stim_num).baseline_global';
baseline_global_X(1,:)=ones(1,size(baseline_global_Y,2));
baseline_global_X(2,:)=2*ones(1,size(baseline_global_Y,2));
E =peaks_stat(stim_num).baseline_global_std(1,:);
        g1=figure;
hold on
line(baseline_global_X,baseline_global_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(baseline_global_X(:,1), mean(baseline_global_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(peaks_stat(stim_num).baseline_global)
hold off
        x1limits = [0.75 2.25]; x1ticks = [1,2];        
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Trace Baseline [mV]', 'FontSize', 28,'fontname', 'arial');
        title(['Global Baseline, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_change_baseline_global)] ,'FontSize', 20,'fontname', 'arial'); 

                 % Plot adaptation amplitude ratio (amp8-10/amp1)
adapt_amp_Y=peaks_stat(stim_num).adapt_amp';
adapt_amp_X(1,:)=ones(1,size(adapt_amp_Y,2));
adapt_amp_X(2,:)=2*ones(1,size(adapt_amp_Y,2));
E =peaks_stat(stim_num).adapt_amp_std;
        g2=figure;
hold on
line(adapt_amp_X,adapt_amp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(adapt_amp_X(:,1), nanmean(adapt_amp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(peaks_stat(stim_num).adapt_amp)
hold off
        x1limits = [0.75 2.25];  x1ticks = [1,2];       
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Adaptation Amplitude-ratio', 'FontSize', 28,'fontname', 'arial');
        title(['Adaptation Amplitude ratio, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_change_adapt_amp)] ,'FontSize', 20,'fontname', 'arial');   

               % Plot adaptation power ratio (F1)
F1_Y=peaks_stat(stim_num).F1';
F1_X(1,:)=ones(1,size(F1_Y,2));
F1_X(2,:)=2*ones(1,size(F1_Y,2));
E =peaks_stat(stim_num).F1_std;
        g3=figure;
hold on
line(F1_X,F1_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(F1_X(:,1), mean(F1_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(peaks_stat(stim_num).F1)

hold off
        x1limits = [0.75 2.25];   x1ticks = [1,2];   
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('F1', 'FontSize', 28,'fontname', 'arial');
        title(['Adaptation power ratio, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_change_F1)] ,'FontSize', 20,'fontname', 'arial');   
        
                         % Plot pre-response STD
Vm_pre_response_STD_Y=peaks_stat(stim_num).pre_response_STD';
Vm_pre_response_STD_X(1,:)=ones(1,size(Vm_pre_response_STD_Y,2));
Vm_pre_response_STD_X(2,:)=2*ones(1,size(Vm_pre_response_STD_Y,2));
E =peaks_stat(stim_num).pre_response_STD_std;
        g4=figure;
hold on
line(Vm_pre_response_STD_X,Vm_pre_response_STD_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Vm_pre_response_STD_X(:,1), mean(Vm_pre_response_STD_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(peaks_stat(stim_num).pre_response_STD)
hold off
        x1limits = [0.75 2.25];        x1ticks = [1,2];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('std [mV]', 'FontSize', 28,'fontname', 'arial');
        title(['pre-train STD, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_change_pre_response_STD)] ,'FontSize', 20,'fontname', 'arial');   

        % Plot evoked train STD
Vm_res_train_STD_Y=peaks_stat(stim_num).Vm_res_train_STD';
Vm_res_train_STD_X(1,:)=ones(1,size(Vm_res_train_STD_Y,2));
Vm_res_train_STD_X(2,:)=2*ones(1,size(Vm_res_train_STD_Y,2));
E =peaks_stat(stim_num).Vm_res_train_STD_std;
        g5=figure;
hold on
line(Vm_res_train_STD_X,Vm_res_train_STD_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Vm_res_train_STD_X(:,1), mean(Vm_res_train_STD_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(peaks_stat(stim_num).Vm_res_train_STD)
hold off
        x1limits = [0.75 2.25];         x1ticks = [1,2];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('std [mV]', 'FontSize', 28,'fontname', 'arial');
        title(['evoked train STD, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_change_Vm_res_train_STD)] ,'FontSize', 20,'fontname', 'arial');   
        
                % Plot post-train STD
Vm_post_train_STD_Y= peaks_stat(stim_num).post_train_STD';
Vm_post_train_STD_X(1,:)=ones(1,size(Vm_post_train_STD_Y,2));
Vm_post_train_STD_X(2,:)=2*ones(1,size(Vm_post_train_STD_Y,2));
E =peaks_stat(stim_num).post_train_STD_std;
        g6=figure;
hold on
line(Vm_post_train_STD_X,Vm_post_train_STD_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Vm_post_train_STD_X(:,1), mean(Vm_post_train_STD_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(peaks_stat(stim_num).post_train_STD)
hold off
        x1limits = [0.75 2.25];   x1ticks = [1,2];    
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('std [mV]', 'FontSize', 28,'fontname', 'arial');
        title(['post-train STD, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_change_post_train_STD)] ,'FontSize', 20,'fontname', 'arial');   

    %% Scatter plots for first peak parameters
stim_num=1;
% if print_flag==1;
min_line=floor(min(min(peaks_stat(stim_num).amp)));
max_line=ceil(max(max(peaks_stat(stim_num).amp)));
 g1=figure;
 hold on
scatter(peaks_stat(stim_num).amp(:,1),peaks_stat(stim_num).amp(:,2),180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30);  ylabel('NB+','FontSize', 30); 
%         title(['Mean Response Amplitude [mV], stim 1, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_amp)] ,'FontSize', 20);   
        set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line],'FontSize', 30);
         title(['Response Amplitude, stim ',num2str(stim_num)] ,'FontSize', 26);   
 cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'
saveas(g1,['Mean Response Amplitude_stim ',num2str(stim_num),'.fig']) 
print(g1,['Mean Response Amplitude_stim ',num2str(stim_num)],'-dpng','-r600','-opengl')

% %
% min_line=floor(min(min(peaks_stat(stim_num).amp_abs)));
% max_line=ceil(max(max(peaks_stat(stim_num).amp_abs)));
%  g2=figure;
%  hold on
% scatter(peaks_stat(stim_num).amp_abs(:,1),peaks_stat(stim_num).amp_abs(:,2),180,'markerfacecolor','k')
% line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
% hold off
%         xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30); 
%         title(['Mean Response Global Amplitude [mV], stim 1,n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_amp_abs)] ,'FontSize', 20);   
%        set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30); 
% 
% %
% min_line=floor(min(min(peaks_stat(stim_num).half_width.*1000)));
% max_line=ceil(max(max(peaks_stat(stim_num).half_width.*1000)));
%  g3=figure;
%  hold on
% scatter(peaks_stat(stim_num).half_width(:,1).*1000,peaks_stat(stim_num).half_width(:,2).*1000,180,'markerfacecolor','k')
% line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
% hold off
%         xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30); 
%         title(['Mean Half width [mV], stim 1, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_half_width)] ,'FontSize', 20);   
%         set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);
% 
%         %
% min_line=floor(min(min(peaks_stat(stim_num).per10_lat.*1000)));
% max_line=ceil(max(max(peaks_stat(stim_num).per10_lat.*1000)));
%  g4=figure;
%  hold on
% scatter(peaks_stat(stim_num).per10_lat(:,1).*1000,peaks_stat(stim_num).per10_lat(:,2).*1000,180,'markerfacecolor','k')
% line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
% hold off
%         xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30); 
%         title(['Mean Onset Latency [mV], stim 1,n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_per10_lat)] ,'FontSize', 20);   
%        set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30); 
%   
% %
% min_line=floor(min(min(peaks_stat(stim_num).latency.*1000)));
% max_line=ceil(max(max(peaks_stat(stim_num).latency.*1000)));
%  g5=figure;
%  hold on
% scatter(peaks_stat(stim_num).latency(:,1).*1000,peaks_stat(stim_num).latency(:,2).*1000,180,'markerfacecolor','k')
% line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
% hold off
%         xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30); 
%         title(['Mean peak latency [mV], stim 1, n=' num2str(length(files_to_analyze)), ', p=' num2str(peaks_stat(stim_num).ttest_p_latency)] ,'FontSize', 20);   
%         set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);  
%        
%  %      
% min_line=floor(min(min(peaks_stat(stim_num).F1)));
% max_line=ceil(max(max(peaks_stat(stim_num).F1)));
%         g7=figure;
%  hold on
% scatter(peaks_stat(stim_num).F1(:,1),peaks_stat(stim_num).F1(:,2),180,'markerfacecolor','k')
% line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
% hold off
%         xlabel('NB-' ,'FontSize', 30);  ylabel('NB+','FontSize', 30);  
%         title(['Adaptation power ratio, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_F1)] ,'FontSize', 20);   
%          set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);
%    
% %
% min_line=floor(min(min(peaks_stat(stim_num).baseline_global)));
% max_line=ceil(max(max(peaks_stat(stim_num).baseline_global)));
%  g8=figure;
%  hold on
% scatter(peaks_stat(stim_num).baseline_global(:,1),peaks_stat(stim_num).baseline_global(:,2),180,'markerfacecolor','k')
% line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
% hold off
%         xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30);  
%         title(['Global Baseline, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_baseline_global)] ,'FontSize', 20); 
%         set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);
% 
% %
% min_line=floor(min(min(peaks_stat(stim_num).val)));
% max_line=ceil(max(max(peaks_stat(stim_num).val)));
%  g9=figure;
%  hold on
% scatter(peaks_stat(stim_num).val(:,1),peaks_stat(stim_num).val(:,2),180,'markerfacecolor','k')
% line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
% hold off
%         xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30);  
%         title(['Peak Value, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_val)] ,'FontSize', 20); 
%         set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);
% end
 
%% rmANOVA with built-in matlab function - Vm STD 
pre_response_STD_noNB(:,1)= peaks_stat(1).pre_response_STD(:,1);
pre_response_STD_NB(:,1)= peaks_stat(1).pre_response_STD(:,2);
Vm_res_STD_meanstim_noNB(:,1)=peaks_stat(1).Vm_res_STD_meanstim(:,1);
Vm_res_STD_meanstim_NB(:,1)=peaks_stat(1).Vm_res_STD_meanstim(:,2);
post_train_STD_noNB(:,1)=peaks_stat(1).post_train_STD(:,1);
post_train_STD_NB(:,1)=peaks_stat(1).post_train_STD(:,2);

clear pre_response_STD Vm_res_STD post_train_STD

pre_response_STD(:,1)=[pre_response_STD_noNB;pre_response_STD_NB];
Vm_res_STD(:,1)=[Vm_res_STD_meanstim_noNB;Vm_res_STD_meanstim_NB];
post_train_STD(:,1)=[post_train_STD_noNB;post_train_STD_NB];
 S_1(:,1)=1:length(files_to_analyze); %22;
ta_vector_names={'pre_response_STD_noNB','pre_response_STD_NB','Vm_res_STD_meanstim_noNB','Vm_res_STD_meanstim_NB','post_train_STD_noNB','post_train_STD_NB'};
clear VmSTD ranovatbl ta
ta=table(S_1,pre_response_STD_noNB,pre_response_STD_NB,Vm_res_STD_meanstim_noNB,Vm_res_STD_meanstim_NB,post_train_STD_noNB,post_train_STD_NB,...
    'variablenames', {'cells','Y1','Y2','Y3','Y4','Y5','Y6'});
factorNames = {'Sensory_stim','NB'};
within = table({'Sp';'Sp';'Ev';'Ev';'post';'post'},{'N';'Y';'N';'Y';'N';'Y'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% fit the repeated measures model
rm_VmSTD = fitrm(ta,'Y1-Y6~1','WithinDesign',within);
%test for sphericity
mauchly_VmSTD=mauchly(rm_VmSTD);
% run my repeated measures anova here
[ranovatbl_VmSTD] = ranova(rm_VmSTD, 'WithinModel','Sensory_stim*NB');
%multiple comparisons:
NBbySS_VmSTD = multcompare(rm_VmSTD,'NB','By','Sensory_stim');
SSbyNB_VmSTD = multcompare(rm_VmSTD,'Sensory_stim','By','NB');
VmSTD_rmANOVA_NB_effect_p= ranovatbl_VmSTD(5,5);
VmSTD_rmANOVA_sensory_stim_effect_p= ranovatbl_VmSTD(3,5);


eps_VmSTD = epsilon(rm_VmSTD);
eps=table2cell(eps_VmSTD(1,2));
%in case where the sphericity assumption is violated (mauchly test p<0.05):
%if eps<0.75 - use the Greenhouse-Geisser correction.
%if eps>0.75 - use the Huynh-Feldt correction
if eps{1}<0.75
    VmSTD_p_before=table2cell(NBbySS_VmSTD(1,6));
    VmSTD_p_during=table2cell(NBbySS_VmSTD(3,6));
    VmSTD_p_after=table2cell(NBbySS_VmSTD(5,6));
else
     VmSTD_p_before=table2cell(NBbySS_VmSTD(1,7));
     VmSTD_p_during=table2cell(NBbySS_VmSTD(3,7));
     VmSTD_p_after=table2cell(NBbySS_VmSTD(5,7));
end

VmSTD_stat.table=ta;
VmSTD_stat.table_data_vecs=ta_vector_names;
VmSTD_stat.within_design=within;
VmSTD_stat.rm=rm_VmSTD;
VmSTD_stat.mauchly=mauchly_VmSTD;
VmSTD_stat.eps=eps_VmSTD;
VmSTD_stat.ANOVA=ranovatbl_VmSTD;
VmSTD_stat.multcomp_NBbySS=NBbySS_VmSTD;
VmSTD_stat.multcomp_SSbyNB=SSbyNB_VmSTD;
VmSTD_stat.p_before=VmSTD_p_before;
VmSTD_stat.p_during=VmSTD_p_during;
VmSTD_stat.p_after=VmSTD_p_after;
%% rmANOVA with matlab function - Vm M 
pre_response_M_noNB(:,1)= peaks_stat(1).pre_response_M(:,1);
pre_response_M_NB(:,1)= peaks_stat(1).pre_response_M(:,2);
Vm_res_M_meanstim_noNB(:,1)=peaks_stat(1).Vm_res_M_meanstim(:,1);
Vm_res_M_meanstim_NB(:,1)=peaks_stat(1).Vm_res_M_meanstim(:,2);
post_train_M_noNB(:,1)=peaks_stat(1).post_train_M(:,1);
post_train_M_NB(:,1)=peaks_stat(1).post_train_M(:,2);

clear pre_response_M Vm_res_M post_train_M

pre_response_M(:,1)=[pre_response_M_noNB;pre_response_M_NB];
Vm_res_M(:,1)=[Vm_res_M_meanstim_noNB;Vm_res_M_meanstim_NB];
post_train_M(:,1)=[post_train_M_noNB;post_train_M_NB];
 S_1(:,1)=1:length(files_to_analyze); %22;
ta_vector_names={'pre_response_M_noNB','pre_response_M_NB','Vm_res_M_meanstim_noNB','Vm_res_M_meanstim_NB','post_train_M_noNB','post_train_M_NB'};
clear VmM ranovatbl ta
ta=table(S_1,pre_response_M_noNB,pre_response_M_NB,Vm_res_M_meanstim_noNB,Vm_res_M_meanstim_NB,post_train_M_noNB,post_train_M_NB,...
    'variablenames', {'cells','Y1','Y2','Y3','Y4','Y5','Y6'});
factorNames = {'Sensory_stim','NB'};
within = table({'Sp';'Sp';'Ev';'Ev';'post';'post'},{'N';'Y';'N';'Y';'N';'Y'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% fit the repeated measures model
rm_VmM = fitrm(ta,'Y1-Y6~1','WithinDesign',within);
%test for sphericity
mauchly_VmM=mauchly(rm_VmM);
% run my repeated measures anova here
[ranovatbl_VmM] = ranova(rm_VmM, 'WithinModel','Sensory_stim*NB');
%multiple comparisons:
NBbySS_VmM = multcompare(rm_VmM,'NB','By','Sensory_stim');
SSbyNB_VmM = multcompare(rm_VmM,'Sensory_stim','By','NB');
VmM_rmANOVA_NB_effect_p= ranovatbl_VmM(5,5);
VmM_rmANOVA_sensory_stim_effect_p= ranovatbl_VmM(3,5);

eps_VmM = epsilon(rm_VmM);
eps=table2cell(eps_VmM(1,2));
%in case where the sphericity assumption is violated (mauchly test p<0.05):
%if eps<0.75 - use the Greenhouse-Geisser correction.
%if eps>0.75 - use the Huynh-Feldt correction
if eps{1}<0.75
    VmM_p_before=table2cell(NBbySS_VmM(1,6));
    VmM_p_during=table2cell(NBbySS_VmM(3,6));
    VmM_p_after=table2cell(NBbySS_VmM(5,6));
else
     VmM_p_before=table2cell(NBbySS_VmM(1,7));
    VmM_p_during=table2cell(NBbySS_VmM(3,7));
    VmM_p_after=table2cell(NBbySS_VmM(5,7));
end

VmM_stat.table=ta;
VmM_stat.table_data_vecs=ta_vector_names;
VmM_stat.within_design=within;
VmM_stat.rm=rm_VmM;
VmM_stat.mauchly=mauchly_VmM;
VmM_stat.eps=eps_VmM;
VmM_stat.ANOVA=ranovatbl_VmM;
VmM_stat.multcomp_NBbySS=NBbySS_VmM;
VmM_stat.multcomp_SSbyNB=SSbyNB_VmM;
VmM_stat.rmANOVA_NB_effect_p=VmM_rmANOVA_NB_effect_p;
VmM_stat.rmANOVA_sensory_stim_effect_p=VmM_rmANOVA_sensory_stim_effect_p;
VmM_stat.p_before=VmM_p_before;
VmM_stat.p_during=VmM_p_during;
VmM_stat.p_after=VmM_p_after;
%% rmANOVA with matlab function -amplitude. 
peak_amp_noNB(:,:)= peak_amp{2}(:,:);
peak_amp_NB(:,:)= peak_amp{3}(:,:);

clear pre_response_M Vm_res_M post_train_M

peak_amp_all(:,:)=[peak_amp_noNB;peak_amp_NB];
 S_1(:,1)=1:length(files_to_analyze); %22;

clear VmM ranovatbl ta
ta_vector_names={'peak_amp_noNB(:,1:10)','peak_amp_NB(:,1:10)'};
ta=table(S_1,peak_amp_noNB(:,1),peak_amp_noNB(:,2),peak_amp_noNB(:,3),peak_amp_noNB(:,4),peak_amp_noNB(:,5),peak_amp_noNB(:,6),...
    peak_amp_noNB(:,7),peak_amp_noNB(:,8),peak_amp_noNB(:,9),peak_amp_noNB(:,10),...
    peak_amp_NB(:,1),peak_amp_NB(:,2),peak_amp_NB(:,3),peak_amp_NB(:,4),peak_amp_NB(:,5),peak_amp_NB(:,6),...
    peak_amp_NB(:,7),peak_amp_NB(:,8),peak_amp_NB(:,9),peak_amp_NB(:,10),...
    'variablenames', {'cells','Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9','Y10','Y11','Y12','Y13','Y14','Y15','Y16','Y17','Y18','Y19','Y20'}); %,,'Y21','Y22'
factorNames = {'NB','Time'};
within = table({'N';'N';'N';'N';'N';'N';'N';'N';'N';'N';'Y';'Y';'Y';'Y';'Y';'Y';'Y';'Y';'Y';'Y'},{'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'1';'2';'3';'4';'5';'6';'7';'8';'9';'10'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% fit the repeated measures model
rm_amp = fitrm(ta,'Y1-Y20~1','WithinDesign',within);
%test for sphericity
mauchly_amp=mauchly(rm_amp);
% run my repeated measures anova here
[ranovatbl_amp] = ranova(rm_amp,'withinmodel','NB*Time');
%multiple comparisons:
NBbyTime_amp = multcompare(rm_amp,'NB','By','Time');
TimebyNB_amp = multcompare(rm_amp,'Time','By','NB');

eps_amp = epsilon(rm_amp);
eps=table2cell(eps_amp(1,2));
if eps{1}<0.75
    amp_rmANOVA_NB_effect_p= table2cell(ranovatbl_amp(3,6));
    amp_rmANOVA_Time_effect_p= table2cell(ranovatbl_amp(5,6));
else
     amp_rmANOVA_NB_effect_p= table2cell(ranovatbl_amp(3,7));
     amp_rmANOVA_Time_effect_p= table2cell(ranovatbl_amp(5,7));
end

amp_stat.table=ta;
amp_stat.table_data_vecs=ta_vector_names;
amp_stat.within_design=within;
amp_stat.rm=rm_amp;
amp_stat.mauchly=mauchly_amp;
amp_stat.eps=eps_amp;
amp_stat.ANOVA=ranovatbl_amp;
amp_stat.multcomp_NBbyTime=NBbyTime_amp;
amp_stat.multcomp_TimebyNB=TimebyNB_amp;
amp_stat.rmANOVA_NB_effect_p= amp_rmANOVA_NB_effect_p;
amp_stat.rmANOVA_Time_effect_p= amp_rmANOVA_Time_effect_p;
%% bar plot of trial to trial variability NB- and NB+ before, during and after sensory stimulation
    stim_num=1;
   j1=figure;  
    barwidth1=0.3;
        hold on
        errbar_h1=errorbar((1:3)-0.15,nanmean([peaks_stat(stim_num).pre_response_STD(:,1),peaks_stat(stim_num).Vm_res_STD_meanstim(:,1),peaks_stat(stim_num).post_train_STD(:,1)],1),zeros(1,3), nanstd([peaks_stat(stim_num).pre_response_STD(:,1),peaks_stat(stim_num).Vm_res_STD_meanstim(:,1),peaks_stat(stim_num).post_train_STD(:,1)],0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         fn_errorbar_capsize(errbar_h1,capincrease,asymmetric)
        errbar_h2=errorbar((1:3)+0.15,nanmean([peaks_stat(stim_num).pre_response_STD(:,2),peaks_stat(stim_num).Vm_res_STD_meanstim(:,2),peaks_stat(stim_num).post_train_STD(:,2)],1),zeros(1,3), nanstd([peaks_stat(stim_num).pre_response_STD(:,2),peaks_stat(stim_num).Vm_res_STD_meanstim(:,2),peaks_stat(stim_num).post_train_STD(:,2)],0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         fn_errorbar_capsize(errbar_h2,capincrease,asymmetric)
        bar1=bar((1:3)-0.15, nanmean([peaks_stat(stim_num).pre_response_STD(:,1),peaks_stat(stim_num).Vm_res_STD_meanstim(:,1),peaks_stat(stim_num).post_train_STD(:,1)],1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:));
        bar2=bar((1:3)+0.15, nanmean([peaks_stat(stim_num).pre_response_STD(:,2),peaks_stat(stim_num).Vm_res_STD_meanstim(:,2),peaks_stat(stim_num).post_train_STD(:,2)],1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:));
%        line(pre_response_STD_X,pre_response_STD_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
%setting the errorbars not to appear in the legend:
hAnnotation = get(errbar_h1,'Annotation');  hLegendEntry = get(hAnnotation,'LegendInformation'); set(hLegendEntry,'IconDisplayStyle','off')
hAnnotation = get(errbar_h2,'Annotation');  hLegendEntry = get(hAnnotation,'LegendInformation'); set(hLegendEntry,'IconDisplayStyle','off')
% [legh,objh,outh,outm] = legend('NB Off','NB On','Location','northeast');
l = legend (legend_string,'fontsize',9,'Location','northeast'); legend('boxoff')
ylim_data=[get(gca,'ylim')]';
  my=max([peaks_stat(1).pre_response_STD_std,peaks_stat(1).Vm_res_STD_meanstim_std,peaks_stat(1).post_train_STD_std])*1.1+max([peaks_stat(1).pre_response_STD_m,peaks_stat(1).Vm_res_STD_meanstim_m,peaks_stat(1).post_train_STD_m]); 
linex=[1-0.15; 2-0.15];
liney=[my-0.1; my-0.1];
asterisk_sensory='*';
        if VmSTD_stat.p_before{1,1} >0.05 
    asterisk_before='n.s.';
else if VmSTD_stat.p_before{1,1}<0.05 && VmSTD_stat.p_before{1,1}>0.01
    asterisk_before='*';
    else if VmSTD_stat.p_before{1,1}<0.01 && VmSTD_stat.p_before{1,1}>0.001
            asterisk_before='**';
    else if VmSTD_stat.p_before{1,1}<0.001
             asterisk_before='***';
        end
        end
    end
end

if VmSTD_p_during{1,1} >0.05 
    asterisk_during='n.s.';
else if VmSTD_stat.p_during{1,1}<0.05 && VmSTD_stat.p_during{1,1}>0.01
    asterisk_during='*';
    else if VmSTD_stat.p_during{1,1}<0.01 && VmSTD_stat.p_during{1,1}>0.001
            asterisk_during='**';
    else if VmSTD_stat.p_during{1,1}<0.001
             asterisk_during='***';
        end
        end
    end
end
if VmSTD_stat.p_after{1,1} >0.05 
    asterisk_after='n.s.';
else if VmSTD_stat.p_after{1,1}<0.05 && VmSTD_stat.p_after{1,1}>0.01
    asterisk_after='*';
    else if VmSTD_stat.p_after{1,1}<0.01 && VmSTD_stat.p_after{1,1}>0.001
            asterisk_after='**';
    else if VmSTD_stat.p_after{1,1}<0.001
             asterisk_after='***';
        end
        end
    end
end
text(1,my,asterisk_before,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
% line([1-0.15;1+0.15],[my;my],'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','b')
text(2,my,asterisk_during,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
% line([2-0.15;2+0.15],[my;my],'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','b')
text(3,my,asterisk_after,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
% line([3-0.15;3+0.15],[my;my],'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','b')
% line(linex,liney,'color',[0 0 1],'linewidth',1.5,'markersize',10,'markerfacecolor','b')
% text(1.5-0.15,liney(1,1),asterisk_sensory,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17,'color',[0 0 1])        
        hold off
        ylim_data=[0 14];
        x1limits = [0.5 3.5];   x1ticks = [1,2,3];    
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',ylim_data,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','FontSize', 16,'xticklabel',{'Before','During','After'} ,'box', 'off'); %'fontweight', 'bold', 
        xlabel('Sensory Stimulation' ,'FontSize', 16);
        ylabel('Vm STD [mV]' ,'FontSize', 16);    
        title(['Vm trial-to-trial variability, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
        %% bar plot of mean Vm NB- and NB+ before, during and after sensory stimulation
    stim_num=1;
   j2=figure;  
    barwidth1=0.3;
        hold on
%         errbar_h1=errorbar((1:3)-0.15,nanmean([peaks_stat(stim_num).pre_response_M(:,1),peaks_stat(stim_num).Vm_res_M_meanstim(:,1),peaks_stat(stim_num).post_train_M(:,1)],1),zeros(1,3), nanstd([peaks_stat(stim_num).pre_response_M(:,1),peaks_stat(stim_num).Vm_res_M_meanstim(:,1),peaks_stat(stim_num).post_train_M(:,1)],0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         fn_errorbar_capsize(errbar_h1,capincrease,asymmetric)
%         errbar_h2=errorbar((1:3)+0.15,nanmean([peaks_stat(stim_num).pre_response_M(:,2),peaks_stat(stim_num).Vm_res_M_meanstim(:,2),peaks_stat(stim_num).post_train_M(:,2)],1),zeros(1,3), nanstd([peaks_stat(stim_num).pre_response_M(:,2),peaks_stat(stim_num).Vm_res_M_meanstim(:,2),peaks_stat(stim_num).post_train_M(:,2)],0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         fn_errorbar_capsize(errbar_h2,capincrease,asymmetric)
        errbar_h1=errorbar((1:3)-0.15,nanmean([peaks_stat(stim_num).pre_response_M(:,1),peaks_stat(stim_num).Vm_res_M_meanstim(:,1),peaks_stat(stim_num).post_train_M(:,1)],1), nanstd([peaks_stat(stim_num).pre_response_M(:,1),peaks_stat(stim_num).Vm_res_M_meanstim(:,1),peaks_stat(stim_num).post_train_M(:,1)],0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h2=errorbar((1:3)+0.15,nanmean([peaks_stat(stim_num).pre_response_M(:,2),peaks_stat(stim_num).Vm_res_M_meanstim(:,2),peaks_stat(stim_num).post_train_M(:,2)],1), nanstd([peaks_stat(stim_num).pre_response_M(:,2),peaks_stat(stim_num).Vm_res_M_meanstim(:,2),peaks_stat(stim_num).post_train_M(:,2)],0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar1=bar((1:3)-0.15, nanmean([peaks_stat(stim_num).pre_response_M(:,1),peaks_stat(stim_num).Vm_res_M_meanstim(:,1),peaks_stat(stim_num).post_train_M(:,1)],1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:));
        bar2=bar((1:3)+0.15, nanmean([peaks_stat(stim_num).pre_response_M(:,2),peaks_stat(stim_num).Vm_res_M_meanstim(:,2),peaks_stat(stim_num).post_train_M(:,2)],1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:));
%        line(pre_response_M_X,pre_response_M_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
%setting the errorbars not to appear in the legend:
hAnnotation = get(errbar_h1,'Annotation');  hLegendEntry = get(hAnnotation,'LegendInformation'); set(hLegendEntry,'IconDisplayStyle','off')
hAnnotation = get(errbar_h2,'Annotation');  hLegendEntry = get(hAnnotation,'LegendInformation'); set(hLegendEntry,'IconDisplayStyle','off')
% [legh,objh,outh,outm] = legend('NB Off','NB On','Location','northeast');
l = legend (legend_string,'fontsize',9,'Location','northeast'); legend('boxoff')
ylim_data=[get(gca,'ylim')]';
my=min(ylim_data)+5;
%   my=-1*max([peaks_stat(1).pre_response_M_std,peaks_stat(1).Vm_res_M_meanstim_std,peaks_stat(1).post_train_M_std])*1.1+max([peaks_stat(1).pre_response_M_m,peaks_stat(1).Vm_res_M_meanstim_m,peaks_stat(1).post_train_M_m]); 
linex=[1-0.15; 2-0.15];
liney=[my-0.1; my-0.1];
asterisk_sensory='*';
        if VmM_stat.p_before{1,1} >0.05 
    asterisk_before='n.s.';
else if VmM_stat.p_before{1,1}<0.05 && VmM_stat.p_before{1,1}>0.01
    asterisk_before='*';
    else if VmM_stat.p_before{1,1}<0.01 && VmM_stat.p_before{1,1}>0.001
            asterisk_before='**';
    else if VmM_stat.p_before{1,1}<0.001
             asterisk_before='***';
        end
        end
    end
end

if VmM_stat.p_during{1,1} >0.05 
    asterisk_during='n.s.';
else if VmM_stat.p_during{1,1}<0.05 && VmM_stat.p_during{1,1}>0.01
    asterisk_during='*';
    else if VmM_stat.p_during{1,1}<0.01 && VmM_stat.p_during{1,1}>0.001
            asterisk_during='**';
    else if VmM_stat.p_during{1,1}<0.001
             asterisk_during='***';
        end
        end
    end
end
if VmM_stat.p_after{1,1} >0.05 
    asterisk_after='n.s.';
else if VmM_stat.p_after{1,1}<0.05 && VmM_stat.p_after{1,1}>0.01
    asterisk_after='*';
    else if VmM_stat.p_after{1,1}<0.01 && VmM_stat.p_after{1,1}>0.001
            asterisk_after='**';
    else if VmM_stat.p_after{1,1}<0.001
             asterisk_after='***';
        end
        end
    end
end
text(1,my,asterisk_before,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
% line([1-0.15;1+0.15],[my;my],'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','b')
text(2,my,asterisk_during,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
% line([2-0.15;2+0.15],[my;my],'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','b')
text(3,my,asterisk_after,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
% line([3-0.15;3+0.15],[my;my],'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','b')
% set(gca,'ydir','reverse')
        hold off
        ylim_data=[-70 -45]; %get(gca,'ylim'); %[0 8.5];
        x1limits = [0.5 3.5];   x1ticks = [1,2,3];    
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',ylim_data,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','FontSize', 16,'xticklabel',{'Before','During','After'} ,'box', 'off'); %'fontweight', 'bold', 
        xlabel('Sensory Stimulation' ,'FontSize', 16);
        ylabel('Mean Vm [mV]' ,'FontSize', 16);    
        title(['Mean Vm, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
        %% Plot parameters along the train stim - version 3:bars+error bars of the mean values  
        %under construction, and not used because I am using the single
        %trial analysis for that
%     j3=figure;     
%     barwidth1=0.3;
%     nstim=1:11;
%         hold on
%         errbar_h=errorbar([1:size(peak_amp{2}(:,nstim),2)]-0.3,nanmean(peak_amp{2}(:,nstim),1),zeros(1,size(peak_amp{2}(:,nstim),2)), nanstd(peak_amp{2}(:,nstim),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         errbar_h=errorbar([1:size(peak_amp{3}(:,nstim),2)],nanmean(peak_amp{3}(:,nstim),1),zeros(1,size(peak_amp{3}(:,nstim),2)), nanstd(peak_amp{3}(:,nstim),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         bar([1:size(peak_amp{2}(:,nstim),2)]-0.3, nanmean(peak_amp{2}(:,nstim),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(peak_amp{3}(:,nstim),2)], nanmean(peak_amp{3}(:,nstim),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%   set(gca,'xlim',[0, length(nstim)+1],'xtick',[nstim], 'fontsize',16);
%         ylim_data=[get(gca,'ylim')]';
%         p_amp=max(ylim_data).*ones(size(s_amp));
%         plot(s_amp-0.15,p_amp,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Amplitude [mV]' ,'FontSize', 16);    
%         title(['Mean Response Amplitude - local baseline,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%        
%%        
        % Plot SNR1 (noise is from the Vm during the train)
stim_num=1;
SNR1_Y=peaks_stat(stim_num).SNR1';
SNR1_X(1,:)=ones(1,size(SNR1_Y,2));
SNR1_X(2,:)=2*ones(1,size(SNR1_Y,2));
E =peaks_stat(stim_num).SNR1_std;
linex=[1;2];
my=max(max(SNR1_Y))*1.1; 
liney=[my;my];
if peaks_stat(stim_num).wilcoxon_p_SNR1>0.05 
    asterisk='n.s.';
else if peaks_stat(stim_num).wilcoxon_p_SNR1<0.05 && peaks_stat(stim_num).wilcoxon_p_SNR1>0.01
    asterisk='*';
    else if peaks_stat(stim_num).wilcoxon_p_SNR1<0.01 && peaks_stat(stim_num).wilcoxon_p_SNR1>0.001
            asterisk='**';
    else if peaks_stat(stim_num).wilcoxon_p_SNR1<0.001
             asterisk='***';
        end
        end
    end
end
        k16=figure;
hold on
line(SNR1_X,SNR1_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(SNR1_X(:,1), mean(SNR1_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','fontsize',17) %,'verticalAlignment','bottom'
% boxplot(peaks_stat(stim_num).SNR1);
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [-1.5 11];
        y1ticks = [0,5,10];
        set( gca, 'xlim', x1limits,'xtick', x1ticks, 'ylim', y1limits,'ytick',y1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('SNR', 'FontSize', 28,'fontname', 'arial');
        title(['SNR1, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_SNR1)] ,'FontSize', 20,'fontname', 'arial');   

%plot SNR2
        SNR2_Y=peaks_stat(stim_num).SNR2';
SNR2_X(1,:)=ones(1,size(SNR2_Y,2));
SNR2_X(2,:)=2*ones(1,size(SNR2_Y,2));
E =peaks_stat(stim_num).SNR2_std;
        k17=figure;
hold on
line(SNR2_X,SNR2_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(SNR2_X(:,1), mean(SNR2_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(peaks_stat(stim_num).SNR2);
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('SNR', 'FontSize', 28,'fontname', 'arial');
        title(['SNR2, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_SNR2)] ,'FontSize', 20,'fontname', 'arial');   

        % Plot Amplitude signal 
Amplitude_signal_Y=peaks_stat(stim_num).Amplitude_signal';
Amplitude_signal_X(1,:)=ones(1,size(Amplitude_signal_Y,2));
Amplitude_signal_X(2,:)=2*ones(1,size(Amplitude_signal_Y,2));
E =peaks_stat(stim_num).Amplitude_signal_std;
linex=[1;2];
my=max(max(Amplitude_signal_Y))*1.1+0.3; 
liney=[my;my];
if peaks_stat(stim_num).wilcoxon_p_Amplitude_signal>0.05 
    asterisk='n.s.';
else if peaks_stat(stim_num).wilcoxon_p_Amplitude_signal<0.05 && peaks_stat(stim_num).wilcoxon_p_Amplitude_signal>0.01
    asterisk='*';
    else if peaks_stat(stim_num).wilcoxon_p_Amplitude_signal<0.01 && peaks_stat(stim_num).wilcoxon_p_Amplitude_signal>0.001
            asterisk='**';
    else if peaks_stat(stim_num).wilcoxon_p_Amplitude_signal<0.001
             asterisk='***';
        end
        end
    end
end
        k18=figure;
hold on
line(Amplitude_signal_X,Amplitude_signal_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Amplitude_signal_X(:,1), mean(Amplitude_signal_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
% boxplot(peaks_stat(stim_num).Amplitude_signal);
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 7];
        y1ticks = [0,3,6];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim', y1limits,'ytick', y1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Signal Amplitude', 'FontSize', 28,'fontname', 'arial');
        title(['Amplitude signal, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_Amplitude_signal)] ,'FontSize', 20,'fontname', 'arial');   

        % Plot Amplitude noise1 (noise is from the Vm during the train)
Amplitude_noise1_Y=peaks_stat(stim_num).Amplitude_noise1';
Amplitude_noise1_X(1,:)=ones(1,size(Amplitude_noise1_Y,2));
Amplitude_noise1_X(2,:)=2*ones(1,size(Amplitude_noise1_Y,2));
E =peaks_stat(stim_num).Amplitude_noise1_std;
linex=[1;2];
my=max(max(Amplitude_noise1_Y))*1.3; 
liney=[my;my];
if peaks_stat(stim_num).wilcoxon_p_Amplitude_noise1>0.05 
    asterisk='n.s.';
else if peaks_stat(stim_num).wilcoxon_p_Amplitude_noise1<0.05 && peaks_stat(stim_num).wilcoxon_p_Amplitude_noise1>0.01
    asterisk='*';
    else if peaks_stat(stim_num).wilcoxon_p_Amplitude_noise1<0.01 && peaks_stat(stim_num).wilcoxon_p_Amplitude_noise1>0.001
            asterisk='**';
    else if peaks_stat(stim_num).wilcoxon_p_Amplitude_noise1<0.001
             asterisk='***';
        end
        end
    end
end

        k19=figure;
hold on
line(Amplitude_noise1_X,Amplitude_noise1_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Amplitude_noise1_X(:,1), mean(Amplitude_noise1_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','fontsize',17) %,'verticalAlignment','bottom'
% boxplot(peaks_stat(stim_num).Amplitude_noise1);
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 7];
        y1ticks = [0,3,6];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim', y1limits,'ytick', y1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Noise Amplitude', 'FontSize', 28,'fontname', 'arial');
        title(['Amplitude noise1, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_Amplitude_noise1)] ,'FontSize', 20,'fontname', 'arial');   

                % Plot Amplitude noise2 (noise is from the Vm before the train)
Amplitude_noise2_Y=peaks_stat(stim_num).Amplitude_noise2';
Amplitude_noise2_X(1,:)=ones(1,size(Amplitude_noise2_Y,2));
Amplitude_noise2_X(2,:)=2*ones(1,size(Amplitude_noise2_Y,2));
E =peaks_stat(stim_num).Amplitude_noise2_std;
        k20=figure;
hold on
line(Amplitude_noise2_X,Amplitude_noise2_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Amplitude_noise2_X(:,1), mean(Amplitude_noise2_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(peaks_stat(stim_num).Amplitude_noise2);
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Noise Amplitude', 'FontSize', 28,'fontname', 'arial');
        title(['Amplitude noise2, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).wilcoxon_p_Amplitude_noise2)] ,'FontSize', 20,'fontname', 'arial');   
        %% save figures
if save_flag==1
%     switch exp_type
%         case 1
%             cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'
%         case 2
%             cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Traces+std+mean+summary'
%     end
 cd(path_output)          
% saveas(h1,'Train Amplitude Local', 'fig') 
% print(h1,'Train Amplitude Local','-dpng','-r600','-opengl')
saveas(h2,'Train Amplitude Local', 'fig') 
print(h2,'Train Amplitude Local','-dpng','-r600','-opengl')
saveas(h3,'Train Peak Latency', 'fig') 
print(h3,'Train Peak Latency','-dpng','-r600','-opengl')
saveas(h4,'Train Onset Latency', 'fig') 
print(h4,'Train Onset Latency','-dpng','-r600','-opengl')
% saveas(h5,'Train Half-Width', 'fig') 
% print(h5,'Train Half-Width','-dpng','-r600','-opengl')
% saveas(h6,'Train Peak_Value', 'fig') 
% print(h6,'Train Peak_Value','-dpng','-r600','-opengl')
% saveas(h7,'Train Local Baseline', 'fig') 
% print(h7,'Train Local Baseline','-dpng','-r600','-opengl')
% saveas(h8,'Train Response Mean Vm', 'fig') 
% print(h8,'Train Response Mean Vm','-dpng','-r600','-opengl')
saveas(h9,'Train Response Mean Vm STD', 'fig') 
print(h9,'Train Response Mean Vm STD','-dpng','-r600','-opengl')
% % 
% % 
saveas(g1,'Global Baseline', 'fig')         
print(g1,'Global Baseline','-dpng','-r600','-opengl') 
saveas(g2,'Adaptation Amplitude Ratio', 'fig') 
print(g2,'Adaptation Amplitude Ratio','-dpng','-r600','-opengl') 
saveas(g3,'Adaptation Power Ratio', 'fig') 
print(g3,'Adaptation Power Ratio','-dpng','-r600','-opengl') 
% saveas(g4,'Vm STD_Pre-train', 'fig')
% print(g4,'Vm STD_Pre-train','-dpng','-r600','-opengl') 
% saveas(g5,'Vm STD_Train', 'fig') 
% print(g5,'Vm STD_Train','-dpng','-r600','-opengl') 
% saveas(g6,'Vm STD_Post-train', 'fig') 
% print(g6,'Vm STD_Post-train','-dpng','-r600','-opengl') 

saveas(k16,'SNR1', 'fig') 
print(k16,'SNR1','-dpng','-r600','-opengl') 
saveas(k17,'SNR2', 'fig') 
print(k17,'SNR2','-dpng','-r600','-opengl') 
saveas(k18,'Amplitude_Signal_v2', 'fig') 
print(k18,'Amplitude_Signal_v2','-dpng','-r600','-opengl') 
saveas(k19,'Amplitude_Noise1_v2', 'fig') 
print(k19,'Amplitude_Noise1_v2','-dpng','-r600','-opengl') 
saveas(k20,'Amplitude_Noise2', 'fig') 
print(k20,'Amplitude_Noise2','-dpng','-r600','-opengl') 

saveas(j1,'Vm STD_Before During After Sensory stim_v2', 'fig') 
print(j1,'Vm STD_Before During After Sensory stim_v2','-dpng','-r600','-opengl') 
saveas(j2,'Vm M_Before During After Sensory stim_v4', 'fig') 
print(j2,'Vm M_Before During After Sensory stim_v4','-dpng','-r600','-opengl') 
end


 %%
% filename =  'SNR_plot_12_cells'; %'F62P6X6+9_zero_traces+std+mean_zoom-in';
% % cd 'd:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations';
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH';
% print (6, '-deps', filename)
cd(path_output)
filename='response_parameters'; 
save(filename, 'files_to_analyze', 'peaks', 'peaks_stat','VmSTD_stat','VmM_stat','amp_stat')
