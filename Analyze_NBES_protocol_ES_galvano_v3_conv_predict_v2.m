%% Analyze NBES protocol ES+galvano v3 Ldeconvs - recent version.
% This file was created on 26/04/2016 
%This file is used for the analysis of files created with extract_NBES_Data_v2

%       protocol 5 (NBES+galvnao (5 x-values))
%       protocol 6 (ChAT intracellular current injections light+galvnao (15 x-values))
%       protocol 8 - ChAT intracellular: light only (3 x-values): 1) train
%       20 Hz, 10ms ondur, 50 pulses. 2) single pulse 5 sec ondur. 3)single pulse 10 ms
%       protocol 9 - IC: I-V curve
%       protocol 10 (NBES+galvnao train+test (3 x-values))
%                 11 - NBES: 3 currents+galvano train+test (6 x-values)

%look at fft of signals and kernels. check why the kernel is so noisy. perhaps decide on a certain
%filter and find its parameters. Wiener deconvolution (wiki). need to know the SNR in each
%frequency. 
%kernel fft and remove high frequencies. filter kernel 1-160.
%% for opening workspace saved 
clear all
close all
        %%
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP+Vm filtered\Deconv ExpVar\Cross-Validation Off'
% load('data_for_deconv_analysis');

 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
files_to_analyze =[44,46,48,50,52,56,58,62,72,75];
save_flag=0;
cross_validation_flag=1;
kernel_dur=50; %[ms]           %if kernel_dur is negative, then kernelsize will be the size of the whole interval taken 
w=0.05;
type=1;     %1 for spont., 2 for evoked
type_name{1}='spont'; type_name{2}='evoked';
 kernel_peak_loc=zeros(length(files_to_analyze),2);
 kernel_peak_val=zeros(length(files_to_analyze),2);
 kernel_peak_lag=zeros(length(files_to_analyze),2);
 
    for fileind=1:length(files_to_analyze) ;         
    clearvars -except files_to_analyze fileind files save_flag type type_name cross_validation_flag R R_shuff R_M R_STD R_shuff_M R_shuff_STD G G_shuff...
        kernel_dur kernelc kernelc_shuff interval kernelc_M kernelc_STD kernelc_shuff_M kernelc_shuff_STD kernel_peak_loc kernel_peak_val...
        sig1 sig1_shuff sig2 sig2_shuff sig2_predict sig2_shuff_predict cc_sig1sig2 cc_sig1sig2_predict cc_sig1sig2_shuff cc_sig1sig2_shuff_predict...
        cc_sig2 cc_sig2_predict cc_sig2_shuff cc_sig2_shuff_predict w kernel_peak_lag MSE MSE_M MSE_STD MSE_shuff MSE_shuff_M MSE_shuff_STD...
        Ch2_data_filt data_no_spikes_filt
    
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
 %%   
                sf{1} = Param.sf_Vm;
                sf{2} = Param.sf_I1;
                sf{3} = Param.sf_V2;
                sf{4} = Param.sf_I2;
                dt=1/sf{channel};
                             
                sf_airpuff = Param.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                sf_galvano = Param.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;   
% raw_data{3} = raw_data{3}./20; %dividing by the LFP gain

trials = size(raw_data{3},2);

   %% cross-correlation Vm-LFP - spontaneous activity
       
       % Dividing LFP signal by the amplifier gain:
raw_data{3}=raw_data{3}./20;
        %bandpass filtering to remove noise from LFP channel
        bp_filt=files(files_to_analyze(fileind)).V2_filter;
        for xx=1:3
            for trace = 1:trials;    
                    jl=raw_data{3}(:,trace,xx);
                    Ch2_data_filt{fileind}(:,trace,xx)=bandPass_fft_IL_NEW2016(jl,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
                    jm=Ch2_data_filt{fileind}(:,trace,xx);
                    Ch2_data_filt{fileind}(:,trace,xx)=bandPass_fft_IL_NEW2016(jm,dt,bp_filt(1),bp_filt(2),0,0); %filtering LFP 1-160Hz
                    kl=data_no_spikes{channel}(:,trace,xx);
                    data_no_spikes_filt{fileind}(:,trace,xx)=bandPass_fft_IL_NEW2016(kl,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
                     km=data_no_spikes_filt{fileind}(:,trace,xx);
%                      data_no_spikes_filt{fileind}(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,bp_filt(1),bp_filt(2),0,0); %filtering Vm same as LFP
            end %temp
        end
    end
%                 cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP+Vm filtered\Deconv ExpVar\Cross-Validation Off'
% save('data_for_deconv_analysis');
%%
 for fileind=1:length(files_to_analyze) ;  
             for trace_type= type; %1:2;  %1 for spont., 2 for evoked
        clear  start_time duration start_time start_sample end_sample interval interval_mat x y patch_xdata patch_ydata yex ylim_data
        coeffs=[]; 
        % U=unique(x_value);
        % trace_type=length(U); %for trace_type=1 plot the xcorr before and after the ES according to the two start times. for trace_type=2 plot the xcorr from the same time point in the 2 different x-values.
        switch trace_type
            case 1
                 start_time = [0.4,5]; %[sec] %[0,5]
                 duration = 2.5; %[sec] 
                    for t=1:length(start_time);
                     start_sample(:,t) = ceil(start_time(t).*sf{1});
                        if start_time(t)==0
                            start_sample(:,t) = 1;
                        end
                      end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
                      interval(:,t) = start_sample(:,t):end_sample(:,t);
                      x_value=[1,1];
                    end 
            case 2
                x_value=[2:3];
                 start_sample =stim2_X{x_value(1)}(1,1);  %start with sensory stim
        %          duration = 1; %[sec]
                galvano_nstim = Param.facade(6);
                galvano_freq = Param.facade(7);
                duration = galvano_nstim./galvano_freq+0.05;
                end_sample = start_sample+duration.*sf{1}-1;
                interval(:,1) = round(start_sample:end_sample);
                interval(:,2) = interval(:,1);                
        end
        
        kernelsize=(kernel_dur./1000*sf{1}+1)*2;
             if kernel_dur<0
                 kernelsize=size(interval,1);
             end
clear data_Vm data_LFP data_Vm_noDC data_LFP_noDC data_Vm_filt data_Vm_filt_noDC
            for t=1:2
            %% Subtract mean from the interval (fragment of trace)
%             data_Vm{t}=data_no_spikes{1}(interval(:,t),:,x_value(t));
            data_LFP{t}=Ch2_data_filt{fileind}(interval(:,t),:,x_value(t));
%             data_Vm_noDC{t} = fn_Subtract_Mean(data_Vm{t});
            data_LFP_noDC{t} = fn_Subtract_Mean(data_LFP{t});
            % Subtract mean from the interval (fragment of trace) for filtered data
            data_Vm_filt{t}=data_no_spikes_filt{fileind}(interval(:,t),:,x_value(t));
            data_Vm_filt_noDC{t} = fn_Subtract_Mean(data_Vm_filt{t});
 %% Bootstrap to create shuffled Vm-LFP trace couples
%   xcorr_L=length(interval(:,t))-1;
 xcorr_L=size(data_LFP{t},1)-1; %0.5*sf{1}; %lag of 0.5sec
   CI_std=2;
    iterations=200;
    DC_sub=0; %put 1 if the mean trace was not subtracted from the data.
    for it=1:iterations; 
[data_LFP_shuff{t}, data_Vm_filt_shuff{t}]=fn_bootstrap_shuff_trace_couples(data_LFP{t}(:,:),data_Vm_filt{t}, 1, DC_sub);
% [data_LFP_shuff{t}, data_Vm_filt_shuff{t}]=fn_bootstrap_shuff_trace_couples(data_LFP_noDC{t}(:,:),data_Vm_filt_noDC{t}, 1, DC_sub);
 %deconvolution of the shuffled trace pairs:
         shuff_traces=size(data_LFP_shuff{t},2);
                    sig1_shuff{fileind,t}=data_LFP_shuff{t};
                    sig2_shuff{fileind,t}=data_Vm_filt_shuff{t};
                   
                        for trace = 1:shuff_traces;
                            G_shuff{fileind,t}(:,trace) = Ldeconvs(sig1_shuff{fileind,t}(:,trace),sig2_shuff{fileind,t}(:,trace),kernelsize,w,1);
                        end
                        
                             for trace = 1:shuff_traces;
                                if cross_validation_flag==1
                                    kernel=G_shuff{fileind,t};
                                    kernel(:,trace)=[];
                                    kernelc_shuff{fileind,t}(:,trace)=mean(kernel,2);
                                else
                                    kernelc_shuff{fileind,t}(:,trace) = G_shuff{fileind,t}(:,trace);
                                end
                               sig2_shuff_predict_tmp(:,trace)=conv(kernelc_shuff{fileind,t}(:,trace),sig1_shuff{fileind,t}(:,trace))+mean(sig2_shuff{fileind,t}(:,trace));
                                sig2_shuff_predict{fileind,t}(:,trace)=sig2_shuff_predict_tmp(kernelsize/2:end-kernelsize/2,trace);
                                [cc_sig2_shuff_predict{fileind,t}(:,trace), cc_sig2_shuff_predict_lag{fileind,t}(:,trace),cc_sig2_shuff_predict_ci{fileind,t}(:,trace)] =...
                                    crosscorr(sig2_shuff_predict{fileind,t}(:,trace), sig2_shuff{fileind,t}(:,trace),xcorr_L);          %,xcorr_L,CI_std           
                                [r, p] = corr(sig2_shuff_predict{fileind,t}(:,trace),sig2_shuff{fileind,t}(:,trace));
                                pearson_r_shuff_temp(trace,1)=r;
                                pearson_p_shuff_temp(trace,1)=p;
                                R_shuff_temp(trace,1) =  pearson_r_shuff_temp(trace,1)^2;
                                mse_shuff(trace,1) = fn_MSE(sig2_shuff_predict{fileind,t}(:,trace), sig2_shuff{fileind,t}(:,trace),1);
                                 [cc_sig2_shuff{fileind,t}(:,trace), cc_sig2_shuff_lag{fileind,t}(:,trace), cc_sig2_shuff_ci{fileind,t}(:,trace)] =...
                                     crosscorr(sig2_shuff{fileind,t}(:,trace), sig2_shuff{fileind,t}(:,trace),xcorr_L); %,xcorr_L,CI_std
                                 [cc_sig1sig2_shuff{fileind,t}(:,trace), cc_sig1sig2_shuff_lag{fileind,t}(:,trace), cc_sig1sig2_shuff_ci{fileind,t}(:,trace)] = ...
                                     crosscorr(sig1_shuff{fileind,t}(:,trace), sig2_shuff{fileind,t}(:,trace),xcorr_L); %,xcorr_L,CI_std
                                 [cc_sig1sig2_shuff_predict{fileind,t}(:,trace), cc_sig1sig2_shuff_predict_lag{fileind,t}(:,trace), cc_sig1sig2_shuff_predict_ci{fileind,t}(:,trace)] = ...
                                     crosscorr(sig1_shuff{fileind,t}(:,trace), sig2_shuff_predict{fileind,t}(:,trace),xcorr_L); %,xcorr_L,CI_std                                
                             end
                           pearson_r_shuff{fileind,t}(:,it) = pearson_r_shuff_temp;
                           pearson_p_shuff{fileind,t}(:,it) = pearson_p_shuff_temp;
                           R_shuff{fileind,t}(:,it) = R_shuff_temp;
                           MSE_shuff{fileind,t}(:,it)=mse_shuff;
                           kernelc_shuff_m{fileind,t}(:,it)=mean(kernelc_shuff{fileind,t}(:,:),2); %mean across traces
    end
    R_shuff_M{fileind,t}=mean(mean(R_shuff{fileind,t}(:,:),1));
    R_shuff_STD{fileind,t}=std(mean(R_shuff{fileind,t}(:,:),1)); %std across iterations, after averaging over the trials.                      
     MSE_shuff_M{fileind,t}=mean(mean(MSE_shuff{fileind,t}(:,1),1));
     MSE_shuff_STD{fileind,t}=std(mean(MSE_shuff{fileind,t}(:,1),1));
     kernelc_shuff_M{fileind,t}= mean(kernelc_shuff_m{fileind,t}(:,:),2);
     kernelc_shuff_STD{fileind,t}= std(kernelc_shuff_m{fileind,t}(:,:),0,2);
%             end %temporary end of t loop
%% Deconvolution of the real data          
%deconvolution of sig_1 and sig_2, convolving the kernel with sig_1 to create sig_2_predict.   
                    sig1{fileind,t}=data_LFP{t};
                    sig2{fileind,t}=data_Vm_filt{t};
%                     sig1{fileind,t}=data_LFP_noDC{t};
%                     sig2{fileind,t}=data_Vm_filt_noDC{t};
                for trace = 1:size(sig1{fileind,t},2);
                    G{fileind,t}(:,trace) = Ldeconvs(sig1{fileind,t}(:,trace),sig2{fileind,t}(:,trace),kernelsize,w,1);
                end
                
                %cross-validation of the kernel - averaging kernels from n-1 trials and testing on the last trial:
                 for trace = 1:1:size(sig1{fileind,t},2);
                    if cross_validation_flag==1    
                         kernel=G{fileind,t};
                         kernel(:,trace)=[];
                         kernelc{fileind,t}(:,trace)=mean(kernel,2);
                    else
                        kernelc{fileind,t}(:,trace) = G{fileind,t}(:,trace);
                    end
                                sig2_predict_tmp(:,trace)=conv(kernelc{fileind,t}(:,trace),sig1{fileind,t}(:,trace))+mean(sig2{fileind,t}(:,trace));
                                sig2_predict{fileind,t}(:,trace)=sig2_predict_tmp(kernelsize/2:end-kernelsize/2,trace);
                                [cc_sig2_predict{fileind,t}(:,trace), cc_sig2_predict_lag{fileind,t}(:,trace), cc_sig2_predict_ci{fileind,t}(:,trace)] =...
                                    crosscorr(sig2_predict{fileind,t}(:,trace), sig2{fileind,t}(:,trace),xcorr_L);      %,xcorr_L,CI_std
                                [r, p] = corr(sig2_predict{fileind,t}(:,trace),sig2{fileind,t}(:,trace));
                                pearson_r_temp(trace,1)=r;
                                pearson_p_temp(trace,1)=p;
                                R_temp(trace,1) =  pearson_r_temp(trace,1)^2;
                                mse(trace,1) = fn_MSE(sig2_predict{fileind,t}(:,trace), sig2{fileind,t}(:,trace),1);
                                 [cc_sig2{fileind,t}(:,trace), cc_sig2_lag{fileind,t}(:,trace), cc_sig2_ci{fileind,t}(:,trace)] =...
                                     crosscorr(sig2{fileind,t}(:,trace), sig2{fileind,t}(:,trace),xcorr_L); %,xcorr_L,CI_std
                                 [cc_sig1sig2{fileind,t}(:,trace), cc_sig1sig2_lag{fileind,t}(:,trace), cc_sig1sig2_ci{fileind,t}(:,trace)] =...
                                     crosscorr(sig1{fileind,t}(:,trace), sig2{fileind,t}(:,trace),xcorr_L); %,xcorr_L,CI_std
                                 [cc_sig1sig2_predict{fileind,t}(:,trace), cc_sig1sig2_predict_lag{fileind,t}(:,trace), cc_sig1sig2_predict_ci{fileind,t}(:,trace)] =...
                                     crosscorr(sig1{fileind,t}(:,trace), sig2_predict{fileind,t}(:,trace),xcorr_L); %,xcorr_L,CI_std
                 end
                           pearson_r{fileind,t}(:,1) = pearson_r_temp;
                           pearson_p{fileind,t}(:,1) = pearson_p_temp;
                           R{fileind,t}(:,1) = R_temp;     
                           R_M{fileind,t}=mean(R{fileind,t}(:,:),1);
                           R_STD{fileind,t}=std(R{fileind,t}(:,:),1); 
                           MSE{fileind,t}(:,1)=mse(:);
                           MSE_M{fileind,t}=mean(MSE{fileind,t}(:,1));
                           MSE_STD{fileind,t}=std(MSE{fileind,t}(:,1));
                           kernelc_M{fileind,t}=mean(kernelc{fileind,t}(:,:),2); %mean across traces
                           kernelc_STD{fileind,t}=std(kernelc{fileind,t}(:,:),0,2);    %mean across traces   
                           [kernel_peak_val(fileind,t), kernel_peak_loc(fileind,t)]=max(abs(kernelc_M{fileind,t}(:)));
                           kernel_peak_lag(fileind,t)=(kernel_peak_loc(fileind,t)-kernelsize/2)*dt*1000;
                          
            end
        end
    end  
    kernel_population(:,:,1)=cell2mat(kernelc_M(:,1)');
    kernel_population(:,:,2)=cell2mat(kernelc_M(:,2)');
    kernel_population_M(:,1)=mean(kernel_population(:,:,1),2);
    kernel_population_M(:,2)=mean(kernel_population(:,:,2),2);
    kernel_population_STD(:,1)=std(kernel_population(:,:,1),0,2);
    kernel_population_STD(:,2)=std(kernel_population(:,:,2),0,2);

    kernel_population_shuff(:,:,1)=cell2mat(kernelc_shuff_M(:,1)');
    kernel_population_shuff(:,:,2)=cell2mat(kernelc_shuff_M(:,2)');
    kernel_population_shuff_M(:,1)=mean(kernel_population_shuff(:,:,1),2);
    kernel_population_shuff_M(:,2)=mean(kernel_population_shuff(:,:,2),2);
    kernel_population_shuff_STD(:,1)=std(kernel_population_shuff(:,:,1),0,2);
    kernel_population_shuff_STD(:,2)=std(kernel_population_shuff(:,:,2),0,2);

            %% Statistics
            %MSE - absolute values + relative change
        deconv_stat.MSE=cell2mat(MSE_M);
        deconv_stat.MSE(8,:)=[];
        deconv_stat.MSE_m=nanmean(deconv_stat.MSE,1);
        deconv_stat.MSE_std=nanstd(deconv_stat.MSE,0,1);
        deconv_stat.change_MSE=[(deconv_stat.MSE(:,2)-deconv_stat.MSE(:,1))./abs(deconv_stat.MSE(:,2))].*100; %percent change
        deconv_stat.change_MSE_m=nanmean(deconv_stat.change_MSE,1);
        deconv_stat.change_MSE_std=nanstd(deconv_stat.change_MSE,0,1);      
        %testing for normal distribution       
        [deconv_stat.lillietest_h_MSE, deconv_stat.lillietest_p_MSE] = lillietest(deconv_stat.MSE(:,2)- deconv_stat.MSE(:,1));
        [deconv_stat.lillietest_h_change_MSE, deconv_stat.lillietest_p_change_MSE] = lillietest( deconv_stat.change_MSE);
        %paired ttest 
        [deconv_stat.ttest_h_MSE, deconv_stat.ttest_p_MSE]= ttest(deconv_stat.MSE(:,1),deconv_stat.MSE(:,2));   
        [deconv_stat.wilcoxon_p_MSE, deconv_stat.wilcoxon_h_MSE]= signrank(deconv_stat.MSE(:,1),deconv_stat.MSE(:,2));
        [deconv_stat.ttest_h_change_MSE, deconv_stat.ttest_p_change_MSE]= ttest(deconv_stat.change_MSE);
        [deconv_stat.wilcoxon_p_change_MSE, deconv_stat.wilcoxon_h_change_MSE]= signrank(deconv_stat.change_MSE);
      
             %MSE shuffled - absolute values + relative change
        deconv_stat.MSE_shuff=cell2mat(MSE_shuff_M);
        deconv_stat.MSE_shuff(8,:)=[];
        deconv_stat.MSE_shuff_m=nanmean(deconv_stat.MSE_shuff,1);
        deconv_stat.MSE_shuff_std=nanstd(deconv_stat.MSE_shuff,0,1);
        deconv_stat.change_MSE_shuff=[(deconv_stat.MSE_shuff(:,2)-deconv_stat.MSE_shuff(:,1))./abs(deconv_stat.MSE_shuff(:,2))].*100; %percent change
        deconv_stat.change_MSE_shuff_m=nanmean(deconv_stat.change_MSE_shuff,1);
        deconv_stat.change_MSE_shuff_std=nanstd(deconv_stat.change_MSE_shuff,0,1);      
        %testing for normal distribution       
        [deconv_stat.lillietest_h_MSE_shuff, deconv_stat.lillietest_p_MSE_shuff] = lillietest(deconv_stat.MSE_shuff(:,2)- deconv_stat.MSE_shuff(:,1));
        [deconv_stat.lillietest_h_change_MSE_shuff, deconv_stat.lillietest_p_change_MSE_shuff] = lillietest( deconv_stat.change_MSE_shuff);
        %paired ttest 
        [deconv_stat.ttest_h_MSE_shuff, deconv_stat.ttest_p_MSE_shuff]= ttest(deconv_stat.MSE_shuff(:,1),deconv_stat.MSE_shuff(:,2));   
        [deconv_stat.wilcoxon_p_MSE_shuff, deconv_stat.wilcoxon_h_MSE_shuff]= signrank(deconv_stat.MSE_shuff(:,1),deconv_stat.MSE_shuff(:,2));
        [deconv_stat.ttest_h_change_MSE_shuff, deconv_stat.ttest_p_change_MSE_shuff]= ttest(deconv_stat.change_MSE_shuff);
        [deconv_stat.wilcoxon_p_change_MSE_shuff, deconv_stat.wilcoxon_h_change_MSE_shuff]= signrank(deconv_stat.change_MSE_shuff);
      
            %Explained Variance - absolute values + relative change
        deconv_stat.R=cell2mat(R_M);
        deconv_stat.R_m=nanmean(deconv_stat.R,1);
        deconv_stat.R_std=nanstd(deconv_stat.R,0,1);
        deconv_stat.change_R=[(deconv_stat.R(:,2)-deconv_stat.R(:,1))./abs(deconv_stat.R(:,2))].*100; %percent change
        deconv_stat.change_R_m=nanmean(deconv_stat.change_R,1);
        deconv_stat.change_R_std=nanstd(deconv_stat.change_R,0,1);      
        %testing for normal distribution       
        [deconv_stat.lillietest_h_R, deconv_stat.lillietest_p_R] = lillietest(deconv_stat.R(:,2)- deconv_stat.R(:,1));
        [deconv_stat.lillietest_h_change_R, deconv_stat.lillietest_p_change_R] = lillietest( deconv_stat.change_R);
        %paired ttest 
        [deconv_stat.ttest_h_R, deconv_stat.ttest_p_R]= ttest(deconv_stat.R(:,1),deconv_stat.R(:,2));   
        [deconv_stat.wilcoxon_p_R, deconv_stat.wilcoxon_h_R]= signrank(deconv_stat.R(:,1),deconv_stat.R(:,2));
        [deconv_stat.ttest_h_change_R, deconv_stat.ttest_p_change_R]= ttest(deconv_stat.change_R);
        [deconv_stat.wilcoxon_p_change_R, deconv_stat.wilcoxon_h_change_R]= signrank(deconv_stat.change_R);
        
            %Explained Variance shuffled - absolute values + relative change
        deconv_stat.R_shuff=cell2mat(R_shuff_M);
        deconv_stat.R_shuff_m=nanmean(deconv_stat.R_shuff,1);
        deconv_stat.R_shuff_std=nanstd(deconv_stat.R_shuff,0,1);
        deconv_stat.change_R_shuff=[(deconv_stat.R_shuff(:,2)-deconv_stat.R_shuff(:,1))./abs(deconv_stat.R_shuff(:,2))].*100; %peR_shuffcent change
        deconv_stat.change_R_shuff_m=nanmean(deconv_stat.change_R_shuff,1);
        deconv_stat.change_R_shuff_std=nanstd(deconv_stat.change_R_shuff,0,1);      
        %testing for normal distribution       
        [deconv_stat.lillietest_h_R_shuff, deconv_stat.lillietest_p_R_shuff] = lillietest(deconv_stat.R_shuff(:,2)- deconv_stat.R_shuff(:,1));
        [deconv_stat.lillietest_h_change_R_shuff, deconv_stat.lillietest_p_change_R_shuff] = lillietest( deconv_stat.change_R_shuff);
        %paired ttest 
        [deconv_stat.ttest_h_R_shuff, deconv_stat.ttest_p_R_shuff]= ttest(deconv_stat.R_shuff(:,1),deconv_stat.R_shuff(:,2));   
        [deconv_stat.wilcoxon_p_R_shuff, deconv_stat.wilcoxon_h_R_shuff]= signrank(deconv_stat.R_shuff(:,1),deconv_stat.R_shuff(:,2));
        [deconv_stat.ttest_h_change_R_shuff, deconv_stat.ttest_p_change_R_shuff]= ttest(deconv_stat.change_R_shuff);
        [deconv_stat.wilcoxon_p_change_R_shuff, deconv_stat.wilcoxon_h_change_R_shuff]= signrank(deconv_stat.change_R_shuff);
        
            %Coupling Lag - absolute values + relative change
        deconv_stat.coupling_lag=kernel_peak_lag;
        deconv_stat.coupling_lag_m=nanmean(deconv_stat.coupling_lag,1);
        deconv_stat.coupling_lag_std=nanstd(deconv_stat.coupling_lag,0,1);
        deconv_stat.change_coupling_lag=[(deconv_stat.coupling_lag(:,2)-deconv_stat.coupling_lag(:,1))./abs(deconv_stat.coupling_lag(:,2))].*100; %percent change
        deconv_stat.change_coupling_lag_m=nanmean(deconv_stat.change_coupling_lag,1);
        deconv_stat.change_coupling_lag_std=nanstd(deconv_stat.change_coupling_lag,0,1);      
        %testing for normal distribution       
        [deconv_stat.lillietest_h_coupling_lag, deconv_stat.lillietest_p_coupling_lag] = lillietest(deconv_stat.coupling_lag(:,2)- deconv_stat.coupling_lag(:,1));
        [deconv_stat.lillietest_h_change_coupling_lag, deconv_stat.lillietest_p_change_coupling_lag] = lillietest( deconv_stat.change_coupling_lag);
        %paired ttest 
        [deconv_stat.ttest_h_coupling_lag, deconv_stat.ttest_p_coupling_lag]= ttest(deconv_stat.coupling_lag(:,1),deconv_stat.coupling_lag(:,2));   
        [deconv_stat.wilcoxon_p_coupling_lag, deconv_stat.wilcoxon_h_coupling_lag]= signrank(deconv_stat.coupling_lag(:,1),deconv_stat.coupling_lag(:,2));
        [deconv_stat.ttest_h_change_coupling_lag, deconv_stat.ttest_p_change_coupling_lag]= ttest(deconv_stat.change_coupling_lag);
        [deconv_stat.wilcoxon_p_change_coupling_lag, deconv_stat.wilcoxon_h_change_coupling_lag]= signrank(deconv_stat.change_coupling_lag); 
        
         %Coupling Lag shuffled - absolute values + relative change
        deconv_stat. coupling_lag_shuff=kernel_peak_lag;
        deconv_stat. coupling_lag_shuff_m=nanmean(deconv_stat.coupling_lag_shuff,1);
        deconv_stat. coupling_lag_shuff_std=nanstd(deconv_stat.coupling_lag_shuff,0,1);
        deconv_stat.change_coupling_lag_shuff=[(deconv_stat.coupling_lag_shuff(:,2)-deconv_stat.coupling_lag_shuff(:,1))./abs(deconv_stat.coupling_lag_shuff(:,2))].*100; %percent change
        deconv_stat.change_coupling_lag_shuff_m=nanmean(deconv_stat.change_coupling_lag_shuff,1);
        deconv_stat.change_coupling_lag_shuff_std=nanstd(deconv_stat.change_coupling_lag_shuff,0,1);      
        %testing for normal distribution       
        [deconv_stat.lillietest_h_coupling_lag_shuff, deconv_stat.lillietest_p_coupling_lag_shuff] = lillietest(deconv_stat.coupling_lag_shuff(:,2)- deconv_stat.coupling_lag_shuff(:,1));
        [deconv_stat.lillietest_h_change_coupling_lag_shuff, deconv_stat.lillietest_p_change_coupling_lag_shuff] = lillietest( deconv_stat.change_coupling_lag_shuff);
        %paired ttest 
        [deconv_stat.ttest_h_coupling_lag_shuff, deconv_stat.ttest_p_coupling_lag_shuff]= ttest(deconv_stat.coupling_lag_shuff(:,1),deconv_stat.coupling_lag_shuff(:,2));   
        [deconv_stat.wilcoxon_p_coupling_lag_shuff, deconv_stat.wilcoxon_h_coupling_lag_shuff]= signrank(deconv_stat.coupling_lag_shuff(:,1),deconv_stat.coupling_lag_shuff(:,2));
        [deconv_stat.ttest_h_change_coupling_lag_shuff, deconv_stat.ttest_p_change_coupling_lag_shuff]= ttest(deconv_stat.change_coupling_lag_shuff);
        [deconv_stat.wilcoxon_p_change_coupling_lag_shuff, deconv_stat.wilcoxon_h_change_coupling_lag_shuff]= signrank(deconv_stat.change_coupling_lag_shuff); 

            %%
    if save_flag==1;
        clear data data_LFP data_LFP_noDC data_LFP_shuff data_Vm_filt data_Vm_filt_noDC data_Vm_filt_shuff data_no_spikes data_no_spikes_filt Ch2_data_filt
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP+Vm filtered\Deconv ExpVar\Cross-Validation On'
%         cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP+Vm filtered\Deconv ExpVar\Cross-Validation Off'
        if trace_type==1;
            filename=['deconv_expvar_spont_',num2str(length(files_to_analyze)),'_cells_w', num2str(w*100)];
        end
         if trace_type==2;
            filename=['deconv_expvar_evoked_',num2str(length(files_to_analyze)),'_cells_w', num2str(w*100)];
         end
        save(filename)
%              'files_to_analyze', 'R','R_shuff', 'R_M', 'R_STD', 'R_shuff_M', 'R_shuff_STD', 'G', 'G_shuff',...
%                 'kernel_dur', 'kernelc', 'kernelc_shuff', 'interval','kernelc_M', 'kernelc_STD', 'kernelc_shuff_M',...
%                 'kernelc_shuff_STD', 'kernel_peak_val', 'kernel_peak_loc', 'kernel_peak_lag', 'w', 'MSE', 'MSE_M', 'MSE_STD','MSE_shuff', 'MSE_shuff_M', 'MSE_shuff_STD',...
%                 'sig1','sig2', 'sig1_shuff', 'sig2_shuff','sig2_predict_shuff', 'sig2_predict' )
    end
%% loading workspace:

% clear all
% type_name{1}='spont'; type_name{2}='evoked';
% load deconv_expvar_spont_10_cells_w20 
% load deconv_expvar_evoked_10_cells_w20
% load deconv_expvar_spont_10_cells_w5
% load deconv_expvar_evoked_10_cells_w5

% saveas(1,['kernel_w',num2str(w*100),'percent_file',num2str(fileind),'_',type_name{type},'.fig'])
% saveas(2,['Vm+Vmpredict_w',num2str(w*100),'percent_file',num2str(fileind),'_',type_name{type},'.fig'])


% % saveas(1,'kernel_w50percent_file1_spont.fig')
% % saveas(2,'Vm+Vmpredict_w50percent_file1_t2_spont.fig')
% % saveas(1,'kernel_w50percent_file1_evoked.fig')
% % saveas(2,'Vm+Vmpredict_w50percent_file1_t2_evoked.fig')


                     %% making figures - trace-specific kernel
%                      color_table=[ 0 0 0; [136 137 138]./256; [13 49 133]./256; [33 118 209]./256];   NB+ in blue
            color_table=[0 0 0; [216 22 22]/256; [136 137 138]/256; [255 153 153]/256]; %NB+ in red
                     fileind=2;
                     trace=3;
                  for t=1:2
                      if t==1
                          f1=figure;
                      else
                  f11= figure;  
                      end
                     hold on    
                      %kernel of specific trace
%                             h1=fn_shadedErrorBar([-1*kernelsize/2:kernelsize/2-1].*dt.*1000,kernelc{fileind,t}(:,trace),std(kernelc{fileind,t},0,2),{'LineWidth',1.5,'color', color_table(t,:)},1);
%                              h2=fn_shadedErrorBar([-1*kernelsize/2:kernelsize/2-1].*dt.*1000,kernelc_shuff{fileind,t}(:,trace),std(kernelc{fileind,t},0,2),{'LineWidth',1.5,'color', color_table(t+2,:)},1);  
                    %mean kernek
                            h1=fn_shadedErrorBar([-1*kernelsize/2:kernelsize/2-1].*dt.*1000,mean(kernelc{fileind,t},2),std(kernelc{fileind,t},0,2),{'LineWidth',1.5,'color', color_table(t,:)},1);
                             h2=fn_shadedErrorBar([-1*kernelsize/2:kernelsize/2-1].*dt.*1000,mean(kernelc_shuff{fileind,t},2),std(kernelc{fileind,t},0,2),{'LineWidth',1.5,'color', color_table(t+2,:)},1);  
                              ylim1=-0.6; %max(abs(mean(kernelc{fileind,1},2)));
                             ylim2=0.4; %-1*ylim1;
                             if trace_type==2
                             ylim1=-1; 
                             ylim2=0.2; 
                             end
                             xlim1=[-1*kernelsize/2,kernelsize/2-1].*dt.*1000;
                             set(gca,'ylim',[ylim1 ylim2],'xlim',xlim1,'linewidth',1, 'ticklength', [0.010 0.010],'fontsize',18,'fontname', 'arial' ,'box', 'off')
                            set(gca,'YDir','Reverse')
                             legend([h1.mainLine h2.mainLine],{'Actual Data','Shuffled Data'},'fontsize',12,'Location','northeast','box', 'off');
% getting from the legend the handles of the objects and then turning off the legend display of
% those that I don't want. An alternative is to access the handles of the vector objh and deleting them
% [legh,objh,outh,outm] = legend(gca,'show');
% hAnnotation = get(outh([1:3,5:7],1),'Annotation');  hLegendEntry = get(cell2mat(hAnnotation),'LegendInformation'); set(cell2mat(hLegendEntry),'IconDisplayStyle','off')
% legend(gca,'off')
% l = legend ('Data','Shuffled','fontsize',10,'Location','northeast','linewidth',2); legend('boxoff')
                            ylim_data=[get(gca,'ylim')]';
                            hold off                            
                  end
%% population kernel
                 for t=1:2
                  if t==1
                          f2=figure;
                      else
                  f21= figure;  
                      end
                     hold on    
                     if t==1
                            title(['Mean Population Kernel, NB Off'],'fontsize',18,'fontname', 'arial');
                     end
                    if t==2
                            title(['Mean Population Kernel, NB On'],'fontsize',18,'fontname', 'arial');
                    end
      %  CI using bootstrap:
      [prcntile1, prcntile2]=fn_get_CI_w_bootstrap(kernel_population(:,:,t),0);
 k=size(kernel_population(:,:,t),2);
                            h1=fn_shadedErrorBar([-1*kernelsize/2:kernelsize/2-1].*dt.*1000,mean(kernel_population(:,:,t),2),std(kernel_population(:,:,t),0,2)./sqrt(k),{'LineWidth',1.5,'color', color_table(t,:)},1);
                             h2=fn_shadedErrorBar([-1*kernelsize/2:kernelsize/2-1].*dt.*1000,mean(kernel_population_shuff(:,:,t),2),std(kernel_population_shuff(:,:,t)./sqrt(k),0,2),{'LineWidth',1.5,'color', color_table(t+2,:)},1);  
%                              plot([-1*kernelsize/2:kernelsize/2-1].*dt.*1000,prcntile1,'k-.')
%                              plot([-1*kernelsize/2:kernelsize/2-1].*dt.*1000,prcntile2,'k-.')
                             ylim1=-0.2; %max(abs(mean(kernelc{fileind,1},2)));
                             ylim2=0.1; %-1*ylim1;
                             if trace_type==2;
                               ylim1=-0.6; %max(abs(mean(kernelc{fileind,1},2)));
                                ylim2=0.15; %-1*ylim1;
                             end
                             xlim1=[-1*kernelsize/2,kernelsize/2-1].*dt.*1000;
                             set(gca,'ylim',[ylim1 ylim2],'xlim',xlim1,'linewidth',1, 'ticklength', [0.010 0.010],'fontsize',18,'fontname', 'arial' ,'box', 'off')
                             set(gca,'YDir','Reverse')
                             legend([h1.mainLine h2.mainLine],{'Actual Data','Shuffled Data'},'fontsize',12,'Location','northeast', 'box', 'off');
                            ylim_data=[get(gca,'ylim')]';
                            hold off                            
                         end
         %% Vm-Vm predicted - actual data
            color_table=[0 0 0; [216 22 22]/256; [136 137 138]/256; [255 153 153]/256]; %NB+ in red
fileind=2;
      trace = 3;                         
  leg_text{1,1}=['Vm ES off'];leg_text{2,1}=['Vm predicted ES off']; leg_text{3,1}=['Vm ES on']; leg_text{4,1}=['Vm predicted ES on'];        
                 for t=1:2;
                     x_axis=[];       
                     x_axis=[0:length(interval(:,t))-1].*dt.*1000;
%                      x_axis_s=[0:length(interval(:,t))-1];
                    if t==1
                          f3=figure;
                      else
                  f31= figure;  
                    end 
%                       set(gcf,'Units','centimeters');
                        sp1=subplot(2,1,1);
%                         ax1_pos=[1,8,5,3];
%                         ax2_pos=[ax1_pos(1), 1,ax1_pos(3),5];
%                         set(sp1,'Units','centimeters')%,'Position',ax1_pos); 
%                         sp1.Position=ax1_pos;
                        hold on     
                            plot(x_axis,sig2{fileind,t}(:,trace),'LineWidth',2,'color', color_table(t+2,:)); 
%                             plot(x_axis,sig1{fileind,t}(:,trace)-60,'LineWidth',2,'color',[0 0 1]); 
                            plot(x_axis,sig2_predict{fileind,t}(:,trace),'LineWidth',2,'color', color_table(t,:));   
%                              title(['f', num2str(files_to_analyze(fileind)),'Vm,Vm predicted - actual data']); 
                            text(x_axis(1),sig2_predict{fileind,t}(1,trace),[num2str(floor(sig2_predict{fileind,t}(1,trace))), ' mV '],'HorizontalAlignment','right','fontsize',13,'fontname','arial')
                             axis tight
                              ylim_data=[get(gca,'ylim')]';
                              xlim_data=[get(gca,'xlim')]';
                             
                             % plotting stim_2:
                     if trace_type==2
                        patch_xdata=[stim2_X{2}-interval(1,t); flipud(stim2_X{2}-interval(1,t))];
                        yex=wextend('1D','sym',ylim_data,1);
                        l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
                        temp_y=wextend('ac','sym',yex,l);
                        patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
                        patch_cdata=ones(size(patch_xdata));
                        patch(patch_xdata.*dt*1000,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
                    end
                    %plotting scale bar
                    yline_start=ylim_data(1)-2; yline_end=yline_start+5;
                    xline_start=x_axis(1)-100-350;
                    if t==2;
                     yline_start=ylim_data(1); yline_end=yline_start+2; 
                 end
                    if trace_type==2
                       yline_start=ylim_data(1)-2; yline_end=yline_start+5; 
                       xline_start=x_axis(1)-200;
                    end
                     xline_end=xline_start+200;                   
                    xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
                    stringh=[num2str(xline(2,1)-xline(1,1)), ' mS'];
                    stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
                    hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
                    hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
                    htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',13); %'color',[0 0 0]
                    htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',13);
            hold off
                    set(gca,'ylim',[yline_start-0.2, ylim_data(2)+0.25*(ylim_data(2)- ylim_data(1))],'xlim',[xline_start-0.2, xlim_data(2)]);                  
%                      legend({leg_text{2*t-1},leg_text{2*t}},'fontsize',12,'Location','northeast','box', 'off');
                    set(gca, 'visible', 'off') ;
                    
                    sp2=subplot(2,1,2);
%                     set(sp2,'Units','centimeters','Position',ax2_pos); 
                    hold on
                   plot(x_axis,sig1{fileind,t}(:,trace)-60,'LineWidth',2,'color',[0 0 1]);  
                   text(x_axis(1),sig1{fileind,t}(1,trace),[num2str(floor(sig1{fileind,t}(1,trace))), ' mV '],'HorizontalAlignment','right','fontsize',13,'fontname','arial')
                             axis tight
                              ylim_data=[get(gca,'ylim')]';
                              xlim_data=[get(gca,'xlim')]';
                              %plotting scale bar
                    yline_start=ylim_data(1)+0.2; yline_end=yline_start+0.2;
                    xline_start=x_axis(1)-100-350;
                    if t==2;
                     yline_start=ylim_data(1); yline_end=yline_start+0.2; 
                 end
                    if trace_type==2
                       yline_start=ylim_data(1)+0.2; yline_end=yline_start+0.5; 
                       xline_start=x_axis(1)-200;
                    end
                     xline_end=xline_start+200;                   
                    xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
                    stringh=[num2str(xline(2,1)-xline(1,1)), ' mS'];
                    stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
%                     hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
                    hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
%                     htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom','fontsize',13); %'color',[0 0 0]
                    htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',13);
            hold off
                    set(gca,'ylim',[yline_start-0.2, ylim_data(2)+0.25*(ylim_data(2)- ylim_data(1))],'xlim',[xline_start-0.2, xlim_data(2)]);    
                    set(gca,'YDir','Reverse')
                    set(gca, 'visible', 'off') ;
                    

   end                       
 %%             Vm-Vm predicted - shuffled data
           
                         fileind=2;
      trace = 3;                 
               
  leg_text{1,1}=['Vm ES off'];leg_text{2,1}=['Vm predicted ES off']; leg_text{3,1}=['Vm ES on']; leg_text{4,1}=['Vm predicted ES on'];                
                 for t=1:2;
                     x_axis=[];       
                     x_axis=[0:length(interval(:,t))-1].*dt.*1000;
                   if t==1
                          f4=figure;
                      else
                  f41= figure;  
                   end
                      sp1=subplot(2,1,1);
                        hold on                        
                            plot(x_axis,sig2_shuff{fileind,t}(:,trace),'LineWidth',2,'color', color_table(t+2,:)); 
                            plot(x_axis,sig2_shuff_predict{fileind,t}(:,trace),'LineWidth',2,'color', color_table(t,:));  
%                              title(['f', num2str(files_to_analyze(fileind)),'Vm,Vm predicted - shuffled data']);       
                             text(x_axis(1),sig2_shuff_predict{fileind,t}(1,trace),[num2str(floor(sig2_shuff_predict{fileind,t}(1,trace))), ' mV '],'HorizontalAlignment','right','fontsize',13,'fontname','arial')
                             axis tight
                             ylim_data=[get(gca,'ylim')]';
                              xlim_data=[get(gca,'xlim')]';
                             
                             % plotting stim_2:
                     if trace_type==2
                        patch_xdata=[stim2_X{2}-interval(1,t); flipud(stim2_X{2}-interval(1,t))];
                        yex=wextend('1D','sym',ylim_data,1);
                        l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
                        temp_y=wextend('ac','sym',yex,l);
                        patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
                        patch_cdata=ones(size(patch_xdata));
                        patch(patch_xdata.*dt*1000,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
                    end
                    %plotting scale bar
                     yline_start=ylim_data(1)-2; yline_end=yline_start+5;
                    xline_start=x_axis(1)-100-350;
                    if t==2;
                     yline_start=ylim_data(1); yline_end=yline_start+2; 
                 end
                    if trace_type==2
                       yline_start=ylim_data(1)-2; yline_end=yline_start+5; 
                       xline_start=x_axis(1)-200;
                    end
                    
                    xline_end=xline_start+200;                   
                    xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
                    stringh=[num2str(xline(2,1)-xline(1,1)), ' mS'];
                    stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
                    hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
                    hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
                    htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',13); %'color',[0 0 0]
                    htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',13);
            hold off
                    set(gca,'ylim',[yline_start-0.2, ylim_data(2)+0.25*(ylim_data(2)- ylim_data(1))],'xlim',[xline_start-0.2, xlim_data(2)]);
                     legend({leg_text{2*t-1},leg_text{2*t}},'fontsize',12,'Location','northeast','box', 'off');
                    set(gca, 'visible', 'off') ;
                    
                    sp2=subplot(2,1,2);
%                     set(sp2,'Units','centimeters','Position',ax2_pos); 
                    hold on
                   plot(x_axis,sig1_shuff{fileind,t}(:,trace)-60,'LineWidth',2,'color',[0 0 1]);  
                   text(x_axis(1),sig1{fileind,t}(1,trace),[num2str(floor(sig1{fileind,t}(1,trace))), ' mV '],'HorizontalAlignment','right','fontsize',13,'fontname','arial')
                             axis tight
                              ylim_data=[get(gca,'ylim')]';
                              xlim_data=[get(gca,'xlim')]';
                              %plotting scale bar
                    yline_start=ylim_data(1)+0.2; yline_end=yline_start+0.2;
                    xline_start=x_axis(1)-100-350;
                    if t==2;
                     yline_start=ylim_data(1); yline_end=yline_start+0.2; 
                 end
                    if trace_type==2
                       yline_start=ylim_data(1)+0.2; yline_end=yline_start+0.5; 
                       xline_start=x_axis(1)-200;
                    end
                     xline_end=xline_start+200;                   
                    xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
                    stringh=[num2str(xline(2,1)-xline(1,1)), ' mS'];
                    stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
%                     hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
                    hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
%                     htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom','fontsize',13); %'color',[0 0 0]
                    htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',13);
            hold off
                    set(gca,'ylim',[yline_start-0.2, ylim_data(2)+0.25*(ylim_data(2)- ylim_data(1))],'xlim',[xline_start-0.2, xlim_data(2)]);    
                    set(gca,'YDir','Reverse')
                    set(gca, 'visible', 'off') ;
                                       
                 end   
                 
                        

 %%              cross-correlation Vm-Vm predicted - actual data
               fileind=2; 
                leg_text{1,1}=['ccVm-Vm'];leg_text{2,1}=['ccVm-Vm predicted']; leg_text{3,1}=['ccVm-Vm']; leg_text{4,1}=['ccVm-Vm predicted'];            
                    for t=1:2;  
                        x_s=-1*xcorr_L:xcorr_L;
                     x_axis=[];       
                     x_axis=[-1*xcorr_L:xcorr_L].*dt*1000;
                      if t==1
                          f5=figure;
                      else
                  f51= figure;  
                      end
                        hold on
                        plot(x_axis,mean(cc_sig2{fileind,t}(:,:),2),'LineWidth',2,'color', color_table(t+2,:));
                        plot(x_axis,mean(cc_sig2_predict{fileind,t}(:,:),2),'LineWidth',2,'color', color_table(t,:));
                        set(gca,'xlim',[-1*xcorr_L,xcorr_L].*dt*1000);
%                         title(['f', num2str(files_to_analyze(fileind)),' mean ccVm-Vm predicted - actual data'])
                         legend({leg_text{2*t-1},leg_text{2*t}},'fontsize',12,'Location','northeast','box', 'off');
%                          axis tight
                y1limits=[-0.3, 1.1];
                if trace_type==2
                    y1limits=[-0.4, 1.1];
                end
                         set(gca, 'linewidth',1, 'ticklength', [0.010 0.010],'fontsize',18,'fontname', 'arial',...
                            'ylim',y1limits )
                          xlabel('Lag [mS]','fontsize',18,'fontname', 'arial');
                          ylabel('Correlation', 'fontsize',18,'fontname', 'arial')
                         hold off
                        end
  %%             cross-correlation Vm-Vm predicted - shuffled data                               
                         leg_text{1,1}=['ccVm-Vm'];leg_text{2,1}=['ccVm-Vm predicted']; leg_text{3,1}=['ccVm-Vm']; leg_text{4,1}=['ccVmVm predicted'];        
                        for t=1:2;
                             if t==1
                          f6=figure;
                      else
                  f61= figure;  
                      end
                        hold on
                        plot([-1*xcorr_L:xcorr_L].*dt*1000,mean(cc_sig2_shuff{fileind,t}(:,:),2),'LineWidth',2,'color', color_table(t+2,:));
                        plot([-1*xcorr_L:xcorr_L].*dt*1000,mean(cc_sig2_shuff_predict{fileind,t}(:,:),2),'LineWidth',2,'color', color_table(t,:));
                        set(gca,'xlim',[-1*xcorr_L,xcorr_L].*dt*1000);
%                         title(['f', num2str(files_to_analyze(fileind)),'mean ccVm-Vm predicted - shuffled data'])
                         legend({leg_text{2*t-1},leg_text{2*t}},'fontsize',12,'Location','northeast', 'box', 'off');
%                          axis tight
            y1limits=[-0.3, 1.1];
                if trace_type==2
                    y1limits=[-0.4, 1.1];
                end
                          set(gca, 'linewidth',1, 'ticklength', [0.010 0.010],'fontsize',18,'fontname', 'arial',...
                              'ylim',y1limits )
                          xlabel('Lag [mS]','fontsize',18,'fontname', 'arial');
                          ylabel('Correlation', 'fontsize',18,'fontname', 'arial')
                         hold off
                        end                     
     %% Explained variance - per single cell
     
% tmp_Y= [R{1,1}(:,:),R{1,2}(:,:)]';
% tmp_X(1,:)=ones(1,size(tmp_Y,2));
% tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
% E = [R_STD{1}, R_STD{2}];
% g1=figure;
% hold on
% line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
%         y1limits = [-0.1 0.7];
%         y1ticks = [-0.1,0,0.1,0.2,0.3,0.4,0.5];
%         set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1, 'ylim', y1limits,'ytick', y1ticks,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Explained Variance', 'FontSize', 28,'fontname', 'arial');
% %         title(['Spontaneous event onset value,  p=' num2str(event_ongoing_stat.wilcoxon_p_onVal_m)] ,'FontSize', 20,'fontname', 'arial');   
     %% Mean Squared Error - all cells    
tmp_Y= deconv_stat.MSE';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(deconv_stat.MSE);
linex=[1;2];
my=max(max(tmp_Y))*1.1; 
liney=[my;my];
if deconv_stat.wilcoxon_p_MSE>0.05 
    asterisk='n.s.';
else if deconv_stat.wilcoxon_p_MSE<0.05 && deconv_stat.wilcoxon_p_MSE>0.01
    asterisk='*';
    else if deconv_stat.wilcoxon_p_MSE<0.01 && deconv_stat.wilcoxon_p_MSE>0.001
            asterisk='**';
    else if deconv_stat.wilcoxon_p_MSE<0.001
             asterisk='***';
        end
        end
    end
end

g1=figure;
hold on
% line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
boxplot(deconv_stat.MSE)
line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 60];
        y1ticks = [0, 10, 20, 30, 40, 50];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1, 'ylim',y1limits, 'ytick',y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Mean Squared Error', 'FontSize', 28,'fontname', 'arial');
%         title(['Mean Squared Error ',type_name{type},  ', p=' num2str(deconv_stat.wilcoxon_p_MSE)] ,'FontSize', 20,'fontname', 'arial');   
clear tmp_Y tmp_X E

  % Mean Squared Error shuffled data - all cells    
tmp_Y= deconv_stat.MSE_shuff';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(deconv_stat.MSE_shuff);
linex=[1;2];
my=max(max(tmp_Y))*1.1; 
liney=[my;my];
if deconv_stat.wilcoxon_p_MSE_shuff>0.05 
    asterisk_shuff='n.s.';
else if deconv_stat.wilcoxon_p_MSE_shuff<0.05 && deconv_stat.wilcoxon_p_MSE_shuff>0.01
    asterisk_shuff='*';
    else if deconv_stat.wilcoxon_p_MSE_shuff<0.01 && deconv_stat.wilcoxon_p_MSE_shuff>0.001
            asterisk_shuff='**';
    else if deconv_stat.wilcoxon_p_MSE_shuff<0.001
             asterisk_shuff='***';
        end
        end
    end
end

g11=figure;
hold on
% line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
boxplot(deconv_stat.MSE_shuff)
line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_shuff,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 60];
        y1ticks = [0, 10, 20, 30, 40, 50];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1, 'ylim',y1limits, 'ytick',y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Mean Squared Error', 'FontSize', 28,'fontname', 'arial');
%         title(['Mean Squared Error ',type_name{type},  ', p=' num2str(deconv_stat.wilcoxon_p_MSE_shuff)] ,'FontSize', 20,'fontname', 'arial');   

g111=figure;
linex=[1 3;2 4];
my=max(max([deconv_stat.MSE, deconv_stat.MSE_shuff]))*1.1; 
liney=[my my;my my];
x1limits = [0.75 4.25];
x1ticks = [1,2,3,4];
hold on
% boxplot([R_mat(:,1), R_shuff_mat(:,1),R_mat(:,2), R_shuff_mat(:,2)])
boxplot([deconv_stat.MSE, deconv_stat.MSE_shuff])
line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(3.5,my,asterisk_shuff,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(1.5,my+5,'actual data','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(3.5,my+5,'shuffled data','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1, 'ylim', y1limits,'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',{'NB-','NB+','NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
    ylabel('Mean Squared Error', 'FontSize', 28,'fontname', 'arial');
clear tmp_Y tmp_X E
%% Explained Variance - all cells
R_mat=deconv_stat.R;
tmp_Y= R_mat';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(R_mat);
linex=[1;2];
my=max(max(tmp_Y))*1.1; 
liney=[my;my];
if deconv_stat.wilcoxon_p_R >0.05 
    asterisk='n.s.';
else if deconv_stat.wilcoxon_p_R<0.05 && deconv_stat.wilcoxon_p_R>0.01
    asterisk='*';
    else if deconv_stat.wilcoxon_p_R<0.01 && deconv_stat.wilcoxon_p_R>0.001
            asterisk='**';
    else if deconv_stat.wilcoxon_p_R<0.001
             asterisk='***';
        end
        end
    end
end
g2=figure;
hold on
% line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
boxplot(R_mat)
line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 0.5];
        if trace_type==2;
            y1limits = [0 0.7];
        end
        y1ticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1, 'ylim', y1limits,'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Explained Variance', 'FontSize', 28,'fontname', 'arial');
%  title(['Explained Variance ',type_name{type},  ', p=' num2str(deconv_stat.wilcoxon_p_R)] ,'FontSize', 20,'fontname', 'arial');   
clear tmp_Y tmp_X E

% Explained Variance shuffled - all cells
R_shuff_mat=deconv_stat.R_shuff;
tmp_Y= R_shuff_mat';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(R_shuff_mat);
linex=[1;2];
my=max(max(tmp_Y))*1.1; 
liney=[my;my];
if deconv_stat.wilcoxon_p_R_shuff >0.05 
    asterisk_shuff='n.s.';
else if deconv_stat.wilcoxon_p_R_shuff<0.05 && deconv_stat.wilcoxon_p_R_shuff>0.01
    asterisk_shuff='*';
    else if deconv_stat.wilcoxon_p_R_shuff<0.01 && deconv_stat.wilcoxon_p_R_shuff>0.001
            asterisk_shuff='**';
    else if deconv_stat.wilcoxon_p_R_shuff<0.001
             asterisk_shuff='***';
        end
        end
    end
end
g21=figure;
hold on
% line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
boxplot(R_shuff_mat)
line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_shuff,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)        
x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 0.5];
        if trace_type==2;
            y1limits = [0 0.7];
        end
        y1ticks = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1, 'ylim', y1limits,'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Explained Variance', 'FontSize', 28,'fontname', 'arial');
%  title(['Explained Variance ',type_name{type},  ', p=' num2str(deconv_stat.wilcoxon_p_R_shuff)] ,'FontSize', 20,'fontname', 'arial');   

g211=figure;
linex=[1 3;2 4];
my=max(max([R_mat, R_shuff_mat]))*1.1; 
liney=[my my;my my];
x1limits = [0.75 4.25];
x1ticks = [1,2,3,4];
hold on
% boxplot([R_mat(:,1), R_shuff_mat(:,1),R_mat(:,2), R_shuff_mat(:,2)])
boxplot([R_mat, R_shuff_mat])
line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(3.5,my,asterisk_shuff,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(1.5,my+0.05,'actual data','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
text(3.5,my+0.05,'shuffled data','HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',13)
set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1, 'ylim', y1limits,'ytick', y1ticks,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',{'NB-','NB+','NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
    ylabel('Explained Variance', 'FontSize', 28,'fontname', 'arial');
clear tmp_Y tmp_X E
%% Coupling lag (lag of kernel peak)
tmp_Y= kernel_peak_lag';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(kernel_peak_lag);
g3=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(kernel_peak_lag);
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [-0.1 0.5];
        y1ticks = [-0.1,0,0.1,0.2,0.3,0.4,0.5];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.0001 0.0001],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('coupling lag [mS]', 'FontSize', 28,'fontname', 'arial');
%  title(['Coupling Lag ',type_name{type},  ', p=' num2str(deconv_stat.wilcoxon_p_coupling_lag)] ,'FontSize', 20,'fontname', 'arial');   
clear tmp_Y tmp_X E
  %%                    
%                         figure
%                          for t=1:2;
%                     subplot(1,2,t); hold on
%                         plot([-1*xcorr_L:xcorr_L].*dt*1000,cc_sig1sig2{fileind,t}(:,trace),'LineWidth',1,'color', color_table(2*t-1,:));
%                         plot([-1*xcorr_L:xcorr_L].*dt*1000,cc_sig1sig2_predict{fileind,t}(:,trace),'LineWidth',1,'color', color_table(2*t,:));
%                         set(gca,'xlim',[-500,500]);
%                         hold off
%                         title('ccLFP-Vm, ccLFP-Vm predicted')                    
%                          end
%%
if save_flag==1;
saveas(f1,['f', num2str(files_to_analyze(fileind)),'_Mean Kernel_NB Off_',type_name{type},'.fig'])         
print(f1,['f', num2str(files_to_analyze(fileind)),'_Mean Kernel_NB Off_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f11,['f', num2str(files_to_analyze(fileind)),'_Mean Kernel_NB On_',type_name{type},'.fig'])
print(f11,['f', num2str(files_to_analyze(fileind)),'_Mean Kernel_NB On_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f2,['Mean Population Kernel_NB Off_',type_name{type},'.fig'])  
print(f2,['Mean Population Kernel_NB Off_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f21,['Mean Population Kernel_NB On_',type_name{type},'.fig']) 
print(f21,['Mean Population Kernel_NB On_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f3,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_VmVmpredict_actual data_NB Off_',type_name{type},'.fig'])     
print(f3,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_VmVmpredict_actual data_NB Off_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f31,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_VmVmpredict_actual data_NB On_',type_name{type},'.fig'])  
print(f31,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_VmVmpredict_actual data_NB On_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f4,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_VmVmpredict_shuffled data_NB Off_',type_name{type},'.fig'])   
print(f4,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_VmVmpredict_shuffled data_NB Off_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f41,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_VmVmpredict_shuffled data_NB On_',type_name{type},'.fig'])   
print(f41,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_VmVmpredict_shuffled data_NB On_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f5,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_ccVmVmpredict_actual data_NB Off_',type_name{type},'.fig'])    
print(f5,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_ccVmVmpredict_actual data_NB Off_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f51,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_ccVmVmpredict_actual data_NB On_',type_name{type},'.fig'])    
print(f51,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_ccVmVmpredict_actual data_NB On_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f6,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_ccVmVmpredict_shuffled data_NB Off_',type_name{type},'.fig'])  
print(f6,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_ccVmVmpredict_shuffled data_NB Off_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(f61,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_ccVmVmpredict_shuffled data_NB On_',type_name{type},'.fig'])  
print(f61,['f', num2str(files_to_analyze(fileind)),'_t',num2str(trace),'_ccVmVmpredict_shuffled data_NB On_',type_name{type}],'-dpng','-r600','-opengl') 
saveas(g1,['Mean Squared Error_', type_name{type},'.fig'])         
print(g1,['Mean Squared Error_', type_name{type}],'-dpng','-r600','-opengl') 
saveas(g11,['Mean Squared Error_shuffled_', type_name{type},'.fig'])         
print(g11,['Mean Squared Error_shuffled_', type_name{type}],'-dpng','-r600','-opengl') 
saveas(g111,['Mean Squared Error_actual+shuffled_', type_name{type},'.fig'])         
print(g111,['Mean Squared Error_actual_shuffled_', type_name{type}],'-dpng','-r600','-opengl') 
saveas(g2,['Explained Variance_', type_name{type},'.fig']) 
print(g2,['Explained Variance_', type_name{type}],'-dpng','-r600','-opengl') 
saveas(g21,['Explained Variance_shuffled_', type_name{type},'.fig']) 
print(g21,['Explained Variance_shuffled_', type_name{type}],'-dpng','-r600','-opengl') 
saveas(g211,['Explained Variance_actual+shuffled_', type_name{type},'.fig']) 
print(g211,['Explained Variance_actual+shuffled_', type_name{type}],'-dpng','-r600','-opengl') 
saveas(g3,['Coupling Lag_', type_name{type},'.fig']) 
print(g3,['Coupling Lag_', type_name{type}],'-dpng','-r600','-opengl') 
% saveas(g13,['Coupling Lag_shuffled_', type_name{type},'.fig']) 
% print(g13,['Coupling Lag_shuffled_', type_name{type}],'-dpng','-r600','-opengl') 
end
%% rmANOVA with matlab function - R
clear all
load deconv_expvar_spont_10_cells_w5.mat
R_mat=cell2mat(R_M);
R_spont_noNB(:,1)=R_mat(:,1);
R_spont_NB(:,1)=R_mat(:,2);

clear R_mat
load deconv_expvar_evoked_10_cells_w5.mat
R_mat=cell2mat(R_M);
R_evoked_noNB(:,1)=R_mat(:,1);
R_evoked_NB(:,1)=R_mat(:,2);
R_spont(:,1)=[R_spont_noNB;R_spont_NB];
R_evoked(:,1)=[R_evoked_noNB;R_evoked_NB];
 S_1(:,1)=1:10;

ta=table(S_1,R_spont_noNB,R_spont_NB,R_evoked_noNB,R_evoked_NB,...
    'variablenames', {'cells','Y1','Y2','Y3','Y4'});
factorNames = {'Sensory_stim','NB'};
within = table({'Sp';'Sp';'Ev';'Ev'},{'N';'Y';'N';'Y'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% fit the repeated measures model
rm = fitrm(ta,'Y1-Y4~1','WithinDesign',within);
%test for sphericity
mauchly(rm)
% run my repeated measures anova here
[ranovatbl] = ranova(rm, 'WithinModel','Sensory_stim*NB');
%multiple comparisons:
c = multcompare(rm,'NB','By','Sensory_stim');
%% rmANOVA with matlab function - R shuffled
clear all
load deconv_expvar_spont_10_cells_w5.mat
R_shuff_mat=cell2mat(R_shuff_M);
R_shuff_spont_noNB(:,1)=R_shuff_mat(:,1);
R_shuff_spont_NB(:,1)=R_shuff_mat(:,2);

clear R_shuff_mat
load deconv_expvar_evoked_10_cells_w5.mat
R_shuff_mat=cell2mat(R_shuff_M);
R_shuff_evoked_noNB(:,1)=R_shuff_mat(:,1);
R_shuff_evoked_NB(:,1)=R_shuff_mat(:,2);
R_shuff_spont(:,1)=[R_shuff_spont_noNB;R_shuff_spont_NB];
R_shuff_evoked(:,1)=[R_shuff_evoked_noNB;R_shuff_evoked_NB];
 S_1(:,1)=1:10;

ta=table(S_1,R_shuff_spont_noNB,R_shuff_spont_NB,R_shuff_evoked_noNB,R_shuff_evoked_NB,...
    'variablenames', {'cells','Y1','Y2','Y3','Y4'});
factorNames = {'Sensory_stim','NB'};
within = table({'Sp';'Sp';'Ev';'Ev'},{'N';'Y';'N';'Y'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% fit the repeated measures model
rm = fitrm(ta,'Y1-Y4~1','WithinDesign',within);
%test for sphericity
mauchly(rm)
% run my repeated measures anova here
[ranovatbl_shuff] = ranova(rm, 'WithinModel','Sensory_stim*NB');
%multiple comparisons:
c_shuff = multcompare(rm,'NB','By','Sensory_stim');
%% variables for repeated-measures ANOVA
% clear all
% load deconv_expvar_spont_10_cells_w5.mat
% R_mat=cell2mat(R_M);
% R_spont_noNB(:,1)=R_mat(:,1);
% R_spont_NB(:,1)=R_mat(:,2);
% F2_spont_noNB(:,1)=zeros(size(R_spont_noNB)); %zeroes for NB-
% F2_spont_NB(:,1)=ones(size(R_spont_NB)); %ones for NB+
% R_spont(:,1)=[R_spont_noNB;R_spont_NB];
% F1_spont(:,1)=zeros(size(R_spont));
% F2_spont(:,1)=[F2_spont_noNB;F2_spont_NB];
% 
% MSE_mat=cell2mat(MSE_M);
% MSE_spont_noNB(:,1)=MSE_mat(:,1);
% MSE_spont_NB(:,1)=MSE_mat(:,2);
% MSE_spont(:,1)=[MSE_spont_noNB;MSE_spont_NB];
% 
% clearvars -except R_spont MSE_spont F1_spont F2_spont
% load deconv_expvar_evoked_10_cells_w5.mat
% R_mat=cell2mat(R_M);
% R_evoked_noNB(:,1)=R_mat(:,1);
% R_evoked_NB(:,1)=R_mat(:,2);
% F2_evoked_noNB(:,1)=zeros(size(R_evoked_noNB)); %zeroes for NB-
% F2_evoked_NB(:,1)=ones(size(R_evoked_NB)); %ones for NB+
% R_evoked(:,1)=[R_evoked_noNB;R_evoked_NB];
% F1_evoked(:,1)=ones(size(R_evoked));
% F2_evoked(:,1)=[F2_evoked_noNB;F2_evoked_NB];
% MSE_mat=cell2mat(MSE_M);
% MSE_evoked_noNB(:,1)=MSE_mat(:,1);
% MSE_evoked_NB(:,1)=MSE_mat(:,2);
% MSE_evoked(:,1)=[MSE_evoked_noNB;MSE_evoked_NB];
%  S_1(:,1)=1:10;
% % all together
% Y_R=[R_spont;R_evoked];
% Y_MSE=[MSE_spont;MSE_evoked];
% F1=[F1_spont;F1_evoked];
% F2=[F2_spont;F2_evoked];
% S(:,1)=[S_1;S_1;S_1;S_1];
% FactNames={'spont/evoked','Nb-/NB+'};
%  filename='rmANOVA_R_actual_vars';
% save(filename,'Y_R','Y_MSE', 'F1', 'F2','S','FactNames')
% 
% stats_R = rm_anova2(Y_R,S,F1,F2,FactNames);
% stats_MSE = rm_anova2(Y_MSE,S,F1,F2,FactNames);
%  