%% Analyze NBES protocol ES+galvano v3 Ldeconvs - recent version.
% This file was created on 26/04/2016 
%This file is used for the analysis of files created with extract_NBES_Data_v2
%old version - before introducing the cross validation method

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
cc_stat=[]; cc_spont=[]; cc_evoked=[]; cc=[]; cc_shuffled=[]; cc_shuf_sub=[];
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
cc_spont_for_xls_mean=[]; cc_evoked_for_xls_mean=[];
files_to_analyze =44; %[44,46,48,50,52,56,58,62,72,75];
save_flag=0;
    for fileind=1:length(files_to_analyze) ;         
    clearvars -except files_to_analyze fileind files save_flag cc_GMsig2_sig2 cc_Gsig2_sig2 cc_Gsig2test_sig2test
    
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
        for trace_type= 1; %1:2;  %1 for spont., 2 for evoked
        if trace_type==1 
            x_value=[1,1];
        end
        if trace_type==2 
        x_value=[2:3]; %2:3; %[1,1];
        end

        %bandpass filtering to remove noise from LFP channel
        bp_filt=files(files_to_analyze(fileind)).V2_filter;
        for xx=1:3
            for trace = 1:trials;    
                    jl=raw_data{3}(:,trace,xx);
                    Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(jl,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
                    jm=data_no_spikes{channel}(:,trace,xx);
                    data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jm,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
                     km=data_no_spikes_filt{channel}(:,trace,xx);
                     data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,bp_filt(1),bp_filt(2),0,0); %filtering Vm same as LFP
            end
        end

        %%
        clear  start_time duration start_time start_sample end_sample interval interval_mat x y patch_xdata patch_ydata yex ylim_data
        coeffs=[]; 
        % U=unique(x_value);
        % trace_type=length(U); %for trace_type=1 plot the xcorr before and after the ES according to the two start times. for trace_type=2 plot the xcorr from the same time point in the 2 different x-values.
        switch trace_type
            case 1
                 start_time = [0.4,5.6]; %[sec] %[0,5]
                 duration = 2.5; %[sec] 
                    for t=1:length(start_time);
                     start_sample(:,t) = ceil(start_time(t).*sf{1});
                        if start_time(t)==0
                            start_sample(:,t) = 1;
                        end
                      end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
                      interval(:,t) = start_sample(:,t):end_sample(:,t);
                    end 
            case 2
                 start_sample =stim2_X{x_value(1)}(1,1);  %start with sensory stim
        %          duration = 1; %[sec]
                galvano_nstim = Param.facade(6);
                galvano_freq = Param.facade(7);
                duration = galvano_nstim./galvano_freq+0.05;
                end_sample = start_sample+duration.*sf{1}-1;
                interval(:,1) = round(start_sample:end_sample);
                interval(:,2) = interval(:,1);
        end
clear data_Vm data_LFP data_Vm_noDC data_LFP_noDC data_Vm_filt data_Vm_filt_noDC sig1 sig2 
            for t=1:2
            %% Subtract mean from the interval (fragment of trace)
            data_Vm{t}=data_no_spikes{1}(interval(:,t),:,x_value(t));
            data_LFP{t}=Ch2_data_filt(interval(:,t),:,x_value(t));
            data_Vm_noDC{t} = fn_Subtract_Mean(data_Vm{t});
            data_LFP_noDC{t} = fn_Subtract_Mean(data_LFP{t});
            % Subtract mean from the interval (fragment of trace) for filtered data
            data_Vm_filt{t}=data_no_spikes_filt{1}(interval(:,t),:,x_value(t));
            data_Vm_filt_noDC{t} = fn_Subtract_Mean(data_Vm_filt{t});
            %% Ldeconvs
            %deconvolution of sig_1 and sig_2, convolving the kernel with sig_1 to create sig_2_predict.
            xcorr_L=length(interval(:,t))-1;
            CI_std=2;
            
                for trace = 1:trials;
                    sig1=data_LFP_noDC{t}(:,trace);
                    sig2=data_Vm_filt_noDC{t}(:,trace);
                    G{fileind,t}(:,trace) = Ldeconvs(sig1,sig2);
                    sig2_predict{fileind,t}(:,trace) = conv(G{fileind,t}(:,trace),sig1,'same');
                    cc_sig2_predict{fileind,t}(:,trace) = crosscorr(sig2_predict{fileind,t}(:,trace), sig2,xcorr_L,CI_std);
                    cc_max{fileind}(trace,t) = max(cc_sig2_predict{fileind,t}(:,trace));
                     cc_sig2{fileind,t}(:,trace) = crosscorr(sig2, sig2,xcorr_L,CI_std);
                     cc_sig1sig2{fileind,t}(:,trace) = crosscorr(sig1, sig2,xcorr_L,CI_std);
                     cc_sig1sig2_predict{fileind,t}(:,trace) = crosscorr(sig1, sig2_predict{fileind,t}(:,trace),xcorr_L,CI_std);
                    figure, subplot(1,4,1); plot([0:length(interval(:,t))-1].*dt,G{fileind,t}(:,trace)); title('kernel');
                    subplot(1,4,2); hold on                    
                    plot([0:length(interval(:,t))-1].*dt,sig2,'k'); 
                    plot([0:length(interval(:,t))-1].*dt,sig2_predict{fileind,t}(:,trace),'r');
                    hold off
                    title('Vm,Vm predicted')
                    subplot(1,4,3); 
%                     corr_x_axis=(size(cc_Vm{fileind,t},1)/2-0.5./dt:size(cc_Vm{fileind,t},1)/2+0.5./dt-1).*dt;
                    hold on
                    plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig2{fileind,t}(:,trace),'k');
                    plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig2_predict{fileind,t}(:,trace),'r');
                    set(gca,'xlim',[-0.5,0.5]);
                    hold off
                    title('ccVm,ccVm predicted')
                    subplot(1,4,4); hold on
                    plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig1sig2{fileind,t}(:,trace),'k'); 
                    plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig1sig2_predict{fileind,t}(:,trace),'r'); 
                     set(gca,'xlim',[-0.5,0.5]);
                    hold off
                    title('ccLFP-Vm, ccLFP-Vm predicted')
                    
                     pause
                    close(gcf);
                    
                    %training on the first half of each trace and test on the second half
                    sig1train=data_LFP_noDC{t}(1:end./2,trace);
                    sig2train=data_Vm_filt_noDC{t}(1:end./2,trace);
                    G_train{fileind,t}(:,trace) = Ldeconvs(sig1train,sig2train);
                end
                
                G_mean{fileind,t}=mean(G{fileind,t}(:,:),2);
                xcorr_Ltest=length(sig1train)-1;
                %Repeat the Ldeconvs with the mean kernel
                
                 for trace = 1:trials;  
                    sig1=data_LFP_noDC{t}(:,trace);
                    sig2=data_Vm_filt_noDC{t}(:,trace);
                    sig2_predict{fileind,t}(:,trace) = conv(G_mean{fileind,t}(:),sig1,'same');
                    cc_sig2_predict{fileind,t}(:,trace) = crosscorr(sig2_predict{fileind,t}(:,trace), sig2,xcorr_L,CI_std);
                    cc_max_GM{fileind}(trace,t) = max(cc_sig2_predict{fileind,t}(:,trace));
                     cc_sig2{fileind,t}(:,trace) = crosscorr(sig2, sig2,xcorr_L,CI_std);
                     cc_sig1sig2{fileind,t}(:,trace) = crosscorr(sig1, sig2,xcorr_L,CI_std);
                     cc_sig1sig2_predict{fileind,t}(:,trace) = crosscorr(sig1, sig2_predict{fileind,t}(:,trace),xcorr_L,CI_std);
                    figure; title('convolution with mean kernal'); subplot(1,4,1); plot([0:length(interval(:,t))-1].*dt,G_mean{fileind,t}); title('mean  kernel');
                    subplot(1,4,2); hold on                    
                    plot([0:length(interval(:,t))-1].*dt,sig2,'k'); 
                    plot([0:length(interval(:,t))-1].*dt,sig2_predict{fileind,t}(:,trace),'r');
                    hold off
                    title('Vm,Vm predicted')
                    subplot(1,4,3); 
%                     corr_x_axis=(size(cc_Vm{fileind,t},1)/2-0.5./dt:size(cc_Vm{fileind,t},1)/2+0.5./dt-1).*dt;
                    hold on
                    plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig2{fileind,t}(:,trace),'k');
                    plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig2_predict{fileind,t}(:,trace),'r');
                    set(gca,'xlim',[-0.5,0.5]);
                    hold off
                    title('ccVm,ccVm predicted')
                    subplot(1,4,4); hold on
                    plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig1sig2{fileind,t}(:,trace),'k'); 
                    plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig1sig2_predict{fileind,t}(:,trace),'r'); 
                     set(gca,'xlim',[-0.5,0.5]);
                    hold off
                    title('ccLFP-Vm, ccLFP-Vm predicted')
                    
                     pause
                    close(gcf);
                    
                    %test on the second half after training on the first half of each trace
                     sig1test=data_LFP_noDC{t}(end./2+1:end,trace);
                    sig2test=data_Vm_filt_noDC{t}(end./2+1:end,trace);
                    sig2test_predict{fileind,t}(:,trace) = conv(G_train{fileind,t}(:,trace),sig1test,'same');
                    cc_sig2test_predict{fileind,t}(:,trace) = crosscorr(sig2test_predict{fileind,t}(:,trace), sig2test,xcorr_Ltest,CI_std);
                    cc_max_test{fileind}(trace,t) = max(cc_sig2test_predict{fileind,t}(:,trace));
                     cc_sig2test{fileind,t}(:,trace) = crosscorr(sig2test, sig2test,xcorr_Ltest,CI_std);
                     cc_sig1sig2test{fileind,t}(:,trace) = crosscorr(sig1test, sig2test,xcorr_Ltest,CI_std);
                     cc_sig1sig2test_predict{fileind,t}(:,trace) = crosscorr(sig1test, sig2test_predict{fileind,t}(:,trace),xcorr_Ltest,CI_std);
                    figure, subplot(1,4,1); plot([0:length(sig1test)-1].*dt,G_train{fileind,t}(:,trace)); title('training kernel');
                    subplot(1,4,2); hold on                    
                    plot([0:length(sig1test)-1].*dt,sig2test,'k'); 
                    plot([0:length(sig1test)-1].*dt,sig2test_predict{fileind,t}(:,trace),'r');
                    hold off
                    title('Vm-test,Vm-test predicted')
                    subplot(1,4,3); 
%                     corr_x_axis=(size(cc_Vm{fileind,t},1)/2-0.5./dt:size(cc_Vm{fileind,t},1)/2+0.5./dt-1).*dt;
                    hold on
                    plot([-1.*(length(sig1test))+1:(length(sig1test))-1]*dt,cc_sig2test{fileind,t}(:,trace),'k');
                    plot([-1.*(length(sig1test))+1:(length(sig1test))-1]*dt,cc_sig2test_predict{fileind,t}(:,trace),'r');
                    set(gca,'xlim',[-0.5,0.5]);
                    hold off
                    title('ccVm-test,ccVm-test predicted')
                    subplot(1,4,4); hold on
                    plot([-1.*(length(sig1test))+1:(length(sig1test))-1]*dt,cc_sig1sig2test{fileind,t}(:,trace),'k'); 
                    plot([-1.*(length(sig1test))+1:(length(sig1test))-1]*dt,cc_sig1sig2test_predict{fileind,t}(:,trace),'r'); 
                     set(gca,'xlim',[-0.5,0.5]);
                    hold off
                    title('ccLFP-Vm-test, ccLFP-Vm-test predicted')
                    
                     pause
                    close(gcf); 
                 end                
            end
            cc_Gsig2_sig2{fileind}=mean(cc_max{fileind},1);
            cc_GMsig2_sig2{fileind}=mean(cc_max_GM{fileind},1);
            cc_Gsig2test_sig2test{fileind}=mean(cc_max_test{fileind},1);
        end
    end