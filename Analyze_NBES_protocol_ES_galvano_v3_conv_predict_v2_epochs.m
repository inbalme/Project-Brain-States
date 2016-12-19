%% Analyze NBES protocol ES+galvano v3 Ldeconvs - recent version.
% This file was created on 26/04/2016, and later modified to add cross-validation method
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
cc_stat=[]; cc_spont=[]; cc_evoked=[]; cc=[]; cc_shuffled=[]; cc_shuf_sub=[];
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
cc_spont_for_xls_mean=[]; cc_evoked_for_xls_mean=[];
files_to_analyze =44; %[44,46,48,50,52,56,58,62,72,75];
save_flag=0;
    for fileind=1:length(files_to_analyze) ;         
    clearvars -except files_to_analyze fileind files save_flag cc_GMsig2_sig2 cc_Gsig2_sig2 cc_Gsig2test_sig2test R R_M R_STD R_shuff R_shuff_M R_shuff_STD
    
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
                galvano_nstim = Param.facade(6);
                galvano_freq = Param.facade(7);
% raw_data{3} = raw_data{3}./20; %dividing by the LFP gain

trials = size(raw_data{3},2);
   %% cross-correlation Vm-LFP - spontaneous activity
%         for trace_type= 3; %1:2;  %1 for spont., 2 for evoked, 3 for whole-trace-epochs spont+evoked
        %bandpass filtering to remove noise from LFP channel
        bp_filt=files(files_to_analyze(fileind)).V2_filter;
        for xx=1:3
            for trace = 1:trials;    
                    jl=raw_data{3}(:,trace,xx);
                    Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(jl,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
                    jm=Ch2_data_filt(:,trace,xx);
                    Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(jm,dt,bp_filt(1),bp_filt(2),0,0); %filtering LFP 1-160Hz
                    kl=data_no_spikes{channel}(:,trace,xx);
                    data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(kl,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
                     km=data_no_spikes_filt{channel}(:,trace,xx);
                     data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,bp_filt(1),bp_filt(2),0,0); %filtering Vm same as LFP
            end
        end

        %%
        clear  start_time duration start_time start_sample end_sample interval interval_mat x y patch_xdata patch_ydata yex ylim_data
        coeffs=[]; 
        
        % U=unique(x_value);
        % trace_type=length(U); %for trace_type=1 plot the xcorr before and after the ES according to the two start times. for trace_type=2 plot the xcorr from the same time point in the 2 different x-values.
%         switch trace_type                         
%             case 3
                 
                 kernel_dur=50; %[ms]
                 interval_dur=500; %[ms]
                 kernelsize=(kernel_dur./1000*sf{1}+1)*2;
                % kernelsize=1002; % 1002 means 50 ms back and forth
                intervalsize=interval_dur/1000*sf{1}; %5000;
                totallength=size(raw_data{3},1); %200000;
                type='ev'; % 'ev' for explained variance (R^2), 'mse' for mean-squared error
                alltimes=cell(1,totallength/intervalsize-1);

                for i=1:totallength/intervalsize-1 
                    alltimes{1,i}=[(i)*intervalsize/10000 (i+1)*intervalsize/10000];  
                end
                for x_value=1

%                 for ba=alltimes
%                 ba=ba{1};
%                 for n=1:size(alltimes,2)
%                 if isequal(ba,alltimes{n})
%                 t=n;
%                 end
%                 end
for t=1:length(alltimes);
%     ba=alltimes{t};
          interval_begin = alltimes{1,t}(1)*sf{1} ; %which is actually ba(1)
           interval_end = alltimes{1,t}(2)*sf{1}-1; %which is actually ba(2)
           interval(:,t)=interval_begin:interval_end;
      
%%
clear data_Vm data_LFP data_Vm_noDC data_LFP_noDC data_Vm_filt data_Vm_filt_noDC sig1 sig2 
%             for t=1:size(interval,2)
            %% Subtract mean from the interval (fragment of trace)
%             data_Vm{t}=data_no_spikes{1}(interval(:,t),:,x_value);
            data_LFP{t}=Ch2_data_filt(interval(:,t),:,x_value)./20;
%             data_Vm_noDC{t} = fn_Subtract_Mean(data_Vm{t});
            data_LFP_noDC{t} = fn_Subtract_Mean(data_LFP{t});
            % Subtract mean from the interval (fragment of trace) for filtered data
            data_Vm_filt{t}=data_no_spikes_filt{1}(interval(:,t),:,x_value);
            data_Vm_filt_noDC{t} = fn_Subtract_Mean(data_Vm_filt{t});
 %% Bootstrap to create shuffled Vm-LFP trace couples
  xcorr_L=length(interval(:,t))-1;
   CI_std=2;
    iterations=20;
    DC_sub=1; %put 1 if the mean trace was not subtracted from the data.
    for it=1:iterations;
[data_LFP_shuff{t}, data_Vm_filt_shuff{t}]=fn_bootstrap_shuff_trace_couples(data_LFP{t},data_Vm_filt{t}, 1, DC_sub);
 %deconvolution of the shuffled trace pairs:
         shuff_traces=size(data_LFP_shuff{t},2);
                    sig1_shuff{t}=data_LFP_shuff{t};
                    sig2_shuff{t}=data_Vm_filt_shuff{t};
%                     sig1{t}=data_LFP_noDC{t};
%                     sig2{t}=data_Vm_filt_noDC{t};
                        for trace = 1:shuff_traces;
                            G_shuff{fileind,t}(:,trace) = Ldeconvs(sig1_shuff{t}(:,trace),sig2_shuff{t}(:,trace),kernelsize,0.5,1);
                        end
                             for trace = 1:shuff_traces;
                                kernel=G_shuff{fileind,t};
                                kernel(:,trace)=[];
                                kernelc_shuff{fileind,t}(:,trace)=mean(kernel,2);
                                sig2_shuff_predict_tmp(:,trace)=conv(kernelc_shuff{fileind,t}(:,trace),sig1_shuff{t}(:,trace));
                                sig2_shuff_predict{t}(:,trace)=sig2_shuff_predict_tmp(size(kernel,1)/2:end-size(kernel,1)/2,trace);
                                cc_sig2_shuff_predict{t}(:,trace) = crosscorr(sig2_shuff_predict{t}(:,trace), sig2_shuff{t}(:,trace),xcorr_L,CI_std);                                
                                [r, p] = corr(sig2_shuff_predict{t}(:,trace),sig2_shuff{t}(:,trace));
                                pearson_r_shuff_temp(trace,1)=r;
                                pearson_p_shuff_temp(trace,1)=p;
                                R_shuff_temp(trace,1) =  pearson_r_shuff_temp(trace,1)^2;
                                 cc_sig2_shuff{fileind,t}(:,trace) = crosscorr(sig2_shuff{t}(:,trace), sig2_shuff{t}(:,trace),xcorr_L,CI_std);
                                 cc_sig1sig2_shuff(:,trace) = crosscorr(sig1_shuff{t}(:,trace), sig2_shuff{t}(:,trace),xcorr_L,CI_std);
                                 cc_sig1sig2_shuff_predict{fileind,t}(:,trace) = crosscorr(sig1_shuff{t}(:,trace), sig2_shuff_predict{t}(:,trace),xcorr_L,CI_std);
                             end
                           pearson_r_shuff{fileind,t}(:,it) = pearson_r_shuff_temp;
                           pearson_p_shuff{fileind,t}(:,it) = pearson_p_shuff_temp;
                           R_shuff{fileind,t}(:,it) = R_shuff_temp;
                           
    end
    R_shuff_M{fileind,t}=mean(R_shuff{fileind,t}(:,:),1);
    R_shuff_STD{fileind,t}=std(R_shuff{fileind,t}(:,:),1);
%             end %temporary end of t loop
%% Deconvolution of the real data          
%deconvolution of sig_1 and sig_2, convolving the kernel with sig_1 to create sig_2_predict.   
                    sig1{t}=data_LFP_noDC{t};
                    sig2{t}=data_Vm_filt_noDC{t};
                for trace = 1:trials;
                    G{fileind,t}(:,trace) = Ldeconvs(sig1{t}(:,trace),sig2{t}(:,trace),kernelsize,0.5,1);
                end
                
                %cross-validation of the kernel - averaging kernels from n-1 trials and testing on the last trial:
                 for trace = 1:trials;
                                kernel=G{fileind,t};
                                kernel(:,trace)=[];
                                kernelc{fileind,t}(:,trace)=mean(kernel,2);
                                sig2_predict_tmp(:,trace)=conv(kernelc{fileind,t}(:,trace),sig1{t}(:,trace));
                                sig2_predict{t}(:,trace)=sig2_predict_tmp(size(kernel,1)/2:end-size(kernel,1)/2,trace);
                                cc_sig2_predict{t}(:,trace) = crosscorr(sig2_predict{t}(:,trace), sig2{t}(:,trace),xcorr_L,CI_std);                                
                                [r, p] = corr(sig2_predict{t}(:,trace),sig2{t}(:,trace));
                                pearson_r_temp(trace,1)=r;
                                pearson_p_temp(trace,1)=p;
                                R_temp(trace,1) =  pearson_r_temp(trace,1)^2;
                                 cc_sig2{fileind,t}(:,trace) = crosscorr(sig2{t}(:,trace), sig2{t}(:,trace),xcorr_L,CI_std);
                                 cc_sig1sig2(:,trace) = crosscorr(sig1{t}(:,trace), sig2{t}(:,trace),xcorr_L,CI_std);
                                 cc_sig1sig2_predict{fileind,t}(:,trace) = crosscorr(sig1{t}(:,trace), sig2_predict{t}(:,trace),xcorr_L,CI_std);
                 end
                           pearson_r{fileind,t}(:,1) = pearson_r_temp;
                           pearson_p{fileind,t}(:,1) = pearson_p_temp;
                           R{fileind,t}(:,1) = R_temp;     
                           R_M{fileind,t}=mean(R{fileind,t}(:,:),1);
                           R_STD{fileind,t}=std(R{fileind,t}(:,:),1);
                end  
                end   
    end         %temporary end to file loop  
    if save_flag==1;
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP filtered 49-51Hz\Deconv ExpVar'
            filename=['deconv_expvar_epochs_',num2str(fileind),'_cells'];
            save(filename, 'files_to_analyze', 'R','R_shuff', 'R_M', 'R_STD', 'R_shuff_M', 'R_shuff_STD', 'G', 'G_shuff',...
                'kernel_dur','interval_dur','alltimes', 'kernelc', 'kernelc_shuff', 'interval','sig1','sig1_shuff','sig2','sig2_shuff','sig2_predict','sig2_shuff_predict')
    end
                     %% making figures - trace-specific kernel
                     color_table=[ 0 0 0; [136 137 138]./256; [13 49 133]./256; [33 118 209]./256];
                     fileind=1;
                 figure; 
                    subplot(1,2,1); 
                     hold on
                         for t=1:length(alltimes);
                            title('kernel data');
                            fn_shadedErrorBar([0:length(interval(:,t))-1].*dt,mean(kernelc{fileind,t},2),std(kernelc{fileind,t},0,2),{'LineWidth',1,'color', color_table(2*t-1,:)},1);
                            hold off
                         end
                     subplot(1,2,2); 
                         t=1:length(alltimes);
                            hold on
                            title('kernel shuffled');
                            fn_shadedErrorBar([0:kernelsize-1].*dt,mean(kernelc_shuff{fileind,t},2),std(kernelc{fileind,t},0,2),{'LineWidth',1,'color', color_table(2*t,:)},1);                           
                        hold off
%                         plot([0:length(interval(:,t))-1].*dt,mean(kernelc{fileind,t},2),'color',color_table(2*t-1),:)); 
%                         plot([0:length(interval(:,t))-1].*dt,mean(kernelc_shuff{fileind,t},2),'color',color_table(2*t),:)); 
%%
      trace = 2;                 
                figure
  leg_text{1,1}=['Vm ES off'];leg_text{2,1}=['Vm predicted ES off']; leg_text{3,1}=['Vm ES on']; leg_text{4,1}=['Vm predicted ES on'];        

                 for t=1:length(alltimes);
                    subplot(1,2,t); 
                        hold on                        
                            plot([0:length(interval(:,t))-1].*dt,sig2{t}(:,trace),'LineWidth',1,'color', color_table(2*t,:)); 
                            plot([0:length(interval(:,t))-1].*dt,sig2_predict{fileind,t}(:,trace),'LineWidth',1,'color', color_table(2*t-1,:));   
                             title('Vm,Vm predicted'); legend(leg_text{2*t-1},leg_text{2*t});
                 end   
                        hold off
                                               
                         figure
                         leg_text{1,1}=['Vm ES off'];leg_text{2,1}=['Vm predicted ES off']; leg_text{3,1}=['Vm ES on']; leg_text{4,1}=['Vm predicted ES on'];        
                 for t=1:length(alltimes);
                    subplot(1,2,t); 
                        hold on                        
                            plot([0:length(interval(:,t))-1].*dt,sig2_shuff{t}(:,trace),'LineWidth',1,'color', color_table(2*t,:)); 
                            plot([0:length(interval(:,t))-1].*dt,sig2_shuff_predict{fileind,t}(:,trace),'LineWidth',1,'color', color_table(2*t-1,:));  
                             legend(leg_text{2*t-1},leg_text{2*t}); title('Vm,Vm predicted - shuffled'); 
                 end   
                        hold off                  
                       
%                         legend('Vm ES off', 'Vm predicted ES off', 'Vm ES on','Vm predicted ES on')
                        
                     figure
                leg_text{1,1}=['ccVmVm ES off'];leg_text{2,1}=['ccVmVm" ES off']; leg_text{3,1}=['ccVmVm ES on']; leg_text{4,1}=['ccVmVm" ES on'];                            
                for t=1:length(alltimes);
                            subplot(1,2,t); 
                        hold on
                        plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig2{fileind,t}(:,trace),'LineWidth',1,'color', color_table(2*t-1,:));
                        plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig2_predict{fileind,t}(:,trace),'LineWidth',1,'color', color_table(2*t,:));
                        set(gca,'xlim',[-0.5,0.5]);
                        title('ccVm,ccVm predicted')
                         legend(leg_text{2*t-1},leg_text{2*t});
                         hold off
                        end
                        
                         figure
                         leg_text{1,1}=['ccVmVm ES off'];leg_text{2,1}=['ccVmVm" ES off']; leg_text{3,1}=['ccVmVm ES on']; leg_text{4,1}=['ccVmVm" ES on'];        
                        for t=1:length(alltimes);
                            subplot(1,2,t); 
                        hold on
                        plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig2_shuff{fileind,t}(:,trace),'LineWidth',1,'color', color_table(2*t-1,:));
                        plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig2_shuff_predict{fileind,t}(:,trace),'LineWidth',1,'color', color_table(2*t,:));
                        set(gca,'xlim',[-0.5,0.5]);
                        title('ccVm,ccVm predicted - shuffled')
                         legend(leg_text{2*t-1},leg_text{2*t});
                         hold off
                        end
                        
%                     subplot(1,4,4); hold on
%                         plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig1sig2{fileind,t}(:,trace),'k'); 
%                         plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig1sig2_predict{fileind,t}(:,trace),'r'); 
%                         plot([-1.*(length(interval(:,t)))+1:(length(interval(:,t)))-1]*dt,cc_sig1sig2M_predict{fileind,t}(:,trace),':r'); 
%                          set(gca,'xlim',[-0.5,0.5]);
%                         hold off
%                         title('ccLFP-Vm, ccLFP-Vm predicted')                    
                     pause
                    close(gcf);          
                    


%         end
    