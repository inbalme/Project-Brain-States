%Analyze_VmLFPcoherence
% This file was created on 10/2/2017
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
exp_type=1; %1-NBES, 2-ChAT
trace_type_input=[3,2]; %[3,2] for exp_type=1; %for exp_type=2 or 3 use [1,2]
analyze_time_before_train=0;
analyze_train_only_flag=1;
add_to_plot=0.15; %seconds from each side of the trace
save_flag=0;
print_flag=0;
norm_flag=0;
clamp_flag=[]; %[]; %3; %clamp_flag=1 for hyperpolarization traces, clamp_flag=2 for depolarization traces and clamp_flag=3 for no current traces (only clamp to resting Vm)
BP50HzLFP_flag=1; %removing 50Hz noise from LFP signal
BP50HzVm_flag=1; %removing 50Hz noise from Vm signal
BPLFP_flag=1; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=[0.1,200];%[0.1,200]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=0; %filtering LFP and Vm same as LFP was filtered in the multiclamp
bp_manual_Vm=[30,50];%[0,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
    
switch exp_type
    case 1
        files_to_analyze =82; %[44,46,48,52,56,58,62,72,75,82,84]; %[44,46,48,50,52,56,58,62,72,75]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'}; y_ax_label={'Vm'}; y_ax_units={'mV'};   
        legend_string_shuff={'NB+ shuffled', 'NB- shuffled'};       

    case 2
        files_to_analyze =[76,77,80,82,84,87];
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light On', 'Light Off'}; y_ax_label={'Vm'}; y_ax_units={'mV'};   
        legend_string_shuff={'Light On shuffled', 'Light Off shuffled'};
        
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
        legend_string legend_string_shuff analyze_time_before_train analyze_train_only_flag clamp_flag y_ax_label y_ax_units add_to_plot
   
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
   if ~isempty(current_data_filt)
     current_data=current_data_filt;
 end
%% set the path for saving figures and variables
if BP50HzLFP_flag==1 && BP50HzVm_flag==1 && BPVm_flag==1 && BPLFP_flag==1
    path_output=['LFP_50Hz+BP Vm_ 50Hz+BP\BP',num2str(bp_manual_Vm(1)),'-',num2str(bp_manual_Vm(2))];
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

clear color_table    
    whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;

switch exp_type
    case 1
        color_table=[0 0 0; [30,75,14]/256; [136 137 138]/256; [112,172,90]/256; [216 22 22]/256; [255 153 153]/256];
        path_output=['D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP coherence\',path_output];   
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        cd(path_output)       
    case 2
        color_table=[0 0 0; [30,75,14]/256; [136 137 138]/256; [112,172,90]/256; [0 0 204]/256; [255 153 153]/256];
        path_output=['D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm-LFP coherence\',path_output]; 
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
%         interval=interval(1:2000,:);
%%
 fVm=[]; fLFP=[]; f=[]; n=0; %n is length of segment in case we don't want the default length(y)
 fVmConj=[]; fLFPConj=[]; Cohxy=[]; coh=[]; fc=[];
 VmLFPCrossSpectrum=[]; VmVmAutoSpectrum=[]; LFPVmCrossSpectrum=[]; LFPLFPAutoSpectrum=[];
 
for t=1:2
%% Subtract mean from the interval (fragment of trace)
data_LFP{t}=Ch2_data_filt(interval(:,t),:,x_value(t));
data_LFP_noDC{t} = fn_Subtract_Mean(data_LFP{t});
data_Vm{t}=current_data(interval(:,t),:,x_value(t));   
data_Vm_noDC{t} = fn_Subtract_Mean(data_Vm{t});
data_Vm_filt{t}=current_data(interval(:,t),:,x_value(t));
data_Vm_filt_noDC{t} = fn_Subtract_Mean(data_Vm_filt{t});
% for plots:
     interval_plot(:,t)=[interval(1,t)-add_to_plot*sf{1}:interval(end,t)+add_to_plot*sf{1}]; 
    data_Vm_plot{t}=current_data(interval_plot(:,t),:,x_value(t));
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

%Vm fft: 
DC_Vm=[]; spec_mat_Vm_noDC=[]; spec_mat_Vm = []; 
spec_mat_Vm = data_Vm_filt_noDC{t};
DC_Vm= mean(spec_mat_Vm,1);
spec_mat_Vm_noDC=bsxfun(@minus,spec_mat_Vm,DC_Vm);
% l=size(spec_mat,1);
% l=l/2-1;
% DC=wextend('addrow', 'per', DC, l);
% spec_mat_noDC = spec_mat-DC;
 L = size(spec_mat_Vm_noDC,1);
if  n~=0;
    NFFT = 2^nextpow2(n); % Next power of 2 from length of n;
else
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
end
% Y(:,:,t) = fft(spec_mat_noDC,NFFT)/L;
fVm(:,:,t) = fft(spec_mat_Vm_noDC,NFFT);
f = sf{1}/2*linspace(0,1,NFFT/2+1);

%LFP fft: 
DC_LFP=[]; spec_mat_LFP_noDC=[]; spec_mat_LFP = []; 
spec_mat_LFP = data_LFP_noDC{t};
DC_LFP= mean(spec_mat_LFP,1);
spec_mat_LFP_noDC=bsxfun(@minus,spec_mat_LFP,DC_LFP);
% l=size(spec_mat,1);
% l=l/2-1;
% DC=wextend('addrow', 'per', DC, l);
% spec_mat_noDC = spec_mat-DC;
 L = size(spec_mat_LFP_noDC,1);
if  n~=0;
    NFFT = 2^nextpow2(n); % Next power of 2 from length of n;
else
    NFFT = 2^nextpow2(L); % Next power of 2 from length of y
end
fLFP(:,:,t) = fft(spec_mat_LFP_noDC,NFFT);

 fVmConj(:,:,t) = conj(fVm(:,:,t)); %complex conjugate of the average LFP spectrum %average?
  fLFPConj(:,:,t) = conj(fLFP(:,:,t)); %complex conjugate of the average LFP spectrum
  
   VmLFPCrossSpectrum(:,:,t) = fVmConj(:,:,t).* fLFP(:,:,t); %cross spectrum Vm*LFP
    VmVmAutoSpectrum(:,:,t)  = fVmConj(:,:,t).* fVm(:,:,t); %auto spectrum Vm*Vm
    
    LFPVmCrossSpectrum(:,:,t) = fLFPConj(:,:,t).* fVm(:,:,t); %cross spectrum R*S
    LFPLFPAutoSpectrum(:,:,t)  = fLFPConj(:,:,t).* fLFP(:,:,t); %auto spectrum R*R
%     coh{t} = (VmLFPCrossSpectrum(:,:,t) .* LFPVmCrossSpectrum(:,:,t)) ./ (VmVmAutoSpectrum(:,:,t) .* LFPLFPAutoSpectrum(:,:,t));
    coh{t} = (abs(VmLFPCrossSpectrum(:,:,t)).^2) ./ (VmVmAutoSpectrum(:,:,t) .* LFPLFPAutoSpectrum(:,:,t));
    [Cohxy{t}, fc] = mscohere(data_Vm_filt_noDC{t},data_LFP_noDC{t},[],[],NFFT,sf{1}); %[Cxy,F] = mscohere(x,y,window,noverlap,nfft,fs)
end
figure
plot(fc,mean(Cohxy{1},2),'r','DisplayName','NB-');
hold on;
plot(fc,mean(Cohxy{2},2),'k','DisplayName','NB+');
ylim([0 1])
xlim([0 50])
legend('show')
title ('single-sided coherence between Vm and LFP')
xlabel('frequency [Hz]')
ylabel('\gamma^2')
end
    end