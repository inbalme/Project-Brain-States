% data_preprocessing
%this script sets the preprocessing of the data (filtering to remove noise,
%etc...) prior to the analysis

                sf{1} = Param.sf_Vm;
                sf{2} = Param.sf_I1;
                sf{3} = Param.sf_V2;
                sf{4} = Param.sf_I2;
                dt=1/sf{channel};
                             
                sf_airpuff = Param.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                sf_galvano = Param.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;   
    
    current_data_filt=[]; Ch2_data_filt=[];

%% bandpass filtering to remove 50Hz noise from LFP and Vm.
if BP50HzLFP_flag==1;
    if isempty(Ch2_data_filt)
        tmpMat=Ch2_data;
    else
        tmpMat=Ch2_data_filt;
    end
        for xx=1:3
            for trace= 1:size(tmpMat,2)    
                    jl=tmpMat(:,trace,xx);
                    Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(jl,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
            end
        end
end
tmpMat=[];

if BP50HzVm_flag==1;
    if isempty(current_data_filt)
        tmpMat=current_data;
    else
        tmpMat=current_data_filt;
    end
        for xx=1:3
            for trace= 1:size(current_data,2)    
                jm=tmpMat(:,trace,xx);
                current_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(jm,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
            end
        end
end
 tmpMat=[];
bp_filt_LFP=files(files_to_analyze(fileind)).V2_filter;
if BPLFP_flag==1;  %filtering both LFP and VM same as LFP was filtered during the experiment via multiclamp
    if ~isempty(bp_manual_LFP)
        bp_filt_LFP=bp_manual_LFP;
    end
        if isempty(Ch2_data_filt)
        tmpMat=Ch2_data;
    else
        tmpMat=Ch2_data_filt;
        end
       
        for xx=1:3
            for trace= 1:size(tmpMat,2)   
                kl=tmpMat(:,trace,xx);
                if bp_filt_LFP(1)==0 %Low-Pass Filtering
                    Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(kl,dt,-1,bp_filt_LFP(2),0,0); 
                else if bp_filt_LFP(2)==0 %High-Pass Filtering
                        Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(kl,dt,bp_filt_LFP(1),-1,0,0); 
                    else
                        Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(kl,dt,bp_filt_LFP(1),bp_filt_LFP(2),0,0); 
                    end
                end
            end
        end
end
 tmpMat=[];
 bp_filt_Vm=files(files_to_analyze(fileind)).V2_filter; 
if BPVm_flag==1;  %filtering both LFP and VM same as LFP was filtered during the experiment via multiclamp
    if ~isempty(bp_manual_Vm)
        bp_filt_Vm=bp_manual_Vm;
    end      
        if isempty(current_data_filt)
            tmpMat=current_data;
        else
            tmpMat=current_data_filt;
        end
    for xx=1:3
        for trace= 1:size(tmpMat,2)   
               km=tmpMat(:,trace,xx);
                     if bp_filt_Vm(1)==0 %Low-Pass Filtering out above bp_filt_Vm(2)
                         current_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,-1,bp_filt_Vm(2),0,0); 
                        else if bp_filt_Vm(2)==0 %High-Pass Filtering out below bp_filt_Vm(1)
                                 current_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,bp_filt_Vm(1),-1,0,0); 
                            else %band-pass filtering out below bp_filt_Vm(1) and above bp_filt_Vm(2)
                                current_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,bp_filt_Vm(1),bp_filt_Vm(2),0,0); %filtering Vm same as LFP
                            end
                     end
        end
    end
end
