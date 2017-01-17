%% NBES Rin

clear all
close all

global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data Ch2_data stim2_X stim1_X 
 global exp_type
 channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
exp_type=1; %1-NBES, 2-ChAT
% trace_type_input=[3,2]; %1:3
save_flag= 0;
print_flag=0;
norm_flag=0;
BP50HzLFP_flag=0; %removing 50Hz noise from LFP signal
BP50HzVm_flag=0; %removing 50Hz noise from Vm signal
BPLFP_flag=0; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=[0.1,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=0; %filtering LFP and Vm same as LFP was filtered in the multiclamp
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
        files_to_analyze =84; %[62,72,75,82,84]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};       

    case 2
        files_to_analyze =87; %[74,76,80,82,84,87];
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light On', 'Light Off'};
end
        
% for fileind=1;
    for fileind=1:length(files_to_analyze) ;
    close all
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
   
    if ~isempty(current_data_filt)
     current_data=current_data_filt;
    end

%% 
x_value=2:3;
        curr_inj_1=Param.facade(20); %in sec
        curr_inj_dur=Param.facade(22); %in sec
        curr_inj_2=Param.facade(20)+curr_inj_dur+Param.facade(23); %in sec        
        curr_inj_1_end=curr_inj_1+curr_inj_dur;
        curr_inj_2_end=curr_inj_2+curr_inj_dur;
        

curr_inj_interval_before=[round((curr_inj_1+0.025)*sf{1}):round((curr_inj_1_end-0.025)*sf{1})];
curr_inj_interval_after=[round((curr_inj_2+0.025)*sf{1}):round((curr_inj_2_end-0.025)*sf{1})];
Vm_x2_curr_inj(1,1)=mean(mean(data_no_spikes{1}(curr_inj_interval_before,:,2)));
Vm_x2_curr_inj(1,2)=mean(mean(data_no_spikes{1}(curr_inj_interval_after,:,2)));
Rin_x2_curr_inj(1,1)=Vm_x2_curr_inj(1,1)./mean(mean(raw_data{2}(curr_inj_interval_before,:,2)));
Rin_x2_curr_inj(1,2)=Vm_x2_curr_inj(1,2)./mean(mean(raw_data{2}(curr_inj_interval_after,:,2)));
Vm_x3_curr_inj(1,1)=mean(mean(data_no_spikes{1}(curr_inj_interval_before,:,3)));
Vm_x3_curr_inj(1,2)=mean(mean(data_no_spikes{1}(curr_inj_interval_after,:,3)));
Rin_x3_curr_inj(1,1)=Vm_x3_curr_inj(1,1)./mean(mean(raw_data{2}(curr_inj_interval_before,:,3)));
Rin_x3_curr_inj(1,2)=Vm_x3_curr_inj(1,2)./mean(mean(raw_data{2}(curr_inj_interval_after,:,3)));
%here I didn't subtract the value of Vm at beginning of current inj. from
%the value it reached during the curr. inj.
%Need also to get the Rin from each trace and test for each cell if it is
%different with ttest.

curr_inj_interval_before=[round((curr_inj_1)*sf{1}):round((curr_inj_1_end)*sf{1})];
curr_inj_interval_after=[round((curr_inj_2)*sf{1}):round((curr_inj_2_end)*sf{1})];
for trace=1:7;
Iinj=Param.facade(21)/2.5;
curr_inj_del=[Param.facade(20), Param.facade(20)+Param.facade(22)+Param.facade(23)];
curr_inj_duration=Param.facade(22);
voltages1=raw_data{1}(curr_inj_del(1)*sf{1}:(curr_inj_del(1)+curr_inj_duration)*sf{1},trace,3)';
voltages2=raw_data{1}(curr_inj_del(2)*sf{1}:(curr_inj_del(2)+curr_inj_duration)*sf{1},trace,3)';
before_ES(trace,:)  = fn_CalculateCapacitResist(Iinj, voltages1, sf{1}/1000, 1,Vm_x3_curr_inj(1,1));
after_ES(trace,:)  = fn_CalculateCapacitResist(Iinj, voltages2, sf{1}/1000, 1,Vm_x3_curr_inj(1,2));
    end
cell_props(fileind).fname=files_to_analyze(fileind);
cell_props(fileind).trace=trace;
cell_props(fileind).R_electrode=[before_ES(:,1), after_ES(:,1)];
cell_props(fileind).C_electrode=[before_ES(:,2), after_ES(:,2)];
cell_props(fileind).R_cell=[before_ES(:,3), after_ES(:,3)];
cell_props(fileind).R_cell_m=nanmean(cell_props(fileind).R_cell);
cell_props(fileind).R_cell_std=nanstd(cell_props(fileind).R_cell);
cell_props(fileind).C_cell=[before_ES(:,4), after_ES(:,4)];
cell_props(fileind).Error1=[before_ES(:,5), after_ES(:,5)];
cell_props(fileind).Error2=[before_ES(:,6), after_ES(:,6)];
cell_props(fileind).R_total=[before_ES(:,1)+before_ES(:,3), after_ES(:,1)+after_ES(:,3)];

    end
    