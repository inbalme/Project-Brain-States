%% NBES Rin

clear all
close all

 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
save_flag= 1;
print_flag=1;
type_strng = {'ongoing', 'evoked'};
files_to_analyze =[62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,44,46,48,50,52,56,58,62,72,75]; 
% for fileind=1;
    for fileind=1:length(files_to_analyze) ;
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    clear data_no_spike_no_DC
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
    color_table=[0 0 0;color_table(1:6,:)];
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
                
%% low-pass filtering below 300Hz                
                for xx=2:3
    for trace= 1:size(data_no_spikes{channel},2)    
            jj=data_no_spikes{channel}(:,trace,xx);
            data_no_spikes{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jj,dt,-1,300,0,0); 
    end
end
%% 
x_value=2:3;
curr_inj_interval_before=[0.125*sf{1}:0.175*sf{1}];
curr_inj_interval_after=[5.225*sf{1}:5.275*sf{1}];
% curr_inj_interval_before=[0.05*sf{1}:0.225*sf{1}];
% curr_inj_interval_after=[5.150*sf{1}:5.325*sf{1}];
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
for trace=1:7;
Iinj=Param.facade(21)/2.5;
curr_inj_del=[Param.facade(20), Param.facade(20)+Param.facade(22)+Param.facade(23)];
curr_inj_duration=Param.facade(22);
voltages1=raw_data{1}(curr_inj_del(1)*sf{1}:(curr_inj_del(1)+curr_inj_duration)*sf{1},trace,3)';
voltages2=raw_data{1}(curr_inj_del(2)*sf{1}:(curr_inj_del(2)+curr_inj_duration)*sf{1},trace,3)';
before_ES(trace,:)  = CalculateCapacitResist(Iinj, voltages1, sf{1}/1000, 1,Vm_x3_curr_inj(1,1));
after_ES(trace,:)  = CalculateCapacitResist(Iinj, voltages2, sf{1}/1000, 1,Vm_x3_curr_inj(1,2));
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
    