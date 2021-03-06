%% Analyze_ChAT
%calls the functions: 

clear all
close all

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
load ChAT_Files

for fileind = 3;
    
    fname = files(fileind).extracted_name;        
    path = files(fileind).extracted_path;
    cd(path)
    load(fname) 
    
    cd  'D:\Inbal M.Sc\MATLAB\Project Brain States';
    
        for header = 1:length(files(1, fileind).headers);
            stim_type = Param.header(header).stim.type;
            %%  Airpuff/Galvano (no protocol) stimulus:
            if stim_type== 1 || stim_type== 2 || stim_type== 3 || stim_type==4 || stim_type==5 ||stim_type==6

                sf_Vm = Param.header(header).stim.sf; %[1/sec]
                dt_Vm=1/sf_Vm; %[sec]
                time_axis_Vm=[dt_Vm:dt_Vm:size(data.header(header).Vm,1)*dt_Vm];
                sf_airpuff = Param.header(header).stim.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                time_axis_airpuff = [dt_airpuff:dt_airpuff:size(data.header(header).airpuff,1)*dt_airpuff];
                sf_galvano = Param.header(header).stim.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;
                time_axis_galvano = [dt_galvano:dt_galvano:size(data.header(header).galvano,1)*dt_galvano];
     
                  data.header(header).Vm_mean = mean(data.header(header).Vm,2);
                  data.header(header).Vm_std = std(data.header(header).Vm');
                  data.header(header).Ext1_mean = mean(data.header(header).Ext1,2);
                  data.header(header).Ext1_std = std(data.header(header).Ext1');
            end
            
            %%  Airpuff/Galvano protocol stimulus:
            if stim_type== 7 || stim_type== 8
                sf_Vm = Param.header(header).stim.sf; %[1/sec]
                dt_Vm=1/sf_Vm; %[sec]
                time_axis_Vm=[dt_Vm:dt_Vm:size(data.header(header).x_value(1).Vm,1)*dt_Vm];
                sf_airpuff = Param.header(header).stim.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                time_axis_airpuff = [dt_airpuff:dt_airpuff:size(data.header(header).x_value(1).airpuff,1)*dt_airpuff];
                sf_galvano = Param.header(header).stim.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;
                time_axis_galvano = [dt_galvano:dt_galvano:size(data.header(header).x_value(1).galvano,1)*dt_galvano];
                sf_laser = Param.header(header).stim.sf_laser; %[1/sec]
                dt_laser = 1/sf_laser;
                time_axis_laser = [dt_laser:dt_laser:size(data.header(header).x_value(2).laser,1)*dt_laser];

                for x_value = 1:size(data.header(header).x_value,2)
                    
                      data.header(header).x_value(x_value).Vm_mean = mean(data.header(header).x_value(x_value).Vm,2);
                      data.header(header).x_value(x_value).Vm_std = std(data.header(header).x_value(x_value).Vm');
                      data.header(header).x_value(x_value).Ext1_mean = mean(data.header(header).x_value(x_value).Ext1,2);
                      data.header(header).x_value(x_value).Ext1_std = std(data.header(header).x_value(x_value).Ext1');
                end
            end
        end
end
    
%% file 3:
% depth 433um: header 1 - airpuff protocol. airpuff delay 2sec.
%              header 2 - airpuff protocol. airpuff delay 1.2sec. 
% depth 754um: header 3 - airpuff protocol. airpuff delay 1.2sec.


%% Depth 433 header 1

for header = 1:3;
    time_axis_airpuff = [dt_airpuff:dt_airpuff:size(data.header(header).x_value(1).airpuff,1)*dt_airpuff];    
    airpuff_vec = data.header(header).x_value(1).airpuff(:,1)./max(data.header(header).x_value(1).airpuff(:,1));
    sf_laser = Param.header(header).stim.sf_laser; %[1/sec]
    dt_laser = 1/sf_laser;
    time_axis_laser = [dt_laser:dt_laser:size(data.header(header).x_value(2).laser,1)*dt_laser];
    laser_vec = zeros(length(data.header(header).x_value(2).laser(:,1)));
    laser_vec(data.header(header).x_value(2).laser > 2)=1;
    
    LFP_mean_puff_only = data.header(header).x_value(1).Ext1_mean;
    LFP_mean_puff_light = data.header(header).x_value(2).Ext1_mean;

 
figure(header)
        set(gcf,'color','w')
        subplot(2,1,1)
            hold on 
            plot(LFP_mean_puff_only, '-b', 'linewidth', 2)
%             plot(LFP_mean_puff_only - std(LFP_mean_puff_only),'k');
            plot(LFP_mean_puff_light, '-r', 'linewidth', 1)
            hold off
        
        subplot(2,1,2)
            plot(time_axis_airpuff, airpuff_vec, '-k','linewidth', 2);

end