%% Analyze_ChAT
%calls the functions: 

clear all
close all

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
load ChAT_Files

for fileind = 2;
    
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

                for x_value = 1:size(data.header(header).x_value,2)
                    
                      data.header(header).x_value(x_value).Vm_mean = mean(data.header(header).x_value(x_value).Vm,2);
                      data.header(header).x_value(x_value).Vm_std = std(data.header(header).x_value(x_value).Vm');
                      data.header(header).x_value(x_value).Ext1_mean = mean(data.header(header).x_value(x_value).Ext1,2);
                      data.header(header).x_value(x_value).Ext1_std = std(data.header(header).x_value(x_value).Ext1');
                end
            end
        end
end
    
%% file 1:
% depth 200um: header 4 - airpuff protocol
%depth 300um: header 5 - airpuff protocol
%depth 550um: header 7,9,12,14 - airpuff alone
%             header 8,10,11,13 - airpuff + light
%             header 15,16 - airpuff protocol
%depth 800um: header 17,18 - airpuff protocol
%             header 19,21,23,25 - airpuff alone
%             header 20,22,24,26 - airpuff + light


airpuff_vec = data.header(header).airpuff(:,1)./max(data.header(header).airpuff(:,1));

%         figure
%         set(gcf,'color','w')
%         subplot(2,1,1)
%             hold on
%             errorbar(data.header(7).Ext1_mean,data.header(7).Ext1_std, 'b')
%             errorbar(data.header(9).Ext1_mean,data.header(9).Ext1_std, 'b')
%             errorbar(data.header(12).Ext1_mean,data.header(12).Ext1_std, 'b')
%             errorbar(data.header(14).Ext1_mean,data.header(14).Ext1_std, 'b')
%             errorbar(data.header(8).Ext1_mean,data.header(8).Ext1_std, 'r')
%             errorbar(data.header(10).Ext1_mean,data.header(10).Ext1_std, 'r')
%             errorbar(data.header(11).Ext1_mean,data.header(11).Ext1_std, 'r')
%             errorbar(data.header(13).Ext1_mean,data.header(13).Ext1_std, 'r')    
%             hold off
%             
%         subplot(2,1,2)
%             plot(time_axis_airpuff, airpuff_vec, '-k','linewidth', 2);
            
        figure
        set(gcf,'color','w')
        subplot(2,1,1)
            hold on 
            plot(data.header(7).Ext1_mean, '-b', 'linewidth', 2)
%             plot(data.header(7).Ext1_mean - data.header(7).Ext1_std','k');
            plot(data.header(8).Ext1_mean, '-r', 'linewidth', 1)
            hold off
        
        subplot(2,1,2)
            plot(time_axis_airpuff, airpuff_vec, '-k','linewidth', 2);
%% depth 550um
airpuff_vec = data.header(header).airpuff(:,1)./max(data.header(header).airpuff(:,1));

pulled_puff_alone = [data.header(7).Ext1 data.header(9).Ext1 data.header(11).Ext1 data.header(13).Ext1];   
pulled_puff_alone_mean = mean(pulled_puff_alone, 2);
pulled_puff_alone_std = (std(pulled_puff_alone')');
pulled_puff_light = [data.header(8).Ext1 data.header(10).Ext1 data.header(12).Ext1 data.header(14).Ext1];   
pulled_puff_light_mean = mean(pulled_puff_light, 2);
pulled_puff_light_std = (std(pulled_puff_light')');
 
figure
        set(gcf,'color','w')
        subplot(2,1,1)
            hold on 
            plot(pulled_puff_alone_mean, '-b', 'linewidth', 2)
%             plot(pulled_puff_alone_mean - pulled_puff_alone_std,'k');
            plot(pulled_puff_light_mean, '-r', 'linewidth', 2)
            hold off
        
        subplot(2,1,2)
            plot(time_axis_airpuff, airpuff_vec, '-k','linewidth', 2);
% light on 200-700 ms            
%% Depth 800
airpuff_vec = data.header(header).airpuff(:,1)./max(data.header(header).airpuff(:,1));

pulled_puff_alone_d2 = [data.header(19).Ext1 data.header(21).Ext1 data.header(23).Ext1 data.header(25).Ext1];   
pulled_puff_alone_mean_d2 = mean(pulled_puff_alone_d2, 2);
pulled_puff_alone_std_d2 = (std(pulled_puff_alone_d2')');
pulled_puff_light_d2 = [data.header(20).Ext1 data.header(22).Ext1 data.header(24).Ext1 data.header(26).Ext1];   
pulled_puff_light_mean_d2 = mean(pulled_puff_light_d2, 2);
pulled_puff_light_std_d2 = (std(pulled_puff_light_d2')');
 
figure
        set(gcf,'color','w')
        subplot(2,1,1)
            hold on 
            plot(pulled_puff_alone_mean_d2, '-b', 'linewidth', 2)
%             plot(pulled_puff_alone_mean_d2 - pulled_puff_alone_std_d2,'k');
            plot(pulled_puff_light_mean_d2, '-r', 'linewidth', 2)
            hold off
        
        subplot(2,1,2)
            plot(time_axis_airpuff, airpuff_vec, '-k','linewidth', 2);
% light on 200-700 ms            
