%% Analyze_ChAT
%calls the functions: 

clear all
close all

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
load ChAT_Files

for fileind = 1
   
    files(fileind).extracted_name = 'Set2-2012-11-20-001';
    files(fileind).extracted_path = 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
    fname = files(fileind).extracted_name;        
    path = files(fileind).extracted_path;
    cd(path)
    load(fname) 
    
    cd  'D:\Inbal M.Sc\MATLAB\Project Brain States';
    
        for header = 1:1 %length(files(1, fileind).headers);
            stim_type = Param.header(header).stim.type;
            %%  Airpuff/Galvano (no protocol) stimulus:
            if stim_type== 1 || stim_type== 2 || stim_type== 3 || stim_type==4 || stim_type==5 ||stim_type==6

                sf_Vm = Param.header(header).stim.sf; %[1/sec]
                dt_Vm=1/sf_Vm; %[sec]
                time_axis_Vm=(1:size(data.header(header).Vm,1))*dt_Vm;
                sf_airpuff = Param.header(header).stim.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                time_axis_airpuff = (1:size(data.header(header).airpuff,1))*dt_airpuff;
                sf_galvano = Param.header(header).stim.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;
                time_axis_galvano = (1:size(data.header(header).galvano,1))*dt_galvano;
                sf_Ext1 = sf_Vm/3; %[1/sec]
                dt_Ext1=1/sf_Ext1; %[sec]
                time_axis_Ext1=(1:size(data.header(header).Ext1,1))*dt_Ext1;
     
                  data.header(header).Vm_mean = mean(data.header(header).Vm,2);
                  data.header(header).Vm_std = std(data.header(header).Vm');
                  data.header(header).Ext1_mean = mean(data.header(header).Ext1,2);
                  data.header(header).Ext1_std = std(data.header(header).Ext1');
            end
            
            %%  Airpuff/Galvano protocol stimulus:
            if stim_type== 7 || stim_type== 8
                sf_Vm = Param.header(header).stim.sf; %[1/sec]
                dt_Vm=1/sf_Vm; %[sec]
                time_axis_Vm=(1:size(data.header(header).x_value(1).Vm,1))*dt_Vm;
                sf_airpuff = Param.header(header).stim.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                time_axis_airpuff = (1:size(data.header(header).x_value(1).airpuff,1))*dt_airpuff;
                sf_galvano = Param.header(header).stim.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;
                time_axis_galvano = (1:size(data.header(header).x_value(1).galvano,1))*dt_galvano;
                sf_Ext1 = sf_Vm/3; %[1/sec]
                dt_Ext1=1/sf_Ext1; %[sec]
                time_axis_Ext1=(1:size(data.header(header).x_value(1).Ext1,1))*dt_Ext1;
                
                laser_begin_time = 1; %[sec]
                laser_end_time = 2; %[sec]  
                laser_begin_sample = floor(laser_begin_time*sf_Ext1); 
                laser_end_sample = floor(laser_end_time*sf_Ext1); 

                laser_off_begin_sample = 1;
                laser_off_end_sample = floor(laser_begin_time*sf_Ext1);
                
                laser_off_control_begin_time = 10; %[sec]
                laser_off_control_end_time = 11; %[sec]
                laser_off_control_begin_sample = floor(laser_off_control_begin_time*sf_Ext1);
                laser_off_control_end_sample = floor(laser_off_control_end_time*sf_Ext1);
                
                galvano_begin = 2; %[sec]
                galvano_end = 4; %[sec]

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
% depth 200um:header 1-3 - galvano protocol 
%             header 4 - airpuff protocol
%depth 300um: header 5 - airpuff protocol
%             header 6 - galvano protocol
%depth 550um: header 7,9,12,14 - airpuff alone
%             header 8,10,11,13 - airpuff + light
%             header 15,16 - airpuff protocol
%depth 800um: header 17,18 - airpuff protocol
%             header 19,21,23,25 - airpuff alone
%             header 20,22,24,26 - airpuff + light
%% Depth 200 - header 1-3 galvano

header = 1;
high_pass_freq = 0.5; %[Hz]
low_pass_freq = 1; %[Hz]

% dividing each trace to two parts: laser off and laser on:
data_off = data.header(header).x_value(2).Ext1(laser_off_begin_sample:laser_off_end_sample,:);
data_off_control = data.header(header).x_value(2).Ext1(laser_off_control_begin_sample+1:laser_off_control_end_sample,:);
data_on = data.header(header).x_value(2).Ext1(laser_begin_sample+1:laser_end_sample,:);

%bandpass filtering 0.5-1Hz and plotting
[data_HP] = fn_High_Pass (data.header(header).x_value(2).Ext1, sf_Ext1, high_pass_freq);
[data_BP] = fn_Low_Pass (data_HP, sf_Ext1, low_pass_freq);

for trace = 5
    
    %power spectrum:
    [DATA_OFF,freq_off] = fn_Amp_Spectrum (data_off(:,trace), sf_Ext1);
    [DATA_ON,freq_on] = fn_Amp_Spectrum (data_on(:,trace), sf_Ext1);
end

galvano_vec = data.header(header).x_value(2).galvano(:,1)./max(data.header(header).x_value(2).galvano(:,1));

        figure
        set(gcf,'color','w')
        subplot(2,1,1)
            hold on 
            plot(time_axis_Ext1,data.header(header).x_value(2).Ext1_mean, '-k', 'linewidth', 2)
            plot(time_axis_Ext1,data.header(header).x_value(2).Ext1_mean, '-c', 'linewidth', 1)
            hold off
        
        subplot(2,1,2)
            plot(time_axis_galvano, galvano_vec, '-k','linewidth', 2);
            
        pos_galvano(1,:)= [0.15 , 0.35 , 0.7 , 0.55];
        pos_galvano(2,:)= [0.15 , 0.15 , 0.7 , 0.15];
        
%         galvano_vec = galvano_vec-galvano_vec(1);
%         galvano_vec = galvano_vec./max(galvano_vec);
            
        figure
            
        set(gcf,'color','w')
x1limits = [min(time_axis_Ext1) max(time_axis_Ext1)];
x1ticks = [];
y1limits = [-0.6 0.4];
y1ticks = [-0.6:0.2:0.4];

x2limits = [0 15];
x2ticks = [0:1:15];
y2limits = [-0.5 3.5];
y2ticks = [0:1:3];

    axes('position',pos_galvano(1,:) ); %[left, bottom,width, height]
    plot(time_axis_Ext1,data_BP(:,trace),'-k','linewidth', 1)
         set( gca, 'xlim', x1limits, 'xtick', x1ticks, 'ylim', y1limits, 'ytick', y1ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010], 'fontname', 'helvetica', 'fontweight', 'bold','box', 'off' );
    title(['Cell ',num2str(fileind)],'FontSize', 14,'fontname', 'helvetica', 'fontweight', 'bold');
    ylabel('Vm [mV]', 'FontSize', 12,'fontname', 'helvetica', 'fontweight', 'bold');
    
    axes('position',pos_galvano(2,:) ); %[left, bottom,width, height]
    plot(time_axis_galvano, galvano_vec, '-k','linewidth', 2);
    set( gca, 'xlim', x2limits, 'xtick', x2ticks, 'xminortick', 'on', 'ylim', y2limits, 'ytick', y2ticks,...
        'ticklength', [0.010 0.010], 'fontname', 'helvetica', 'fontweight', 'bold','box', 'off' );
    xlabel('Time [sec]' ,'FontSize', 12, 'fontname', 'helvetica', 'fontweight', 'bold');
            
%%  Depth 550 - header 7+8 

airpuff_vec = data.header(7).airpuff(:,1)./max(data.header(7).airpuff(:,1));            
        figure
        set(gcf,'color','w')
        subplot(2,1,1)
            hold on 
            plot(time_axis_Ext1,data.header(7).Ext1_mean, '-k', 'linewidth', 2)
%             plot(data.header(7).Ext1_mean - data.header(7).Ext1_std','k');
            plot(time_axis_Ext1,data.header(8).Ext1_mean, '-c', 'linewidth', 1)
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
            plot(time_axis_Ext1,pulled_puff_alone_mean, '-kb', 'linewidth', 2)
%             plot(pulled_puff_alone_mean - pulled_puff_alone_std,'k');
            plot(time_axis_Ext1,pulled_puff_light_mean, '-c', 'linewidth', 2)
            hold off
        
        subplot(2,1,2)
            plot(time_axis_airpuff, airpuff_vec, '-k','linewidth', 2);
% light on 200-700 ms            
%% Depth 800
airpuff_vec = data.header(19).airpuff(:,1)./max(data.header(19).airpuff(:,1));

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
            plot(time_axis_Ext1,pulled_puff_alone_mean_d2, '-k', 'linewidth', 2)
%             plot(pulled_puff_alone_mean_d2 - pulled_puff_alone_std_d2,'k');
            plot(time_axis_Ext1,pulled_puff_light_mean_d2, '-c', 'linewidth', 2)
            hold off
        
        subplot(2,1,2)
            plot(time_axis_airpuff, airpuff_vec, '-k','linewidth', 2);
% light on 200-700 ms            
