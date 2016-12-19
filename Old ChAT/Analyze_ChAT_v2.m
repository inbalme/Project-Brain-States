%% Analize Chat V2
% this file was used for analysis for the data in my phd proposal
clear all
close all

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
load ChAT_Files

%file 7 (2013-02-13-002) headers 1-3 (2,3,4) LFP S1 with fiber above brain
%file 19 (2013-02-27-002) headers 1,2 - spont+laser in BF 
%file 21 (2013-03-20-001) headers 1-6 (2-7) M1

for fileind = 31 %1:length(files)
    
    fname = files(fileind).extracted_name;        
    path = files(fileind).extracted_path;
    cd(path)
    load(fname) 

cd 'C:\Users\Inbal\Dropbox\Inbal M.Sc\MATLAB\Project Brain States'
%     cd  'D:\Inbal M.Sc\MATLAB\Project Brain States';
    
    color_table = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

    data_laser = [];
    
    for header = 2; %1:length(files(1, fileind).headers);
        
exp_type = Param.header(header).general.exp_type;
exp_type_logical = strcmp('Standby', exp_type);
if exp_type_logical
    switch_type = 2;
else
    switch_type = 2;
end

switch switch_type

    case 2
        
%     header = 6;
%     x_value = 1;
    laser_flag = zeros(1,size(data.header(header).x_value,2));
    for x_value = 1:size(data.header(header).x_value,2)
    
        channel = 1; % V1 - 1, Ext1 - 2

    switch channel

    case 1
        raw_data(:,:,x_value) = data.header(header).x_value(x_value).Vm;
        sf = Param.header(1,header).stim.sf;

    case 2
        raw_data(:,:,x_value) = -1.*data.header(header).x_value(x_value).Ext1; %for extracellular recording, spikes are upside-down
        sf = Param.header(1,header).stim.sf_Ext1;

    end

        dt=1/sf; %[sec]
        time_axis(:,x_value) = (1:size(raw_data,1))*dt;
        sf_airpuff = Param.header(header).stim.sf_airpuff; %[1/sec]
        dt_airpuff = 1/sf_airpuff;
        sf_galvano = Param.header(header).stim.sf_galvano; %[1/sec]
        time_axis_airpuff(:,x_value) = size(data.header(header).x_value(x_value).airpuff,1)*dt_airpuff;    
        airpuff_vec(:,x_value) = data.header(header).x_value(x_value).airpuff(:,1)./...
            max(data.header(header).x_value(x_value).airpuff(:,1));
        dt_galvano = 1/sf_galvano;
        time_axis_galvano = size(data.header(header).x_value(x_value).galvano,1)*dt_galvano;
        sf_laser = Param.header(header).stim.sf_laser; %[1/sec]
        dt_laser = 1/sf_laser;
        time_axis_laser(:,x_value) = (1:size(data.header(header).x_value(x_value).laser,1))*dt_laser; 
        laser_vec(:,x_value) = zeros(length(data.header(header).x_value(x_value).laser(:,1)),1);
        laser_vec(data.header(header).x_value(x_value).laser(:,1) > 2, x_value)=1;   %turning the laser trace into binary
        locations_x_laser = [];
        
        if isempty(find(laser_vec(:,x_value)==1, 1))~=1%isempty(laser_begin(:,x_value))~=1 %condition will be fulfilled if there was laser activation.
            laser_vec_shifted(:,x_value) = laser_vec(2:end,x_value)-laser_vec(1:end-1, x_value);
            laser_begin(:,x_value) = find(laser_vec_shifted(:,x_value)==1); %find the locations where laser pulse starts
            laser_begin(:,x_value) = laser_begin(:,x_value)+1; %correction for the shift
            laser_end(:,x_value) = find(laser_vec_shifted(:,x_value)==-1); %find the locations where laser pulse ends
        
        
            laser_flag(x_value) = 1; %for each x_value: flag=0 for laser off, flag=1 for laser on
            locations_x_laser(1,:) = laser_begin(:,x_value); %arranging the laser begin and end locations in one variable for plotting
            locations_x_laser(2,:) = laser_end(:,x_value);
        end

    high_pass_freq = 0.2; %[Hz]
    low_pass_freq = 150; %[Hz]
    
    [data_HP(:,:,x_value)] = fn_High_Pass (raw_data(:,:,x_value), sf, high_pass_freq);
    [data_BP(:,:,x_value)] = fn_Low_Pass (data_HP(:,:,x_value), sf, low_pass_freq);
    
    data_mean(:,x_value) = mean(data_BP(:,:,x_value),2);
    data_std(:,x_value) = std(data_BP(:,:,x_value),0,2);
        
    
%     for x_value = 1:size(data.header(header).x_value,2)
%                     
%                       data.header(header).x_value(x_value).Vm_mean = mean(data.header(header).x_value(x_value).Vm,2);
%                       data.header(header).x_value(x_value).Vm_std = std(data.header(header).x_value(x_value).Vm');
%                       data.header(header).x_value(x_value).Ext1_mean = mean(data.header(header).x_value(x_value).Ext1,2);
%                       data.header(header).x_value(x_value).Ext1_std = std(data.header(header).x_value(x_value).Ext1');
%                 end
            end
        end
       data_laser = [data_laser(:,:,1) data_BP(:,:,laser_flag==1)]; %data_laser is the data in the traces in which there was a laser stimulus
    end
    data_laser_mean = mean(data_laser,2);
    data_laser_std = std(data_laser,0,2);
%% Plotting mean LFP of whole trace
header_to_plot = header;
trace_to_plot = [1:size(data_mean,2)];
stim_X = locations_x_laser;

[Fig,h] = fn_Plot_Trace(data_mean, header_to_plot, trace_to_plot, dt, dt_laser, stim_X);
Leg_Label{laser_flag==0} = 'laser off';
Leg_Label{laser_flag==1} = 'laser on';
legend(h,Leg_Label)

%% Plotting mean LFP of whole trace, from several headers merged, only laser traces
header_to_plot = header;
trace_to_plot = [1:size(data_laser_mean,2)];
stim_X = locations_x_laser;

[Fig,h] = fn_Plot_Trace_std(data_laser_mean, data_laser_std, trace_to_plot, dt, dt_laser, stim_X);

%% Comparing mean LFP within laser trace
% header_to_plot = header;
% trace_to_plot = [1 2];
% laser_delay = Param.header(1,header).stim.facade(11); %[sec]
% data_off = data_BP(1:laser_delay*sf,:,laser_flag==1);
% data_on = data_BP(laser_delay*sf+1:2.*laser_delay*sf,:,laser_flag==1);
% mean_data_off = data_mean(1:laser_delay*sf,laser_flag==1);
% mean_data_on = data_mean(laser_delay*sf+1:2.*laser_delay*sf,laser_flag==1);
% laser_trace = [mean_data_off mean_data_on];
% stim_X_devided = locations_x_laser-laser_delay*sf_laser;
% 
% [Fig2,h2] = fn_Plot_Trace(laser_trace, header_to_plot, trace_to_plot, dt, dt_laser, stim_X_devided);
% Leg_Label{1} = 'laser off';
% Leg_Label{2} = 'laser on';
% legend(h2,Leg_Label)

%% Plot single LFP traces
for trace_to_plot = [2 6 12 17] %1:size(data_BP(:,:,laser_flag==1),2);
[Fig3,h3] = fn_Plot_Trace(data_BP(:,:,laser_flag==1), header_to_plot, trace_to_plot, dt, dt_laser, stim_X);
end
end
%%
% figure
%         set(gcf,'color','w')
%           hold on 
%             for x_value = 1:size(data.header(header).x_value,2)
%                 plot(time_axis, data_mean(:,x_value), 'color', color_table(x_value,:), 'linewidth', 2)
%             end
%             
%           hold off
% 
