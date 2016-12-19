%% Analize Chat protocol light only v2
% this file was created on 4/5/2014.
%This file is used for the analysis of files created with extract_ChAT_Data_v2

clear all
% close all

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
load ChAT_Files_v2


% for fileind = 1:length(files)
    fileind = 71;
    channel = 1; % V1 - 1, V2 - 3
    
    fname = files(fileind).extracted_name;        
    path = files(fileind).extracted_path;
    cd(path)
    load(fname) 
    
    original_header = 2;
    trace = []; %1:5; %[]; % if trace is empty the defaulat is to take all traces into analysis

cd 'C:\Users\Inbal\Dropbox\Inbal M.Sc\MATLAB\Project Brain States'
%     cd  'D:\Inbal M.Sc\MATLAB\Project Brain States';
    
    color_table = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

    data_BP_laser = [];
    
%     for header = 1:length(files(1, fileind).headers);

 laser_flag = zeros(1,size(data.x_value,2)); %for each x-value, laser_flag is1 if there was light activation
 
protocol = Param.protocol;
switch protocol

    case 1
        laser_flag(1,:) = [1]; 
    case 2
        laser_flag(1,:) = [1]; 
    case 3
        laser_flag(1,:) = [1]; 
end
    
    for x_value = 1:size(data.x_value,2) %matching END in line 104
      
        if isempty(trace)
            trace = 1:size(data.x_value(x_value).Vm,2);
        end
        
         switch channel
             case 1
                 raw_data(:,trace,x_value) = data.x_value(x_value).Vm(:,trace);
                 sf = Param.sf_Vm;

             case 2
                  raw_data(:,trace,x_value) = data.x_value(x_value).V2(:,trace); 
                 sf = Param.sf_V2;
         end

        dt=1/sf; %[sec]
        time_axis(:,x_value) = (1:size(raw_data,1))*dt;
        sf_laser = Param.sf_laser; %[1/sec]
        dt_laser = 1/sf_laser;
        
 
  locations_x_laser = [];
  
        if isfield(data.x_value(x_value), 'laser') ==1
            time_axis_laser(:,x_value) = (1:size(data.x_value(x_value).laser,1))*dt_laser; 
            laser_vec(:,x_value) = zeros(length(data.x_value(x_value).laser(:,1)),1);
            laser_vec(data.x_value(x_value).laser(:,1) > 2, x_value)=1;   %turning the laser trace into binary
           

            if isempty(find(laser_vec(:,x_value)==1, 1))~=1 %condition will be fulfilled if there was laser activation.
                laser_vec_shifted(:,x_value) = laser_vec(2:end,x_value)-laser_vec(1:end-1, x_value);
                laser_begin(:,x_value) = find(laser_vec_shifted(:,x_value)==1); %find the locations where laser pulse starts
                laser_begin(:,x_value) = laser_begin(:,x_value)+1; %correction for the shift
                laser_end(:,x_value) = find(laser_vec_shifted(:,x_value)==-1); %find the locations where laser pulse ends
                laser_flag(1,x_value) = 1;

                locations_x_laser(1,:,x_value) = laser_begin(:,x_value); %arranging the laser begin and end locations in one variable for plotting
                locations_x_laser(2,:,x_value) = laser_end(:,x_value);
            end
        end
       
    high_pass_freq = 1; %[Hz]
    low_pass_freq = 150; %[Hz]
    
    [data_HP(:,:,x_value)] = fn_High_Pass (raw_data(:,:,x_value), sf, high_pass_freq);
    [data_BP(:,:,x_value)] = fn_Low_Pass (data_HP(:,:,x_value), sf, low_pass_freq);
    
    data_BP_mean(:,x_value) = mean(data_BP(:,:,x_value),2);
    data_BP_std(:,x_value) = std(data_BP(:,:,x_value),0,2);
        
                    
                      data_mean(:,x_value) = mean(raw_data(:,:,x_value),2);
                      data_std(:,x_value) = std(raw_data(:,:,x_value),0,2);
                      
                      
                      %Smoothing the data:
    smooth_data(:,:,x_value) = sgolayfilt(raw_data(:,:,x_value), 1, 29);
    smooth_data_mean(:,x_value) = mean(smooth_data(:,:,x_value),2);
    smooth_data_std(:,x_value) = std(smooth_data(:,:,x_value),0,2);
            end
        
%%  Parameters for plots
 x_value = 1;
 if isempty(locations_x_laser)
   stim1_X=0;
 else
    stim1_X = locations_x_laser(:,:,x_value);
 end
%%  Plotting mean unfiltered trace without polygon std

[Fig,h] = fn_Plot_Trace_v2(data_mean(:,x_value), dt, dt_laser, stim1_X);
%% Plot single unfiltered traces on the same axes

trace_ind = 1:5; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
trace_to_plot = raw_data(:,trace_ind, x_value) ;

    [Fig1,h1] = fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim1_X);
    %% Plot single unfiltered traces on different plots
trace_ind = 1:5; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
trace_to_plot = raw_data(:,trace_ind, x_value) ;
for i=1:length(trace_ind)
    [Fig2(i),h2(i)] = fn_Plot_Trace_v2(trace_to_plot(:,i), dt, dt_laser, stim1_X);
end
%% Plotting mean LFP of whole trace
trace_to_plot = data_BP_mean(:,x_value);

[Fig,h] = fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim1_X);
% Leg_Label{laser_flag==0} = 'light off';
% Leg_Label{laser_flag==1} = 'light on';
% legend(h,Leg_Label)
