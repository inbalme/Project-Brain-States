%% Analize Chat V3
% this file was created on 7/4/2014.
clear all
close all

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
load ChAT_Files

% file 34 (2013-10-28-001.exp2) header 1 depth 460 LFP fiber above cortex mouse ChAT-ChR2 (+/-). no effect...
%                                                  header 2 depth 650
% file 35 (2013-10-28-002.exp2) header 1 depth 250 LFP fiber in NB. no effect...
%                                                  header 2 depth 450
%                                                  header 3 depth 650
%                                                  header 4 depth 650
% file 36 (2013-10-28-003.exp2) header 2 light only depth 346 fiber in NB
%                                                  header 3 puff only 
%                                                  header 4 light+puff 
% file 37 (Set2-2014-03-12-001.exp2) header 2, trace all, LFP light amp 6 depth 230
%                                                          header 3, trace all, LFP light amp 4
% file 38 (Set2-2014-03-12-002.exp2) header 4, trace all, LFP light amp 6 depth 450
%                                                          header 3, trace all, LFP light amp 4                                                       
% file 39 (Set2-2014-03-12-003.exp2) header 2, trace all, LFP light amp 6 depth 700
%                                                          header 3, trace all, LFP light amp 4                                                    
% file 40 (Set2-2014-03-12-004.exp2) header 2 trace 1-3, whole cell.
% file 41 (Set2-2014-03-12-005.exp2) header 2 trace 1-9, cell attached
% file 42 (Set2-2014-03-12-006.exp2) header 2 trace 9, whole cell
% file 45 (Set2-2014-03-17-001.exp2) header 1 trace 1-21 intracellular in Barrel ctx with fiber above brain. cell responds to light
% file 54 (Set2-2014-04-03-001.exp2) header 4 
% file 55 (Set2-2014-04-03-002.exp2) header 3 trace 1-27 

% for fileind = 1:length(files)
    fileind = 36;
   
    fname = files(fileind).extracted_name;        
    path = files(fileind).extracted_path;
    cd(path)
    load(fname) 
    
    original_header = 2;
    header = find(Param.orig_headers==original_header);
    trace = []; % if trace is empty the defaulat is to take all traces into analysis

cd 'C:\Users\Inbal\Dropbox\Inbal M.Sc\MATLAB\Project Brain States'
%     cd  'D:\Inbal M.Sc\MATLAB\Project Brain States';
    
    color_table = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

    data_BP_laser = [];
    
%     for header = 1:length(files(1, fileind).headers);
exp_type = Param.header(header).general.exp_type;
exp_type_logical = strcmp('Standby', exp_type);
if exp_type_logical
    switch_type = 2;
else
    switch_type = 2;
end

switch switch_type

    case 2
        
    laser_flag = zeros(1,size(data.header(header).x_value,2));
    
    for x_value = 1:size(data.header(header).x_value,2)
        if isempty(trace)
            trace = 1:size(data.header(header).x_value(x_value).Vm,2);
        end
        
        channel = 1; % V1 - 1, LFP - 2

         switch channel
             case 1
                 raw_data(:,trace,x_value) = data.header(header).x_value(x_value).Vm(:,trace);
                 sf = Param.header(1,header).stim.sf;

             case 2
                 raw_data(:,trace,x_value) = -1.*data.header(header).x_value(x_value).LFP(:,trace); %for extracellular recording, spikes are upside-down
                 sf = Param.header(1,header).stim.sf_LFP;
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

    high_pass_freq = 1; %[Hz]
    low_pass_freq = 150; %[Hz]
    
    [data_HP(:,:,x_value)] = fn_High_Pass (raw_data(:,:,x_value), sf, high_pass_freq);
    [data_BP(:,:,x_value)] = fn_Low_Pass (data_HP(:,:,x_value), sf, low_pass_freq);
    
    data_BP_mean(:,x_value) = mean(data_BP(:,:,x_value),2);
    data_BP_std(:,x_value) = std(data_BP(:,:,x_value),0,2);
        
                    
                      data_mean(:,x_value) = mean(raw_data,2);
                      data_std(:,x_value) = std(raw_data,0,2);
                      data_LFP_mean(:,x_value) = mean(data.header(header).x_value(x_value).LFP,2);
                      data_LFP_std(:,x_value) = std(data.header(header).x_value(x_value).LFP,0,2);

%Smoothing the data:
    smooth_data = sgolayfilt(raw_data, 1, 29);
    smooth_data_mean(:,x_value) = mean(smooth_data,2);
    smooth_data_std(:,x_value) = std(smooth_data,0,2);
            end
        end
       data_BP_laser = [data_BP_laser(:,:,1) data_BP(:,:,laser_flag==1)]; %data_laser is the data in the traces in which there was a laser stimulus
%     end
    data_BP_laser_mean = mean(data_BP_laser,2);
    data_BP_laser_std = std(data_BP_laser,0,2);
%%  Parameters for plots
 header_to_plot = header;
 x_value = 1;
stim_X = locations_x_laser(:,:,x_value);
trace_to_plot = [1:size(data_mean,2)];

%% Plotting mean unfiltered trace with std

% [Fig,h] = fn_Plot_Trace_std(data_mean, data_std, trace_to_plot, dt, dt_laser, stim_X);
 
%%  Plotting mean unfiltered trace

[Fig,h] = fn_Plot_Trace_v2(data_mean, dt, dt_laser, stim_X);
 
 %%  Plotting mean filtered trace

 [Fig,h] = fn_Plot_Trace_v2(smooth_data_mean, dt, dt_laser, stim_X);
%% Plot single unfiltered traces on the same axes
trace_ind = 1:5; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
trace_to_plot = raw_data(:,trace_ind, x_value) ;

    [Fig1,h1] = fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim_X);
    
    %% Plot single filtered traces on the same axes
trace_ind = 1:5; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
trace_to_plot = smooth_data(:,trace_ind, x_value) ;

    [Fig1,h1] = fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim_X);
    %% Plot single unfiltered traces on different plots
trace_ind = 1:5; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
trace_to_plot = raw_data(:,trace_ind, x_value) ;
for i=1:length(trace_ind)
    [Fig2(i),h2(i)] = fn_Plot_Trace_v2(trace_to_plot(:,i), dt, dt_laser, stim_X);
end

%% Plot single filtered traces on different plots
trace_ind = 1:5; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
trace_to_plot = smooth_data(:,trace_ind, x_value) ;
for i=1:length(trace_ind)
    [Fig2(i),h2(i)] = fn_Plot_Trace_v2(trace_to_plot(:,i), dt, dt_laser, stim_X);
end
%% Plotting mean LFP of whole trace
trace_to_plot = data_BP_mean(:,x_value);

[Fig,h] = fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim_X);
Leg_Label{laser_flag==0} = 'light off';
Leg_Label{laser_flag==1} = 'light on';
legend(h,Leg_Label)
% %% binning the mean trace and comparing to the baseline using two-tailed t-test
% bin_win = 0.1; %[sec]
% bin_size = bin_win./dt; %#data points
% bin_num = size(raw_data,1)./bin_size;
% baseline_win = 0.100; %[sec]
% baseline_size = baseline_win./dt; %#data points
% mean_baseline_vec = mean(raw_data(1:baseline_size,:), 1);
% mean_baseline = mean(mean_baseline_vec);
% reshape_mean = reshape(data_mean,bin_size, bin_num);
% for i=1:bin_num
% [h(i),p(i)]=ttest(reshape_mean(:,i),data_mean(1:bin_size,1), 'Alpha', 0.01);
% end
%% binning the matrix and comparing mean to baseline
% bin_win = 0.01; %[sec]
% 
% [bin_mean bin_std] = fn_binning_data(raw_data,dt, bin_win );
% 
% baseline_win = 0.100; %[sec]
% baseline_size = baseline_win./dt; %#data points
% mean_baseline_vec = mean(raw_data(1:baseline_size,:), 1);
% std_baseline_vec = std(raw_data(1:baseline_size,:),0,1);
% mean_baseline = mean(mean_baseline_vec);
% mean_std_baseline = mean(std_baseline_vec);
% hold on
% errorbar(mean(bin_mean,2),mean(bin_std,2))
% errorbar(mean_baseline.*ones(size(bin_mean,1)),mean_std_baseline.*ones(size(bin_mean,1)), 'r-')
% hold off
%% Plotting mean LFP of whole trace, from several headers merged, only laser traces
% header_to_plot = header;
% trace_to_plot = [1:size(data_BP_laser_mean,2)];
% stim_X = locations_x_laser;
% 
% [Fig1,h1] = fn_Plot_Trace_std(data_BP_laser_mean, data_BP_laser_std, trace_to_plot, dt, dt_laser, stim_X);

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

% %% Plot single LFP traces
% for trace_to_plot = [2 6 12 17] %1:size(data_BP(:,:,laser_flag==1),2);
% [Fig3,h3] = fn_Plot_Trace(data_BP(:,:,laser_flag==1), header_to_plot, trace_to_plot, dt, dt_laser, stim_X);
% end
% % end
