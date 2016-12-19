%% Analyze Chat protocol light+galvano v2
% This file was created on 7/5/2014.
%This file is used for the analysis of files created with extract_ChAT_Data_v2
% with protocol 4 (ChAT intracellular light+galvnao (3 x-values))
%       protocol 5 (ChAT intracellular light+galvnao (5 x-values))
%       protocol 6 (ChAT intracellular current injections light+galvnao (15 x-values))
%       protocol 7 - ChAT: Stim+light (yonatan's protocol). 4 x-values:
%                 light train, single galvo, light train+galvo+depo, single galvo+depo
%       protocol 8 - ChAT intracellular: light only (3 x-values): 1) train
%       20 Hz, 10ms ondur, 50 pulses. 2) single pulse 5 sec ondur. 3)single pulse 10 ms
%       protocol 9 - IC: I-V curve

clear all
% close all

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
load ChAT_Files_v2


% for fileind = 1:length(files)
    fileind =19;
   channel = 1; % V1 - 1, V2 - 3
   
    fname = files(fileind).extracted_name;        
    path = files(fileind).extracted_path;
    cd(path)
    load(fname) 
    
    original_header = 2;
    trace =[]; %230:260; %find(Param.orig_headers==original_header); %[]; % if trace is empty the defaulat is to take all traces into analysis

cd 'C:\Users\Inbal\Dropbox\Inbal M.Sc\MATLAB\Project Brain States'
%     cd  'D:\Inbal M.Sc\MATLAB\Project Brain States';
    
    color_table = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

    data_BP_laser = [];
    
%     for header = 1:length(files(1, fileind).headers);
% exp_type = Param.exp_type;
% exp_type_logical = strcmp('Standby', exp_type);
% if exp_type_logical
%     switch_type = 2; %for standby
% else
%     switch_type = 2; %for protocol
% end

laser_flag = zeros(1,size(data.x_value,2)); %for each x-value, laser_flag is1 if there was light activation
galvano_flag = zeros(1,size(data.x_value,2)); %for each x-value, galvano_flag is1 if there was galvano activation
 
protocol = Param.protocol;

if isempty(protocol)
    
else
    switch protocol
        case 1
            laser_flag(1,:) = [1]; 
        case 2
            laser_flag(1,:) = [1]; 
        case 3
            laser_flag(1,:) = [1]; 
        case 4
            laser_flag(1,:) = [1 0 1];
            galvano_flag(1,:) = [0 1 1];     
        case 5
             laser_flag(1,:) = [1 0 1 0 1];
             galvano_flag(1,:) = [0 1 1 1 1];
        case 6
            if length(laser_flag(1,:))==9
                laser_flag(1,:) = [1 1 1 0 0 0 1 1 1];
                galvano_flag(1,:) = [0 0 0 1 1 1 1 1 1];
            else
             laser_flag(1,:) = [1 1 1 0 0 0 1 1 1 0 0 0 1 1 1];
             galvano_flag(1,:) = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1];   
            end
        case 8
             laser_flag(1,:) = [1 1 1];
             galvano_flag(1,:) = [0 0 0];  
        case 9
             laser_flag(1,:) = [0 0 0 0 0 0 0];
             galvano_flag(1,:) = [0 0 0 0 0 0 0];  
        case 10
            laser_flag(1,:) = [0 0 0];
            galvano_flag(1,:) = [0 0 0 ];  
    end
end
         
    locations_x_galvano = [];
    locations_x_laser = [];

    for x_value = 1:size(data.x_value,2) 
      
                if isempty(trace)
                    trace = 1:size(data.x_value(x_value).Vm,2);
                end

                 switch channel
                     case 1
                         raw_data(:,1:length(trace),x_value) = data.x_value(x_value).Vm(:,trace);
                         sf = Param.sf_Vm;

                     case 2
                         raw_data(:,1:length(trace),x_value) = data.x_value(x_value).V2(:,trace); 
                         sf = Param.sf_V2;

                 end

                dt=1/sf; %[sec]
                time_axis(:,x_value) = (1:size(raw_data,1))*dt;
                sf_airpuff = Param.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                sf_galvano = Param.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;
                sf_laser = Param.sf_laser; %[1/sec]
                dt_laser = 1/sf_laser;
                sub_baseline_time = 10; %[ms]
                sub_baseline_time = sub_baseline_time/1000; %[sec]
                sub_baseline_interval = sub_baseline_time*sf;

                if Param.sf_airpuff==0
                else
                time_axis_airpuff(:,x_value) = size(data.x_value(x_value).airpuff,1)*dt_airpuff;    
                airpuff_vec(:,x_value) = data.x_value(x_value).airpuff(:,1)./...
                    max(data.x_value(x_value).airpuff(:,1));
                end

                 if Param.sf_galvano==0
                else
                time_axis_galvano = size(data.x_value(x_value).galvano,1)*dt_galvano;
                galvano_vec(:,x_value) = zeros(length(data.x_value(x_value).galvano(:,1)),1);
                galvano_threshold = abs(mean(data.x_value(x_value).galvano(:,1))) + 3.*abs(std(data.x_value(x_value).galvano(:,1)));
                galvano_vec(abs(data.x_value(x_value).galvano(:,1)) > galvano_threshold, x_value)=1;   %turning the laser trace into binary
                 end

                time_axis_laser(:,x_value) = (1:size(data.x_value(x_value).laser,1))*dt_laser; 
                laser_vec(:,x_value) = zeros(length(data.x_value(x_value).laser(:,1)),1);
                laser_vec(data.x_value(x_value).laser(:,1) > 2, x_value)=1;   %turning the laser trace into binary

                if isempty(find(laser_vec(:,x_value)==1, 1))~=1 %condition will be fulfilled if there was laser activation.
                    if laser_flag(1,x_value)==1 %condition will be fulfilled if there was laser activation.
                        laser_vec_shifted(:,x_value) = laser_vec(2:end,x_value)-laser_vec(1:end-1, x_value);
                        laser_begin(:,x_value) = find(laser_vec_shifted(:,x_value)==1); %find the locations where laser pulse starts
                        laser_begin(:,x_value) = laser_begin(:,x_value)+1; %correction for the shift
                        laser_end(:,x_value) = find(laser_vec_shifted(:,x_value)==-1); %find the locations where laser pulse ends
                        locations_x_laser(1,:,x_value) = laser_begin(:,x_value); %arranging the laser begin and end locations in one variable for plotting
                        locations_x_laser(2,:,x_value) = laser_end(:,x_value);
                    end
                end
                    if galvano_flag(1,x_value) ==1; %condition will be fulfilled if there was galvano activation.
        %          if isempty(find(galvano_vec(:,x_value)==1, 1))~=1 %condition will be fulfilled if there was galvano activation.
                    galvano_vec_shifted(:,x_value) = galvano_vec(2:end,x_value)-galvano_vec(1:end-1, x_value);
                    galvano_begin{x_value}(:,1) = find(galvano_vec_shifted(:,x_value)==1); %find the locations where galvano pulse starts
                    galvano_begin{x_value}(:,1) = galvano_begin{x_value}(:,1)+1; %correction for the shift
                    galvano_end{x_value}(:,1) = find(galvano_vec_shifted(2:end,x_value)==-1); %find the locations where galvano pulse ends
                        if length(galvano_begin{x_value}(:,1))>length(galvano_end{x_value}(:,1)) %if galvano_begin is larger than galvano_end, take only the points in galvano_begin which has a match in galvano_end
                            galvano_begin_trunc{x_value}(:,1)=galvano_begin{x_value}(1:length(galvano_end{x_value}(:,1)),1);
                            galvano_begin{x_value}(:,1)=[];
                            galvano_begin{x_value}=galvano_begin_trunc{x_value}(:,1);
                        end
                    locations_x_galvano{x_value}(1,:) = galvano_begin{x_value}(:,1); %arranging the laser begin and end locations in one variable for plotting
                    locations_x_galvano{x_value}(2,:) = galvano_end{x_value}(:,1);            
        %             galvano_begin(:,x_value) = find(galvano_vec_shifted(:,x_value)==1); %find the locations where galvano pulse starts
        %             galvano_begin(:,x_value) = galvano_begin(:,x_value)+1; %correction for the shift
        %             galvano_end(:,x_value) = find(galvano_vec_shifted(:,x_value)==-1); %find the locations where galvano pulse ends
        %             locations_x_galvano(1,:,x_value) = galvano_begin(:,x_value); %arranging the laser begin and end locations in one variable for plotting
        %             locations_x_galvano(2,:,x_value) = galvano_end(:,x_value);
                end

                high_pass_freq = 1; %[Hz]
                low_pass_freq = 150; %[Hz]

                [data_HP(:,:,x_value)] = fn_High_Pass (raw_data(:,:,x_value), sf, high_pass_freq);
                [data_BP(:,:,x_value)] = fn_Low_Pass (data_HP(:,:,x_value), sf, low_pass_freq);

                data_BP_mean(:,x_value) = mean(data_BP(:,:,x_value),2);
                data_BP_std(:,x_value) = std(data_BP(:,:,x_value),0,2);


                                  data_mean(:,x_value) = mean(raw_data(:,:,x_value),2);
                                  data_std(:,x_value) = std(raw_data(:,:,x_value),0,2);
                                  data_var(:,x_value) = var(raw_data(:,:,x_value),0,2);
                                  data_fano_factor(:,x_value) = (-1).*data_var(:,x_value)./data_mean(:,x_value); %this is incorrect. need to subtract baseline!

                                  %Smoothing the data:
                smooth_data(:,:,x_value) = sgolayfilt(raw_data(:,:,x_value), 1, 29);
                smooth_data_mean(:,x_value) = mean(smooth_data(:,:,x_value),2);
                smooth_data_std(:,x_value) = std(smooth_data(:,:,x_value),0,2);
                smooth_data_var(:,x_value) = var(smooth_data(:,:,x_value),0,2);
    end
%%    for file 16 and 17 only
%  med_trace = [];
%  laser_latency =  90/1000*sf; %locations_x_laser(1,1,x_value);          
%  med_trace(:,:,x_value) = median(smooth_data(laser_latency-sub_baseline_interval:laser_latency,:,x_value),1);    
% for i=1:size(raw_data,2)
%  data_sub_median(:,i,x_value) = smooth_data(:,i,x_value)- med_trace(:,i,x_value);
% end
%%  Stim vectors
 for x_value = 1:size(data.x_value,2) 
         if isempty(locations_x_laser)
             stim1_X{x_value} = [];
         else
            stim1_X{x_value} = locations_x_laser(:,:,x_value);
         end

        if isempty(locations_x_galvano)
            stim2_X{x_value} = [];
        else
            stim2_X{x_value} = locations_x_galvano{x_value}(:,:);
        end

%         stim2_X{x_value} = [12500 13500 14500 15500 16500 17500 18500 19500 20500 21500; 12700 13700 14700 15700 16700 17700 18700 19700 20700 21700];%cell #38
%         stim2_X{x_value} = [30500 31500 32500 33500 34500 35500 36500 37500 38500 39500; 30700 31700 32700 33700 34700 35700 36700 37700 38700 39700]; %cell #37
%         dt_galvano=dt;
 end
 %% Variance plots - two x_values on each plot
x_value = 4:5;

% plots for first x-value:
       trace_ind = 2:6; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
       l=size(raw_data(:,trace_ind, x_value(1)),1);
        l=l/2-1;
        DC=10.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        trace_to_plot = raw_data(:,trace_ind, x_value(1))+DC ;
        [Fig1,h1]= fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim1_X{x_value(1)}, dt_galvano, stim2_X{x_value(1)});
        trace_to_plot = raw_data(:,trace_ind, x_value(2))+DC ;
        [Fig2,h2]= fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim1_X{x_value(2)}, dt_galvano, stim2_X{x_value(2)});
        
        [Fig3,h3]= fn_Plot_Trace_v2(data_std(:,x_value), dt, dt_laser, stim1_X{x_value(2)}, dt_galvano, stim2_X{x_value(2)});
        [Fig4,h4] = fn_Plot_Trace_v2(data_mean(:,x_value), dt, dt_laser, stim1_X{x_value(2)}, dt_galvano, stim2_X{x_value(2)});

     
        ax1 = get(Fig1, 'children');
        pos1 = [0.08 , 0.79 , 0.8 , 0.2];
        top1 = pos1(1,2)+pos1(1,4);

        ax2 = get(Fig2, 'children');
        pos2 = [0.08 , 0.56 , 0.8 , 0.2];
        top2 = pos2(1,2)+pos2(1,4);

        ax3 = get(Fig3, 'children');
        pos3 = [0.08 , 0.33 , 0.8 , 0.2];
        top3 = pos3(1,2)+pos3(1,4);

        ax4 = get(Fig4, 'children');
        pos4 = [0.08 , 0.1 , 0.8 , 0.2];
        top4 = pos4(1,2)+pos4(1,4);
        
        F = figure;
        set(gcf,'color','w');
        set(gcf,'DefaultAxesFontSize',12);
        set(gcf,'DefaultAxesFontName','helvetica');
        set(gcf, 'PaperType', 'A4');
        set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

        ax_copy1 = copyobj(ax1,F); % ax1 to new fig
        set(ax_copy1(1),'position',pos1(1,:)) % Set its position  

        ax_copy2 = copyobj(ax2,F); % ax2 to new fig
        set(ax_copy2(1),'position',pos2(1,:)) % Set its position  

        ax_copy3 = copyobj(ax3,F); % ax3 to new fig
        set(ax_copy3(1),'position',pos3(1,:)) % Set its position  
        
        ax_copy4 = copyobj(ax4,F); % ax3 to new fig
        set(ax_copy4(1),'position',pos4(1,:)) % Set its position  
        
        close Figure 1 Figure 2 Figure 3 Figure 4
%% Variance plot - one x_value at a time
       x_value = 1;
       trace_ind = 2:6; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
       l=size(raw_data(:,trace_ind, x_value),1);
        l=l/2-1;
        DC=10.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        trace_to_plot = raw_data(:,trace_ind, x_value)+DC ;
        [Fig1,h1]= fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim1_X{x_value}, dt_galvano, stim2_X{x_value});
        [Fig2,h2]= fn_Plot_Trace_v2(data_std(:,x_value), dt, dt_laser, stim1_X{x_value}, dt_galvano, stim2_X{x_value});
        [Fig3,h3] = fn_Plot_Trace_v2(data_mean(:,x_value), dt, dt_laser, stim1_X{x_value}, dt_galvano, stim2_X{x_value});



        ax1 = get(Fig1, 'children');
        pos1 = [0.08 , 0.7 , 0.8 , 0.2];
        top1 = pos1(1,2)+pos1(1,4);

        ax2 = get(Fig2, 'children');
        pos2 = [0.08 , 0.4 , 0.8 , 0.2];
        top2 = pos2(1,2)+pos2(1,4);

        ax3 = get(Fig3, 'children');
        pos3 = [0.08 , 0.1 , 0.8 , 0.2];
        top3 = pos3(1,2)+pos3(1,4);

        F = figure;
        set(gcf,'color','w');
        set(gcf,'DefaultAxesFontSize',12);
        set(gcf,'DefaultAxesFontName','helvetica');
        set(gcf, 'PaperType', 'A4');
        set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

        ax_copy1 = copyobj(ax1,F); % ax1 to new fig
        set(ax_copy1(1),'position',pos1(1,:)) % Set its position  

        ax_copy2 = copyobj(ax2,F); % ax2 to new fig
        set(ax_copy2(1),'position',pos2(1,:)) % Set its position  

        ax_copy3 = copyobj(ax3,F); % ax3 to new fig
        set(ax_copy3(1),'position',pos3(1,:)) % Set its position  
        close Figure 1 Figure 2 Figure 3
%%  Plotting mean unfiltered trace without polygon std
x_value = 5;
% [Fig,h] = fn_Plot_Trace_v2(data_mean(:,x_value), dt, dt_laser, stim1_X{x_value}, dt_galvano, stim2_X{x_value});

data_mat=[];
for i=[1,4,5]
data_mat =[data_mat data_mean(:,i)-mean(data_mean(:,i))];
end
data_mat=[data_mat data_mat(:,1)+ data_mat(:,2)];

[Fig,h] = fn_Plot_Trace_v2(data_mat, dt, dt_laser, stim1_X{x_value}, dt_galvano, stim2_X{x_value});

%% Plot power spectrum
start_time = 0; %[sec]
duration = 10; %[sec]
x_value = 1;

DC=[]; spec_mat_noDC=[]; ;spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = []; interval = []; Y = []; f = [];

start_sample = start_time.*sf;
if start_time==0
    start_sample = 1;
end
end_sample = start_sample+duration.*sf-1;
interval = start_sample:end_sample;
spec_mat = raw_data(interval,:,x_value);
DC= mean(spec_mat,1);
l=size(spec_mat,1);
l=l/2-1;
DC=wextend('addrow', 'per', DC, l);
spec_mat_noDC = spec_mat-DC;
spec_mat_mean = mean(spec_mat,2);

[Y,f] = fn_Amp_Spectrum(spec_mat_noDC,sf);
% 
% for i=1:size(spec_mat_noDC,2)
%     [Y(i),f] = fn_Amp_Spectrum(spec_mat_noDC(:,i),sf);

%% Plot single unfiltered traces on the same axes

trace_ind = 1:5; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
trace_to_plot = raw_data(:,trace_ind, x_value) ;
% trace_to_plot = data_BP(:,trace_ind, x_value); %for a single trace of LFP band pass filtered

    [Fig1,h1] = fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim1_X{x_value}, dt_galvano, stim2_X{x_value});
    %% Plot single unfiltered traces on different plots
trace_ind = 1:5; %[1:size(raw_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
trace_to_plot = raw_data(:,trace_ind, x_value) ;
for i=1:length(trace_ind)
    [Fig2(i),h2(i)] = fn_Plot_Trace_v2(trace_to_plot(:,i), dt, dt_laser, stim1_X{x_value}, dt_galvano, stim2_X{x_value});
end
%% Plotting mean LFP of whole trace
trace_to_plot = data_BP_mean(:,x_value);

[Fig,h] = fn_Plot_Trace_v2(trace_to_plot, dt, dt_laser, stim1_X{x_value}, dt_galvano, stim2_X{x_value});
% Leg_Label{laser_flag==0} = 'light off';
% Leg_Label{laser_flag==1} = 'light on';
% legend(h,Leg_Label)

%% 
filename = 'F20P5X1+4+5+expected'; % 'F72P8X2_LFP_power spectrum_0-5Hz'; %'F62P6X6+9_zero_traces+std+mean_zoom-in';
% cd 'c:\Users\Inbal\Dropbox\Inbal M.Sc\Data PhD\ChAT homozygous\Current inj+light+galvano\For Ai';
cd 'c:\Users\Inbal\Dropbox\Inbal M.Sc\Data PhD\ChAT homozygous\For Ai';
% cd 'c:\Users\Inbal\Dropbox\Inbal M.Sc\Data PhD\Figures for ISF 2014'
% cd 'c:\Users\Inbal\Dropbox\Inbal M.Sc\Data PhD\ChAT homozygous\LFP Figures\For Ai';
print (8, '-depsc2', filename)
