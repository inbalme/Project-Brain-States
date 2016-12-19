%% Analyze NBES protocol ES+galvano v2
% This file was created on 7/5/2014. updated on 30/3/15
%This file is used for the analysis of files created with extract_NBES_Data_v2

%       protocol 5 (NBES+galvnao (5 x-values))
%       protocol 6 (ChAT intracellular current injections light+galvnao (15 x-values))
%       protocol 8 - ChAT intracellular: light only (3 x-values): 1) train
%       20 Hz, 10ms ondur, 50 pulses. 2) single pulse 5 sec ondur. 3)single pulse 10 ms
%       protocol 9 - IC: I-V curve
%       protocol 10 (NBES+galvnao train+test (3 x-values))
%                 11 - NBES: 3 currents+galvano train+test (6 x-values)
%% for opening workspace saved 
clear all
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
%  cd 'D:\Inbal M.Sc\Data PhD\DoubleNB\Extracted Data';
% load f1_workspace
% load f8_workspace
load f22_workspace

%%
%     name = 'f22_workspace';
%     cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
% %     cd 'D:\Inbal M.Sc\Data PhD\DoubleNB\Extracted Data';
%     save( name)
%% 
% close all
clear all
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2


% for fileind = 1:length(files)
    fileind =40;
   channel = 1; % V1 - 1, I1 - 2, V2 - 3, I2 - 4
   
    fname = files(fileind).extracted_name;        
    path = files(fileind).extracted_path;
    cd(path)
    load(fname) 
    
    original_header = 2;
    trace =[]; %22:29; %[]; %230:260; %find(Param.orig_headers==original_header); %[]; % if trace is empty the defaulat is to take all traces into analysis

cd 'D:\Dropbox\Inbal M.Sc\MATLAB\Project Brain States'
    
    color_table = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

    data_BP_laser = [];
    
ES_flag = zeros(1,size(data.x_value,2)); %for each x-value, ES_flag is1 if there was ES activation
galvano_flag = zeros(1,size(data.x_value,2)); %for each x-value, galvano_flag is1 if there was galvano activation
 
protocol = Param.protocol;

if isempty(protocol)
    
else
    switch protocol
            
        case 5
             ES_flag(1,:) = [1 0 1 0 1];
             galvano_flag(1,:) = [0 1 1 1 1];
        case 6
            if length(ES_flag(1,:))==9
                ES_flag(1,:) = [1 1 1 0 0 0 1 1 1];
                galvano_flag(1,:) = [0 0 0 1 1 1 1 1 1];
            else
             ES_flag(1,:) = [1 1 1 0 0 0 1 1 1 0 0 0 1 1 1];
             galvano_flag(1,:) = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1];   
            end
        case 8
             ES_flag(1,:) = [1 1 1];
             galvano_flag(1,:) = [0 0 0];  
        case 9
             ES_flag(1,:) = [0 0 0 0 0 0 0];
             galvano_flag(1,:) = [0 0 0 0 0 0 0];  
        case 10
            ES_flag(1,:) = [1 0 1];
            galvano_flag(1,:) = [0 1 1 ];  
        case 11
            ES_flag(1,:) = [1 1 1 1 1 1];
            galvano_flag(1,:) = [0 0 0 1 1 1 ];  
    end
end
         
    locations_x_galvano = [];
    locations_x_ES = [];

    for x_value = 1:size(data.x_value,2) 
      
                if isempty(trace)
                    trace = 1:size(data.x_value(x_value).Vm,2);
                end

                 switch channel
                     case 1
                         raw_data(:,1:length(trace),x_value) = data.x_value(x_value).Vm(:,trace);
                         sf = Param.sf_Vm;
                     case 2
                         raw_data(:,1:length(trace),x_value) = data.x_value(x_value).I1(:,trace); 
                         sf = Param.sf_I1;                       
                      case 3
                         raw_data(:,1:length(trace),x_value) = data.x_value(x_value).V2(:,trace); 
                         sf = Param.sf_V2;
                      case 4
                         raw_data(:,1:length(trace),x_value) = data.x_value(x_value).I2(:,trace); 
                         sf = Param.sf_I2;
                 end

                    if Param.EEG_flag==1;
                    data_EEG(:,1:length(trace),x_value) = data.x_value(x_value).EEG(:,trace);
                         sf_EEG = Param.sf_EEG;
                         dt_EEG = 1/sf_EEG;
                    end
                                            
                
                dt=1/sf; %[sec]
                time_axis(:,x_value) = (1:size(raw_data,1))*dt;
                sf_airpuff = Param.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                sf_galvano = Param.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;


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
                    galvano_sub_mean = data.x_value(x_value).galvano(:,1) - mean(data.x_value(x_value).galvano(:,1));
                    galvano_threshold = abs(mean(galvano_sub_mean)) + 3.*abs(std(galvano_sub_mean));
                    galvano_vec(abs(galvano_sub_mean) > galvano_threshold, x_value)=1;   %turning the galvano trace into binary
                 end

             ES_begin_time = Param. ES_param_val(2);
             ES_duration_time = Param. ES_param_val(5)*(1/Param. ES_param_val(4)); %number of pulses times 1/f
             ES_end_time = ES_begin_time+ES_duration_time;
             ES_vec(:,1) = zeros(length(time_axis),1);
             ES_vec(ES_begin_time*sf:ES_end_time*sf,1)=1;
             
             locations_x_ES(1,:) = ES_begin_time*sf; %arranging the galvano begin and end locations in one variable for plotting
             locations_x_ES(2,:) = ES_end_time*sf; 
             
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
                    locations_x_galvano{x_value}(1,:) = galvano_begin{x_value}(:,1); %arranging the galvano begin and end locations in one variable for plotting
                    locations_x_galvano{x_value}(2,:) = galvano_end{x_value}(:,1);            

                    end

                                  data_mean(:,x_value) = mean(raw_data(:,:,x_value),2);
                                  data_std(:,x_value) = std(raw_data(:,:,x_value),0,2);
                                  data_var(:,x_value) = var(raw_data(:,:,x_value),0,2);
                                  data_CV(:,x_value) = (-1).*data_std(:,x_value)./data_mean(:,x_value); 

                                  %Smoothing the data:
                smooth_data(:,:,x_value) = sgolayfilt(raw_data(:,:,x_value), 1, 29);
                smooth_data_mean(:,x_value) = mean(smooth_data(:,:,x_value),2);
                smooth_data_std(:,x_value) = std(smooth_data(:,:,x_value),0,2);
                smooth_data_var(:,x_value) = var(smooth_data(:,:,x_value),0,2);
                   
                
%                 %removing spikes
%                 for traces=1:size(raw_data(:,:,x_value),2)
%                     [code(:,traces,x_value) data_no_spikes(:,traces,x_value)] = fn_RemoveSpikesMO2(raw_data(:,traces,x_value)', dt);
%                 end
%                data_no_spikes_mean(:,x_value) = mean(data_no_spikes(:,:,x_value),2);
%                data_no_spikes_std(:,x_value) = std(data_no_spikes(:,:,x_value),0,2);
%                data_no_spikes_var(:,x_value) = var(data_no_spikes(:,:,x_value),0,2);
%                data_no_spikes_CV(:,x_value) = (-1).*data_no_spikes_std(:,x_value)./data_no_spikes_mean(:,x_value);
    
    end
    
     %  Stim vectors
 for x_value = 1:size(data.x_value,2) 
         if isempty(locations_x_ES)
             stim1_X = [];
         else
            stim1_X = locations_x_ES(:,:);
         end

        if isempty(locations_x_galvano)
            stim2_X{x_value} = [];
        else
            stim2_X{x_value} = locations_x_galvano{x_value}(:,:);
        end
 end
%% Subtrace mean trace from data without spikes
meansubtract_start_time = 0; %[sec]
meansubtract_duration = 3; %[sec]

 meansubtract_start_sample = meansubtract_start_time.*sf;
if meansubtract_start_time==0
    meansubtract_start_sample = 1;
end
meansubtract_end_sample = meansubtract_start_sample+meansubtract_duration.*sf-1;
meansubtract_interval = round(meansubtract_start_sample:meansubtract_end_sample);

  for x_value = 1:size(data.x_value,2) 
    data_no_spike_no_DC(:,:,x_value) = fn_Subtract_Mean(data_no_spikes(:,:,x_value),meansubtract_interval);
    data_no_DC(:,:,x_value) = fn_Subtract_Mean(raw_data(:,:,x_value),meansubtract_interval);
  end
 
 %% calculating meanVm and stdVm according to Golshani
 % spontaneous trial-to-trial variability
 x_value = 1;
 start_time =[0.5 4.5]; %[0.5 4.5]; %[4:0.1:4.9]; %[5:0.05:5.45]; %[sec]
 duration = 1; %1; %0.05; %[sec]

data_mat=[]; mean_mat=[]; std_mat=[]; std_mat_2=[]; interval_mat=[]; average_mean=[]; average_std=[]; average_std_2=[];

data_mat=data_no_spike_no_DC;

for t=1:length(start_time);
DC=[]; start_sample = []; end_sample = []; interval = []; 

start_sample = start_time(t).*sf;
if start_time(t)==0
    start_sample = 1;
end
end_sample = start_sample+duration.*sf-1;
interval = round(start_sample:end_sample);
interval_mat(:,:,t) = data_mat(interval,:,x_value);
mean_mat(t,:) = mean(interval_mat(:,:,t),1);
std_mat(t,:) = std(interval_mat(:,:,t),0,1);
std_mat_2(t,:) = std(interval_mat(:,:,t),0,2);
end

average_mean = mean(mean_mat,2);
average_std = mean(std_mat,2);
average_std_2 = mean(std_mat_2,2);

 %%
 x_value = 2:3;
figure
  for t=1:length(x_value);
     figure(gcf)
      hold on
    shadedErrorBar((1:size(raw_data,1))'.*dt,mean(raw_data(:,:,x_value(t)),2),std(raw_data(:,:,x_value(t)),0,2),{'LineWidth',1,'color', color_table(t,:)},1);
%     xlim([0 35]); ylim([-0.1 1.5])
  end

%% Variance plots - two x_values on each plot, for flexible data
x_value = 2:3;

plot_data=raw_data; %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
for i=1:length(x_value)
if ES_flag(x_value(i))==1
    plot_data(29990:35010,:,x_value(i)) =nan; %Ignoring the ES artifact
end
end
plot_data_mean = mean(plot_data,2);
plot_data_mean(29990:35010,:,:)=nan;
plot_data_std =  std(plot_data,0,2);
plot_data_std(29990:35010,:,:)=nan;
plot_data_CV = (-1).*plot_data_std./plot_data_mean; 
plot_data_CV(29990:35010,:,:)=nan;

% plots for first x-value:
       trace_ind = [1:size(plot_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
       l=size(plot_data(:,trace_ind, x_value(1)),1);
        l=l/2-1;
        DC=20.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        trace_to_plot = plot_data(:,trace_ind, x_value(1))+DC ; %DC is added just to make space between the traces
        [Fig1,h1]= fn_Plot_Trace_v2(trace_to_plot, dt, dt, stim1_X, dt_galvano, stim2_X{x_value(1)});
        trace_to_plot = plot_data(:,trace_ind, x_value(2))+DC ;
        [Fig2,h2]= fn_Plot_Trace_v2(trace_to_plot, dt, dt, stim1_X, dt_galvano, stim2_X{x_value(2)});
        
        [Fig3,h3]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X, dt_galvano, stim2_X{x_value(2)});
        [Fig4,h4] = fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt, stim1_X, dt_galvano, stim2_X{x_value(2)});
%         [Fig5,h5]= fn_Plot_Trace_v2(data_no_spikes_CV(:,x_value), dt, dt, stim1_X, dt_galvano, stim2_X{x_value(2)});

     
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
        
%         close Figure 1 Figure 2 Figure 3 Figure 4
%% Variance plot - one x_value at a time, for flexible data
       x_value = 2;
       plot_data=raw_data; %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
        plot_data_mean = mean(plot_data,2);
        plot_data_std =  std(plot_data,0,2);
        plot_data_CV = (-1).*plot_data_std./plot_data_mean; 

       trace_ind = [1:size(plot_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
       l=size(plot_data(:,trace_ind, x_value),1);
        l=l/2-1;
        DC=10.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        trace_to_plot = plot_data(:,trace_ind, x_value)+DC ;
        [Fig1,h1]= fn_Plot_Trace_v2(trace_to_plot, dt, dt, stim1_X, dt_galvano, stim2_X{x_value});
        [Fig2,h2]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X, dt_galvano, stim2_X{x_value});
        [Fig3,h3] = fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt, stim1_X, dt_galvano, stim2_X{x_value});



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
%         close Figure 1 Figure 2 Figure 3

%% Plot power spectrum
start_time = [0.5,4.5]; %[sec]
duration = 2; %[sec]
x_value = 1;
Y_abs = []; f = [];

for t=1:length(start_time);
DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = []; interval = []; 

start_sample = start_time(t).*sf;
if start_time(t)==0
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

[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf);
end
% Plot power spectrum with shaded error bar
% figure
%  for t=1:length(start_time);
%       hold on
%     fn_shadedErrorBar(f(:,t),mean(Y_abs(:,:,t),2),std(Y_abs(:,:,t),0,2),{'LineWidth',1,'color', color_table(t,:)},1);
%     hold off
%     xlim([0 35]); ylim([-0.01 0.04])
%  end

 % Plot power spectrum without shaded error bar
figure
 for t=1:length(start_time);
      hold on
        plot(f(:,t),mean(Y_abs(:,:,t),2),'color', color_table(t,:)) 
      hold off
        xlim([0 35]); ylim([-0.01 0.04])
%         title('')
        set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2)
        xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
        ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
 end


%% 
filename =  '2015-07-09-004_x3'; %'F62P6X6+9_zero_traces+std+mean_zoom-in';

cd 'd:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH';

print (4, '-depsc2', filename)
