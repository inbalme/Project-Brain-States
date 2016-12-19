%% Analyze_adaptation_ongoing
% This file was created on 3/4/2016 
%This file is used for the analysis of files created with extract_NBES_Data_v2

%       protocol 5 (NBES+galvnao (5 x-values))
%       protocol 6 (ChAT intracellular current injections light+galvnao (15 x-values))
%       protocol 8 - ChAT intracellular: light only (3 x-values): 1) train
%       20 Hz, 10ms ondur, 50 pulses. 2) single pulse 5 sec ondur. 3)single pulse 10 ms
%       protocol 9 - IC: I-V curve
%       protocol 10 (NBES+galvnao train+test (3 x-values))
%                 11 - NBES: 3 currents+galvano train+test (6 x-values)
%% for opening workspace saved 
clear all; close all
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
peaks_for_xls=[]; peak_for_xls_mean=[]; 
save_flag=1;
files_to_analyze =19; %[8,10,11,12,14,15,16,18,19,20,22,36,37,44,46,48,52,56,58,62,72,75];%[1,44,46,48,50,52]; %[1,8,16,22,44,46,48,50,52];
% for fileind=1;
    for fileind=1:length(files_to_analyze) ;
            clearvars -except cc_stat cc_spont cc_evoked  files_to_analyze fileind files cc_spont_for_xls_mean cc_evoked_for_xls_mean cc lags cc_shuffled cc_mean cc_shuff_sub_mean save_flag
    close all
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
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
                
                w_stim_freq(fileind) = Param.facade(7);
               %% Subtract mean trace from data with or without spikes
                meansubtract_start_time = 0; %[sec]
                meansubtract_duration = 3; %[sec]
% for x_value = 1:size(data.x_value,2)      
        x_value = 2;
                 meansubtract_start_sample = meansubtract_start_time.*sf{channel};
                if meansubtract_start_time==0
                    meansubtract_start_sample = 1;
                end
                meansubtract_end_sample = meansubtract_start_sample+meansubtract_duration.*sf{channel}-1;
                meansubtract_interval = round(meansubtract_start_sample:meansubtract_end_sample);
               
                data_no_spike_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(data_no_spikes{channel}(:,:,x_value),meansubtract_interval);
%                 data_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(raw_data{channel}(:,:,x_value),meansubtract_interval);    
% end   %x_value loop
%%
clear  start_time duration start_time start_sample end_sample interval interval_mat x y patch_xdata patch_ydata yex ylim_data
coeffs=[];
trace_type=1;
% U=unique(x_value);
% trace_type=length(U); %for trace_type=1 plot the xcorr before and after the ES according to the two start times. for trace_type=2 plot the xcorr from the same time point in the 2 different x-values.
         start_time = [0.5,8]; %[sec] %[0,5]
         start_time(2) = stim2_X{x_value}(1,end).*dt+0.1; %start 0.5 sec after the last sensory stim
         duration = 2; %[sec] 
            for t=1:length(start_time);
             start_sample(:,t) = ceil(start_time(t).*sf{1});
                if start_time(t)==0
                    start_sample(:,t) = 1;
                end
              end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
              interval(:,t) = start_sample(:,t):end_sample(:,t);
            end 
 %% Variance plots - two x_values on each plot, for flexible data
x_value = 2;
for t=1:2;
plot_data(:,:,t)=data_no_spikes{channel}( interval(:,t),:,x_value); %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
plot_data_no_DC(:,:,t)=data_no_spike_no_DC{channel}( interval(:,t),:,x_value);
end
plot_data_mean = mean(plot_data,2);
plot_data_std =  std(plot_data_no_DC,0,2);
plot_data_var=var(plot_data_no_DC,0,2);

% plots for first interval:
       trace_ind = 1:size(plot_data,2); %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
       l=size(plot_data,1);
        l=l/2-1;
        DC=20.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        trace_to_plot = plot_data(:,trace_ind, 1)+DC ; %DC is added just to make space between the traces        
        [Fig1,h1]= fn_Plot_Trace_v2_no_stim(trace_to_plot, dt, dt, stim1_X{channel}, dt, stim2_X{x_value});       
        hline = findobj(gcf, 'type', 'line'); set(hline(1:end),'color',[0 0 0]);
        y1lim=[get(gca,'ylim')]';
        trace_to_plot = plot_data(:,trace_ind, 2)+DC ;
        [Fig2,h2]= fn_Plot_Trace_v2_no_stim(trace_to_plot, dt, dt, stim1_X{channel}, dt, stim2_X{x_value});
         hline = findobj(gcf, 'type', 'line'); set(hline(1:end),'color',[0 0 1]);
        set(gca,'ylim',y1lim);
        [Fig3,h3]= fn_Plot_Trace_v2_no_stim(plot_data_std(:,:), dt, dt, stim1_X{channel}, dt, stim2_X{x_value});
        ylabel('std [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
%         hline = findobj(gcf, 'type', 'line'); set(hline(13),'color',[0 0 0]);  set(hline(14),'color',[0 0 1]);
        [Fig4,h4] = fn_Plot_Trace_v2_no_stim(plot_data_mean(:,:), dt, dt, stim1_X{channel}, dt, stim2_X{x_value});
        ylabel('mean Vm [mV]', 'FontSize', 16);          
        [Fig5,h5]= fn_Plot_Trace_v2_no_stim(plot_data_std(:,:), dt, dt, stim1_X{channel}, dt, stim2_X{x_value});        
%         [Fig5,h5]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
        ylabel('std [mV]', 'FontSize', 16);  legend('Before adaptation', 'After adaptation'); legend('boxoff','Location','northeast'); 
%         hline = findobj(gcf, 'type', 'line'); set(hline(13),'color',[0 0 0]);  set(hline(14),'color',[0 0 1]);
        [Fig6,h6] = fn_Plot_Trace_v2_no_stim(plot_data_mean(:,:), dt, dt, stim1_X{channel}, dt, stim2_X{x_value});
        ylabel('mean Vm [mV]', 'FontSize', 16);  legend('Before adaptation', 'After adaptation'); legend('boxoff','Location','northeast');

                Fig7=figure;
                plot_mean_std_a=plot_data_mean+plot_data_std;
                plot_mean_std_b=plot_data_mean-plot_data_std;
        hold on
        fn_shadedErrorBar([1:size(plot_data,1)].*dt,plot_data_mean(:,:,1),plot_data_std(:,:,1),{'LineWidth',1,'color', color_table(7,:)},1);        
            max_time_axis = ceil(size(plot_data,1).*dt);
            max_data = max(max(plot_mean_std_a));
            min_data = min(min(plot_mean_std_b));
            five_percent = (max_data-min_data).*0.02;
            min_data_wide = min_data-2.*five_percent;
            max_data_wide = max_data+2.*five_percent;
        fn_shadedErrorBar([1:size(plot_data,1)].*dt,plot_data_mean(:,:,2),plot_data_std(:,:,2),{'LineWidth',1,'color', color_table(1,:)},1);

 x1limits = [0 max_time_axis];
 x1ticks = [0 0.5.*max_time_axis max_time_axis];
 if abs(min_data_wide) < 1
     min_data_wide = (floor(min_data_wide.*10))./10;
 else
     min_data_wide = floor(min_data_wide);
 end
 if abs(max_data_wide) < 1 
      max_data_wide = (ceil(max_data_wide.*10))./10;
 else
      max_data_wide = ceil(max_data_wide);
 end
 y1limits = [min_data_wide max_data_wide]; %[-0.3 0.5]; [0.5 1.5]
 y1ticks =  [min_data_wide max_data_wide]; %[-0.3 0 0.5]; %y1limits; [0.5 1 1.5]
 set( gca, 'xlim', x1limits, 'ylim', y1limits);

         hold off
        
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
        
        ax5 = get(Fig5, 'children'); ax6 = get(Fig6, 'children'); ax7 = get(Fig7, 'children'); 
        
        Fig8 = figure;
        set(gcf,'color','w');
        set(gcf,'DefaultAxesFontSize',12);
        set(gcf,'DefaultAxesFontName','helvetica');
        set(gcf, 'PaperType', 'A4');
        set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

        ax_copy1 = copyobj(ax1,Fig8); % ax1 to new fig
        set(ax_copy1(1),'position',pos1(1,:)) % Set its position  

        ax_copy2 = copyobj(ax2,Fig8); % ax2 to new fig
        set(ax_copy2(1),'position',pos2(1,:)) % Set its position  

        ax_copy3 = copyobj(ax3,Fig8); % ax3 to new fig
        set(ax_copy3(1),'position',pos3(1,:)) % Set its position  
        
        ax_copy4 = copyobj(ax4,Fig8); % ax3 to new fig
        set(ax_copy4(1),'position',pos4(1,:)) % Set its position  
%         
        close(Fig3); close(Fig4);  
        %% Power spectrum before and after adaptation
       Y_abs = []; f = []; Y_abs_exp=[];  f_exp = [];
       for t=1:2;
       DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = [];
 spec_mat =plot_data(:,:,t);
 [spec_mat_noDC] = fn_subtract_mean_trace(spec_mat);
 [Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf{channel});
       end 
 figure(10) %(fileind)
      hold on
        plot(f(:,1),mean(Y_abs(:,:,1),2),'k','linewidth',1.5) 
        plot(f(:,2),mean(Y_abs(:,:,2),2),'r','linewidth',1.5) 

      hold off
       xlim([0 50]); %ylim([-0.01 25])

        set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2)
        xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
        ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
        title(['file ', num2str(files_to_analyze(fileind)),' ',num2str(w_stim_freq),'Hz'],'FontSize', 12); 
        leg=legend('Before adaptation','After adaptation');
        set(leg,'box','off', 'FontSize', 12);
       
                %% analysis of Vm-LFP crosscorrelation
%             if files_to_analyze(fileind)>36
% 
%                         %bandpass filtering to remove noise from LFP channel
%                 for xx=1:3
%                     for trace= 1:size(raw_data{3},2)    
%                             jj=raw_data{3}(:,trace,xx);
%                             raw_data{3}(:,trace,xx)=bandPass_fft_IL_NEW2016(jj,dt,49,51,1,0);
%                     end
%                 end
%             end
    end