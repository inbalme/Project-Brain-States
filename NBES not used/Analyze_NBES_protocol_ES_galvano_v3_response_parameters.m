%% Analyze NBES protocol ES+galvano v3_response_parameters
% This file was created on 5/10/2015 and updated on 19/1/2016
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
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
peaks_for_xls=[]; peak_for_xls_mean=[]; 
save_flag= 0;
print_flag=0;
files_to_analyze =[1,44,46,48,50,52,56,58,62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,16,22,44,46,48,50,52,56,58,62,72,75];%[1,44,46,48,50,52]; %[1,8,16,22,44,46,48,50,52];
% for fileind=1;
    for fileind=1:length(files_to_analyze) ;
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
%% low-pass filtering below 300Hz                
                for xx=2:3
    for trace= 1:size(data_no_spikes{channel},2)    
            jj=data_no_spikes{channel}(:,trace,xx);
            data_no_spikes{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jj,dt,-1,300,0,0); 
    end
end
                %% Subtract mean trace from data with or without spikes
                meansubtract_start_time = 0; %[sec]
                meansubtract_duration = 3; %[sec]
for x_value = 1:size(data.x_value,2) 
     
                 meansubtract_start_sample = meansubtract_start_time.*sf{channel};
                if meansubtract_start_time==0
                    meansubtract_start_sample = 1;
                end
                meansubtract_end_sample = meansubtract_start_sample+meansubtract_duration.*sf{channel}-1;
                meansubtract_interval = round(meansubtract_start_sample:meansubtract_end_sample);
               
                data_no_spike_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(data_no_spikes{channel}(:,:,x_value),meansubtract_interval);
%                 data_no_DC{channel}(:,:,x_value) = fn_Subtract_Mean(raw_data{channel}(:,:,x_value),meansubtract_interval);    
end
%% Vm histogram - only for spontaneous activity
%  start_time = [0.5,5]; %[sec] %[0,5]
%          duration = 2.5; %[sec] 
%             for t=1:length(start_time);
%              start_sample(:,t) = ceil(start_time(t).*sf{1});
%                 if start_time(t)==0
%                     start_sample(:,t) = 1;
%                 end
%               end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
%               interval(:,t) = start_sample(:,t):end_sample(:,t);             
%               data_vec(:,t)=reshape(data_no_spikes{1}(interval(:,t),:,1),numel(data_no_spikes{1}(interval(:,t),:,1)),1); %take only data from x_value=1;
%               nbin=ceil(range(data_vec(:,1))); %set the range according to the interval before ES
%             end
%             figure
%             hist(data_vec,nbin)
%% Variance plots - two x_values on each plot, for flexible data
x_value = 2:3;

plot_data=data_no_spikes{channel}; %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
% for i=1:length(x_value)
% if ES_flag(x_value(i))==1
% %     plot_data(29910:35010,:,x_value(i)) =nan; %Ignoring the ES artifact
% end
% end
plot_data_mean = mean(plot_data,2);
% plot_data_mean(29910:35010,:,:)=nan;
plot_data_std =  std(data_no_spike_no_DC{channel},0,2);
% plot_data_var=var(plot_data,0,2);
plot_data_var=var(data_no_spike_no_DC{channel},0,2);
plot_data_ff = (-1).*plot_data_var./plot_data_mean; 
% plot_data_CV(29910:35010,:,:)=nan;

% plots for first x-value:
       trace_ind = 1:size(plot_data,2); %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
       l=size(plot_data(:,trace_ind, x_value(1)),1);
        l=l/2-1;
        DC=20.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        trace_to_plot = plot_data(:,trace_ind, x_value(1))+DC ; %DC is added just to make space between the traces
        x2lim=[stim2_X{x_value(1)}(1,1).*dt-0.5,stim2_X{x_value(1)}(1,1).*dt+2];
        
        [Fig1,h1]= fn_Plot_Trace_v2(trace_to_plot, dt, dt, stim1_X{channel}, dt, stim2_X{x_value(1)});       
        hline = findobj(gcf, 'type', 'line'); set(hline(1:end),'color',[0 0 0]);
        y1lim=[get(gca,'ylim')]';
        trace_to_plot = plot_data(:,trace_ind, x_value(2))+DC ;
        [Fig2,h2]= fn_Plot_Trace_v2(trace_to_plot, dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
         hline = findobj(gcf, 'type', 'line'); set(hline(1:end),'color',[0 0 1]);
        set(gca,'ylim',y1lim);
        [Fig3,h3]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
        ylabel('std [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
%         hline = findobj(gcf, 'type', 'line'); set(hline(13),'color',[0 0 0]);  set(hline(14),'color',[0 0 1]);
        [Fig4,h4] = fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
        ylabel('mean Vm [mV]', 'FontSize', 16);          
        [Fig5,h5]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});        
%         [Fig5,h5]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
        ylabel('std [mV]', 'FontSize', 16);  legend('NB-', 'NB+'); legend('boxoff','Location','northeast'); 
%         hline = findobj(gcf, 'type', 'line'); set(hline(13),'color',[0 0 0]);  set(hline(14),'color',[0 0 1]);
        [Fig6,h6] = fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
        ylabel('mean Vm [mV]', 'FontSize', 16);  legend('NB-', 'NB+'); legend('boxoff','Location','northeast');
        [Fig9,h9]= fn_Plot_Trace_v2(plot_data_var(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
        ylabel('Variance [mV^2]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);

                Fig7=figure;
                plot_mean_std_a=plot_data_mean+plot_data_std;
                plot_mean_std_b=plot_data_mean-plot_data_std;
        hold on
        fn_shadedErrorBar([1:size(plot_data,1)].*dt,plot_data_mean(:,:,2),plot_data_std(:,:,2),{'LineWidth',1,'color', color_table(7,:)},1);        
            max_time_axis = ceil(size(plot_data,1).*dt);
            max_data = max(max(plot_mean_std_a(floor(stim2_X{x_value(1)}(1,1)):end,:,:)));
            min_data = min(min(plot_mean_std_b(floor(stim2_X{x_value(1)}(1,1)):end,:,:)));
            five_percent = (max_data-min_data).*0.02;
            min_data_wide = min_data-2.*five_percent;
            max_data_wide = max_data+2.*five_percent;
        fn_shadedErrorBar([1:size(plot_data,1)].*dt,plot_data_mean(:,:,3),plot_data_std(:,:,3),{'LineWidth',1,'color', color_table(1,:)},1);

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
        %plot ES stim
        stim1_Y = ones(size(stim1_X{channel})).*(max_data+five_percent); 
        line(stim1_X{channel}.*dt,stim1_Y,'LineWidth',6,'Color','b');
        %plot whisker stim
        ylim_data=[get(gca,'ylim')]';
        patch_xdata=[stim2_X{x_value(2)}; flipud(stim2_X{x_value(2)})];
        yex=wextend('1D','sym',ylim_data,1);
        l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
        temp_y=wextend('ac','sym',yex,l);
        patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
        patch_cdata=ones(size(patch_xdata));
        patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
        set(gca,'linewidth',1.2)
                ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);

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
        
        ax5 = get(Fig5, 'children'); ax6 = get(Fig6, 'children'); ax7 = get(Fig7, 'children'); ax9 = get(Fig9, 'children'); 
        
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
%% saving figures
if save_flag==1;
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10Hz'  
        saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x2.fig']) 
        print(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x2'],'-dpng','-r600','-opengl') 
        saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x3.fig']) 
        print(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x3'],'-dpng','-r600','-opengl') 
        saveas(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt.fig']) 
        print(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt'],'-dpng','-r600','-opengl') 
        saveas(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3.fig']) 
        print(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3'],'-dpng','-r600','-opengl') 
        saveas(Fig7,['f' num2str(files_to_analyze(fileind)) '_mean+std_x2+3.fig']) 
        print(Fig7,['f' num2str(files_to_analyze(fileind)) '_mean+std_x2+3'],'-dpng','-r600','-opengl') 
        saveas(Fig8,['f' num2str(files_to_analyze(fileind)) '_summary.fig']) 
        print(Fig8,['f' num2str(files_to_analyze(fileind)) '_summary'],'-dpng','-r600','-opengl')     
        saveas(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x2+3_mean-subt.fig']) 
        print(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x2+3_mean-subt'],'-dpng','-r600','-opengl')     
        
        set(ax1,'xlim',x2lim); set(ax2,'xlim',x2lim); set(ax5,'xlim',x2lim); set(ax6,'xlim',x2lim); set(ax7,'xlim',x2lim); set(ax9,'xlim',x2lim); 
        
        saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x2_zoom-in.fig']) 
        print(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x2_zoom-in'],'-dpng','-r600','-opengl') 
        saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x3_zoom-in.fig']) 
        print(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x3_zoom-in'],'-dpng','-r600','-opengl') 
        saveas(Fig5,['f' num2str(files_to_analyze(fileind)) '_Vm std_x2+3_mean-subt_zoom-in.fig']) 
        print(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x2+3_mean-subt_zoom-in'],'-dpng','-r600','-opengl') 
        saveas(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_zoom-in.fig']) 
        print(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x2+3_zoom-in'],'-dpng','-r600','-opengl') 
        saveas(Fig7,['f' num2str(files_to_analyze(fileind)) '_mean+std_x2+3_zoom-in.fig']) 
        print(Fig7,['f' num2str(files_to_analyze(fileind)) '_mean+std_x2+3_zoom-in'],'-dpng','-r600','-opengl') 
        saveas(Fig8,['f' num2str(files_to_analyze(fileind)) '_summary.fig']) 
        print(Fig8,['f' num2str(files_to_analyze(fileind)) '_summary'],'-dpng','-r600','-opengl') 
        saveas(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x2+3_mean-subt_zoom-in.fig']) 
        print(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x2+3_mean-subt_zoom-in'],'-dpng','-r600','-opengl') 
        
end        
close all
% close(Fig1, Fig2, Fig5, Fig6, Fig7, Fig8,Fig9)
%% Response parameters - amplitude, latency and more...
%find peak values and locations in interval following each whisker stim.

int_time=50; %[ms]

    clear data_mean data_std data_mean_smooth mean_peak_val_all mean_peak_loc_all mean_peak_val mean_peak_loc baseline_interval baseline_mean_local baseline_trace...
        mean_peak_half_peak_rise mean_peak_half_peak_fall mean_peak_half_peak_width temp temp_mean half_peak_val response_rise_interval...
        response_fall_interval baseline_mean_global mean_10per_val mean_10per_loc mean_10per_val mean_10per_loc...
        baseline_std_local baseline_std_global baseline_trace_global
t=0;
Y_abs = []; f = []; Y_abs_exp=[];  f_exp = [];
for x_value = 2:3;
    t=t+1; %index for the power spectrum
    data_mean = mean(data_no_spikes{channel}(:,:,x_value),2); %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC 
    data_mean_smooth=sgolayfilt(data_mean, 1, 29);
    data_deriv1=diff(data_mean_smooth);
    data_deriv2=diff(data_mean_smooth,2);

            for stim_num=1:size(stim2_X{x_value},2);
                min_res_del=4; %[ms]
                response_begin=stim2_X{x_value}(1,stim_num)+(min_res_del./1000).*sf{1}; %time onset of whisker stim.+ minimum response latency.
                response_end=response_begin+int_time.*sf{1}./1000;  
                baseline_time=3; %[ms] %local baseline before each peak
                baseline_time_total=100; %[ms] %global baseline for the whole trace
%                 baseline_interval(:,1)=[stim2_X{x_value}(1,stim_num):stim2_X{x_value}(1,stim_num)+(baseline_time./1000).*sf{1}-1];%includes
%                 the stim artifact...
evoked_STD(:,stim_num) = mean(std(data_no_spikes{channel}(response_begin:response_end,:,x_value),0,2)); %mean std across traces (trial-to-trial)
                baseline_interval(:,1)=[stim2_X{x_value}(1,stim_num)-(baseline_time./1000).*sf{1}:stim2_X{x_value}(1,stim_num)-1];
                baseline_std_local(stim_num,:) = std(data_no_spikes{channel}(baseline_interval,:,x_value),0,1); %std along time, within the trace
                baseline_trace(1,:) = mean(data_no_spikes{channel}(baseline_interval,:,x_value),1); %mean along time, within the trace
                baseline_mean_local(stim_num,1)=mean(baseline_trace); %mean across traces
                baseline_trace_global(1,:) = mean(data_no_spikes{channel}(stim2_X{x_value}(1,1)-(baseline_time_total./1000).*sf{1}:stim2_X{x_value}(1,1)-1,:,x_value),1); %mean along time, within the trace
%                 baseline_mean_global(stim_num,1)=mean(data_mean(stim2_X{x_value}(1,1)-(baseline_time_total./1000).*sf{1}:stim2_X{x_value}(1,1)-1)); %mean over 100ms prior to the first stim in the train
                baseline_mean_global(stim_num,1) = mean(baseline_trace_global); %mean across traces
                baseline_std_global(:,1) = std(data_no_spikes{channel}(stim2_X{x_value}(1,1)-(baseline_time_total./1000).*sf{1}:stim2_X{x_value}(1,1)-1,:,x_value),0,2); %std across traces (trial-to-trial)
                baseline_stdM_global(stim_num,1) = mean(baseline_std_global); %mean trial-to-trial variability along time
                minpeak=mean(data_mean_smooth(baseline_interval))+2.*mean(baseline_std_local(stim_num,1)); %threshold criterion for minimal peak height
                [mean_peak_val_all{x_value,stim_num}, mean_peak_loc_all{x_value,stim_num}] =  findpeaks(data_mean_smooth(response_begin:response_end,1),'sortstr','descend');
%                 [mean_peak_val_all{x_value,stim_num}, mean_peak_loc_all{x_value,stim_num}] =  findpeaks(data_mean_smooth(response_begin:response_end,1),'minpeakheight',minpeak);               
%                 mean_peak_val{x_value}(stim_num,1)=mean_peak_val_all{x_value,stim_num}(1);
                mean_peak_loc{x_value}(stim_num,1)=mean_peak_loc_all{x_value,stim_num}(1);
                mean_peak_loc{x_value}(stim_num,1)=mean_peak_loc{x_value}(stim_num,1)+response_begin;
                mean_peak_val{x_value}(stim_num,1)=data_mean(mean_peak_loc{x_value}(stim_num,1));
                 mean_peak_amp_abs{x_value}(stim_num,1)=mean_peak_val{x_value}(stim_num,1)-baseline_mean_global(stim_num,1);
                mean_peak_amp{x_value}(stim_num,1)=mean_peak_val{x_value}(stim_num,1)-baseline_mean_local(stim_num,1);
                
%                 max_deriv2_loc(:,stim_num)=find(data_deriv2(response_begin:response_end,1)==0);
                if mean_peak_amp{x_value}(stim_num,1)<0
%                      mean_peak_amp{x_value}(stim_num,1)=0;
                    mean_peak_amp{x_value}(stim_num,1)=nan;
                    mean_peak_rise_loc{x_value}(stim_num,1)=nan;
                    half_peak_val(stim_num,1)=nan;
                    mean_10per_loc{x_value}(stim_num,1) = nan;  
                    mean_90per_loc{x_value}(stim_num,1) = nan; 
                    mean_peak_half_peak_rise{x_value}(stim_num,1)=nan;
                    mean_peak_half_peak_fall{x_value}(stim_num,1)=nan;
                    mean_peak_half_peak_width{x_value}(stim_num)=nan;
                else
                                               
                %finding half-peak width and 10%peak value and latency from the mean trace:
                half_peak_val(stim_num,1)=baseline_mean_local(stim_num,1)+mean_peak_amp{1,x_value}(stim_num)./2;
                mean_10per_val{x_value}(stim_num,1)=baseline_mean_local(stim_num,1)+mean_peak_amp{1,x_value}(stim_num).*0.1;
                mean_90per_val{x_value}(stim_num,1)=baseline_mean_local(stim_num,1)+mean_peak_amp{1,x_value}(stim_num).*0.9;
                mean_peak_rise_loc{x_value}(stim_num,1)=find(data_deriv1(response_begin:response_end,1)==max(data_deriv1(response_begin:response_end,1)));
                mean_peak_rise_loc{x_value}(stim_num,1)=mean_peak_rise_loc{x_value}(stim_num,1)+response_begin;
                mean_peak_rise_val{x_value}(stim_num,1)=data_mean(mean_peak_rise_loc{x_value}(stim_num,1),1);
                response_rise_interval=[response_begin:mean_peak_loc{x_value}(stim_num,1)];
                response_fall_interval=[mean_peak_loc{x_value}(stim_num,1):stim2_X{x_value}(1,stim_num)+(100./1000).*sf{1}];        
%                  if isempty(find(data_mean(response_rise_interval)>= mean_10per_val{x_value}(stim_num,1),1,'first'))
%                     mean_10per_loc{x_value}(stim_num,1) = nan;                   
%                  else                    
                    mean_10per_loc{x_value}(stim_num,1)=response_begin+find(data_mean(response_rise_interval)>= mean_10per_val{x_value}(stim_num,1),1,'first');
%                  end
%                 if isempty(find(data_mean(response_rise_interval)>= mean_10per_val{x_value}(stim_num,1),1,'first'))
%                     mean_10per_loc{x_value}(stim_num,1) = nan;                   
%                 else                     
                    mean_90per_loc{x_value}(stim_num,1)=response_begin+find(data_mean(response_rise_interval)>= mean_90per_val{x_value}(stim_num,1),1,'first');
%                 end
                mean_peak_half_peak_rise{x_value}(stim_num,1)=response_begin+find(data_mean(response_rise_interval)>=half_peak_val(stim_num,1),1,'first');
                    if isempty(find(data_mean(response_fall_interval)<=half_peak_val(stim_num,1),1,'first'))
                        mean_peak_half_peak_fall{x_value}(stim_num)=mean_peak_loc{x_value}(stim_num)+find(data_mean(response_fall_interval)==min(data_mean(response_fall_interval)),1,'first');
                    else                    
                        mean_peak_half_peak_fall{x_value}(stim_num,1)=mean_peak_loc{x_value}(stim_num,1)+find(data_mean(response_fall_interval)<=half_peak_val(stim_num,1),1,'first');
                    end
                mean_peak_half_peak_width{x_value}(stim_num)= (mean_peak_half_peak_fall{x_value}(stim_num)-mean_peak_half_peak_rise{x_value}(stim_num)).*dt; 
                end
 % parameters from single traces:               
                %try keti's way of finding max response in single traces.
%                loop on all traces, findpeaks on smoothed data without spikes
 %               [mean_peak_val_all{x_value,stim_num}, mean_peak_loc_all{x_value,stim_num}] =  findpeaks(data_mean_smooth(response_begin:response_end,1),'sortstr','descend');                        
%taking the value in the raw data traces corresponding to the location of the peak from the mean trace
                peaks(fileind).stim_num(stim_num).val(:,x_value) = [data_no_spikes{channel}(mean_peak_loc{x_value}(stim_num),:,x_value)]';
                peaks(fileind).stim_num(stim_num).amp(:,x_value) = peaks(fileind).stim_num(stim_num).val(:,x_value)-baseline_trace';
    %parameters for single traces            
                temp(:,1)=fileind(ones(size(data_no_spikes{1},2),1));
                temp(:,2)=stim_num(ones(size(data_no_spikes{1},2),1));
                temp(:,3)=x_value(ones(size(data_no_spikes{1},2),1));
                temp(:,4)=peaks(fileind).stim_num(stim_num).val(:,x_value);
                temp(:,5)=peaks(fileind).stim_num(stim_num).amp(:,x_value);
                peaks_for_xls=[peaks_for_xls; temp];
            end             
%%    % verification of detected peaks in the raw data:
%      [Fig6,h6] = fn_Plot_Trace_v2([data_mean_smooth], dt, dt, stim1_X{channel}, dt, stim2_X{x_value}); %add to plot: raw_data{channel}(:,:,x_value)
%         ylabel('mean Vm [mV]', 'FontSize', 16);  legend('NB-', 'NB+'); legend('boxoff','Location','northeast');
%                     hold on 
% %                         for trace = 1:size(data_no_spikes{1},2)         
% %                              h1=plot([1:size(data_no_spikes{channel},1)].*dt, data_no_spikes{channel}(:,trace,x_value),'b');
% %                              h2=scatter(mean_peak_loc{x_value}.*dt,data_no_spikes{channel}(mean_peak_loc{x_value},trace,x_value),'r','fill'); %mark the mean peaks on the single traces
%                              h2=scatter(mean_peak_loc{x_value}.*dt,data_mean(mean_peak_loc{x_value}),'r','fill'); %mark the mean peaks on the mean traces
% %                              h2=scatter(mean_peak_loc{x_value}.*dt,data_mean_smooth(mean_peak_loc{x_value}),'r', 'fill'); %mark the peaks on the smoothed mean trace
% %                              h3=scatter(mean_peak_half_peak_rise{x_value}.*dt,data_no_spikes{channel}(mean_peak_half_peak_rise{x_value},trace,x_value),'g','fill');
% %                              h4=scatter(mean_peak_half_peak_fall{x_value}.*dt,data_no_spikes{channel}(mean_peak_half_peak_fall{x_value},trace,x_value),'g','fill');                       
%                              h3=scatter(mean_peak_half_peak_rise{x_value}.*dt,data_mean(mean_peak_half_peak_rise{x_value}),'g','fill');
%                              h4=scatter(mean_peak_half_peak_fall{x_value}.*dt,data_mean(mean_peak_half_peak_fall{x_value}),'g','fill');                       
%                                 h5=scatter(stim2_X{x_value}(1,:).*dt,baseline_mean_local,'c','fill');
%                                  h6=scatter(stim2_X{x_value}(1,:).*dt,baseline_mean_global,'y','fill');
% %                                  h7=scatter(mean_peak_rise_loc{x_value}.*dt,data_mean_smooth(mean_peak_rise_loc{x_value}),'r', 'fill'); %mark the peaks on the mean trace
%                                 h7=scatter(mean_10per_loc{x_value}.*dt,data_mean(mean_10per_loc{x_value}),'r','fill'); %mark the mean peaks on the mean traces
%                                  pause
% %                              delete(h1)
%                              delete(h2)
%                              delete(h3)
%                              delete(h4)
%                              delete(h5)
%                              delete(h6)
%                              delete(h7)
%                              figure(Fig6) 
% %                          end
%                       hold off
%                         close (gcf)
%% power spectrum of sensory response for adaptation F1
fstim=Param.facade(7); %train stim frequency
nstim=Param.facade(6); %number of stim in the train
stim1_duration=1000./fstim./1000; %1000ms/galvano frequency converted to [sec]
duration = stim1_duration.*nstim; %stim duration multiplied by the number of stim
DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = []; interval = []; expected_mat=[];

start_sample = stim2_X{x_value}(1,1); %[sec]
end_sample = start_sample+duration.*sf{channel}-1;
interval = start_sample:end_sample; %taking the interval of the sensory stim train
spec_mat = data_no_spikes{channel}(interval,:,x_value);
expected_mat=repmat(spec_mat(1:stim1_duration.*sf{channel},:),nstim,1); %replicating the response for the 1st stim in the train
 [spec_mat_noDC] = fn_subtract_mean_trace(spec_mat);
 [expected_mat_noDC] = fn_subtract_mean_trace(expected_mat);
% DC= mean(spec_mat,1);
% l=size(spec_mat,1);
% l=l/2-1;
% DC=wextend('addrow', 'per', DC, l);
% spec_mat_noDC = spec_mat-DC;
[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf{channel});
[Y_abs_exp(:,:,t),f_exp(:,t)] = fn_Amp_Spectrum(expected_mat_noDC,sf{channel});
max_res(1,t)=max(mean(Y_abs(f(:,t)>(fstim-1) &f(:,t)<(fstim+1),:,t),2));
max_res_exp(1,t)=max(mean(Y_abs_exp(f_exp(:,t)>(fstim-1) &f_exp(:,t)<(fstim+1) ,:,t),2));
F1(1,t)=max_res(1,t)./max_res_exp(1,t);

figure(10) %(fileind)
      hold on
        plot(f(:,t),mean(Y_abs(:,:,t),2),'color', color_table(t,:)) 
        plot(f_exp(:,t),mean(Y_abs_exp(:,:,t),2),'color', color_table(t+2,:)) 
      hold off
        xlim([0 50]); ylim([-0.01 25])

        set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2)
        xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
        ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
        title(['Response Power file ', num2str(files_to_analyze(fileind))],'FontSize', 12); 
        leg=legend('NB-','expected NB-','NB+','expected NB+');
        set(leg,'box','off', 'FontSize', 12);
%%
%Parameters from mean trace                    peak_10per_lat    
peaks(fileind).fname = fname; 
peaks(fileind).baseline_time_local=baseline_time; %baseline before each whisker stim.
peaks(fileind).baseline_time_global=baseline_time_total; %global baseline for the whole mean trace
peaks(fileind).baseline_local(:,x_value)=baseline_mean_local; 
peaks(fileind).baseline_global(:,x_value)=baseline_mean_global; 
peaks(fileind).baseline_global_std(:,x_value)=baseline_stdM_global;
peaks(fileind).evoked_STD(:,x_value)=evoked_STD;
peaks(fileind).val(:,x_value) = mean_peak_val{1,x_value};
peaks(fileind).amp(:,x_value) = mean_peak_amp{1,x_value};
peaks(fileind).amp_abs(:,x_value) = mean_peak_amp_abs{1,x_value};
peaks(fileind).loc(:,x_value) = mean_peak_loc{1,x_value};
peaks(fileind).latency(:,x_value)=(mean_peak_loc{1,x_value}-(stim2_X{x_value}(1,:))').*dt;
peaks(fileind).per10_val(:,x_value) = mean_10per_val{x_value}(:,:);
peaks(fileind).per10_lat(:,x_value)= (mean_10per_loc{x_value}-(stim2_X{x_value}(1,:))').*dt;
peaks(fileind).per90_val(:,x_value) = mean_90per_val{x_value}(:,:);
peaks(fileind).per90_lat(:,x_value)= (mean_90per_loc{x_value}-(stim2_X{x_value}(1,:))').*dt;
peaks(fileind).half_width(:,x_value) = mean_peak_half_peak_width{x_value};
peaks(fileind).adapt_amp(:,x_value)=mean(mean_peak_amp{1,x_value}(8:10))./mean_peak_amp{1,x_value}(1);
peaks(fileind).fft_max_res(:,x_value)=max_res(1,t);
peaks(fileind).fft_max_res_exp(:,x_value)=max_res_exp(1,t);
peaks(fileind).F1(:,x_value)=F1(1,t);

    
                temp_mean(:,1)=fileind(ones(size(peaks(fileind).val(:,x_value),1),1));
                temp_mean(:,2)=[1:size(peaks(fileind).val(:,x_value),1)];
                temp_mean(:,3)=x_value(ones(size(peaks(fileind).val(:,x_value),1),1));
                temp_mean(:,4)=peaks(fileind).val(:,x_value);
                temp_mean(:,5)=peaks(fileind).amp(:,x_value);
                temp_mean(:,6)=peaks(fileind).amp_abs(:,x_value);
                temp_mean(:,7)=peaks(fileind).latency(:,x_value);
                temp_mean(:,8)=peaks(fileind).half_width(:,x_value);
                temp_mean(:,9)=peaks(fileind).per10_lat(:,x_value);
                temp_mean(:,10)=peaks(fileind).per10_lat(:,x_value);
                temp_mean(:,11)=peaks(fileind).baseline_local(:,x_value);   
                temp_mean(:,12)=peaks(fileind).evoked_STD(:,x_value); 
                
                peak_for_xls_mean=[peak_for_xls_mean; temp_mean];
                
%                 peak_adapt_amp{x_value}(fileind,1)=peaks(fileind).adapt_amp(:,x_value);
%                 peak_F1{x_value}(fileind,1)=peaks(fileind).F1(:,x_value);
%                 peak_baseline_global{x_value}(fileind,1)=peaks(fileind).baseline_global(1,x_value);
%                 peak_baseline_global_std{x_value}(fileind,1)=peaks(fileind).baseline_global_std(1,x_value);
         for stim_num=1:10;
            peak_val{x_value}(fileind,stim_num)=peaks(fileind).val(stim_num,x_value);
            peak_amp{x_value}(fileind,stim_num)=peaks(fileind).amp(stim_num,x_value);
            peak_amp_abs{x_value}(fileind,stim_num)=peaks(fileind).amp_abs(stim_num,x_value);
            peak_lat{x_value}(fileind,stim_num)=peaks(fileind).latency(stim_num,x_value);
            peak_half_width{x_value}(fileind,stim_num)=peaks(fileind).half_width(stim_num,x_value);
            peak_per10_lat{x_value}(fileind,stim_num)=peaks(fileind).per10_lat(stim_num,x_value);
            peak_per90_lat{x_value}(fileind,stim_num)=peaks(fileind).per90_lat(stim_num,x_value);
            peak_baseline_local{x_value}(fileind,stim_num)=peaks(fileind).baseline_local(stim_num,x_value);
            peak_evoked_STD{x_value}(fileind,stim_num)=peaks(fileind).evoked_STD(stim_num,x_value);
            peak_adapt_amp{x_value}(fileind,stim_num)=peaks(fileind).adapt_amp(1,x_value);
            peak_F1{x_value}(fileind,stim_num)=peaks(fileind).F1(1,x_value);
            peak_baseline_global{x_value}(fileind,stim_num)=peaks(fileind).baseline_global(stim_num,x_value);
            peak_baseline_global_std{x_value}(fileind,stim_num)=peaks(fileind).baseline_global_std(stim_num,x_value);
           
         end
         
end
close all
clear data_no_spike_no_DC
    end %end of files loop
 
    %% Statistics for peaks
%t-tests between cells      
    for stim_num=1; %1:10;
        %normalize peak values and perform ttest
%         peak_val_noES=peak_for_xls_mean([find(peak_for_xls_mean(:,2)==stim_num & peak_for_xls_mean(:,3)==2)],4);
%         peak_val_ES=peak_for_xls_mean([find(peak_for_xls_mean(:,2)==stim_num & peak_for_xls_mean(:,3)==3)],4);
        
        peaks_stat(stim_num).val=[peak_val{2}(:,stim_num),  peak_val{3}(:,stim_num)];
        norm_peak_val_noES= peak_val{2}(:,stim_num)./ peak_val{2}(:,stim_num);
        norm_peak_val_ES= peak_val{3}(:,stim_num)./ peak_val{2}(:,stim_num);
        peaks_stat(stim_num).norm_val=[norm_peak_val_noES, norm_peak_val_ES];
        peaks_stat(stim_num).norm_val_m=mean(peaks_stat(stim_num).norm_val,1);
        peaks_stat(stim_num).norm_val_std=std(peaks_stat(stim_num).norm_val,0,1);
        %testing for normal distribution
        diff_val= peak_val{3}(:,stim_num)- peak_val{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_val, peaks_stat(stim_num).lillietest_p_mean_val] = lillietest(diff_val);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_val, peaks_stat(stim_num).ttest_p_norm_val]= ttest(peaks_stat(stim_num).norm_val(:,1),peaks_stat(stim_num).norm_val(:,2));

        clear norm_peak_val_noES norm_peak_val_ES diff_val
        %normalize peak amp and perform ttest
%         peak_amp_noES=peak_for_xls_mean([find(peak_for_xls_mean(:,2)==stim_num & peak_for_xls_mean(:,3)==2)],5);
%         peak_amp_ES=peak_for_xls_mean([find(peak_for_xls_mean(:,2)==stim_num & peak_for_xls_mean(:,3)==3)],5);
        peaks_stat(stim_num).amp=[peak_amp{2}(:,stim_num),  peak_amp{3}(:,stim_num)];
        norm_peak_amp_noES= peak_amp{2}(:,stim_num)./ peak_amp{2}(:,stim_num);
        norm_peak_amp_ES= peak_amp{3}(:,stim_num)./ peak_amp{2}(:,stim_num);
        peaks_stat(stim_num).norm_amp=[norm_peak_amp_noES, norm_peak_amp_ES];
        peaks_stat(stim_num).norm_amp_m=mean(peaks_stat(stim_num).norm_amp,1);
        peaks_stat(stim_num).norm_amp_std=std(peaks_stat(stim_num).norm_amp,0,1);
        %testing for normal distribution
        diff_amp= peak_amp{3}(:,stim_num)- peak_amp{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_amp, peaks_stat(stim_num).lillietest_p_mean_amp] = lillietest(diff_amp);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_amp, peaks_stat(stim_num).ttest_p_norm_amp]= ttest(peaks_stat(stim_num).norm_amp(:,1),peaks_stat(stim_num).norm_amp(:,2));

        clear norm_peak_amp_noES norm_peak_amp_ES diff_amp
        
        %normalize peak amp_abs and perform ttest
        peaks_stat(stim_num).amp_abs=[peak_amp_abs{2}(:,stim_num),  peak_amp_abs{3}(:,stim_num)];
        norm_peak_amp_noES= peak_amp_abs{2}(:,stim_num)./ peak_amp_abs{2}(:,stim_num);
        norm_peak_amp_ES= peak_amp_abs{3}(:,stim_num)./ peak_amp_abs{2}(:,stim_num);
        peaks_stat(stim_num).norm_amp_abs=[norm_peak_amp_noES, norm_peak_amp_ES];
        peaks_stat(stim_num).norm_amp_abs_m=mean(peaks_stat(stim_num).norm_amp_abs,1);
        peaks_stat(stim_num).norm_amp_abs_std=std(peaks_stat(stim_num).norm_amp_abs,0,1);
        %testing for normal distribution
        diff_amp_abs= peak_amp_abs{3}(:,stim_num)- peak_amp_abs{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_amp_abs, peaks_stat(stim_num).lillietest_p_mean_amp_abs] = lillietest(diff_amp_abs);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_amp, peaks_stat(stim_num).ttest_p_norm_amp]= ttest(peaks_stat(stim_num).norm_amp(:,1),peaks_stat(stim_num).norm_amp(:,2));

        clear norm_peak_amp_noES norm_peak_amp_ES diff_amp_abs
        
        %normalize peak latency and perform ttest
        peaks_stat(stim_num).latency=[peak_lat{2}(:,stim_num),  peak_lat{3}(:,stim_num)];
        norm_peak_latency_noES= peak_lat{2}(:,stim_num)./peak_lat{2}(:,stim_num);
        norm_peak_latency_ES= peak_lat{3}(:,stim_num)./peak_lat{2}(:,stim_num);
        peaks_stat(stim_num).norm_latency=[norm_peak_latency_noES, norm_peak_latency_ES];
        peaks_stat(stim_num).norm_latency_m=mean(peaks_stat(stim_num).norm_latency,1);
        peaks_stat(stim_num).norm_latency_std=std(peaks_stat(stim_num).norm_latency,0,1);
        %testing for normal distribution
        diff_latency= peak_lat{3}(:,stim_num)- peak_lat{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_latency, peaks_stat(stim_num).lillietest_p_mean_latency] = lillietest(diff_latency);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_latency, peaks_stat(stim_num).ttest_p_norm_latency]= ttest(peaks_stat(stim_num).norm_latency(:,1),peaks_stat(stim_num).norm_latency(:,2));

        clear norm_peak_latency_noES norm_peak_latency_ES diff_latency
        
        %normalize peak half_width and perform ttest
        peaks_stat(stim_num).half_width=[peak_half_width{2}(:,stim_num),  peak_half_width{3}(:,stim_num)];
        norm_peak_half_width_noES= peak_half_width{2}(:,stim_num)./ peak_half_width{2}(:,stim_num);
        norm_peak_half_width_ES= peak_half_width{3}(:,stim_num)./ peak_half_width{2}(:,stim_num);
        peaks_stat(stim_num).norm_half_width=[norm_peak_half_width_noES, norm_peak_half_width_ES];
        peaks_stat(stim_num).norm_half_width_m=mean(peaks_stat(stim_num).norm_half_width,1);
        peaks_stat(stim_num).norm_half_width_std=std(peaks_stat(stim_num).norm_half_width,0,1);
        %testing for normal distribution
        diff_half_width= peak_half_width{3}(:,stim_num)- peak_half_width{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_half_width, peaks_stat(stim_num).lillietest_p_mean_half_width] = lillietest(diff_half_width);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_half_width, peaks_stat(stim_num).ttest_p_norm_half_width]= ttest(peaks_stat(stim_num).norm_half_width(:,1),peaks_stat(stim_num).norm_half_width(:,2));

clear norm_peak_half_width_noES norm_peak_half_width_ES diff_half_width

%normalize peak 10 percent latency and perform ttest
        peaks_stat(stim_num).per10_lat=[peak_per10_lat{2}(:,stim_num),  peak_per10_lat{3}(:,stim_num)];
        norm_peak_per10_lat_noES= peak_per10_lat{2}(:,stim_num)./ peak_per10_lat{2}(:,stim_num);
        norm_peak_per10_lat_ES= peak_per10_lat{3}(:,stim_num)./ peak_per10_lat{2}(:,stim_num);
        peaks_stat(stim_num).norm_per10_lat=[norm_peak_per10_lat_noES, norm_peak_per10_lat_ES];
        peaks_stat(stim_num).norm_per10_lat_m=mean(peaks_stat(stim_num).norm_per10_lat,1);
        peaks_stat(stim_num).norm_per10_lat_std=std(peaks_stat(stim_num).norm_per10_lat,0,1);
        %testing for normal distribution
        diff_per10_lat= peak_per10_lat{3}(:,stim_num)- peak_per10_lat{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_per10_lat, peaks_stat(stim_num).lillietest_p_mean_per10_lat] = lillietest(diff_per10_lat);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_per10_lat, peaks_stat(stim_num).ttest_p_norm_per10_lat]= ttest(peaks_stat(stim_num).norm_per10_lat(:,1),peaks_stat(stim_num).norm_per10_lat(:,2));

        clear norm_peak_per10_lat_noES norm_peak_per10_lat_ES diff_per10_lat
        
        %normalize peak 90 percent latency and perform ttest
        peaks_stat(stim_num).per90_lat=[peak_per90_lat{2}(:,stim_num),  peak_per90_lat{3}(:,stim_num)];
        norm_peak_per90_lat_noES= peak_per90_lat{2}(:,stim_num)./ peak_per90_lat{2}(:,stim_num);
        norm_peak_per90_lat_ES= peak_per90_lat{3}(:,stim_num)./ peak_per90_lat{2}(:,stim_num);
        peaks_stat(stim_num).norm_per90_lat=[norm_peak_per90_lat_noES, norm_peak_per90_lat_ES];
        peaks_stat(stim_num).norm_per90_lat_m=mean(peaks_stat(stim_num).norm_per90_lat,1);
        peaks_stat(stim_num).norm_per90_lat_std=std(peaks_stat(stim_num).norm_per90_lat,0,1);
        %testing for normal distribution
        diff_per90_lat= peak_per90_lat{3}(:,stim_num)- peak_per90_lat{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_per90_lat, peaks_stat(stim_num).lillietest_p_mean_per90_lat] = lillietest(diff_per90_lat);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_per90_lat, peaks_stat(stim_num).ttest_p_norm_per90_lat]= ttest(peaks_stat(stim_num).norm_per90_lat(:,1),peaks_stat(stim_num).norm_per90_lat(:,2));

        clear norm_peak_per90_lat_noES norm_peak_per90_lat_ES diff_per90_lat
        
                %normalize local baseline and perform ttest
         peaks_stat(stim_num).baseline_local=[peak_baseline_local{2}(:,stim_num),  peak_baseline_local{3}(:,stim_num)];
        norm_peak_baseline_local_noES= peak_baseline_local{2}(:,stim_num)./ peak_baseline_local{2}(:,stim_num);
        norm_peak_baseline_local_ES= peak_baseline_local{3}(:,stim_num)./ peak_baseline_local{2}(:,stim_num);
        peaks_stat(stim_num).norm_baseline_local=[norm_peak_baseline_local_noES, norm_peak_baseline_local_ES];
        peaks_stat(stim_num).norm_baseline_local_m=mean(peaks_stat(stim_num).norm_baseline_local,1);
        peaks_stat(stim_num).norm_baseline_local_std=std(peaks_stat(stim_num).norm_baseline_local,0,1);
        %testing for normal distribution
        diff_baseline_local= peak_baseline_local{3}(:,stim_num)- peak_baseline_local{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_baseline_local, peaks_stat(stim_num).lillietest_p_mean_baseline_local] = lillietest(diff_baseline_local);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_baseline_local, peaks_stat(stim_num).ttest_p_norm_baseline_local]= ttest(peaks_stat(stim_num).norm_baseline_local(:,1),peaks_stat(stim_num).norm_baseline_local(:,2));

        clear norm_peak_baseline_local_noES norm_peak_baseline_local_ES diff_baseline_local
        
         %normalize global baseline and perform ttest
         peaks_stat(stim_num).baseline_global=[peak_baseline_global{2}(:,stim_num),  peak_baseline_global{3}(:,stim_num)];
        norm_peak_baseline_global_noES= peak_baseline_global{2}(:,stim_num)./ peak_baseline_global{2}(:,stim_num);
        norm_peak_baseline_global_ES= peak_baseline_global{3}(:,stim_num)./ peak_baseline_global{2}(:,stim_num);
        peaks_stat(stim_num).norm_baseline_global=[norm_peak_baseline_global_noES, norm_peak_baseline_global_ES];
        peaks_stat(stim_num).norm_baseline_global_m=mean(peaks_stat(stim_num).norm_baseline_global,1);
        peaks_stat(stim_num).norm_baseline_global_std=std(peaks_stat(stim_num).norm_baseline_global,0,1);
        %testing for normal distribution
        diff_baseline_global= peak_baseline_global{3}(:,stim_num)- peak_baseline_global{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_baseline_global, peaks_stat(stim_num).lillietest_p_mean_baseline_global] = lillietest(diff_baseline_global);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_baseline_global, peaks_stat(stim_num).ttest_p_norm_baseline_global]= ttest(peaks_stat(stim_num).norm_baseline_global(:,1),peaks_stat(stim_num).norm_baseline_global(:,2));

        clear norm_peak_baseline_global_noES norm_peak_baseline_global_ES diff_baseline_global
        
            %normalize evoked STD (mean trial-to-trial STD) and perform ttest
        peaks_stat(stim_num).evoked_STD=[peak_evoked_STD{2}(:,stim_num),  peak_evoked_STD{3}(:,stim_num)];
        norm_peak_evoked_STD_noES= peak_evoked_STD{2}(:,stim_num)./ peak_evoked_STD{2}(:,stim_num);
        norm_peak_evoked_STD_ES= peak_evoked_STD{3}(:,stim_num)./ peak_evoked_STD{2}(:,stim_num);
        peaks_stat(stim_num).norm_evoked_STD=[norm_peak_evoked_STD_noES, norm_peak_evoked_STD_ES];
        peaks_stat(stim_num).norm_evoked_STD_m=mean(peaks_stat(stim_num).norm_evoked_STD,1);
        peaks_stat(stim_num).norm_evoked_STD_std=std(peaks_stat(stim_num).norm_evoked_STD,0,1);
        %testing for normal distribution
        diff_evoked_STD= peak_evoked_STD{3}(:,stim_num)- peak_evoked_STD{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_evoked_STD, peaks_stat(stim_num).lillietest_p_mean_evoked_STD] = lillietest(diff_evoked_STD);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_evoked_STD, peaks_stat(stim_num).ttest_p_norm_evoked_STD]= ttest(peaks_stat(stim_num).norm_evoked_STD(:,1),peaks_stat(stim_num).norm_evoked_STD(:,2));

        clear norm_peak_evoked_STD_noES norm_peak_evoked_STD_ES diff_evoked_STD
        
%normalize adapt_amp and perform ttest
        peaks_stat(stim_num).adapt_amp=[peak_adapt_amp{2}(:,stim_num),  peak_adapt_amp{3}(:,stim_num)];
        norm_peak_adapt_amp_noES= peak_adapt_amp{2}(:,stim_num)./ peak_adapt_amp{2}(:,stim_num);
        norm_peak_adapt_amp_ES= peak_adapt_amp{3}(:,stim_num)./ peak_adapt_amp{2}(:,stim_num);
        peaks_stat(stim_num).norm_adapt_amp=[norm_peak_adapt_amp_noES, norm_peak_adapt_amp_ES];
        peaks_stat(stim_num).norm_adapt_amp_m=mean(peaks_stat(stim_num).norm_adapt_amp,1);
        peaks_stat(stim_num).norm_adapt_amp_std=std(peaks_stat(stim_num).norm_adapt_amp,0,1);
        %testing for normal distribution
        diff_adapt_amp= peak_adapt_amp{3}(:,stim_num)- peak_adapt_amp{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_adapt_amp, peaks_stat(stim_num).lillietest_p_mean_adapt_amp] = lillietest(diff_adapt_amp);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_adapt_amp, peaks_stat(stim_num).ttest_p_norm_adapt_amp]= ttest(peaks_stat(stim_num).norm_adapt_amp(:,1),peaks_stat(stim_num).norm_adapt_amp(:,2));
 
        clear norm_peak_adapt_amp_noES norm_peak_adapt_amp_ES diff_adapt_amp
        
        %normalize F1 and perform ttest
        peaks_stat(stim_num).F1=[peak_F1{2}(:,stim_num),  peak_F1{3}(:,stim_num)];
        norm_peak_F1_noES= peak_F1{2}(:,stim_num)./ peak_F1{2}(:,stim_num);
        norm_peak_F1_ES= peak_F1{3}(:,stim_num)./ peak_F1{2}(:,stim_num);
        peaks_stat(stim_num).norm_F1=[norm_peak_F1_noES, norm_peak_F1_ES];
        peaks_stat(stim_num).norm_F1_m=mean(peaks_stat(stim_num).norm_F1,1);
        peaks_stat(stim_num).norm_F1_std=std(peaks_stat(stim_num).norm_F1,0,1);
        %testing for normal distribution
        diff_F1= peak_F1{3}(:,stim_num)- peak_F1{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_F1, peaks_stat(stim_num).lillietest_p_mean_F1] = lillietest(diff_F1);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_F1, peaks_stat(stim_num).ttest_p_norm_F1]= ttest(peaks_stat(stim_num).norm_F1(:,1),peaks_stat(stim_num).norm_F1(:,2));

        clear norm_peak_F1_noES norm_peak_F1_ES diff_F1
         
  %normalize global baseline std and perform ttest
        peaks_stat(stim_num).baseline_global_std=[peak_baseline_global_std{2}(:,stim_num),  peak_baseline_global_std{3}(:,stim_num)];
        norm_peak_baseline_global_std_noES= peak_baseline_global_std{2}(:,stim_num)./ peak_baseline_global_std{2}(:,stim_num);
        norm_peak_baseline_global_std_ES= peak_baseline_global_std{3}(:,stim_num)./ peak_baseline_global_std{2}(:,stim_num);
        peaks_stat(stim_num).norm_baseline_global_std=[norm_peak_baseline_global_std_noES, norm_peak_baseline_global_std_ES];
        peaks_stat(stim_num).norm_baseline_global_std_m=mean(peaks_stat(stim_num).norm_baseline_global_std,1);
        peaks_stat(stim_num).norm_baseline_global_std_std=std(peaks_stat(stim_num).norm_baseline_global_std,0,1);
        %testing for normal distribution
        diff_baseline_global_std= peak_baseline_global_std{3}(:,stim_num)- peak_baseline_global_std{2}(:,stim_num);
        [peaks_stat(stim_num).lillietest_h_mean_baseline_global_std, peaks_stat(stim_num).lillietest_p_mean_baseline_global_std] = lillietest(diff_baseline_global_std);
        %paired ttest 
        [peaks_stat(stim_num).ttest_h_norm_baseline_global_std, peaks_stat(stim_num).ttest_p_norm_baseline_global_std]= ttest(peaks_stat(stim_num).norm_baseline_global_std(:,1),peaks_stat(stim_num).norm_baseline_global_std(:,2));
        
        
        clear norm_peak_baseline_global_std_noES norm_peak_baseline_global_std_ES diff_baseline_global_std
       
    end
    %% Plot response amplitude Vs. local baseline
%     h1=figure;
%     cell_num=1;
    for cell_num=1:length(files_to_analyze)
        linfit_ES_Off=polyfit(peaks(cell_num).baseline_local(:,2),peaks(cell_num).amp(:,2),1);
        linval_ES_Off=polyval(linfit_ES_Off,peaks(cell_num).baseline_local(:,2));
        linfit_ES_On=polyfit(peaks(cell_num).baseline_local(:,3),peaks(cell_num).amp(:,3),1);
        linval_ES_On=polyval(linfit_ES_On,peaks(cell_num).baseline_local(:,3));
       

        h(cell_num)=figure;
        color(cell_num,:)=rand(1,3);
 hold on
scatter(peaks(cell_num).baseline_local(:,2),peaks(cell_num).amp(:,2),30,'markeredgecolor',color(cell_num,:))
scatter(peaks(cell_num).baseline_local(:,3),peaks(cell_num).amp(:,3),30,'markerfacecolor',color(cell_num,:),'markeredgecolor',color(cell_num,:))
plot(peaks(cell_num).baseline_local(:,2),linval_ES_Off,'k-')
plot(peaks(cell_num).baseline_local(:,3),linval_ES_On,'b-')

hold off
        xlabel('Baseline' ,'FontSize', 16);  ylabel('Amplitude','FontSize', 16); 
        title(['Response amplitude Vs. Baseline, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
%         set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line]);
    end   
    
    %% Plot response amplitude Vs. local baseline - all cells on same plot
     h=figure;
    for cell_num=1:length(files_to_analyze)    
        color(cell_num,:)=rand(1,3);
 hold on
scatter(peaks(cell_num).baseline_local(:,2),peaks(cell_num).amp(:,2),30,'markeredgecolor',color(cell_num,:))
scatter(peaks(cell_num).baseline_local(:,3),peaks(cell_num).amp(:,3),30,'markerfacecolor',color(cell_num,:),'markeredgecolor',color(cell_num,:))
hold off
        xlabel('Baseline' ,'FontSize', 16);  ylabel('Amplitude','FontSize', 16); 
        title(['Response amplitude Vs. Baseline, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
%         set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line]);
    end   
      %% Plot response amplitude to first stim Vs. local baseline - all cells on same plot
    baseline_vec=[]; amp_vec=[];
    for cell_num=1:length(files_to_analyze)   
        baseline_vec = [baseline_vec; peaks(cell_num).baseline_local(:,2:3)];
        amp_vec = [amp_vec; peaks(cell_num).amp(:,2:3)];
    end
    if print_flag==1;
         h=figure;
 hold on
scatter(baseline_vec(:,1),amp_vec(:,1),50,'markeredgecolor',color(cell_num,:))
scatter(baseline_vec(:,2),amp_vec(:,2),50,'markerfacecolor',color(cell_num,:),'markeredgecolor',color(cell_num,:))
hold off
        xlabel('Baseline' ,'FontSize', 16);  ylabel('Amplitude','FontSize', 16); 
        title(['Response amplitude Vs. Baseline, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
%         set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line]);
    end     
%% Plot parameters along the train stim

color_table_rand = rand(length(files_to_analyze),3);
    if print_flag==1;
h1=figure;
        figure(h1);     
        hold on
        errorbar([1:size(peak_amp{3},2)], mean(peak_amp{3}./peak_amp{2},1),std(peak_amp{3}./peak_amp{2},0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        line([0;12],[1;1],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean Response Amplitude [mV]' ,'FontSize', 16);    
        title(['Mean Response Amplitude - local baseline,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);     
        
         h2=figure;
        figure(h2);     
        hold on
        errorbar([1:size(peak_amp_abs{3},2)],mean(peak_amp_abs{3}./peak_amp_abs{2},1),std(peak_amp_abs{3}./peak_amp_abs{2},0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        line([0;12],[1;1],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean Response Amplitude [mV]' ,'FontSize', 16);     
        title(['Mean Response Amplitude - global baseline, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);     


         h3=figure;
        figure(h3);     
        hold on
        errorbar([1:size(peak_lat{3},2)],mean(peak_lat{3}./peak_lat{2},1),std(peak_lat{3}./peak_lat{2},0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        line([0;12],[1;1],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean Onset Latency' ,'FontSize', 16);     
        title(['Mean Onset Latency, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);     


        h4=figure;
        figure(h4);     
        hold on
        errorbar([1:size(peak_per10_lat{3},2)],mean(peak_per10_lat{3}./peak_per10_lat{2},1),std(peak_per10_lat{3}./peak_per10_lat{2},0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        line([0;12],[1;1],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean Peak Latency' ,'FontSize', 16);     
        title(['Mean 10% Peak Latency, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);    

         h5=figure;
        figure(h5);     
        hold on
        errorbar([1:size(peak_half_width{3},2)],mean(peak_half_width{3}./peak_half_width{2},1),std(peak_half_width{3}./peak_half_width{2},0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        line([0;12],[1;1],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean half width' ,'FontSize', 16);     
        title(['Mean half width, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);    

         h6=figure;
        figure(h6);     
        hold on
        errorbar([1:size(peak_val{3},2)],mean(peak_val{3}./peak_lat{2},1),std(peak_val{3}./peak_val{2},0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        line([0;12],[1;1],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean peak value' ,'FontSize', 16);     
        title(['Mean peak value,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
 
         h7=figure;
        figure(h7);     
        hold on
        errorbar([1:size(peak_baseline_local{3},2)],mean(peak_baseline_local{3}./peak_baseline_local{2},1),std(peak_baseline_local{3}./peak_baseline_local{2},0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        line([0;12],[1;1],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean local baseline [mV]' ,'FontSize', 16);    
        title(['Mean local baseline, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);         
    %% Scatter plots for first peak parameters
stim_num=1;
%
min_line=floor(min(min(peaks_stat(stim_num).amp)));
max_line=ceil(max(max(peaks_stat(stim_num).amp)));
 g1=figure;
 hold on
scatter(peaks_stat(stim_num).amp(:,1),peaks_stat(stim_num).amp(:,2),180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30);  ylabel('NB+','FontSize', 30); 
        title(['Mean Response Local Amplitude [mV], stim 1, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_amp)] ,'FontSize', 20);   
        set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line],'FontSize', 30);

%
min_line=floor(min(min(peaks_stat(stim_num).amp_abs)));
max_line=ceil(max(max(peaks_stat(stim_num).amp_abs)));
 g2=figure;
 hold on
scatter(peaks_stat(stim_num).amp_abs(:,1),peaks_stat(stim_num).amp_abs(:,2),180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30); 
        title(['Mean Response Global Amplitude [mV], stim 1,n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_amp_abs)] ,'FontSize', 20);   
       set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30); 

%
min_line=floor(min(min(peaks_stat(stim_num).half_width.*1000)));
max_line=ceil(max(max(peaks_stat(stim_num).half_width.*1000)));
 g3=figure;
 hold on
scatter(peaks_stat(stim_num).half_width(:,1).*1000,peaks_stat(stim_num).half_width(:,2).*1000,180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30); 
        title(['Mean Half width [mV], stim 1, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_half_width)] ,'FontSize', 20);   
        set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);

        %
min_line=floor(min(min(peaks_stat(stim_num).per10_lat.*1000)));
max_line=ceil(max(max(peaks_stat(stim_num).per10_lat.*1000)));
 g4=figure;
 hold on
scatter(peaks_stat(stim_num).per10_lat(:,1).*1000,peaks_stat(stim_num).per10_lat(:,2).*1000,180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30); 
        title(['Mean Onset Latency [mV], stim 1,n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_per10_lat)] ,'FontSize', 20);   
       set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30); 
  
%
min_line=floor(min(min(peaks_stat(stim_num).per10_lat.*1000)));
max_line=ceil(max(max(peaks_stat(stim_num).per10_lat.*1000)));
 g5=figure;
 hold on
scatter(peaks_stat(stim_num).per10_lat(:,1).*1000,peaks_stat(stim_num).per10_lat(:,2).*1000,180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30); 
        title(['Mean 10% peak latency [mV], stim 1, n=' num2str(length(files_to_analyze)), ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_per10_lat)] ,'FontSize', 20);   
        set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);
   
       
%      
min_line=floor(min(min(peaks_stat(stim_num).adapt_amp)));
max_line=ceil(max(max(peaks_stat(stim_num).adapt_amp)));
        g6=figure;
 hold on
scatter(peaks_stat(stim_num).adapt_amp(:,1),peaks_stat(stim_num).adapt_amp(:,2),180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30); 
        title(['Adaptation amplitude ratio,n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_adapt_amp)] ,'FontSize', 20);  
        set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);
      
 %      
min_line=floor(min(min(peaks_stat(stim_num).F1)));
max_line=ceil(max(max(peaks_stat(stim_num).F1)));
        g7=figure;
 hold on
scatter(peaks_stat(stim_num).F1(:,1),peaks_stat(stim_num).F1(:,2),180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30);  ylabel('NB+','FontSize', 30);  
        title(['Adaptation power ratio, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_F1)] ,'FontSize', 20);   
         set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);
   
%
min_line=floor(min(min(peaks_stat(stim_num).baseline_global)));
max_line=ceil(max(max(peaks_stat(stim_num).baseline_global)));
 g8=figure;
 hold on
scatter(peaks_stat(stim_num).baseline_global(:,1),peaks_stat(stim_num).baseline_global(:,2),180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30);  
        title(['Global Baseline, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_baseline_global)] ,'FontSize', 20); 
        set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);

%
min_line=floor(min(min(peaks_stat(stim_num).val)));
max_line=ceil(max(max(peaks_stat(stim_num).val)));
 g9=figure;
 hold on
scatter(peaks_stat(stim_num).val(:,1),peaks_stat(stim_num).val(:,2),180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30); ylabel('NB+','FontSize', 30);  
        title(['Peak Value, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_val)] ,'FontSize', 20); 
        set(gca, 'xlim', [min_line max_line], 'ylim', [min_line max_line],'FontSize', 30);
%% save figures
if save_flag==1
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10Hz'
print(h1,'Train_response_amp','-dpng','-r600','-opengl')
saveas(h2,'Train_response_amp_global.fig') 
print(h2,'Train_response_amp_global','-dpng','-r600','-opengl')
saveas(h3,'Train_onset_latency.fig') 
print(h3,'Train_onset_latency','-dpng','-r600','-opengl')
saveas(h4,'Train_10%peak_latency.fig') 
print(h4,'Train_10%peak_latency','-dpng','-r600','-opengl')
saveas(h5,'Train_half-width.fig') 
print(h5,'Train_half-width','-dpng','-r600','-opengl')
saveas(h6,'Train_peak_value.fig') 
print(h6,'Train_peak_value','-dpng','-r600','-opengl')
saveas(h7,'Train_local_baseline.fig') 
print(h7,'Train_local_baseline','-dpng','-r600','-opengl')
% 
% 
saveas(g1,'Stim1_local_amp.fig')         
print(g1,'Stim1_local_amp','-dpng','-r600','-opengl') 
saveas(g2,'Stim1_global_amp.fig') 
print(g2,'Stim1_global_amp','-dpng','-r600','-opengl') 
saveas(g3,'Stim1_half_width.fig') 
print(g3,'Stim1_half_width','-dpng','-r600','-opengl') 
saveas(g4,'Stim1_onset_latency.fig')
print(g4,'Stim1_onset_latency','-dpng','-r600','-opengl') 
saveas(g5,'Stim1_10%peak_latency.fig') 
print(g5,'Stim1_10%peak_latency','-dpng','-r600','-opengl') 
saveas(g6,'adaptation_amp_ratio.fig') 
print(g6,'adaptation_amp_ratio','-dpng','-r600','-opengl') 
saveas(g7,'adaptation_power_ratio.fig') 
print(g7,'adaptation_power_ratio','-dpng','-r600','-opengl') 
saveas(g8,'global_baseline.fig') 
print(g8,'global_baseline','-dpng','-r600','-opengl') 
saveas(g9,'peak_value.fig') 
print(g9,'peak_value','-dpng','-r600','-opengl') 
end
    end
%% t-tests within cell    
%     for fileind=1:length(files_to_analyze) ;%           
%                for stim_num=1:max(peak_for_xls_mean(peak_for_xls_mean(:,1)==fileind,2));
%                     [peaks(fileind).val_ttest_h(stim_num), peaks(fileind).val_ttest_p(stim_num)]=ttest(peaks(fileind).stim_num(stim_num).val(:,2),peaks(fileind).stim_num(stim_num).val(:,3));
%                     [peaks(fileind).amp_ttest_h(stim_num), peaks(fileind).amp_ttest_p(stim_num)]=ttest(peaks(fileind).stim_num(stim_num).amp(:,2),peaks(fileind).stim_num(stim_num).amp(:,3));                    
%                end%                
%     end
%% Plot parameters along the train stim
% color_table_rand = rand(length(files_to_analyze),3);
%  h1=figure;
%  for fileind=1:length(files_to_analyze) ;
%         figure(h1);     
%         hold on
%         plot(1:size(peaks(fileind).amp,1),peaks(fileind).amp(:,2),'--', 'color', color_table(fileind,:),'LineWidth',1.5,'marker','o','MarkerFaceColor', color_table(fileind,:),'MarkerEdgeColor','k','MarkerSize',7)
%         plot(1:size(peaks(fileind).amp,1),peaks(fileind).amp(:,3), 'color', color_table(fileind,:),'LineWidth',1.5,'marker','o','MarkerFaceColor', color_table(fileind,:),'MarkerEdgeColor','k','MarkerSize',7)     
%         hold off
%         xlabel('Stimulus number' ,'FontSize', 16);
%         ylabel('Mean Response Amplitude [mV]' ,'FontSize', 16);     
%  end
    %% Plot power spectrum
% start_time = [0.5,4.5]; %[sec]
% duration = 2; %[sec]
% x_value = 1;
% Y_abs = []; f = [];
% 
% for t=1:length(start_time);
% DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = []; interval = []; 
% 
% start_sample = start_time(t).*sf{channel};
% if start_time(t)==0
%     start_sample = 1;
% end
% end_sample = start_sample+duration.*sf{channel}-1;
% interval = start_sample:end_sample;
% spec_mat = raw_data{channel}(interval,:,x_value);
% DC= mean(spec_mat,1);
% l=size(spec_mat,1);
% l=l/2-1;
% DC=wextend('addrow', 'per', DC, l);
% spec_mat_noDC = spec_mat-DC;
% 
% [Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf{channel});
% end
% % Plot power spectrum with shaded error bar
% % figure
% %  for t=1:length(start_time);
% %       hold on
% %     fn_shadedErrorBar(f(:,t),mean(Y_abs(:,:,t),2),std(Y_abs(:,:,t),0,2),{'LineWidth',1,'color', color_table(t,:)},1);
% %     hold off
% %     xlim([0 35]); ylim([-0.01 0.04])
% %  end
% 
%  % Plot power spectrum without shaded error bar
% figure
%  for t=1:length(start_time);
%       hold on
%         plot(f(:,t),mean(Y_abs(:,:,t),2),'color', color_table(t,:)) 
%       hold off
%         xlim([0 35]); ylim([-0.01 0.04])
% %         title('')
%         set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2)
%         xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
%         ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
%  end
 %%
% filename =  'SNR_plot_12_cells'; %'F62P6X6+9_zero_traces+std+mean_zoom-in';
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\n13';
% % cd 'd:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations';
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH';
% print (6, '-deps', filename)
% filename='response_parameters'; 
% save(filename, 'files_to_analyze', 'peaks', 'peaks_stat', 'peaks_for_xls', 'peak_for_xls_mean')
