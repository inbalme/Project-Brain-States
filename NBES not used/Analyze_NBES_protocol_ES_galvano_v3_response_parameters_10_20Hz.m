%% Analyze NBES protocol ES+galvano v3_response_parameters
% This file was created on 5/10/2015 and updated on 19/1/2016. OLD VERSION. NOT USED!
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
files_to_analyze =[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,16,22,44,46,48,50,52,56,58,62,72,75];%[1,44,46,48,50,52]; %[1,8,16,22,44,46,48,50,52];
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
%  start_time = [0.5,5.5]; %[sec] %[0,5]
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
% %     plot_data(29990:35010,:,x_value(i)) =nan; %Ignoring the ES artifact
% end
% end
plot_data_mean = mean(plot_data,2);
% plot_data_mean(29990:35010,:,:)=nan;
plot_data_std =  std(data_no_spike_no_DC{channel},0,2);
% plot_data_var=var(plot_data,0,2);
plot_data_var=var(data_no_spike_no_DC{channel},0,2);
plot_data_ff = (-1).*plot_data_var./plot_data_mean; 
% plot_data_CV(29990:35010,:,:)=nan;

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
        response_fall_interval baseline_mean_global mean_90per_val mean_90per_loc mean_10per_val mean_10per_loc...
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
                baseline_interval(:,1)=[stim2_X{x_value}(1,stim_num)-(baseline_time./1000).*sf{1}:stim2_X{x_value}(1,stim_num)-1];
                baseline_std_local(stim_num,:) = std(data_no_spikes{channel}(baseline_interval,:,x_value),0,1); %std along time, within the trace
                baseline_trace(1,:) = mean(data_no_spikes{channel}(baseline_interval,:,x_value),1); %mean along time, within the trace
                baseline_mean_local(stim_num,1)=mean(baseline_trace); %mean across traces
                baseline_trace_global(1,:) = mean(data_no_spikes{channel}(stim2_X{x_value}(1,1)-(baseline_time_total./1000).*sf{1}:stim2_X{x_value}(1,1)-1,:,x_value),1); %mean along time, within the trace
%                 baseline_mean_global(stim_num,1)=mean(data_mean(stim2_X{x_value}(1,1)-(baseline_time_total./1000).*sf{1}:stim2_X{x_value}(1,1)-1)); %mean over 100ms prior to the first stim in the train
                baseline_mean_global(stim_num,1) = mean(baseline_trace_global); %mean across traces
                baseline_std_global(:,1) = std(data_no_spikes{channel}(stim2_X{x_value}(1,1)-(baseline_time_total./1000).*sf{1}:stim2_X{x_value}(1,1)-1,:,x_value),0,2); %std acrross traces (trial-to-trial)
                baseline_stdM_global(stim_num,1) = mean(baseline_std_global); %mean trial-to-trial variability along time
                           end             
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
%Parameters from mean trace                  
peaks(fileind).fname = fname; 
peaks(fileind).baseline_time_local=baseline_time; %baseline before each whisker stim.
peaks(fileind).baseline_time_global=baseline_time_total; %global baseline for the whole mean trace
peaks(fileind).baseline_local(:,x_value)=baseline_mean_local; 
peaks(fileind).baseline_global(:,x_value)=baseline_mean_global; 
peaks(fileind).baseline_global_std(:,x_value)=baseline_stdM_global;
peaks(fileind).fft_max_res(:,x_value)=max_res(1,t);
peaks(fileind).fft_max_res_exp(:,x_value)=max_res_exp(1,t);
peaks(fileind).F1(:,x_value)=F1(1,t);

    
                temp_mean(:,1)=fileind(ones(size(peaks(fileind).baseline_local(:,x_value),1),1));
                temp_mean(:,2)=[1:size(peaks(fileind).baseline_local(:,x_value),1)];
                temp_mean(:,3)=x_value(ones(size(peaks(fileind).baseline_local(:,x_value),1),1));
                temp_mean(:,4)=peaks(fileind).baseline_local(:,x_value);              
                
                F1_RI{x_value}(fileind,1)=peaks(fileind).F1(:,x_value);
                baseline_global{x_value}(fileind,1)=peaks(fileind).baseline_global(1,x_value);
                baseline_global_std{x_value}(fileind,1)=peaks(fileind).baseline_global_std(1,x_value);

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
         
              
peak_for_xls_mean=[peak_for_xls_mean; temp_mean];
end
close all
clear data_no_spike_no_DC
    end %end of files loop
    
    %% Statistics for peaks
%t-tests between cells      
    for stim_num=1:10;
      
                %normalize local baseline and perform ttest
        peak_baseline_local_noES=peak_for_xls_mean([find(peak_for_xls_mean(:,2)==stim_num & peak_for_xls_mean(:,3)==2)],4);
        peak_baseline_local_ES=peak_for_xls_mean([find(peak_for_xls_mean(:,2)==stim_num & peak_for_xls_mean(:,3)==3)],4);
        peaks_stat(stim_num).peak_baseline_local=[peak_baseline_local_noES, peak_baseline_local_ES];
        peaks_stat(stim_num).peak_baseline_local_m=mean(peaks_stat(stim_num).peak_baseline_local,1);
        peaks_stat(stim_num).peak_baseline_local_std=std(peaks_stat(stim_num).peak_baseline_local,0,1);
        norm_peak_baseline_local_noES=peak_baseline_local_noES./peak_baseline_local_noES;
        norm_peak_baseline_local_ES=peak_baseline_local_ES./peak_baseline_local_noES;
        %testing for normal distribution
        diff_val=peak_baseline_local_ES-peak_baseline_local_noES;
        [peaks_stat(stim_num).lillietest_h_peak_baseline_local, peaks_stat(stim_num).lillietest_p_peak_baseline_local] = lillietest(diff_val);
        %paired ttest
        peaks_stat(stim_num).norm_peak_baseline_local=[norm_peak_baseline_local_noES, norm_peak_baseline_local_ES];
        mean_baseline_local_std(stim_num,:)=std(peaks_stat(stim_num).norm_peak_baseline_local,0,1);
        mean_baseline_local_m(stim_num,:)=mean(peaks_stat(stim_num).norm_peak_baseline_local,1);
        peaks_stat(stim_num).norm_peak_baseline_local_std=mean_baseline_local_std(stim_num,:);
        [peaks_stat(stim_num).ttest_h_norm_peak_baseline_local, peaks_stat(stim_num).ttest_p_norm_peak_baseline_local]= ttest( peaks_stat(stim_num).norm_peak_baseline_local(:,1), peaks_stat(stim_num).norm_peak_baseline_local(:,2));
        [peaks_stat(stim_num).ttest_h_peak_baseline_local, peaks_stat(stim_num).ttest_p_peak_baseline_local]= ttest( peaks_stat(stim_num).peak_baseline_local(:,1), peaks_stat(stim_num).peak_baseline_local(:,2));

        %normalize F1 and perform ttest
        peak_F1_RI_noES= F1_RI{2};
        peak_F1_RI_ES=F1_RI{3};
        peaks_stat(stim_num).F1_RI=[peak_F1_RI_noES, peak_F1_RI_ES];
        peaks_stat(stim_num).F1_RI_m=mean(peaks_stat(stim_num).F1_RI,1);
        peaks_stat(stim_num).F1_RI_std=std(peaks_stat(stim_num).F1_RI,0,1);
        norm_peak_F1_RI_noES=peak_F1_RI_noES./peak_F1_RI_noES;
        norm_peak_F1_RI_ES=peak_F1_RI_ES./peak_F1_RI_noES;
        %testing for normal distribution
        diff_val=peak_F1_RI_ES-peak_F1_RI_noES;       
        [peaks_stat(stim_num).lillietest_h_F1_RI, peaks_stat(stim_num).lillietest_p_F1_RI] = lillietest(diff_val);
        %paired ttest
        peaks_stat(stim_num).norm_F1_RI=[norm_peak_F1_RI_noES, norm_peak_F1_RI_ES];
        F1_RI_std(stim_num,:)=std(peaks_stat(stim_num).norm_F1_RI,0,1);
        F1_RI_m(stim_num,:)=mean(peaks_stat(stim_num).norm_F1_RI,1);
        peaks_stat(stim_num).norm_F1_RI_std=F1_RI_std(stim_num,:);
        [peaks_stat(stim_num).ttest_h_norm_F1_RI, peaks_stat(stim_num).ttest_p_norm_F1_RI]= ttest(peaks_stat(stim_num).norm_F1_RI(:,1),peaks_stat(stim_num).norm_F1_RI(:,2));
        [peaks_stat(stim_num).ttest_h_F1_RI, peaks_stat(stim_num).ttest_p_F1_RI]= ttest(peaks_stat(stim_num).F1_RI(:,1),peaks_stat(stim_num).F1_RI(:,2));

         %normalize global baseline and perform ttest
        peak_baseline_global_noES= baseline_global{2};
        peak_baseline_global_ES=baseline_global{3};
        peaks_stat(stim_num).baseline_global=[peak_baseline_global_noES, peak_baseline_global_ES];
        peaks_stat(stim_num).baseline_global_m=mean(peaks_stat(stim_num).baseline_global,1);
        peaks_stat(stim_num).baseline_global_std=std(peaks_stat(stim_num).baseline_global,0,1);
        norm_peak_baseline_global_noES=peak_baseline_global_noES./peak_baseline_global_noES;
        norm_peak_baseline_global_ES=peak_baseline_global_ES./peak_baseline_global_noES;
        %testing for normal distribution
        diff_val=peak_baseline_global_ES-peak_baseline_global_noES;       
        [peaks_stat(stim_num).lillietest_h_baseline_global, peaks_stat(stim_num).lillietest_p_baseline_global] = lillietest(diff_val);
        %paired ttest
        peaks_stat(stim_num).norm_baseline_global=[norm_peak_baseline_global_noES, norm_peak_baseline_global_ES];
        mean_baseline_global_std(stim_num,:)=std(peaks_stat(stim_num).norm_baseline_global,0,1);
        mean_baseline_global_m(stim_num,:)=mean(peaks_stat(stim_num).norm_baseline_global,1);
        peaks_stat(stim_num).norm_baseline_global_std=mean_baseline_global_std(stim_num,:);
        [peaks_stat(stim_num).ttest_h_norm_baseline_global, peaks_stat(stim_num).ttest_p_norm_baseline_global]= ttest(peaks_stat(stim_num).norm_baseline_global(:,1),peaks_stat(stim_num).norm_baseline_global(:,2));
        [peaks_stat(stim_num).ttest_h_baseline_global, peaks_stat(stim_num).ttest_p_baseline_global]= ttest(peaks_stat(stim_num).baseline_global(:,1),peaks_stat(stim_num).baseline_global(:,2));

  %normalize global baseline std and perform ttest
        peak_baseline_global_std_noES= baseline_global_std{2};
        peak_baseline_global_std_ES=baseline_global_std{3};
        peaks_stat(stim_num).baseline_global_STD = [peak_baseline_global_std_noES, peak_baseline_global_std_ES];
        norm_peak_baseline_global_std_noES=peak_baseline_global_std_noES./peak_baseline_global_std_noES;
        norm_peak_baseline_global_std_ES=peak_baseline_global_std_ES./peak_baseline_global_std_noES;
        %testing for normal distribution
        diff_val=peak_baseline_global_std_ES-peak_baseline_global_std_noES;       
        [peaks_stat(stim_num).lillietest_h_baseline_global_STD, peaks_stat(stim_num).lillietest_p_baseline_global_STD] = lillietest(diff_val);
        %paired ttest
        peaks_stat(stim_num).norm_baseline_global_STD=[norm_peak_baseline_global_std_noES, norm_peak_baseline_global_std_ES];
        mean_baseline_global_STD_std(stim_num,:)=std(peaks_stat(stim_num).norm_baseline_global_STD,0,1);
        mean_baseline_global_STD_m(stim_num,:)=mean(peaks_stat(stim_num).norm_baseline_global_STD,1);
        peaks_stat(stim_num).norm_baseline_global_STD_std=mean_baseline_global_STD_std(stim_num,:);
        [peaks_stat(stim_num).ttest_h_norm_baseline_global_STD, peaks_stat(stim_num).ttest_p_norm_baseline_global_STD]= ttest(peaks_stat(stim_num).norm_baseline_global_STD(:,1),peaks_stat(stim_num).norm_baseline_global_STD(:,2));
        [peaks_stat(stim_num).ttest_h_baseline_global_STD, peaks_stat(stim_num).ttest_p_baseline_global_STD]= ttest(peaks_stat(stim_num).baseline_global_STD(:,1),peaks_stat(stim_num).baseline_global_STD(:,2));
        
        
        clear norm_peak_F1_RI_noES norm_peak_F1_RI_ES...
            peak_baseline_global_noES peak_baseline_global_ES norm_peak_baseline_global_noES norm_peak_baseline_global_ES...
            peak_baseline_local_noES peak_baseline_local_ES norm_peak_baseline_local_noES norm_peak_baseline_local_ES...
            peak_baseline_global_std_ES peak_baseline_global_std_noES norm_peak_baseline_global_std_ES norm_peak_baseline_global_std_noES...
            baseline_std_local baseline_std_global baseline_mean_local baseline_mean_global baseline_stdM_local
       
    end  
%% Plot parameters along the train stim

color_table_rand = rand(length(files_to_analyze),3);
    if print_flag==1;
 
         h7=figure;
        figure(h7);     
        hold on
        errorbar([1:size(mean_baseline_local_m,1)],mean_baseline_local_m(:,2),mean_baseline_local_std(:,2),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        line([0;12],[1;1],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean local baseline [mV]' ,'FontSize', 16);    
        title(['Mean local baseline, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);         
    %% Scatter plots for first peak parameters
stim_num=1;
min_line=floor(min(min(peaks_stat(stim_num).F1_RI)));
max_line=ceil(max(max(peaks_stat(stim_num).F1_RI)));
        g7=figure;
 hold on
scatter(peaks_stat(stim_num).F1_RI(:,1),peaks_stat(stim_num).F1_RI(:,2),180,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('NB-' ,'FontSize', 30);  ylabel('NB+','FontSize', 30);  
        title(['Adaptation power ratio, n=' num2str(length(files_to_analyze)) ', p=' num2str(peaks_stat(stim_num).ttest_p_norm_F1_RI)] ,'FontSize', 20);   
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

%% save figures
if save_flag==1
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10Hz'
print(h1,'Train_response_amp','-dpng','-r600','-opengl')
saveas(h2,'Train_response_amp_global.fig') 
print(h2,'Train_response_amp_global','-dpng','-r600','-opengl')
saveas(h3,'Train_onset_latency.fig') 
print(h3,'Train_onset_latency','-dpng','-r600','-opengl')
saveas(h4,'Train_90%peak_latency.fig') 
print(h4,'Train_90%peak_latency','-dpng','-r600','-opengl')
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
saveas(g5,'Stim1_90%peak_latency.fig') 
print(g5,'Stim1_90%peak_latency','-dpng','-r600','-opengl') 
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
%         plot(1:size(peaks(fileind).mean_amp,1),peaks(fileind).mean_amp(:,2),'--', 'color', color_table(fileind,:),'LineWidth',1.5,'marker','o','MarkerFaceColor', color_table(fileind,:),'MarkerEdgeColor','k','MarkerSize',7)
%         plot(1:size(peaks(fileind).mean_amp,1),peaks(fileind).mean_amp(:,3), 'color', color_table(fileind,:),'LineWidth',1.5,'marker','o','MarkerFaceColor', color_table(fileind,:),'MarkerEdgeColor','k','MarkerSize',7)     
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
