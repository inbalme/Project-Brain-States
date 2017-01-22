%% Analyze NBES protocol curr inj VC v3_response_parameters
% This file was created on 5/4/2016
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
close all
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
peaks_for_xls=[]; peak_for_xls_mean=[]; 
save_flag=0;
print_flag=1;
print_flag_cc=1;
fontsize1=20;
    
files_to_analyze =51; %[31,38,42,51,61,64,67,69,71,74,77];
% for fileind=1;
    for fileind=1:length(files_to_analyze) ;
%     close all
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
                
                 curr_inj_del = [Param.facade(27) Param.facade(30)] ;
            if length(Param.facade)>32
                curr_inj_del = [Param.facade(27) Param.facade(30) Param.facade(33)];
            end
            curr_inj_depo = Param.facade(1);
            curr_inj_hyper = Param.facade(2);
            curr_inj_dur = Param.facade(7);
if isempty(data_no_spikes)
    data_no_spikes{channel}=raw_data{channel};
end
raw_data{3} = raw_data{3}./20; %dividing by the LFP gain
           %% low-pass filtering below 300Hz                
                for xx=1:6
    for trace= 1:size(data_no_spikes{channel},2)    
            jj=data_no_spikes{channel}(:,trace,xx);
            data_no_spikes{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jj,dt,-1,300,0,0); 
    end
end

%% bandpass filtering to remove noise from LFP channel
bp_filt=files(files_to_analyze(fileind)).V2_filter;
for xx=1:6
    for trace= 1:size(raw_data{3},2)    
            jl=raw_data{3}(:,trace,xx);
            Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(jl,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
            jm=data_no_spikes{channel}(:,trace,xx);
            data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jm,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
             km=data_no_spikes_filt{channel}(:,trace,xx);
%              data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,bp_filt(1),bp_filt(2),0,0); %filtering Vm same as LFP
    end
end
%%
clear  start_time duration start_time start_sample end_sample interval interval_mat x y patch_xdata patch_ydata yex ylim_data
coeffs=[];
         start_time = curr_inj_del+0.1; %[sec] %[0,5]
%          start_time(2) = stim2_X{x_value}(1,end).*dt+0.1; %start 0.5 sec after the last sensory stim
         duration = 1.2; %curr_inj_dur-0.03; %[sec] 
            for t=1:length(start_time);
             start_sample(:,t) = ceil(start_time(t).*sf{1});
              start_sample(mod(start_sample(:,t),2)==1,t)=start_sample(mod(start_sample(:,t),2)==1,t)-1;                 
                if start_time(t)==0
                    start_sample(:,t) = 1;
                end
              end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1});
             interval_temp(:,t)= start_sample(:,t):end_sample(:,t);
              if mod(size(interval_temp(:,t),1),2)
                   interval(:,t) =interval_temp(1:end-1,t);
              else
                   interval(:,t) =interval_temp(:,t);
              end
            end 
            
%% Variance plots - two x_values on each plot, for flexible data
for x_value = [3,6]; %[3,6]; % 1:size(data.x_value,2) 
    clear color_table
    color_table=[0 0 0; [216 22 22]/256; [136 137 138]/256; [255 153 153]/256];
%     close all
for t=1:2;
plot_data(:,:,t)=(-1).*data_no_spikes_filt{channel}( interval(:,t),:,x_value); %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
plot_data_mean(:,t) = mean(plot_data(:,:,t),2);
plot_data_std(:,t) =  std(plot_data(:,:,t),0,2);
plot_data_var(:,t)=var(plot_data(:,:,t),0,2);
end
         if isempty(stim2_X{x_value}); stim2_loc_temp=[]; else
                 stim2_loc_temp=stim2_X{x_value}(1:2,:)-start_sample(1);    
         end
% plots:
if print_flag==1;
       trace_ind = [3,4,5,6]; %1:size(plot_data,2); %trace_ind is the numbers of traces to be plotted
       l=size(plot_data,1);
        l=l/2-1;
        DC=50.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        x2lim=[stim2_X{4}(1,1).*dt-0.5,stim2_X{4}(1,1).*dt+2];
        stim1=[];
        
        trace_to_plot = plot_data(:,trace_ind,1)+DC ; %DC is added just to make space between the traces
        [Fig1,h1]= fn_Plot_Trace_v2_no_stim(trace_to_plot, dt, dt, stim1, dt, stim2_loc_temp);     
         ylabel('Im [pA]', 'FontSize', fontsize1);
        hline = findobj(gcf, 'type', 'line'); set(hline(1:end),'color',color_table(1,:));
        set(gca,'fontsize',fontsize1);
        y1lim=[get(gca,'ylim')]';
       if files_to_analyze(fileind)==77; hpatch = findobj(gcf, 'type', 'patch'); delete(hpatch); end
        trace_to_plot = plot_data(:,trace_ind,2)+DC ;
        [Fig2,h2]= fn_Plot_Trace_v2_no_stim(trace_to_plot, dt, dt, stim1, dt, stim2_loc_temp);
        ylabel('Im [pA]', 'FontSize', fontsize1);
         hline = findobj(gcf, 'type', 'line'); set(hline(1:end),'color',color_table(2,:));
          set(gca,'fontsize',fontsize1);
        set(gca,'ylim',y1lim);
          if files_to_analyze(fileind)==77; hpatch = findobj(gcf, 'type', 'patch'); delete(hpatch);end
        [Fig3,h3]= fn_Plot_Trace_v2_no_stim(plot_data_std, dt, dt, stim1, dt, stim2_loc_temp);
        ylabel('std [pA]', 'FontSize', fontsize1);  xlabel('Time [sec]' ,'FontSize', fontsize1);
        hline = findobj(gcf, 'type', 'line'); set(hline(1),'color',color_table(2,:));  set(hline(2),'color',color_table(1,:));
         set(gca,'fontsize',fontsize1);
        [Fig4,h4] = fn_Plot_Trace_v2_no_stim(plot_data_mean, dt, dt,stim1, dt, stim2_loc_temp);
        ylabel('mean Im [pA]', 'FontSize', fontsize1);   
         hline = findobj(gcf, 'type', 'line'); set(hline(1),'color',color_table(2,:));  set(hline(2),'color',color_table(1,:));
         set(gca,'fontsize',fontsize1);
        [Fig5,h5]= fn_Plot_Trace_v2_no_stim(plot_data_std, dt, dt, stim1, dt, stim2_loc_temp);        
%         [Fig5,h5]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt, stim2_X{x_value(2)});
        ylabel('std [pA]', 'FontSize', fontsize1);  legend('NB-', 'NB+'); legend('boxoff','Location','northeast'); 
        hline = findobj(gcf, 'type', 'line'); set(hline(1),'color',color_table(2,:));  set(hline(2),'color',color_table(1,:));     
         set(gca,'fontsize',fontsize1);
       [Fig6,h6] = fn_Plot_Trace_v2_no_stim(plot_data_mean, dt, dt, stim1, dt, stim2_loc_temp);
        ylabel('mean Im [pA]', 'FontSize', fontsize1);  legend('NB-', 'NB+'); legend('boxoff','Location','northeast');
        hline = findobj(gcf, 'type', 'line'); set(hline(1),'color',color_table(2,:));  set(hline(2),'color',color_table(1,:));
        [Fig9,h9]= fn_Plot_Trace_v2_no_stim(plot_data_var, dt, dt, stim1, dt, stim2_loc_temp);
        ylabel('Variance [pA^2]', 'FontSize', fontsize1);  xlabel('Time [sec]' ,'FontSize',fontsize1);

                Fig7=figure;
                plot_mean_std_a=plot_data_mean+plot_data_std;
                plot_mean_std_b=plot_data_mean-plot_data_std;
        hold on
        fn_shadedErrorBar([1:size(plot_data,1)].*dt,plot_data_mean(:,1),plot_data_std(:,1),{'LineWidth',1,'color', color_table(1,:)},1);        
%             max_time_axis = ceil(size(plot_data,1).*dt);
            max_time_axis = size(plot_data,1).*dt;
            max_data = max(max(plot_mean_std_a(floor(stim2_X{4}(1,1).*dt./dt):end,:,:)));
            min_data = min(min(plot_mean_std_b(floor(stim2_X{4}(1,1).*dt./dt):end,:,:)));
            five_percent = (max_data-min_data).*0.02;
            min_data_wide = min_data-2.*five_percent;
            max_data_wide = max_data+2.*five_percent;
        fn_shadedErrorBar([1:size(plot_data,1)].*dt,plot_data_mean(:,2),plot_data_std(:,2),{'LineWidth',1,'color', color_table(2,:)},1);
            ylabel('Im mean+std [pA]', 'FontSize', fontsize1);  xlabel('Time [sec]' ,'FontSize', fontsize1);
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
%         stim1_Y = ones(size(stim1_X{channel})).*(max_data+five_percent); 
%         line(stim1_X{channel}.*dt,stim1_Y,'LineWidth',6,'Color','b');
        %plot whisker stim
        if x_value>3
            ylim_data=[get(gca,'ylim')]';
            stim2_loc_temp=stim2_X{4}(1:2,:)-start_sample(1);
            patch_xdata=[stim2_loc_temp; flipud(stim2_loc_temp)];
            yex=wextend('1D','sym',ylim_data,1);
            l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
            temp_y=wextend('ac','sym',yex,l);
            patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
            patch_cdata=ones(size(patch_xdata));
            patch(patch_xdata.*dt,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
            set(gca,'linewidth',1.2)
                    
        end
        hold off
        
        ax1 = get(Fig1, 'children');
        pos1 = [0.13 , 0.79 , 0.8 , 0.2];
        top1 = pos1(1,2)+pos1(1,4);

        ax2 = get(Fig2, 'children');
        pos2 = [0.13 , 0.56 , 0.8 , 0.2];
        top2 = pos2(1,2)+pos2(1,4);

        ax3 = get(Fig4, 'children');
        pos3 = [0.13 , 0.33 , 0.8 , 0.2];
        top3 = pos3(1,2)+pos3(1,4);

        ax4 = get(Fig3, 'children');
        pos4 = [0.13 , 0.1 , 0.8 , 0.2];
        top4 = pos4(1,2)+pos4(1,4);
        
        ax5 = get(Fig5, 'children'); ax6 = get(Fig6, 'children'); ax7 = get(Fig7, 'children'); ax9 = get(Fig9, 'children'); 
        
        Fig8 = figure;
        set(gcf,'color','w');
        set(gcf,'DefaultAxesFontSize',18);
        set(gcf,'DefaultAxesFontName','helvetica');
        set(gcf, 'PaperType', 'A4');
        set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

        ax_copy1 = copyobj(ax1,Fig8); % ax1 to new fig
        set(ax_copy1(1),'position',pos1(1,:),'xticklabel',[]); % Set its position  
        set(Fig8, 'currentaxes', ax_copy1(1));  xl=xlabel('');
%          set(F, 'currentaxes', spont_cc_Off_f46_ax_copy(1));  xl=xlabel('');  

        ax_copy2 = copyobj(ax2,Fig8); % ax2 to new fig
        set(ax_copy2(1),'position',pos2(1,:),'xticklabel',[]); % Set its position  
         set(Fig8, 'currentaxes', ax_copy2(1)); xl=xlabel('');

        ax_copy3 = copyobj(ax3,Fig8); % ax3 to new fig
        set(ax_copy3(1),'position',pos3(1,:),'xticklabel',[]); % Set its position  
        set(Fig8, 'currentaxes', ax_copy3(1));  xl=xlabel('');
        
        ax_copy4 = copyobj(ax4,Fig8); % ax3 to new fig
        set(ax_copy4(1),'position',pos4(1,:)) % Set its position
                
        close(Fig3); close(Fig4);      
%% saving figures
if save_flag==1;
    cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\Curr_inj_protocol'  
        saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x' num2str(x_value) '_ES-'],'fig') 
        print(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x' num2str(x_value) '_ES-'],'-dpng','-r600','-opengl') 
        saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x' num2str(x_value) '_ES+'],'fig') 
        print(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x' num2str(x_value) '_ES+'],'-dpng','-r600','-opengl') 
        saveas(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x' num2str(x_value) '_mean-subt'],'fig') 
        print(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x' num2str(x_value) '_mean-subt'],'-dpng','-r600','-opengl') 
        saveas(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x' num2str(x_value)],'fig') 
        print(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x' num2str(x_value)],'-dpng','-r600','-opengl') 
        saveas(Fig7,['f' num2str(files_to_analyze(fileind)) '_mean+std_x' num2str(x_value)],'fig') 
        print(Fig7,['f' num2str(files_to_analyze(fileind)) '_mean+std_x' num2str(x_value)],'-dpng','-r600','-opengl') 
        saveas(Fig8,['f' num2str(files_to_analyze(fileind)) '_summary_x' num2str(x_value)],'fig') 
        print(Fig8,['f' num2str(files_to_analyze(fileind)) '_summary_x' num2str(x_value)],'-dpng','-r600','-opengl')     
        saveas(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x' num2str(x_value) '_mean-subt'],'fig') 
        print(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x' num2str(x_value) '_mean-subt'],'-dpng','-r600','-opengl')     
        
%         set(ax1,'xlim',x2lim); set(ax2,'xlim',x2lim); set(ax5,'xlim',x2lim); set(ax6,'xlim',x2lim); set(ax7,'xlim',x2lim); set(ax9,'xlim',x2lim); 
%         
%        saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x' num2str(x_value) 'ES-_zoom-in.fig']) 
%         print(Fig1,['f' num2str(files_to_analyze(fileind)) '_traces_x' num2str(x_value) 'ES-_zoom-in'],'-dpng','-r600','-opengl') 
%         saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x' num2str(x_value) 'ES+_zoom-in.fig']) 
%         print(Fig2,['f' num2str(files_to_analyze(fileind)) '_traces_x' num2str(x_value) 'ES+_zoom-in'],'-dpng','-r600','-opengl') 
%         saveas(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x' num2str(x_value) '_mean-subt_zoom-in.fig']) 
%         print(Fig5,['f' num2str(files_to_analyze(fileind)) '_std_x' num2str(x_value) '_mean-subt_zoom-in'],'-dpng','-r600','-opengl') 
%         saveas(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x' num2str(x_value) '_zoom-in.fig']) 
%         print(Fig6,['f' num2str(files_to_analyze(fileind)) '_mean_x' num2str(x_value) '_zoom-in'],'-dpng','-r600','-opengl') 
%         saveas(Fig7,['f' num2str(files_to_analyze(fileind)) '_mean+std_x' num2str(x_value) '_zoom-in.fig']) 
%         print(Fig7,['f' num2str(files_to_analyze(fileind)) '_mean+std_x' num2str(x_value) '_zoom-in'],'-dpng','-r600','-opengl') 
%         saveas(Fig8,['f' num2str(files_to_analyze(fileind)) '_summary_x' num2str(x_value) '_zoom-in.fig']) 
%         print(Fig8,['f' num2str(files_to_analyze(fileind)) '_summary_x' num2str(x_value) '_zoom-in'],'-dpng','-r600','-opengl')     
%         saveas(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x' num2str(x_value) '_mean-subt_zoom-in.fig']) 
%         print(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x' num2str(x_value) '_mean-subt_zoom-in'],'-dpng','-r600','-opengl')     
        
end  
end
% close(Fig1, Fig2, Fig5, Fig6, Fig7, Fig8,Fig9)
%% Vm-LFP cross-correlations
%% bootstrap
 % bootstrap to calculate shuffeled cc
 for t=1:2;
     data_Vm{t}=data_no_spikes{1}(interval(:,t),:,x_value);
     data_Vm_filt{t}=data_no_spikes_filt{1}(interval(:,t),:,x_value);
     data_LFP{t}=Ch2_data_filt(interval(:,t),:,x_value);
     data_LFP_noDC{t} = fn_Subtract_Mean(data_LFP{t});
     data_Vm_filt_noDC{t} = fn_Subtract_Mean(data_Vm_filt{t});

    trials=size(data_no_spikes_filt{1},2);
    iterations=60;  
    DC_sub=1; %put 1 if the mean trace was not subtracted from the data.
    for i=1:iterations
        cc_shuffled_it{fileind,t}(:,i)=fn_bootstrap_cc(data_LFP{t}(:,:),data_Vm_filt{t}, 1, DC_sub);
    end
cc_shuffled_mean{fileind}(:,t)=mean(cc_shuffled_it{fileind,t}(:,:),2);
        % Calculate Correlations trace by trace
  for trace = 1:trials
        [cc{fileind,t}(:,trace),lags{fileind,t}(:,1)] = xcorr(data_LFP_noDC{t}(:,trace),data_Vm_filt_noDC{t}(:,trace),'coeff') ; 
%         [cc_bound{fileind,t}(:,trace), lags_bound{fileind,t}(:,1), bounds{fileind,t}(:,trace)] = crosscorr(data_Vm_filt{t}(:,trace),data_LFP{t}(:,trace));
        cc_shuff_sub{fileind,t}(:,trace)=cc{fileind,t}(:,trace)-cc_shuffled_mean{fileind}(:,t);     
  end
  
% %using shuffled-subtracted cc:
cc_shuff_sub_mean{fileind}(:,t) = mean(cc_shuff_sub{fileind,t}(:,:),2);
l_cc=length(cc_shuff_sub_mean{fileind}(:,t));
[c_mean_max_abs_val(t), c_mean_max_r(t)]=max(abs(cc_shuff_sub_mean{fileind}(ceil(l_cc./2)-50.*sf{1}./1000:ceil(l_cc./2)+50.*sf{1}./1000,t))); %look for the maximum in an interval of 50ms around the lag-zero point
% c_mean_max_r(t) = c_mean_max_r(t)+ceil(l_cc./2)-50.*sf{1}./1000;
c_mean_max_val(t)=cc_shuff_sub_mean{fileind}(c_mean_max_r(1)+ceil(l_cc./2)-50.*sf{1}./1000,t);
c_mean_max_val_shuff(t)=cc_shuffled_mean{fileind}(c_mean_max_r(1)+ceil(l_cc./2)-50.*sf{1}./1000,t); %not sure if this is correct to take it
c_mean_max_lag(t)=c_mean_max_r(t)-50.*sf{1}./1000;
cc_mean{fileind}(:,t) = mean(cc{fileind,t}(:,:),2);
%  [c_mean_max_abs_val(t), c_mean_max_r(t)]=max(abs(cc_mean{fileind}(:,t)));
% c_mean_max_val(t)=cc_mean{fileind}(c_mean_max_r(1),t);
c_max(:,t)=cc_shuff_sub{fileind,t}(c_mean_max_r(1),:);
c_lag0(:,t)=cc_shuff_sub{fileind,t}(lags{fileind,1}==0,:);
end
[max_diff_val, max_diff_loc]=max(abs(cc_shuff_sub_mean{fileind}(:,1)-cc_shuff_sub_mean{fileind}(:,2)));
% [max_diff_val, max_diff_loc]=max(abs(cc_mean{fileind}(:,1)-cc_mean{fileind}(:,2)));
c_maxdiff(:,1)=cc_shuff_sub{fileind,1}(max_diff_loc,:);
c_maxdiff(:,2)=cc_shuff_sub{fileind,2}(max_diff_loc,:);
%% Plots
trace_fontsize=12;
scalebar_fontsize=12;
clear color_table
    color_table=[0 0 0; [30,75,14]/256; [136 137 138]/256; [112,172,90]/256; [216 22 22]/256; [255 153 153]/256];
    whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;
    if print_flag_cc==1;
% plotting one trace of data and LFP against each other -
 trace=3;%1:size(data_no_spikes{1},2); %5;
 
% before ES
        
        data_Vm_plot=data_Vm{1}(:,trace);
        data_LFP_plot=data_LFP{1}(:,trace);
        interval_plot(:,1)=interval(:,1);   %[1:size(data_Vm_plot,1)].*dt*1000  ; %interval(:,1);   

        min_dist=min(data_LFP_plot.*20-data_Vm_plot);
        shift_LFP=5-min_dist;
        data_LFP_plot=shift_LFP+data_LFP_plot.*20;
        
Fig{fileind}(x_value)=figure;
    hold on
        p1=plot(interval_plot(:,1).*dt*1000,data_Vm_plot, 'color',color_table(1,:),'LineWidth',1.2);
        p2=plot(interval_plot(:,1).*dt*1000,data_LFP_plot,'color',color_table(2,:),'LineWidth',1.2);
        text(interval_plot(1,1).*dt*1000,data_Vm_plot(1),[num2str(floor(data_Vm_plot(1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
             axis tight
              ylim_data=[get(gca,'ylim')]';
              xlim_data=[get(gca,'xlim')]';

        if x_value>3
            patch_xdata=[stim2_X{x_value}(1:2,:); flipud(stim2_X{x_value}(1:2,:))];
            yex=wextend('1D','sym',ylim_data,1);
            l=(size(patch_xdata,2)-1)/2;
            patch_ydata=wextend('ac','sym',yex,l);
            patch_cdata=ones(size(patch_xdata));
            patch(patch_xdata.*dt*1000,patch_ydata,patch_cdata,'faceColor',whisker_stim_color(1,:),'edgecolor','none'); %'faceAlpha', 0.3
            set(gca,'xlim',xlim)
            set(gca,'children',flipud(get(gca,'children')))
        end
    %plotting scale bar
                    yline_start(1)=ylim_data(1)-2; yline_end(1)=yline_start(1)+5;
                    yline_start(2)=data_LFP_plot(1)-2; yline_end(2)=yline_start(2)+5;
                    xline_start=interval_plot(1,1).*dt*1000-650;
%                     if t==2;
%                      yline_start=ylim_data(1); yline_end=yline_start+2; 
%                  end
                    if x_value>3
                       yline_start(1)=ylim_data(1)-5; yline_end(1)=yline_start(1)+10; 
                       yline_start(2)=data_LFP_plot(1)-2; yline_end(2)=yline_start(2)+10;
                       xline_start=interval_plot(1,1).*dt*1000-350;
                    end
                     xline_end=xline_start+200;                   
                    xline=[xline_start xline_start;xline_end xline_start]; 
                    yline1=[yline_start(1) yline_start(1); yline_start(1) yline_end(1)];
                    yline2=[yline_start(2) yline_start(2); yline_start(2) yline_end(2)];
                    stringh=[num2str(xline(2,1)-xline(1,1)), ' mS'];
                    stringv1=[num2str(yline1(2,2)-yline1(2,1)), ' mV'];
                    stringv2=[num2str((yline1(2,2)-yline1(2,1))./20), ' mV'];
                    hline_h(1)=line(xline(:,1),yline1(1,:),'color',[0 0 0],'linewidth',2); %horizontal scale bar for Vm
                    hline_v(1)=line(xline(:,2),yline1(:,2),'color',[0 0 0],'linewidth',2); %vertical scale bar for Vm
                    hline_v(2)=line(xline(:,2),yline2(:,2),'color',[0 0 0],'linewidth',2); %vertical scale bar for LFP
%                     htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline1(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',10); %'color',[0 0 0]
                    htext_h=text(xline(1,1),yline1(1,1),stringh,'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
%                     htext_v(1)=text(xline(1,1),yline1(2,1)+(yline1(2,2)-yline1(2,1))/2,stringv1,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',10);
                    htext_v(1)=text(xline(1,1),yline1(2,1),stringv1,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
                    htext_v(2)=text(xline(1,1),yline2(2,1)+(yline2(2,2)-yline2(2,1))/2,stringv2,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);

    hold off
title('Vm-LFP single trace ES Off','FontSize', 16); 
% [legh,objh,outh,outm]=legend([p1,p2],{'Vm','LFP'},'fontsize',12,'box','off', 'location','northeast'); 
% set(objh,'linewidth',1.5); 
% set(objh(1:2),'fontsize',11); 
ylabel('Potential [mV]' ,'FontSize', 14); xlabel('Time [S]','FontSize', 14)
set(gca,'ylim',[yline_start(1)-5, ylim_data(2)+0.25*(ylim_data(2)- ylim_data(1))],'xlim',[xline_start-200, xlim_data(2)],'fontsize',14);                  
l=legend([p2,p1],{'LFP','Vm'},'linewidth',1.5,'Location','northeast', 'box', 'off'); 
l.LineWidth=1.5;
set(gca, 'visible', 'off') ;

% pause
% close(gcf);
% end %end of trace for loop

% after ES  
        interval_plot(:,2)=interval(:,2);        
        data_Vm_plot=data_Vm{2}(:,trace);
        data_LFP_plot=data_LFP{2}(:,trace);  
%         interval_plot(:,2)=[1:size(data_Vm_plot,1)].*dt*1000;
   
        min_dist=min(data_LFP_plot.*20-data_Vm_plot);
        shift_LFP=5-min_dist;
        data_LFP_plot=shift_LFP+data_LFP_plot.*20;
        
Fig{fileind}(x_value+1)=figure;
    hold on
        p1=plot(interval_plot(:,1).*dt*1000,data_Vm_plot,'color',color_table(5,:),'LineWidth',1.2);
       p2=plot(interval_plot(:,1).*dt*1000,data_LFP_plot,'color',color_table(4,:),'LineWidth',1.2);
        text(interval_plot(1,1).*dt*1000,data_Vm_plot(1),[num2str(floor(data_Vm_plot(1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
             axis tight
              ylim_data=[get(gca,'ylim')]';
              xlim_data=[get(gca,'xlim')]';
              
    if x_value>3
            patch_xdata=[stim2_X{x_value}(1:2,:); flipud(stim2_X{x_value}(1:2,:))];           
            yex=wextend('1D','sym',ylim_data,1);
            l=(size(patch_xdata,2)-1)/2;
            patch_ydata=wextend('ac','sym',yex,l);
            patch_cdata=ones(size(patch_xdata));
            patch(patch_xdata.*dt*1000,patch_ydata,patch_cdata,'faceColor',whisker_stim_color(1,:),'edgecolor','none'); %'faceAlpha', 0.3
            set(gca,'children',flipud(get(gca,'children')))
        end
    %plotting scale bar
                    yline_start(1)=ylim_data(1)-2; yline_end(1)=yline_start(1)+5;
                    yline_start(2)=data_LFP_plot(1)-2; yline_end(2)=yline_start(2)+5;
                    xline_start=interval_plot(1,1).*dt*1000-650;
%                     if t==2;
%                      yline_start=ylim_data(1); yline_end=yline_start+2; 
%                  end
                    if x_value>3
                       yline_start(1)=ylim_data(1)-5; yline_end(1)=yline_start(1)+10; 
                       yline_start(2)=data_LFP_plot(1)-2; yline_end(2)=yline_start(2)+10;
                       xline_start=interval_plot(1,1).*dt*1000-350;
                    end
                     xline_end=xline_start+200;                   
                    xline=[xline_start xline_start;xline_end xline_start]; 
                    yline1=[yline_start(1) yline_start(1); yline_start(1) yline_end(1)];
                    yline2=[yline_start(2) yline_start(2); yline_start(2) yline_end(2)];
                    stringh=[num2str(xline(2,1)-xline(1,1)), ' mS'];
                    stringv1=[num2str(yline1(2,2)-yline1(2,1)), ' mV'];
                    stringv2=[num2str((yline2(2,2)-yline2(2,1))./20), ' mV'];
                    hline_h(1)=line(xline(:,1),yline1(1,:),'color',[0 0 0],'linewidth',2); %horizontal scale bar for Vm
                    hline_v(1)=line(xline(:,2),yline1(:,2),'color',[0 0 0],'linewidth',2); %vertical scale bar for Vm 
                    hline_v(2)=line(xline(:,2),yline2(:,2),'color',[0 0 0],'linewidth',2); %vertical scale bar for LFP
%                     htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline1(1,1),stringh,'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',10); %'color',[0 0 0]
                    htext_h=text(xline(1,1),yline1(1,1),stringh,'HorizontalAlignment', 'left','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
%                     htext_v(1)=text(xline(1,1),yline1(2,1)+(yline1(2,2)-yline1(2,1))/2,stringv1,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',10);
                    htext_v(1)=text(xline(1,1),yline1(2,1),stringv1,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
                    htext_v(2)=text(xline(1,1),yline2(2,1)+(yline2(2,2)-yline2(2,1))/2,stringv2,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);
    hold off
title('Vm-LFP single trace ES On','FontSize', 16); 
% [legh,objh,outh,outm]=legend([p1,p2],{'Vm','LFP'},'fontsize',12,'box','off', 'location','northeast'); 
% set(objh,'linewidth',1.5); 
% set(objh(1:2),'fontsize',11); 
ylabel('Potential [mV]' ,'FontSize', 14); xlabel('Time [S]','FontSize', 14)
set(gca,'ylim',[yline_start(1)-5, ylim_data(2)+0.25*(ylim_data(2)- ylim_data(1))],'xlim',[xline_start-200, xlim_data(2)], 'fontsize',14);    
l=legend([p2,p1],{'LFP','Vm'},'Location','northeast', 'box', 'off');
l.LineWidth=1.5;
 set(gca, 'visible', 'off') ;
 %% Plot Vm-LFP cross-correlation of actual and shuffled data on separate figures
cc_shuff_sub_sem1=std(cc_shuff_sub{fileind,1}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,1}(:,:),2));
 cc_shuff_sub_sem2=std(cc_shuff_sub{fileind,2}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,2}(:,:),2));
 cc_shuffled_sem1=std(cc_shuffled_it{fileind,1}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,1}(:,:),2));
 cc_shuffled_sem2=std(cc_shuffled_it{fileind,2}(:,:),0,2)./sqrt(size(cc_shuff_sub{fileind,2}(:,:),2));
 
Fig{fileind}(x_value+10)=figure;
 hold on
 h1=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuff_sub_mean{fileind}(:,1),cc_shuff_sub_sem1,{'LineWidth',1.5,'color', color_table(1,:)},1);
h2=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuff_sub_mean{fileind}(:,2),cc_shuff_sub_sem2,{'LineWidth',1.5,'color', color_table(5,:)},1);
% ylim=get(gca,'ylim');
ylim=[-0.4 0.4];
yticks=[-0.4, -0.2 0, 0.2];
% hline_zero=line([0 0],[ylim(1) ylim(2)],'linestyle','-.','color',[136 137 138]./256,'linewidth',1);
% title('Vm-LFP crosscorrelation','FontSize', 16); 
l=legend([h1.mainLine h2.mainLine ],{'NB-','NB+'},'Location','northeast', 'box', 'off');
l.LineWidth=1.5;
l.FontSize=11;
% legend('ES off mean cc','ES on mean cc'); legend('boxoff');%legend('boxoff','fontsize',6,'linewidth',1.5,'Position',[0.39,0.45,0.25,0.1]);
xlabel('Lags [S]' ,'FontSize', 14); ylabel('CC' ,'FontSize', 14);
set(gca,'xlim',[-0.5 0.5],'xtick',[-0.5 0 0.5],'ylim', ylim, 'ytick',yticks,'fontsize',14); %set(gca,'xlim',[-0.2 0.2]); set(gca,'xtick',[-0.2 0 0.2]);
hold off

Fig{fileind}(x_value+11)=figure;
 hold on
h3=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuffled_mean{fileind}(:,1),cc_shuffled_sem1,{'LineWidth',1.5,'linestyle','--','color', color_table(1,:)},1);
h4=fn_shadedErrorBar(lags{fileind,1}.*dt,cc_shuffled_mean{fileind}(:,2),cc_shuffled_sem2,{'LineWidth',1.5,'linestyle','--','color', color_table(5,:)},1);
%  ylim=get(gca,'ylim');
ylim=[-0.2, 0.2];
 yticks=[-0.2,0,0.2];
 if x_value>3
ylim=[-0.8 0.8];
yticks=[-0.8, -0.4, 0, 0.4, 0.8];
 end
% hline_zero=line([0 0],[ylim(1) ylim(2)],'linestyle','-.','color',[136 137 138]./256,'linewidth',1);
% title('Vm-LFP crosscorrelation','FontSize', 16); 
l=legend([h3.mainLine h4.mainLine ],{'NB- shuffled','NB+ shuffled'},'Location','northeast', 'box', 'off');
l.LineWidth=1.5;
l.FontSize=11;
% legend('ES off mean cc','ES on mean cc'); legend('boxoff');%legend('boxoff','fontsize',6,'linewidth',1.5,'Position',[0.39,0.45,0.25,0.1]);
xlabel('Lags [S]' ,'FontSize', 14); ylabel('CC' ,'FontSize', 14);
set(gca,'xlim',[-0.5 0.5],'xtick',[-0.5 0 0.5],'ylim',ylim,'ytick',yticks,'fontsize',14); %set(gca,'xlim',[-0.2 0.2]); set(gca,'xtick',[-0.2 0 0.2]);
hold off
end
end
    end
%% Response parameters - amplitude, latency and more...
%find peak values and locations in interval following each whisker stim.

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
% save(filename, 'peaks', 'peaks_stat', 'peaks_for_xls', 'peak_for_xls_mean')
