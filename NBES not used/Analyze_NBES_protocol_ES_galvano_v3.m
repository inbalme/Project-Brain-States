%% Analyze NBES protocol ES+galvano v3
% This file was created on 5/10/2015
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
 global dt sf dt_galvano data data_no_spikes files Param
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2

    fileind = 44;
    
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    
    fname = files(fileind).extracted_name;        
    path = files(fileind).extracted_path;
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
                
               
        
%% Smooth data
 for x_value = 1:size(data.x_value,2) 
      for ii=1:4
                smooth_data{ii}(:,:,x_value) = sgolayfilt(raw_data{ii}(:,:,x_value), 1, 29);
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

               
                data_no_spike_no_DC{1}(:,:,x_value) = fn_Subtract_Mean(data_no_spikes{1}(:,:,x_value),meansubtract_interval);
                data_no_DC{1}(:,:,x_value) = fn_Subtract_Mean(raw_data{1}(:,:,x_value),meansubtract_interval);
    
end
 %% calculating meanVm and stdVm according to Golshani
 % spontaneous trial-to-trial variability
 x_value = 1;
 start_time =[0.5 4.5]; %[0.5 4.5]; %[4:0.1:4.9]; %[5:0.05:5.45]; %[sec]
 duration = 1; %1; %0.05; %[sec]

data_mat=[]; mean_mat=[]; std_mat=[]; std_mat_2=[]; interval_mat=[]; average_mean=[]; average_std=[]; average_std_2=[];

data_mat=data_no_spike_no_DC{channel};

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
    shadedErrorBar((1:size(raw_data{channel},1))'.*dt,mean(raw_data{channel}(:,:,x_value(t)),2),std(raw_data{channel}(:,:,x_value(t)),0,2),{'LineWidth',1,'color', color_table(t,:)},1);
%     xlim([0 35]); ylim([-0.1 1.5])
  end

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
plot_data_std =  std(plot_data,0,2);
% plot_data_std(29990:35010,:,:)=nan;
plot_data_CV = (-1).*plot_data_std./plot_data_mean; 
% plot_data_CV(29990:35010,:,:)=nan;

% plots for first x-value:
       trace_ind = [1:size(plot_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
       l=size(plot_data(:,trace_ind, x_value(1)),1);
        l=l/2-1;
        DC=20.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        trace_to_plot = plot_data(:,trace_ind, x_value(1))+DC ; %DC is added just to make space between the traces
        [Fig1,h1]= fn_Plot_Trace_v2(trace_to_plot, dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value(1)});
        hline = findobj(gcf, 'type', 'line'); 
        trace_to_plot = plot_data(:,trace_ind, x_value(2))+DC ;
        [Fig2,h2]= fn_Plot_Trace_v2(trace_to_plot, dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value(2)});
        
        [Fig3,h3]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value(2)});
        ylabel('std [mV]', 'FontSize', 16); 
%         hline = findobj(gcf, 'type', 'line'); set(hline(13),'color',[0 0 0]);  set(hline(14),'color',[0 0 1]);
        [Fig4,h4] = fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value(2)});
        ylabel('mean Vm [mV]', 'FontSize', 16);  
        
        [Fig5,h5]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value(2)});
        ylabel('std [mV]', 'FontSize', 16);  legend('ES Off', 'ES On'); legend('boxoff','Location','northeast'); 
%         hline = findobj(gcf, 'type', 'line'); set(hline(13),'color',[0 0 0]);  set(hline(14),'color',[0 0 1]);
        [Fig6,h6] = fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value(2)});
        ylabel('mean Vm [mV]', 'FontSize', 16);  legend('ES Off', 'ES On'); legend('boxoff','Location','northeast');
     
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
        
        close(Fig3); close(Fig4);
%% Response parameters - amplitude, latency and more...
%find peak values and locations in interval following each whisker stim.
int_time=50; %[ms]
x_value = 2;
stim_num=1;
data_mean = mean(data_no_spikes{channel}(:,:,x_value),2); %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
stim_begin=stim2_X{x_value}(1,stim_num);
stim_end=stim_begin+int_time.*sf_galvano./1000;
data_mean_smooth=sgolayfilt(data_mean, 1, 29);
%  [mean_peak_val{x_value,stim_num}, mean_peak_loc{x_value,stim_num}] =  findpeaks(data_mean_smooth(stim_begin:stim_end,1),'npeaks',1);   
[mean_peak_val_all{x_value,stim_num}, mean_peak_loc_all{x_value,stim_num}] =  findpeaks(data_mean_smooth(stim_begin:stim_end,1),'sortstr','descend');
mean_peak_val{x_value,stim_num}=mean_peak_val_all{x_value,stim_num}(1);
mean_peak_loc{x_value,stim_num}=mean_peak_loc_all{x_value,stim_num}(1);
 mean_peak_loc{x_value,stim_num}=mean_peak_loc{x_value,stim_num}+stim_begin;
[Fig6,h6] = fn_Plot_Trace_v2([data_mean_smooth], dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value}); %add to plot: raw_data{channel}(:,:,x_value)
        ylabel('mean Vm [mV]', 'FontSize', 16);  legend('ES Off', 'ES On'); legend('boxoff','Location','northeast');
        hold on
%           scatter(mean_peak_loc{x_value,stim_num}.*dt,mean_peak_val{x_value,stim_num},'r', 'fill')
          

            for trace = 1:size(data_no_spikes{1},2)         
                 h1=plot([1:size(raw_data{channel},1)].*dt, raw_data{channel}(:,trace,x_value),'k');
                 h2=scatter(mean_peak_loc{x_value,stim_num}.*dt,raw_data{channel}(mean_peak_loc{x_value,stim_num},trace,x_value),'r', 'fill');
                 pause
                 delete(h1)
                 delete(h2)
                 figure(Fig6)
            end
        hold off
           
% peaks(x_value,trace,stim_num)        
 a=1:size(raw_data{channel},1);
  %use the peak locations to find the peaks in the single traces
%% Variance plot - one x_value at a time, for flexible data
       x_value = 2;
       plot_data=raw_data{channel}; %data_no_spikes %raw_data %data_no_spike_no_DC %data_no_DC
        plot_data_mean = mean(plot_data,2);
        plot_data_std =  std(plot_data,0,2);
        plot_data_CV = (-1).*plot_data_std./plot_data_mean; 

       trace_ind = [1:size(plot_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
       l=size(plot_data(:,trace_ind, x_value),1);
        l=l/2-1;
        DC=10.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        trace_to_plot = plot_data(:,trace_ind, x_value)+DC ;
        [Fig1,h1]= fn_Plot_Trace_v2(trace_to_plot, dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value});
        [Fig2,h2]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value});
        [Fig3,h3] = fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt, stim1_X{channel}, dt_galvano, stim2_X{x_value});



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

start_sample = start_time(t).*sf{channel};
if start_time(t)==0
    start_sample = 1;
end
end_sample = start_sample+duration.*sf{channel}-1;
interval = start_sample:end_sample;
spec_mat = raw_data{channel}(interval,:,x_value);
DC= mean(spec_mat,1);
l=size(spec_mat,1);
l=l/2-1;
DC=wextend('addrow', 'per', DC, l);
spec_mat_noDC = spec_mat-DC;

[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf{channel});
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

%% cross-correlation Vm-LFP - spontaneous activity
clear c lags c_mean lags_new start_time duration start_time start_sample end_sample interval interval_mat x y patch_xdata patch_ydata yex ylim_data
x_value=2:3; %2:3; %1;
coeffs=[];
trace_type=length(x_value); %for trace_type=1 plot the xcorr before and after the ES according to the two start times. for trace_type=2 plot the xcorr from the same time point in the 2 different x-values.
switch trace_type
    case 1
 start_time = [0.5,5]; %[sec] %[0,5]
 duration = 2.5; %[sec] 
for t=1:length(start_time);
 start_sample(:,t) = start_time(t).*sf{1};
if start_time(t)==0
    start_sample(:,t) = 1;
end
end_sample(:,t) = start_sample(:,t)+duration.*sf{1}-1;
interval(:,t) = round(start_sample(:,t):end_sample(:,t));
% interval_mat1(:,:,t) = data_downsamp(interval(:,t),:,x_value);
% interval_mat2(:,:,t) = data_EEG(interval(:,t),:,x_value);

% trace=1;
  for trace = 1:size(data_no_spikes{1},2)
          x{t}(:,trace)=data_no_spikes{1}(interval(:,t),trace,x_value)-mean(data_no_spikes{1}(interval(:,t),trace,x_value));
          y{t}(:,trace)=raw_data{3}(interval(:,t),trace,x_value)-mean(raw_data{3}(interval(:,t),trace,x_value));
        [c{t}(:,trace),lags{t}(:,1)] = xcorr(x{t}(:,trace),y{t}(:,trace),'coeff') ;
        
        [r,p] = corrcoef(x{t}(:,trace),y{t}(:,trace));
        coeffs{t}(1,trace)= r(1,2);
        pvals{t}(1,trace)= p(1,2);
     end      
c_mean(:,t) = mean(c{t}(:,:),2);
coeffs_mean(1,t)=mean(coeffs{t});
pvals_mean(1,t)=mean(pvals{t});
end

% plotting one trace of data and LFP against each other -
trace=5;
% before ES
figure
subplot(3,1,1)
hold on
plot(interval(:,1).*dt,x{1}(:,trace), 'k-','LineWidth',1.2)
plot(interval(:,1).*dt,y{1}(:,trace),'g-','LineWidth',1.2)
title('Vm-LFP single trace ES Off','FontSize', 16); legend('Vm','LFP'); legend('boxoff', 'location','northeast'); 
ylabel('Potential [mV]' ,'FontSize', 16);
hold off

% after ES
subplot(3,1,2)
hold on
plot(interval(:,2).*dt,x{2}(:,trace), 'k-','LineWidth',1.2)
plot(interval(:,2).*dt,y{2}(:,trace),'g-','LineWidth',1.2)
title('Vm-LFP single trace ES On','FontSize', 16); legend('Vm','LFP'); legend('boxoff', 'location','northeast'); 
ylabel('Potential [mV]' ,'FontSize', 16);
hold off

%plotting the crosscorrelation for a single trace
subplot(3,1,3)
hold on
plot( lags{1,1}.*dt,c{1}(:,trace), 'k-.')
plot( lags{1,1}.*dt,c{2}(:,trace), 'b-.')
plot( lags{1,1}.*dt,c_mean(:,1),'k-', 'LineWidth',1.5)
plot( lags{1,1}.*dt,c_mean(:,2), 'b-','LineWidth',1.5)
title('Vm-LFP crosscorrelation','FontSize', 16); legend('ES off single trace','ES on single trace','ES off average','ES on average'); legend('boxoff');%legend('boxoff','fontsize',6,'linewidth',1.5,'Position',[0.39,0.45,0.25,0.1]);
xlabel('Time [sec]' ,'FontSize', 16); ylabel('Correlation' ,'FontSize', 16);
hold off

    case 2
 start_sample =stim2_X{x_value(1)}(1,1);  %start with sensory stim
 duration = 2; %[sec] 
end_sample = start_sample+duration.*sf{1}-1;
interval(:,1) = round(start_sample:end_sample);
% interval_mat1(:,:,t) = data_downsamp(interval(:,t),:,x_value);
% interval_mat2(:,:,t) = data_EEG(interval(:,t),:,x_value);
for i=1:length(x_value);
     for trace = 1:size(data_no_spikes{1},2) 
% trace=1;
         x{i}(:,trace)=data_no_spikes{1}(interval(:,1),trace,x_value(i))-mean(data_no_spikes{1}(interval(:,1),trace,x_value(i)));
         y{i}(:,trace)=raw_data{3}(interval(:,1),trace,x_value(i))-mean(raw_data{3}(interval(:,1),trace,x_value(i)));
        [c{i}(:,trace),lags{i}(:,1)] = xcorr( x{i}(:,trace), y{i}(:,trace),'coeff') ;
        
        [r,p] = corrcoef( x{i}(:,trace), y{i}(:,trace));
        coeffs{i}(1,trace)= r(1,2);
        pvals{i}(1,trace)= p(1,2);  
end
c_mean(:,i) = mean(c{i}(:,:),2);
coeffs_mean(1,i)=mean(coeffs{i});
pvals_mean(1,i)=mean(pvals{i});
end  

% plotting one trace of data and LFP against each other -
trace=5;
% no ES
figure
subplot(3,1,1)
hold on
plot(interval(:,1).*dt,x{1}(:,trace), 'k-')
plot(interval(:,1).*dt,y{1}(:,trace),'g-')
max_data=max(max(x{1}(:,trace)));
% stim2_Y = ceil(ones(size(stim2_X{2})).*(max_data)); 
% line(stim2_X{2}.*dt_galvano,stim2_Y,'LineWidth',6,'Color','r')
% patch_ydata=[zeros(size(stim2_X{2})); stim2_Y];
title('Vm-LFP single trace ES Off','FontSize', 16); legend('Vm','LFP');  legend('boxoff', 'location','northeast'); 
ylabel('Potential [mV]' ,'FontSize', 16);

patch_xdata=[stim2_X{2}; flipud(stim2_X{2})];
ylim_data=[get(gca,'ylim')]';
yex=wextend('1D','sym',ylim_data,1);
l=(size(patch_xdata,2)-1)/2;
patch_ydata=wextend('ac','sym',yex,l);
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt_galvano,patch_ydata,patch_cdata,'faceColor','r','edgecolor','w','faceAlpha', 0.3);
hold off

clear patch_xdata patch_ydata yex ylim_data
% with ES
subplot(3,1,2)
hold on
plot(interval(:,1).*dt,x{2}(:,trace), 'k-')
plot(interval(:,1).*dt,y{2}(:,trace),'g-')
title('Vm-LFP single trace ES On','FontSize', 16); legend('Vm','LFP'); legend('boxoff', 'location','northeast'); 
ylabel('Potential [mV]' ,'FontSize', 16);
patch_xdata=[stim2_X{2}; flipud(stim2_X{2})];
ylim_data=[get(gca,'ylim')]';
yex=wextend('1D','sym',ylim_data,1);
l=(size(patch_xdata,2)-1)/2;
patch_ydata=wextend('ac','sym',yex,l);
patch_cdata=ones(size(patch_xdata));
patch(patch_xdata.*dt_galvano,patch_ydata,patch_cdata,'faceColor','r','edgecolor','w','faceAlpha', 0.3);
hold off

%plotting the crosscorrelation for a single trace
subplot(3,1,3)
hold on
plot( lags{1,1}.*dt,c{1}(:,trace), 'k-.')
plot( lags{1,1}.*dt,c{2}(:,trace), 'b-.')
plot( lags{1,1}.*dt,c_mean(:,1),'k','LineWidth',1.5)
plot( lags{1,1}.*dt,c_mean(:,2), 'b','LineWidth',1.5)
title('Vm-LFP crosscorrelation','FontSize', 16); legend('ES off single trace','ES on single trace','ES off average','ES on average'); legend('boxoff');%legend('boxoff','fontsize',6,'linewidth',1.5,'Position',[0.39,0.45,0.25,0.1]);
xlabel('Time [sec]' ,'FontSize', 16); ylabel('Correlation' ,'FontSize', 16);
hold off

end

%%

%% 
filename =  'SNR_plot_12_cells'; %'F62P6X6+9_zero_traces+std+mean_zoom-in';

cd 'd:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH';

print (6, '-depsc2', filename)
