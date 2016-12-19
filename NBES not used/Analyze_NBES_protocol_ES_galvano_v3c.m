%% Analyze NBES protocol ES+galvano v3c
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
peaks_for_xls=[]; peak_for_xls_mean=[]; cc_for_xls_mean=[];
files_to_analyze = [44,46,48,50,52];
    for fileind=1:length(files_to_analyze) ;
    
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
   %% cross-correlation Vm-LFP - spontaneous activity
clear c lags c_mean lags_new start_time duration start_time start_sample end_sample interval interval_mat x y patch_xdata patch_ydata yex ylim_data
x_value=1; %2:3; %1;
coeffs=[];
trace_type=length(x_value); %for trace_type=1 plot the xcorr before and after the ES according to the two start times. for trace_type=2 plot the xcorr from the same time point in the 2 different x-values.
switch trace_type
    case 1
 start_time = [0.5,5]; %[sec] %[0,5]
 duration = 2.5; %[sec] 
for t=1:length(start_time);
 start_sample(:,t) = ceil(start_time(t).*sf{1});
if start_time(t)==0
    start_sample(:,t) = 1;
end
end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
interval(:,t) = start_sample(:,t):end_sample(:,t);

% Calculate Correlations trace by trace
%   for trace = 1:size(data_no_spikes{1},2)
%           data_Vm{t}(:,trace)=data_no_spikes{1}(interval(:,t),trace,x_value)-mean(data_no_spikes{1}(interval(:,t),trace,x_value));
%           data_LFP{t}(:,trace)=raw_data{3}(interval(:,t),trace,x_value)-mean(raw_data{3}(interval(:,t),trace,x_value));
%         [c{t}(:,trace),lags{t}(:,1)] = xcorr(data_Vm{t}(:,trace),data_LFP{t}(:,trace),'coeff') ;       
%      end      
% c_mean(:,t) = mean(c{t}(:,:),2);

%
% end

%% bootstrap
 % bootstrap to calculate shuffeled cc
    trials=size(data_no_spikes{1},2);
    data_Vm{t}=data_no_spikes{1}(interval(:,t),:,x_value);
    data_LFP{t}=raw_data{3}(interval(:,t),:,x_value);
    
 %subtract DC from Vm matrix        
        mean_traces_data_Vm= mean(data_Vm{t},1);
        l=size(data_Vm{t},1);
        l=l/2-1;
        DC_Vm=wextend('addrow', 'per', mean_traces_data_Vm, l);
        data_Vm_noDC{t} = data_Vm{t}-DC_Vm;
 %subtract DC from LFP matrix
        mean_traces_data_LFP= mean(data_LFP{t},1);
        l=size(data_LFP{t},1);
        l=l/2-1;
        DC_LFP=wextend('addrow', 'per', mean_traces_data_LFP, l);
        data_LFP_noDC{t} = data_LFP{t}-DC_LFP;
        
        clear mean_traces_data_Vm DC_Vm mean_traces_data_LFP DC_LFP
        
        for n = 1:60
                f = 1;
                while f
                    Indexs1 = randperm(trials);
                    Indexs2 = randperm(trials);
                    f = sum(Indexs1==Indexs2); %the loop will run until two sets of shuffled traces with no overlaps are generated
                end
                for j = 1:length(Indexs1) %running on all traces
                    w1 = Indexs1(j);
                    w2 = Indexs2(j);               
                    cc_shuffled_temp(:,j) = xcorr(data_Vm_noDC{t}(:,w1,x_value),data_LFP_noDC{t}(:,w2,x_value),'coeff');                                
                end
                cc_shuffled_temp2(:,n)=mean(cc_shuffled_temp,2); %average across traces
        end
        cc_shuffled(:,t)=mean(cc_shuffled_temp2,2); %average across bootstrap repetitions
        
        % Calculate Correlations trace by trace
  for trace = 1:size(data_no_spikes{1},2)
        [cc{t}(:,trace),lags{t}(:,1)] = xcorr(data_Vm_noDC{t}(:,trace),data_LFP_noDC{t}(:,trace),'coeff') ; 
        cc_shuf_sub{t}(:,trace)=cc{t}(:,trace)-cc_shuffled(:,t);
     end      
cc_mean(:,t) = mean(cc{t}(:,:),2);
cc_shuff_sub_mean(:,t) = mean(cc_shuf_sub{t}(:,:),2);

end
%% Plots
% plotting one trace of data and LFP against each other -
trace=5;
% before ES
figure
subplot(3,1,1)
hold on
plot(interval(:,1).*dt,data_Vm{1}(:,trace), 'k-','LineWidth',1.2)
plot(interval(:,1).*dt,data_LFP{1}(:,trace),'g-','LineWidth',1.2)
title('Vm-LFP single trace ES Off','FontSize', 16); legend('Vm','LFP'); legend('boxoff', 'location','northeast'); 
ylabel('Potential [mV]' ,'FontSize', 16);
hold off

% after ES
subplot(3,1,2)
hold on
plot(interval(:,2).*dt,data_Vm{2}(:,trace), 'k-','LineWidth',1.2)
plot(interval(:,2).*dt,data_LFP{2}(:,trace),'g-','LineWidth',1.2)
title('Vm-LFP single trace ES On','FontSize', 16); legend('Vm','LFP'); legend('boxoff', 'location','northeast'); 
ylabel('Potential [mV]' ,'FontSize', 16);
hold off

%plotting the crosscorrelation for a single trace+the mean
subplot(3,1,3)
hold on
plot( lags{1,1}.*dt,cc_mean(:,1),'k-', 'LineWidth',1.5)
plot( lags{1,1}.*dt,cc_mean(:,2), 'b-','LineWidth',1.5)
plot( lags{1,1}.*dt,cc_shuff_sub_mean(:,1), 'k-.')
plot( lags{1,1}.*dt,cc_shuff_sub_mean(:,2), 'b-.')
title('Vm-LFP crosscorrelation','FontSize', 16); legend('ES off mean','ES on mean','ES off mean corrected','ES on mean corrected'); legend('boxoff');%legend('boxoff','fontsize',6,'linewidth',1.5,'Position',[0.39,0.45,0.25,0.1]);
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
[c_mean_max_r(i) c_mean_max_c(i) c_mean_max_val(i)]=find(max(abs(c_mean(:,i))));
c_max(:,i)=c{i}(c_mean_max_r(i),:);
c_lag0(:,i)=c{i}(lags{1}==0,:);
coeffs_mean(1,i)=mean(coeffs{i});
pvals_mean(1,i)=mean(pvals{i});
end  
% 
[max_diff_val max_diff_loc]=max(abs(c_mean(:,1)-c_mean(:,2)));
c_maxdiff(:,1)=c{1}(max_diff_loc,:);
c_maxdiff(:,2)=c{2}(max_diff_loc,:);
cc_evoked(fileind).fname = fname;
cc_evoked(fileind).cc_max=c_mean_max_val;
cc_evoked(fileind).cc_lag0=c_mean(lags{1}==0,:);
cc_evoked(fileind).cc_maxdiff=c_mean(max_diff_loc,:);
cc_evoked(fileind).cc_max_trace = c_max;
cc_evoked(fileind).cc_lag0_trace = c_lag0;
cc_evoked(fileind).cc_maxdiff_trace = c_maxdiff;

for i=1:length(x_value);
            xi=x_value(i);
                cc_mean(:,1)=fileind(ones(size(cc_evoked(fileind).cc_max(:,i),1),1));
                cc_mean(:,2)=xi(ones(size(cc_evoked(fileind).cc_max(:,i),1),1));
                cc_mean(:,3)=cc_evoked(fileind).cc_max(:,i);
                cc_mean(:,4)=cc_evoked(fileind).cc_lag0(:,i);
                cc_mean(:,5)=cc_evoked(fileind).cc_maxdiff(:,i);                
            cc_for_xls_mean=[cc_for_xls_mean; cc_mean];
end
%% plotting one trace of data and LFP against each other -
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

%% Cross correlation statistics
% [h_cc_max_trace(fileind,1) p_cc_max_trace(fileind,1) ci_cc_max_trace(fileind,:) stats_cc_max_trace(fileind)]= ttest(cc_evoked(fileind).cc_max_trace(:,1),cc_evoked(fileind).cc_max_trace(:,2));
% [h_cc_lag0_trace(fileind,1) p_cc_lag0_trace(fileind,1) ci_cc_lag0_trace(fileind,:) stats_cc_lag0_trace(fileind)]= ttest(cc_evoked(fileind).cc_lag0_trace(:,1),cc_evoked(fileind).cc_lag0_trace(:,2));
% [h_cc_maxdiff_trace(fileind,1) p_cc_maxdiff_trace(fileind,1) ci_cc_maxdiff_trace(fileind,:) stats_cc_maxdiff_trace(fileind)]= ttest(cc_evoked(fileind).cc_maxdiff_trace(:,1),cc_evoked(fileind).cc_maxdiff_trace(:,2));
% 
%     clear c_mean c_mean_max_r c_mean_max_c c_mean_max_val x y c c_max c_lag0 coeffs_mean pvals_mean lags c_maxdiff cc_mean
    end % end of files loop
t_cc_max(:,1)=cc_for_xls_mean([find(cc_for_xls_mean(:,2)==2)],3);
t_cc_max(:,2)=cc_for_xls_mean([find(cc_for_xls_mean(:,2)==3)],3);
t_cc_lag0(:,1)=cc_for_xls_mean([find(cc_for_xls_mean(:,2)==2)],4);
t_cc_lag0(:,2)=cc_for_xls_mean([find(cc_for_xls_mean(:,2)==3)],4);
t_cc_maxdiff(:,1)=cc_for_xls_mean([find(cc_for_xls_mean(:,2)==2)],5);
t_cc_maxdiff(:,2)=cc_for_xls_mean([find(cc_for_xls_mean(:,2)==3)],5);
[h_cc_max(:,1) p_cc_max(:,1) ci_cc_max(1,:) stats_cc_max]= ttest(t_cc_max(:,1),t_cc_max(:,2));

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
 %%
filename =  'SNR_plot_12_cells'; %'F62P6X6+9_zero_traces+std+mean_zoom-in';

cd 'd:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH';

print (6, '-depsc2', filename)
