%% Analyze NBES protocol ES+galvano v3e - recent version.
% This file was created on 5/10/2015 and updated on 19/1/2016 and 7/2/2016
%This file is used for the analysis of files created with extract_NBES_Data_v2

%       protocol 5 (NBES+galvnao (5 x-values))
%       protocol 6 (ChAT intracellular current injections light+galvnao (15 x-values))
%       protocol 8 - ChAT intracellular: light only (3 x-values): 1) train
%       20 Hz, 10ms ondur, 50 pulses. 2) single pulse 5 sec ondur. 3)single pulse 10 ms
%       protocol 9 - IC: I-V curve
%       protocol 10 (NBES+galvnao train+test (3 x-values))
%                 11 - NBES: 3 currents+galvano train+test (6 x-values)
   
%% for opening workspace saved 
cc_stat=[]; cc_spont=[]; cc_evoked=[];
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
cc_spont_for_xls_mean=[]; cc_evoked_for_xls_mean=[];
files_to_analyze =[44,46,48,50,52];

    for fileind=1:length(files_to_analyze) ;         
    clearvars -except cc_stat cc_spont cc_evoked  files_to_analyze files cc_spont_for_xls_mean cc_evoked_for_xls_mean
    
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
for trace_type=1:2; 
if trace_type==1
    x_value=[1,1];
end
if trace_type==2
x_value=[2:3]; %2:3; %[1,1];
end

clear c lags c_mean lags_new start_time duration start_time start_sample end_sample interval interval_mat x y patch_xdata patch_ydata yex ylim_data
coeffs=[]; 
% U=unique(x_value);
% trace_type=length(U); %for trace_type=1 plot the xcorr before and after the ES according to the two start times. for trace_type=2 plot the xcorr from the same time point in the 2 different x-values.
switch trace_type
    case 1
         start_time = [0.5,5]; %[sec] %[0,5]
         duration = 2; %[sec] 
            for t=1:length(start_time);
             start_sample(:,t) = ceil(start_time(t).*sf{1});
                if start_time(t)==0
                    start_sample(:,t) = 1;
                end
              end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
              interval(:,t) = start_sample(:,t):end_sample(:,t);
            end 
    case 2
         start_sample =stim2_X{x_value(1)}(1,1);  %start with sensory stim
         duration = 2; %[sec] 
        end_sample = start_sample+duration.*sf{1}-1;
        interval(:,1) = round(start_sample:end_sample);
        interval(:,2) = interval(:,1);
end

for t=1:2
%% Subtract mean from the interval (fragment of trace)
%  data_Vm{t}=data_no_spikes{1}(:,:,x_value);
%  data_LFP{t}=raw_data{3}(:,:,x_value);
data_Vm{t}=data_no_spikes{1}(interval(:,t),:,x_value(t));
data_LFP{t}=raw_data{3}(interval(:,t),:,x_value(t));
data_Vm_noDC{t} = fn_subtract_mean_trace(data_Vm{t});
data_LFP_noDC{t} = fn_subtract_mean_trace(data_LFP{t});
%% bootstrap
 % bootstrap to calculate shuffeled cc
    trials=size(data_no_spikes{1},2);
    iterations=60;  
    DC_sub=1; %put 1 if the mean trace was not subtracted from the data.
cc_shuffled(:,t)=fn_bootstrap_cc(data_Vm{t},data_LFP{t}(:,:),iterations,DC_sub);
       if trace_type==2
          cc_shuffled(:,t)=zeros(size(cc_shuffled(:,t))); % for the evoked there is no need to subtract the shuffled correlations from the bootstrap.
      end        
        % Calculate Correlations trace by trace
  for trace = 1:trials
        [cc{t}(:,trace),lags{t}(:,1)] = xcorr(data_Vm_noDC{t}(:,trace),data_LFP_noDC{t}(:,trace),'coeff') ; 
        cc_shuf_sub{t}(:,trace)=cc{t}(:,trace)-cc_shuffled(:,t);     
  end
     
cc_mean(:,t) = mean(cc{t}(:,:),2);
cc_shuff_sub_mean(:,t) = mean(cc_shuf_sub{t}(:,:),2);
[c_mean_max_abs_val(t), c_mean_max_r(t)]=max(abs(cc_shuff_sub_mean(:,t)));
c_mean_max_val(t)=cc_shuff_sub_mean(c_mean_max_r(1),t);
c_max(:,t)=cc{t}(c_mean_max_r(1),:);
c_lag0(:,t)=cc{t}(lags{1}==0,:);
end
[max_diff_val, max_diff_loc]=max(abs(cc_shuff_sub_mean(:,1)-cc_shuff_sub_mean(:,2)));
c_maxdiff(:,1)=cc{1}(max_diff_loc,:);
c_maxdiff(:,2)=cc{2}(max_diff_loc,:);
%% Plots
% plotting one trace of data and LFP against each other -
trace=5;
% before ES
Fig{fileind}(trace_type)=figure;
hold on
plot(interval(:,1).*dt,data_Vm_noDC{1}(:,trace), 'k','LineWidth',1.2)
plot(interval(:,1).*dt,data_LFP_noDC{1}(:,trace),'color',[136 137 138]./256,'LineWidth',1.2)
title('Vm-LFP single trace ES Off','FontSize', 16); legend('Vm','LFP'); legend('boxoff', 'location','northeast'); 
ylabel('Potential [mV]' ,'FontSize', 16);

    if trace_type==2;
        patch_xdata=[stim2_X{2}; flipud(stim2_X{2})];
        ylim_data=[get(gca,'ylim')]';
        yex=wextend('1D','sym',ylim_data,1);
        l=(size(patch_xdata,2)-1)/2;
        patch_ydata=wextend('ac','sym',yex,l);
        patch_cdata=ones(size(patch_xdata));
        patch(patch_xdata.*dt_galvano,patch_ydata,patch_cdata,'faceColor','r','edgecolor','w','faceAlpha', 0.3);
    end
hold off
clear patch_xdata patch_ydata yex ylim_data

        ax{fileind}(trace_type,:) = get(Fig{fileind}(trace_type), 'children');

% after ES
Fig{fileind}(trace_type+2)=figure;
hold on
plot(interval(:,2).*dt,data_Vm_noDC{2}(:,trace),'color',[13 49 133]./256,'LineWidth',1.2)
plot(interval(:,2).*dt,data_LFP_noDC{2}(:,trace),'color',[33 118 209]./256,'LineWidth',1.2)
title('Vm-LFP single trace ES On','FontSize', 16); legend('Vm','LFP'); legend('boxoff', 'location','northeast'); 
ylabel('Potential [mV]' ,'FontSize', 16);
    if trace_type==2
        patch_xdata=[stim2_X{2}; flipud(stim2_X{2})];
        ylim_data=[get(gca,'ylim')]';
        yex=wextend('1D','sym',ylim_data,1);
        l=(size(patch_xdata,2)-1)/2;
        patch_ydata=wextend('ac','sym',yex,l);
        patch_cdata=ones(size(patch_xdata));
        patch(patch_xdata.*dt_galvano,patch_ydata,patch_cdata,'faceColor','r','edgecolor','w','faceAlpha', 0.3);
    end
hold off

        ax{fileind}(trace_type+2,:) = get(Fig{fileind}(trace_type+2), 'children');

%plotting the crosscorrelation for a single trace+the mean
Fig{fileind}(trace_type+4)=figure;
hold on
% plot( lags{1,1}.*dt,cc{1}(:,trace),'k-.', 'LineWidth',1)
% plot( lags{1,1}.*dt,cc{2}(:,trace), 'b-.','LineWidth',1)
plot( lags{1,1}.*dt,cc_shuff_sub_mean(:,1), 'k-', 'LineWidth',1.5)
plot( lags{1,1}.*dt,cc_shuff_sub_mean(:,2), 'color',[13 49 133]./256, 'LineWidth',1.5)
title('Vm-LFP crosscorrelation','FontSize', 16); legend('ES off mean cc','ES on mean cc'); legend('boxoff');%legend('boxoff','fontsize',6,'linewidth',1.5,'Position',[0.39,0.45,0.25,0.1]);
xlabel('Time [sec]' ,'FontSize', 16); ylabel('Correlation' ,'FontSize', 16);
set(gca,'xlim',[-0.5 0.5]); set(gca,'xtick',[-0.5 0 0.5]);
hold off

        ax{fileind}(trace_type+4,:) = get(Fig{fileind}(trace_type+4), 'children');
%%

    switch trace_type
        case 1
        cc_spont(fileind).fname = fname;
        cc_spont(fileind).cc_max=c_mean_max_val; %value of maximal correlation 
        cc_spont(fileind).cc_max_lag=lags{1}(c_mean_max_r);         
        cc_spont(fileind).cc_max_time=lags{1}(c_mean_max_r).*dt;
        cc_spont(fileind).cc_lag0=cc_shuff_sub_mean(lags{1}==0,:); %mean value of correlation at lag zero
        cc_spont(fileind).cc_maxdiff=cc_shuff_sub_mean(max_diff_loc,:); %location (row) of maximal difference in correlation (the lag)
        cc_spont(fileind).cc_max_trace = c_max; %location (row) of maximal correlation (the lag)
        cc_spont(fileind).cc_lag0_trace = c_lag0; %mean values of correlation at lag zero (for all traces)
        cc_spont(fileind).cc_maxdiff_trace = c_maxdiff; %mean values of correlation at largest difference point (for all traces)

            for t=1:2;
                        xi=x_value(t);
                            temp_mean(:,1)=fileind(ones(size(cc_spont(fileind).cc_max(:,t),1),1)); %col #1: file index
                            temp_mean(:,2)=xi(ones(size(cc_spont(fileind).cc_max(:,t),1),1)); %col #2: x_value
                            temp_mean(:,3)=t(ones(size(cc_spont(fileind).cc_max(:,t),1),1)); %col #3: ES off/on
                            temp_mean(:,4)=cc_spont(fileind).cc_max(:,t); %col #4: max value
                            temp_mean(:,5)=cc_spont(fileind).cc_lag0(:,t); %col #5: value at lag=zero
                            temp_mean(:,6)=cc_spont(fileind).cc_maxdiff(:,t); %col #6: value at lag with maximal difference                
                        cc_spont_for_xls_mean=[cc_spont_for_xls_mean; temp_mean];
            end    
        case 2
        cc_evoked(fileind).fname = fname;
        cc_evoked(fileind).cc_max=c_mean_max_val; %value of maximal correlation
        cc_evoked(fileind).cc_max_lag=lags{1}(c_mean_max_r);
        cc_evoked(fileind).cc_max_time=lags{1}(c_mean_max_r).*dt;
        cc_evoked(fileind).cc_lag0=cc_shuff_sub_mean(lags{1}==0,:); %mean value of correlation at lag zero
        cc_evoked(fileind).cc_maxdiff=cc_shuff_sub_mean(max_diff_loc,:); %location (row) of maximal difference in correlation (the lag)
        cc_evoked(fileind).cc_max_trace = c_max; %location (row) of maximal correlation (the lag)
        cc_evoked(fileind).cc_lag0_trace = c_lag0; %mean values of correlation at lag zero (for all traces)
        cc_evoked(fileind).cc_maxdiff_trace = c_maxdiff; %mean values of correlation at largest difference point (for all traces)

            for t=1:2;
                        xi=x_value(t);
                            temp_mean(:,1)=fileind(ones(size(cc_evoked(fileind).cc_max(:,t),1),1)); %col #1: file index
                            temp_mean(:,2)=xi(ones(size(cc_evoked(fileind).cc_max(:,t),1),1)); %col #2: x_value
                            temp_mean(:,3)=t(ones(size(cc_evoked(fileind).cc_max(:,t),1),1)); %col #3: ES off/on
                            temp_mean(:,4)=cc_evoked(fileind).cc_max(:,t); %col #4: max value
                            temp_mean(:,5)=cc_evoked(fileind).cc_lag0(:,t); %col #5: value at lag=zero
                            temp_mean(:,6)=cc_evoked(fileind).cc_maxdiff(:,t); %col #6: value at lag with maximal difference                  
                        cc_evoked_for_xls_mean=[cc_evoked_for_xls_mean; temp_mean];
            end
    end
clear c_max c_lag0 temp_mean cc_shuff_sub_mean c_mean_max_r c_mean_max_c...% c_mean_max_val...
    cc_mean max_diff_val max_diff_loc c_maxdiff

        pos1 = [0.08 , 0.79 , 0.6 , 0.2];
        top1 = pos1(1,2)+pos1(1,4);
        pos2 = [0.08 , 0.56 , 0.6 , 0.2];
        top2 = pos2(1,2)+pos2(1,4);
        pos3 = [0.75 , 0.79 , 0.2 , 0.2];
        top3 = pos3(1,2)+pos3(1,4);  
        pos4 = [0.08 , 0.33 , 0.6 , 0.2];
        top4 = pos4(1,2)+pos4(1,4);
        pos5 = [0.08 , 0.1 , 0.6 , 0.2];
        top5 = pos5(1,2)+pos5(1,4);
        pos6 = [0.08 , 0.33 , 0.2 , 0.2];
        top6 = pos6(1,2)+pos6(1,4);     

        F{fileind} = figure;
        set(gcf,'color','w');
        set(gcf,'DefaultAxesFontSize',12);
        set(gcf,'DefaultAxesFontName','helvetica');
        set(gcf, 'PaperType', 'A4');
        set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

        switch trace_type
            case 1
        ax_copy1 = copyobj( ax{fileind}(1,:),F{fileind}); % ax1 to new fig
        set(ax_copy1(1),'position',pos1) % Set its position 
        ax_copy3 = copyobj( ax{fileind}(3,:),F{fileind}); % ax3 to new fig
        set(ax_copy3(1),'position',pos3) % Set its position         
         ax_copy5 = copyobj( ax{fileind}(5,:),F{fileind}); % ax5 to new fig
        set(ax_copy5(1),'position',pos5) % Set its position  

            case 2
        ax_copy2 = copyobj( ax{fileind}(2,:),F{fileind}); % ax2 to new fig
        set(ax_copy2(1),'position',pos2(1,:)) % Set its position  
         ax_copy4 = copyobj( ax{fileind}(4,:),F{fileind}); % ax4 to new fig
        set(ax_copy4(1),'position',pos4(1,:)) % Set its position          
         ax_copy6 = copyobj( ax{fileind}(6,:),F{fileind}); % ax6 to new fig
        set(ax_copy6(1),'position',pos6(1,:)) % Set its position  
        end
        
%         close(Fig{fileind}(1)); close(Fig{fileind}(2)); close(Fig{fileind}(3));
%         close(Fig{fileind}(4)); close(Fig{fileind}(5)); close(Fig{fileind}(6));
    end  
end
%         ax1 = get(Fig{fileind}(1), 'children');

% 
%         ax2 = get(Fig(2), 'children');
%         pos2 = [0.08 , 0.56 , 0.6 , 0.2];
%         top2 = pos2(1,2)+pos2(1,4);
% 
%         ax3 = get(Fig(3), 'children');
%         pos3 = [0.75 , 0.79 , 0.2 , 0.2];
%         top3 = pos3(1,2)+pos3(1,4);  
% 
%         ax4 = get(Fig(4), 'children');
%         pos4 = [0.08 , 0.33 , 0.6 , 0.2];
%         top4 = pos4(1,2)+pos4(1,4);
% 
%         ax5 = get(Fig(5), 'children');
%         pos5 = [0.08 , 0.1 , 0.6 , 0.2];
%         top5 = pos5(1,2)+pos5(1,4);
% 
%         ax6 = get(Fig(6), 'children');
%         pos6 = [0.08 , 0.33 , 0.2 , 0.2];
%         top6 = pos6(1,2)+pos6(1,4);     
%         
%         F = figure;
%         set(gcf,'color','w');
%         set(gcf,'DefaultAxesFontSize',12);
%         set(gcf,'DefaultAxesFontName','helvetica');
%         set(gcf, 'PaperType', 'A4');
%         set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
%         set(gcf,'PaperOrientation','portrait');
%         set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);
% 
%         ax_copy1 = copyobj(ax1,F); % ax1 to new fig
%         set(ax_copy1(1),'position',pos1(1,:)) % Set its position  
% 
%         ax_copy2 = copyobj(ax2,F); % ax2 to new fig
%         set(ax_copy2(1),'position',pos2(1,:)) % Set its position  
% 
%         ax_copy3 = copyobj(ax3,F); % ax3 to new fig
%         set(ax_copy3(1),'position',pos3(1,:)) % Set its position  
%         
%          ax_copy4 = copyobj(ax4,F); % ax3 to new fig
%         set(ax_copy4(1),'position',pos4(1,:)) % Set its position  
%         
%          ax_copy5 = copyobj(ax5,F); % ax3 to new fig
%         set(ax_copy5(1),'position',pos5(1,:)) % Set its position  
%         
%          ax_copy6 = copyobj(ax6,F); % ax3 to new fig
%         set(ax_copy6(1),'position',pos6(1,:)) % Set its position  
%         
%         close(Fig1); close(Fig2); close(Fig3);
%% Statistics for cross-correlation
%t-tests between cells  
% switch trace_type
%     case 1
%Spontaneous
        %normalize cc max values and perform ttest
        cc_max_val_noES=cc_spont_for_xls_mean([find(cc_spont_for_xls_mean(:,3)==1)],4);
        cc_max_val_ES=cc_spont_for_xls_mean([find(cc_spont_for_xls_mean(:,3)==2)],4);
        cc_stat.spont.max_val=[cc_max_val_noES, cc_max_val_ES];
        norm_cc_max_val_noES=cc_max_val_noES./cc_max_val_noES;
        norm_cc_max_val_ES=cc_max_val_ES./cc_max_val_noES;

        %testing for normal distribution
        diff_val=cc_max_val_ES-cc_max_val_noES;
        [cc_stat.spont.lillietest_h_max_val, cc_stat.spont.lillietest_p_max_val] = lillietest(diff_val);
        %paired ttest
        cc_stat.spont.norm_max_val=[norm_cc_max_val_noES, norm_cc_max_val_ES];
        [cc_stat.spont.ttest_h_norm_max_val, cc_stat.spont.ttest_p_norm_max_val]= ttest(cc_stat.spont.norm_max_val(:,1),cc_stat.spont.norm_max_val(:,2));
        [cc_stat.spont.ttest_h_max_val, cc_stat.spont.ttest_p_max_val]= ttest(cc_stat.spont.max_val(:,1),cc_stat.spont.max_val(:,2));
        [cc_stat.spont.wilcoxon_p_max_val, cc_stat.spont.wilcoxon_h_max_val]= ranksum(cc_stat.spont.max_val(:,1),cc_stat.spont.max_val(:,2));
        [cc_stat.spont.wilcoxon_p_norm_max_val, cc_stat.spont.wilcoxon_h_norm_max_val]= ranksum(cc_stat.spont.norm_max_val(:,1),cc_stat.spont.norm_max_val(:,2));

        
         %normalize cc zero-lag values and perform ttest
        cc_lag0_noES=cc_spont_for_xls_mean([find(cc_spont_for_xls_mean(:,3)==1)],5);
        cc_lag0_ES=cc_spont_for_xls_mean([find(cc_spont_for_xls_mean(:,3)==2)],5);
        cc_stat.spont.lag0=[cc_lag0_noES, cc_lag0_ES];
        norm_cc_lag0_noES=cc_lag0_noES./cc_lag0_noES;
        norm_cc_lag0_ES=cc_lag0_ES./cc_lag0_noES;
        %testing for normal distribution
        diff_val=cc_lag0_ES-cc_lag0_noES;
        [cc_stat.spont.lillietest_h_lag0, cc_stat.spont.lillietest_p_lag0] = lillietest(diff_val);
        %paired ttest
        cc_stat.spont.norm_lag0=[norm_cc_lag0_noES, norm_cc_lag0_ES];
        [cc_stat.spont.ttest_h_norm_lag0, cc_stat.spont.ttest_p_norm_lag0]= ttest(cc_stat.spont.norm_lag0(:,1),cc_stat.spont.norm_lag0(:,2));
        [cc_stat.spont.ttest_h_lag0, cc_stat.spont.ttest_p_lag0]= ttest(cc_stat.spont.lag0(:,1),cc_stat.spont.lag0(:,2));
        [cc_stat.spont.wilcoxon_p_lag0, cc_stat.spont.wilcoxon_h_lag0]= ranksum(cc_stat.spont.lag0(:,1),cc_stat.spont.lag0(:,2));
        [cc_stat.spont.wilcoxon_p_norm_lag0, cc_stat.spont.wilcoxon_h_norm_lag0]= ranksum(cc_stat.spont.norm_lag0(:,1),cc_stat.spont.norm_lag0(:,2));
        %normalize cc maxdiff values and perform ttest
        cc_maxdiff_noES=cc_spont_for_xls_mean([find(cc_spont_for_xls_mean(:,3)==1)],6);
        cc_maxdiff_ES=cc_spont_for_xls_mean([find(cc_spont_for_xls_mean(:,3)==2)],6);
        cc_stat.spont.maxdiff=[cc_maxdiff_noES, cc_maxdiff_ES];
        norm_cc_maxdiff_noES=cc_maxdiff_noES./cc_maxdiff_noES;
        norm_cc_maxdiff_ES=cc_maxdiff_ES./cc_maxdiff_noES;
        %testing for normal distribution
        diff_val=cc_maxdiff_ES-cc_maxdiff_noES;
        [cc_stat.spont.lillietest_h_maxdiff, cc_stat.spont.lillietest_p_maxdiff] = lillietest(diff_val);
        %paired ttest
        cc_stat.spont.norm_maxdiff=[norm_cc_maxdiff_noES, norm_cc_maxdiff_ES];
        [cc_stat.spont.ttest_h_norm_maxdiff, cc_stat.spont.ttest_p_norm_maxdiff]= ttest(cc_stat.spont.norm_maxdiff(:,1),cc_stat.spont.norm_maxdiff(:,2));
        [cc_stat.spont.ttest_h_maxdiff, cc_stat.spont.ttest_p_maxdiff]= ttest(cc_stat.spont.maxdiff(:,1),cc_stat.spont.maxdiff(:,2));
        [cc_stat.spont.wilcoxon_p_maxdiff, cc_stat.spont.wilcoxon_h_maxdiff]= ranksum(cc_stat.spont.maxdiff(:,1),cc_stat.spont.maxdiff(:,2));
        [cc_stat.spont.wilcoxon_p_norm_maxdiff, cc_stat.spont.wilcoxon_h_norm_maxdiff]= ranksum(cc_stat.spont.norm_maxdiff(:,1),cc_stat.spont.norm_maxdiff(:,2));

%   Evoked
        %normalize cc max values and perform ttest
        cc_max_val_noES=cc_evoked_for_xls_mean([find(cc_evoked_for_xls_mean(:,2)==2)],4);
        cc_max_val_ES=cc_evoked_for_xls_mean([find(cc_evoked_for_xls_mean(:,2)==3)],4);
        cc_stat.evoked.max_val=[cc_max_val_noES, cc_max_val_ES];
        norm_cc_max_val_noES=cc_max_val_noES./cc_max_val_noES;
        norm_cc_max_val_ES=cc_max_val_ES./cc_max_val_noES;
        %testing for normal distribution
        diff_val=cc_max_val_ES-cc_max_val_noES;
        [cc_stat.evoked.lillietest_h_max_val, cc_stat.evoked.lillietest_p_max_val] = lillietest(diff_val);
        %paired ttest
        cc_stat.evoked.norm_max_val=[norm_cc_max_val_noES, norm_cc_max_val_ES];
        [cc_stat.evoked.ttest_h_norm_max_val, cc_stat.evoked.ttest_p_norm_max_val]= ttest(cc_stat.evoked.norm_max_val(:,1),cc_stat.evoked.norm_max_val(:,2));
        [cc_stat.evoked.ttest_h_max_val, cc_stat.evoked.ttest_p_max_val]= ttest(cc_stat.evoked.max_val(:,1),cc_stat.evoked.max_val(:,2));
        [cc_stat.evoked.wilcoxon_p_max_val, cc_stat.evoked.wilcoxon_h_max_val]= ranksum(cc_stat.evoked.max_val(:,1),cc_stat.evoked.max_val(:,2));
        
         %normalize cc zero-lag values and perform ttest
        cc_lag0_noES=cc_evoked_for_xls_mean([find(cc_evoked_for_xls_mean(:,2)==2)],5);
        cc_lag0_ES=cc_evoked_for_xls_mean([find(cc_evoked_for_xls_mean(:,2)==3)],5);
        cc_stat.evoked.lag0=[cc_lag0_noES, cc_lag0_ES];
        norm_cc_lag0_noES=cc_lag0_noES./cc_lag0_noES;
        norm_cc_lag0_ES=cc_lag0_ES./cc_lag0_noES;
        %testing for normal distribution
        diff_val=cc_lag0_ES-cc_lag0_noES;
        [cc_stat.evoked.lillietest_h_lag0, cc_stat.evoked.lillietest_p_lag0] = lillietest(diff_val);
        %paired ttest+wilcoxon test for the original and normalized values
        cc_stat.evoked.norm_lag0=[norm_cc_lag0_noES, norm_cc_lag0_ES];
        [cc_stat.evoked.ttest_h_norm_lag0, cc_stat.evoked.ttest_p_norm_lag0]= ttest(cc_stat.evoked.norm_lag0(:,1),cc_stat.evoked.norm_lag0(:,2));
        [cc_stat.evoked.ttest_h_lag0, cc_stat.evoked.ttest_p_lag0]= ttest(cc_stat.evoked.lag0(:,1),cc_stat.evoked.lag0(:,2));
        [cc_stat.evoked.wilcoxon_p_lag0, cc_stat.evoked.wilcoxon_h_lag0]= ranksum(cc_stat.evoked.lag0(:,1),cc_stat.evoked.lag0(:,2));
        [cc_stat.evoked.wilcoxon_p_norm_lag0, cc_stat.evoked.wilcoxon_h_norm_lag0]= ranksum(cc_stat.evoked.norm_lag0(:,1),cc_stat.evoked.norm_lag0(:,2));

        %normalize cc maxdiff values and perform ttest
        cc_maxdiff_noES=cc_evoked_for_xls_mean([find(cc_evoked_for_xls_mean(:,2)==2)],6);
        cc_maxdiff_ES=cc_evoked_for_xls_mean([find(cc_evoked_for_xls_mean(:,2)==3)],6);
        cc_stat.evoked.maxdiff=[cc_maxdiff_noES, cc_maxdiff_ES];
        norm_cc_maxdiff_noES=cc_maxdiff_noES./cc_maxdiff_noES;
        norm_cc_maxdiff_ES=cc_maxdiff_ES./cc_maxdiff_noES;
        %testing for normal distribution
        diff_val=cc_maxdiff_ES-cc_maxdiff_noES;
        [cc_stat.evoked.lillietest_h_maxdiff, cc_stat.evoked.lillietest_p_maxdiff] = lillietest(diff_val);
        %paired ttest
        cc_stat.evoked.norm_maxdiff=[norm_cc_maxdiff_noES, norm_cc_maxdiff_ES];
        [cc_stat.evoked.ttest_h_norm_maxdiff, cc_stat.evoked.ttest_p_norm_maxdiff]= ttest(cc_stat.evoked.norm_maxdiff(:,1),cc_stat.evoked.norm_maxdiff(:,2));
        [cc_stat.evoked.ttest_h_maxdiff, cc_stat.evoked.ttest_p_maxdiff]= ttest(cc_stat.evoked.maxdiff(:,1),cc_stat.evoked.maxdiff(:,2));
        [cc_stat.evoked.wilcoxon_p_maxdiff, cc_stat.evoked.wilcoxon_h_maxdiff]= ranksum(cc_stat.evoked.maxdiff(:,1),cc_stat.evoked.maxdiff(:,2));
        [cc_stat.evoked.wilcoxon_p_norm_maxdiff, cc_stat.evoked.wilcoxon_h_norm_maxdiff]= ranksum(cc_stat.evoked.norm_maxdiff(:,1),cc_stat.evoked.norm_maxdiff(:,2));


%% Plot spont and evoked lag zero cc peak, not normalized. p-values taken from normalized values paired t-test 
%
% min_line=floor(min(min(cc_stat.spont.lag0)));
% max_line=ceil(max(max(cc_stat.spont.lag0)));
min_line=min(min(cc_stat.spont.lag0)).*1.05;
max_line=max(max(cc_stat.spont.lag0)).*1.05;
 g1=figure;
 hold on
scatter(cc_stat.spont.lag0(:,1),cc_stat.spont.lag0(:,2),50,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('ES off' ,'FontSize', 16);  ylabel('ES on','FontSize', 16); 
        title(['Spontaneous Zero Lag Cross-Correlation, n=' num2str(length(files_to_analyze)) ', p=' num2str(cc_stat.spont.ttest_p_norm_lag0)] ,'FontSize', 16);   
        set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line]);
        
min_line=min(min(cc_stat.evoked.lag0)).*1.05;
max_line=max(max(cc_stat.evoked.lag0)).*0.95;
 g2=figure;
 hold on
scatter(cc_stat.evoked.lag0(:,1),cc_stat.evoked.lag0(:,2),50,'markerfacecolor','k')
line([min_line;max_line], [min_line;max_line],'linewidth',2,'color','k')
hold off
        xlabel('ES off' ,'FontSize', 16);  ylabel('ES on','FontSize', 16); 
        title(['Evoked Zero Lag Cross-Correlation, n=' num2str(length(files_to_analyze)) ', p=' num2str(cc_stat.evoked.ttest_p_norm_lag0)] ,'FontSize', 16);   
        set(gca, 'xlim', [min_line, max_line], 'ylim', [min_line, max_line]);      

 %% Save figures       
%  cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations'

% saveas(g1,'Spont_lag_zero_CC.fig') 
% saveas(g1,'Spont_lag_zero_CC.png')   
% saveas(g2,'Evoked_lag_zero_CC.fig') 
% saveas(g2,'Evoked_lag_zero_CC.png')  
 %% Paired plot of spont normalized max peak
%  clear CC_Y_max_val CC_X E
% CC_Y_max_val=cc_stat.spont.norm_max_val';
% CC_X(1,:)=ones(1,size(CC_Y_max_val,2));
% CC_X(2,:)=2*ones(1,size(CC_Y_max_val,2));
% E = std(CC_Y_max_val,0,2);
% figure
% hold on
% line(CC_X,CC_Y_max_val,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(CC_X(:,1), mean(CC_Y_max_val,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.25,0.5];
%         set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'Before ES','After ES'} ,'box', 'off'); %'fontweight', 'bold', 
%          ylabel('CC spont', 'FontSize', 20,'fontname', 'arial');
%          title('normalized max peak', 'FontSize',20,'fontname', 'arial')
         %% Paired plot of spont non-normalized max peak
%  clear CC_Y_max_val CC_X E
% CC_Y_max_val=cc_stat.spont.max_val';
% CC_X(1,:)=ones(1,size(CC_Y_max_val,2));
% CC_X(2,:)=2*ones(1,size(CC_Y_max_val,2));
% E = std(CC_Y_max_val,0,2);
% figure
% hold on
% line(CC_X,CC_Y_max_val,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(CC_X(:,1), mean(CC_Y_max_val,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
%         y1limits = [-0.5 0.5];
%         y1ticks = [-0.5,0,0.5];
%         set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'Before ES','After ES'} ,'box', 'off'); %'fontweight', 'bold', 
%          ylabel('CC spont', 'FontSize', 20,'fontname', 'arial');
%           title('max peak', 'FontSize',20,'fontname', 'arial')
          %% Paired plot of evoked normalized max peak
%  clear CC_Y_max_val CC_X E
% CC_Y_max_val=cc_stat.evoked.norm_max_val';
% % CC_Y_max_val=cc_stat.evoked.max_val';
% CC_X(1,:)=ones(1,size(CC_Y_max_val,2));
% CC_X(2,:)=2*ones(1,size(CC_Y_max_val,2));
% E = std(CC_Y_max_val,0,2);
% figure
% hold on
% line(CC_X,CC_Y_max_val,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(CC_X(:,1), mean(CC_Y_max_val,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
%         y1limits = [0.5 1.5];
%         y1ticks = [0.5,1,1.5];
%         set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'Before ES','After ES'} ,'box', 'off'); %'fontweight', 'bold', 
%          ylabel('CC evoked', 'FontSize', 20,'fontname', 'arial');
%           title('normalized max peak', 'FontSize',20,'fontname', 'arial')
          %% Paired plot of evoked non-normalized max peak
%  clear CC_Y_max_val CC_X E
% CC_Y_max_val=cc_stat.evoked.max_val';
% CC_X(1,:)=ones(1,size(CC_Y_max_val,2));
% CC_X(2,:)=2*ones(1,size(CC_Y_max_val,2));
% E = std(CC_Y_max_val,0,2);
% figure
% hold on
% line(CC_X,CC_Y_max_val,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(CC_X(:,1), mean(CC_Y_max_val,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
%         y1limits = [-0.8 0];
%         y1ticks = [-0.8,0];
%         set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'Before ES','After ES'} ,'box', 'off'); %'fontweight', 'bold', 
%          ylabel('CC evoked', 'FontSize', 20,'fontname', 'arial');
%           title('max peak', 'FontSize',20,'fontname', 'arial')
          %% Paired plot of spont normalized lag0
%  clear CC_Y_lag0 CC_X E
% CC_Y_lag0=cc_stat.spont.norm_lag0';
% CC_X(1,:)=ones(1,size(CC_Y_lag0,2));
% CC_X(2,:)=2*ones(1,size(CC_Y_lag0,2));
% E = std(CC_Y_lag0,0,2);
% figure
% hold on
% line(CC_X,CC_Y_lag0,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(CC_X(:,1), mean(CC_Y_lag0,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
%         y1limits = [-0.75 1.1];
%         y1ticks = [-0.75,-0.5,-0.25,0,0.25,0.5,1];
%         set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'Before ES','After ES'} ,'box', 'off'); %'fontweight', 'bold', 
%          ylabel('CC spont', 'FontSize', 20,'fontname', 'arial');
%          title('normalized lag zero peak', 'FontSize',20,'fontname', 'arial')
         %% Paired plot of spont non-normalized lag0
%  clear CC_Y_lag0 CC_X E
% CC_Y_lag0=cc_stat.spont.lag0';
% CC_X(1,:)=ones(1,size(CC_Y_lag0,2));
% CC_X(2,:)=2*ones(1,size(CC_Y_lag0,2));
% E = std(CC_Y_lag0,0,2);
% figure
% hold on
% line(CC_X,CC_Y_lag0,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(CC_X(:,1), mean(CC_Y_lag0,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
%         y1limits = [-0.5 0.5];
%         y1ticks = [-0.5,0,0.5];
%         set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'Before ES','After ES'} ,'box', 'off'); %'fontweight', 'bold', 
%          ylabel('CC spont', 'FontSize', 20,'fontname', 'arial');
%           title('lag zero peak', 'FontSize',20,'fontname', 'arial')
          %% Paired plot of evoked normalized lag0
%  clear CC_Y_lag0 CC_X E
% CC_Y_lag0=cc_stat.evoked.norm_lag0';
% CC_X(1,:)=ones(1,size(CC_Y_lag0,2));
% CC_X(2,:)=2*ones(1,size(CC_Y_lag0,2));
% E = std(CC_Y_lag0,0,2);
% figure
% hold on
% line(CC_X,CC_Y_lag0,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(CC_X(:,1), mean(CC_Y_lag0,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
%         y1limits = [0.5 1.5];
%         y1ticks = [0.5,1,1.5];
%         set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'Before ES','After ES'} ,'box', 'off'); %'fontweight', 'bold', 
%          ylabel('CC evoked', 'FontSize', 20,'fontname', 'arial');
%           title('normalized lag zero peak', 'FontSize',20,'fontname', 'arial')
          %% Paired plot of evoked non-normalized lag0
%  clear CC_Y_lag0 CC_X E
% CC_Y_lag0=cc_stat.evoked.lag0';
% CC_X(1,:)=ones(1,size(CC_Y_lag0,2));
% CC_X(2,:)=2*ones(1,size(CC_Y_lag0,2));
% E = std(CC_Y_lag0,0,2);
% figure
% hold on
% line(CC_X,CC_Y_lag0,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(CC_X(:,1), mean(CC_Y_lag0,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
%         y1limits = [-0.8 0];
%         y1ticks = [-0.8,0];
%         set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'Before ES','After ES'} ,'box', 'off'); %'fontweight', 'bold', 
%          ylabel('CC evoked', 'FontSize', 20,'fontname', 'arial');
%           title('lag zero peak', 'FontSize',20,'fontname', 'arial')