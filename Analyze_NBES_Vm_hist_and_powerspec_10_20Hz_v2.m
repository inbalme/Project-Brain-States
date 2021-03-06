%% Analyze NBES_response_parameters_10_20Hz_v2
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
save_flag= 1;
print_flag=0;
type_strng = {'ongoing', 'evoked'};
data_vec_all{1} = []; data_vec_all{2} = []; data_vec_residual_all{1} = []; data_vec_residual_all{2} = [];
files_to_analyze =[44,46,48,50,52,56,58,62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,44,46,48,50,52,56,58,62,72,75]; 
% for fileind=1;
    for fileind=1:length(files_to_analyze) ;
    close all
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    clear data_no_spike_no_DC
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
    color_table=[0 0 0;color_table(1:6,:)];
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
%                 stim_temp=stim2_X{2}(:,1); %taking only the location of the first stim in the train
%                 stim2_X=[];
%                 stim2_X{2}=stim_temp; stim2_X{3}=stim_temp;
                
%% low-pass filtering below 300Hz                
                for xx=1:3
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

%% Vm histogram - for spontaneous and evoked activity
 for trace_type= 1:2; %1:2;  %1 for spont., 2 for evoked
     interval=[]; data_vec=[]; data_vec_residual=[]; data_vec_5traces=[]; data_vec_5traces_residual=[]; 
        if trace_type==1 
            x_value=[1,1];
        end
        if trace_type==2 
        x_value=[2:3]; %2:3; %[1,1];
        end
        
       switch trace_type
            case 1
                 start_time = [0.4,5.6]; %[sec] %[0,5]
                 duration = 2.5; %[sec] 
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
        %          duration = 1; %[sec]
                galvano_nstim = Param.facade(6);
                galvano_freq = Param.facade(7);
                duration = galvano_nstim./galvano_freq+0.05;
                end_sample = start_sample+duration.*sf{1}-1;
                interval(:,1) = round(start_sample:end_sample);
                interval(:,2) = interval(:,1);
       end
        %%
        for t=1:size(interval,2);
              data_vec(:,t)=reshape(data_no_spikes{1}(interval(:,t),:,x_value(t)),numel(data_no_spikes{1}(interval(:,t),:,x_value(t))),1); %take only data from x_value=1, and connect all traces into one long vector;
              data_vec_residual(:,t)=data_vec(:,t)-mean(data_vec(:,t));
              data_vec_5traces(:,t)=reshape(data_no_spikes{1}(interval(:,t),1:5,x_value(t)),numel(data_no_spikes{1}(interval(:,t),1:5,x_value(t))),1); %take only data from x_value=1, and connect all traces into one long vector;
              data_vec_5traces_residual(:,t)=data_vec_5traces(:,t)-mean(data_vec_5traces(:,1));
%               nbin=ceil(range(data_vec(:,1))); %set the range according to the interval before ES
            end
             nbin = 30;
             data_vec_median{trace_type}(fileind,:) = median(data_vec);
%              data_vec_residual_median{trace_type}(fileind,:) = median(data_vec_residual);
             data_vec_5prctile{trace_type}(fileind,:)  = prctile(data_vec,5,1);
%              data_vec_mean{trace_type}(fileind,:) = mean(data_vec,1);
%              data_vec_std{trace_type}(fileind,:) = std(data_vec,1);
%              data_vec_CV{trace_type}(fileind,:) =data_vec_std{trace_type}(fileind,:)./ abs(data_vec_mean{trace_type}(fileind,:));
            data_vec_all{trace_type}=[data_vec_all{trace_type}; data_vec_5traces];
            data_vec_residual_all{trace_type}=[data_vec_residual_all{trace_type}; data_vec_5traces_residual];
  switch trace_type
      case 1
        ongoing(fileind).Vm_median=data_vec_median{trace_type}(fileind,:);  
        ongoing(fileind).Vm_5prctile=data_vec_5prctile{trace_type}(fileind,:);   
%         ongoing(fileind).Vm_mean=data_vec_mean{trace_type}(fileind,:);
%         ongoing(fileind).Vm_std=data_vec_std{trace_type}(fileind,:);
%         ongoing(fileind).Vm_CV=data_vec_CV{trace_type}(fileind,:);
        hists.ongoing=ongoing;
      case 2
        evoked(fileind).Vm_median=data_vec_median{trace_type}(fileind,:);  
        evoked(fileind).Vm_5prctile=data_vec_5prctile{trace_type}(fileind,:);   
%         evoked(fileind).Vm_mean=data_vec_mean{trace_type}(fileind,:);
%         evoked(fileind).Vm_std=data_vec_std{trace_type}(fileind,:);
%         evoked(fileind).Vm_CV=data_vec_CV{trace_type}(fileind,:);
        hists.evoked=evoked;
  end
 %% plotting Vm histogram for each cell:
 if print_flag==1;
            Fig1=figure;
            hold on
            [ncounts(1,:),nbins(1,:)]=hist(data_vec(:,1),nbin);
            ncounts(1,:)=ncounts(1,:)./length(data_vec(:,1));
            [ncounts(2,:),nbins(2,:)]=hist(data_vec(:,2),nbin);
            ncounts(2,:)=ncounts(2,:)./length(data_vec(:,2));
            bar(nbins(1,:),ncounts(1,:));
            bar(nbins(2,:),ncounts(2,:));
            
%             obj = gmdistribution.fit(data_vec(:,1),2);
            h= findobj(gca,'Type','patch');
            set(h(1),'FaceColor',[0 0 1],'EdgeColor','w')
            set(h(2),'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
            ylim1=get(gca,'ylim');
            median1(1)=line([data_vec_median{trace_type}(fileind,1) data_vec_median{trace_type}(fileind,1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
            median1(2)=line([data_vec_median{trace_type}(fileind,2) data_vec_median{trace_type}(fileind,2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
            hold off
%             pause
%% plotting mean subtracted-Vm histogram for each cell:
            Fig2=figure;
            hold on
%             hist(data_vec_residual(:,1),nbin)
%             hist(data_vec_residual(:,2),nbin)
            [ncounts_residual(1,:),nbins_residual(1,:)]=hist(data_vec_residual(:,1),nbin);
            ncounts_residual(1,:)=ncounts_residual(1,:)./length(data_vec_residual(:,1));
            [ncounts_residual(2,:),nbins_residual(2,:)]=hist(data_vec_residual(:,2),nbin);
            ncounts_residual(2,:)=ncounts_residual(2,:)./length(data_vec_residual(:,2));
            bar(nbins_residual(1,:),ncounts_residual(1,:));
            bar(nbins_residual(2,:),ncounts_residual(2,:));
%             obj = gmdistribution.fit(data_vec(:,1),2);
            h= findobj(gca,'Type','patch');
            set(h(1),'FaceColor',[0 0 1],'EdgeColor','w')
            set(h(2),'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
            ylim1=get(gca,'ylim');
            median_residual(1)=line([data_vec_residual_median{trace_type}(fileind,1) data_vec_residual_median{trace_type}(fileind,1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
            median_residual(2)=line([data_vec_residual_median{trace_type}(fileind,2) data_vec_residual_median{trace_type}(fileind,2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
            hold off
%             pause
            
       clear data_no_spike_no_DC data_vec data_vec_residual
       %% saving figures
            if save_flag==1;
                cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'  
                    saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_', type_strng{trace_type},'.fig']) 
                    print(Fig1,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_', type_strng{trace_type}],'-dpng','-r600','-opengl') 
                    saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_residuals_',type_strng{trace_type},'.fig']) 
                    print(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_residuals_', type_strng{trace_type}],'-dpng','-r600','-opengl') 
            end
 end
    end %temp for file loop
    end 
             
%% plotting a Vm histogram for all the cells 
if print_flag==1;
for trace_type=1:2;
             data_vec_all_median{trace_type} = median(data_vec_all{trace_type});
             data_vec_residual_all_median{trace_type} = median(data_vec_residual_all{trace_type});
        Fig_all=figure;
        hold on
%             hist(data_vec_all(:,1),nbin) %problem because there are different number of traces from each cell - different weight for each cell.
%             hist(data_vec_all(:,2),nbin)
            [ncounts_all(1,:),nbins_all(1,:)]=hist(data_vec_all{trace_type}(:,1),nbin);
            ncounts_all(1,:)=ncounts_all(1,:)./length(data_vec_all{trace_type}(:,1));
            [ncounts_all(2,:),nbins_all(2,:)]=hist(data_vec_all{trace_type}(:,2),nbin);
            ncounts_all(2,:)=ncounts_all(2,:)./length(data_vec_all{trace_type}(:,2));
            bar(nbins_all(1,:),ncounts_all(1,:));
            bar(nbins_all(2,:),ncounts_all(2,:));
            h= findobj(gca,'Type','patch');
            set(h(1),'FaceColor',[0 0 1],'EdgeColor','w')
            set(h(2),'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
            ylim1=get(gca,'ylim');
            median1=line([data_vec_all_median{trace_type}(1) data_vec_all_median{trace_type}(1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
            median2=line([data_vec_all_median{trace_type}(2) data_vec_all_median{trace_type}(2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
         hold off
         
  % plotting a mean subtracted (residuals) Vm histogram for all the cells             
        Fig_all_residual=figure;
        hold on
%             hist(data_vec_residual_all(:,1),nbin) %problem because there are different number of traces from each cell - different weight for each cell.
%             hist(data_vec_residual_all(:,2),nbin)
            [ncounts_residual_all(1,:),nbins_residual_all(1,:)]=hist(data_vec_residual_all{trace_type}(:,1),nbin);
            ncounts_residual_all(1,:)=ncounts_residual_all(1,:)./length(data_vec_residual_all{trace_type}(:,1));
            [ncounts_residual_all(2,:),nbins_residual_all(2,:)]=hist(data_vec_residual_all{trace_type}(:,2),nbin);
            ncounts_residual_all(2,:)=ncounts_residual_all(2,:)./length(data_vec_residual_all{trace_type}(:,2));
            bar(nbins_residual_all(1,:),ncounts_residual_all(1,:));
            bar(nbins_residual_all(2,:),ncounts_residual_all(2,:));
            h= findobj(gca,'Type','patch');
            set(h(1),'FaceColor',[0 0 1],'EdgeColor','w')
            set(h(2),'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
            ylim1=get(gca,'ylim');
            median1=line([data_vec_residual_all_median{trace_type}(1) data_vec_residual_all_median{trace_type}(1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
            median2=line([data_vec_residual_all_median{trace_type}(2) data_vec_residual_all_median{trace_type}(2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
         hold off
         
          if save_flag==1;
                cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'  
                    saveas(Fig_all, ['Vm_histogram_all_', type_strng{trace_type},'.fig']) 
                    print(Fig_all,['Vm_histogram_all_', type_strng{trace_type}],'-dpng','-r600','-opengl') 
                    saveas(Fig_all_residual,['Vm_histogram_residuals_all_', type_strng{trace_type},'.fig']) 
                    print(Fig_all_residual,['Vm_histogram_residuals_all_',type_strng{trace_type}],'-dpng','-r600','-opengl') 
          end
end
end
%%   statistics       
          %stats for Vm median when x_value=1
          %it is meaningless to take the median or mean or std of the connected traces. the correct thing is to take the median/mean/std of each trace and then average.
  for trace_type=1:2;

       hists_stat(trace_type).Vm_median=data_vec_median{trace_type};
        norm_Vm_median_noES= data_vec_median{trace_type}(:,1)./ data_vec_median{trace_type}(:,1);
        norm_Vm_median_ES= data_vec_median{trace_type}(:,2)./ data_vec_median{trace_type}(:,1);
       hists_stat(trace_type).norm_Vm_median=[norm_Vm_median_noES, norm_Vm_median_ES];
       hists_stat(trace_type).Vm_median_m=mean(hists_stat(trace_type).Vm_median,1);
       hists_stat(trace_type).Vm_median_std=std(hists_stat(trace_type).Vm_median,0,1);
        %testing for normal distribution
        diff_val= data_vec_median{trace_type}(:,2)- data_vec_median{trace_type}(:,1);
        [hists_stat(trace_type).lillietest_h_Vm_median,hists_stat(trace_type).lillietest_p_Vm_median] = lillietest(diff_val);
        %paired ttest 
        [hists_stat(trace_type).ttest_h_Vm_median,hists_stat(trace_type).ttest_p_Vm_median]= ttest(hists_stat(trace_type).Vm_median(:,1),hists_stat(trace_type).Vm_median(:,2));
        [hists_stat(trace_type).ttest_h_norm_Vm_median,hists_stat(trace_type).ttest_p_norm_Vm_median]= ttest(hists_stat(trace_type).norm_Vm_median(:,1),hists_stat(trace_type).norm_Vm_median(:,2));
        [hists_stat(trace_type).wilcoxon_p_Vm_median,hists_stat(trace_type).wilcoxon_h_Vm_median]= signrank(hists_stat(trace_type).Vm_median(:,1),hists_stat(trace_type).Vm_median(:,2));
        clear norm_Vm_median_noES norm_Vm_median_ES diff_val
        
           %stats for Vm 5 percentile when x_value=1
       hists_stat(trace_type).Vm_5prctile=[data_vec_5prctile{trace_type}(:,:)];
        norm_Vm_5prctile_noES= data_vec_5prctile{trace_type}(:,1)./ data_vec_5prctile{trace_type}(:,1);
        norm_Vm_5prctile_ES= data_vec_5prctile{trace_type}(:,2)./ data_vec_5prctile{trace_type}(:,1);
       hists_stat(trace_type).norm_Vm_5prctile=[norm_Vm_5prctile_noES, norm_Vm_5prctile_ES];
       hists_stat(trace_type).Vm_5prctile_m=mean(hists_stat(trace_type).Vm_5prctile,1);
       hists_stat(trace_type).Vm_5prctile_std=std(hists_stat(trace_type).Vm_5prctile,0,1);
        %testing for normal distribution
        diff_val= data_vec_5prctile{trace_type}(:,2)- data_vec_5prctile{trace_type}(:,1);
        [hists_stat(trace_type).lillietest_h_Vm_5prctile,hists_stat(trace_type).lillietest_p_Vm_5prctile] = lillietest(diff_val);
        %paired ttest 
        [hists_stat(trace_type).ttest_h_Vm_5prctile,hists_stat(trace_type).ttest_p_Vm_5prctile]= ttest(hists_stat(trace_type).Vm_5prctile(:,1),hists_stat(trace_type).Vm_5prctile(:,2));
        [hists_stat(trace_type).ttest_h_norm_Vm_5prctile,hists_stat(trace_type).ttest_p_norm_Vm_5prctile]= ttest(hists_stat(trace_type).norm_Vm_5prctile(:,1),hists_stat(trace_type).norm_Vm_5prctile(:,2));
        [hists_stat(trace_type).wilcoxon_p_Vm_5prctile,hists_stat(trace_type).wilcoxon_h_Vm_5prctile]= signrank(hists_stat(trace_type).Vm_5prctile(:,1),hists_stat(trace_type).Vm_5prctile(:,2));
        clear norm_Vm_5prctile_noES norm_Vm_5prctile_ES diff_val
%%        
%         %stats for Vm mean when x_value=1
%        hists_stat(trace_type).Vm_mean=[data_vec_mean{trace_type}(:,:)];
%         norm_Vm_mean_noES= data_vec_mean{trace_type}(:,1)./ data_vec_mean{trace_type}(:,1);
%         norm_Vm_mean_ES= data_vec_mean{trace_type}(:,2)./ data_vec_mean{trace_type}(:,1);
%        hists_stat(trace_type).norm_Vm_mean=[norm_Vm_mean_noES, norm_Vm_mean_ES];
%        hists_stat(trace_type).Vm_mean_m=mean(hists_stat(trace_type).Vm_mean,1);
%        hists_stat(trace_type).Vm_mean_std=std(hists_stat(trace_type).Vm_mean,0,1);
%         %testing for normal distribution
%         diff_val= data_vec_mean{trace_type}(:,2)- data_vec_mean{trace_type}(:,1);
%         [hists_stat(trace_type).lillietest_h_Vm_mean,hists_stat(trace_type).lillietest_p_Vm_mean] = lillietest(diff_val);
%         %paired ttest 
%         [hists_stat(trace_type).ttest_h_Vm_mean,hists_stat(trace_type).ttest_p_Vm_mean]= ttest(hists_stat(trace_type).Vm_mean(:,1),hists_stat(trace_type).Vm_mean(:,2));
%         [hists_stat(trace_type).ttest_h_norm_Vm_mean,hists_stat(trace_type).ttest_p_norm_Vm_mean]= ttest(hists_stat(trace_type).norm_Vm_mean(:,1),hists_stat(trace_type).norm_Vm_mean(:,2));
%         [hists_stat(trace_type).wilcoxon_p_Vm_mean,hists_stat(trace_type).wilcoxon_h_Vm_mean]= signrank(hists_stat(trace_type).Vm_mean(:,1),hists_stat(trace_type).Vm_mean(:,2));
%         clear norm_Vm_mean_noES norm_Vm_mean_ES diff_val
%         
%          %stats for Vm std when x_value=1
%        hists_stat(trace_type).Vm_std=[data_vec_std{trace_type}(:,:)];
%         norm_Vm_std_noES= data_vec_std{trace_type}(:,1)./ data_vec_std{trace_type}(:,1);
%         norm_Vm_std_ES= data_vec_std{trace_type}(:,2)./ data_vec_std{trace_type}(:,1);
%        hists_stat(trace_type).norm_Vm_std=[norm_Vm_std_noES, norm_Vm_std_ES];
%        hists_stat(trace_type).Vm_std_m=std(hists_stat(trace_type).Vm_std,1);
%        hists_stat(trace_type).Vm_std_std=std(hists_stat(trace_type).Vm_std,0,1);
%         %testing for normal distribution
%         diff_val= data_vec_std{trace_type}(:,2)- data_vec_std{trace_type}(:,1);
%         [hists_stat(trace_type).lillietest_h_Vm_std,hists_stat(trace_type).lillietest_p_Vm_std] = lillietest(diff_val);
%         %paired ttest 
%         [hists_stat(trace_type).ttest_h_Vm_std,hists_stat(trace_type).ttest_p_Vm_std]= ttest(hists_stat(trace_type).Vm_std(:,1),hists_stat(trace_type).Vm_std(:,2));
%         [hists_stat(trace_type).ttest_h_norm_Vm_std,hists_stat(trace_type).ttest_p_norm_Vm_std]= ttest(hists_stat(trace_type).norm_Vm_std(:,1),hists_stat(trace_type).norm_Vm_std(:,2));
%         [hists_stat(trace_type).wilcoxon_p_Vm_std,hists_stat(trace_type).wilcoxon_h_Vm_std]= signrank(hists_stat(trace_type).Vm_std(:,1),hists_stat(trace_type).Vm_std(:,2));
%         clear norm_Vm_std_noES norm_Vm_std_ES diff_val
%         
%          %stats for Vm CV when x_value=1
%        hists_stat(trace_type).Vm_CV=data_vec_CV{trace_type};
%         norm_Vm_CV_noES= data_vec_CV{trace_type}(:,1)./ data_vec_CV{trace_type}(:,1);
%         norm_Vm_CV_ES= data_vec_CV{trace_type}(:,2)./ data_vec_CV{trace_type}(:,1);
%        hists_stat(trace_type).norm_Vm_CV=[norm_Vm_CV_noES, norm_Vm_CV_ES];
%        hists_stat(trace_type).Vm_CV_m=mean(hists_stat(trace_type).Vm_CV,1);
%        hists_stat(trace_type).Vm_CV_std=std(hists_stat(trace_type).Vm_CV,0,1);
%         %testing for normal distribution
%         diff_val= data_vec_CV{trace_type}(:,2)- data_vec_CV{trace_type}(:,1);
%         [hists_stat(trace_type).lillietest_h_Vm_CV,hists_stat(trace_type).lillietest_p_Vm_CV] = lillietest(diff_val);
%         %paired ttest 
%         [hists_stat(trace_type).ttest_h_Vm_CV,hists_stat(trace_type).ttest_p_Vm_CV]= ttest(hists_stat(trace_type).Vm_CV(:,1),hists_stat(trace_type).Vm_CV(:,2));
%         [hists_stat(trace_type).ttest_h_norm_Vm_CV,hists_stat(trace_type).ttest_p_norm_Vm_CV]= ttest(hists_stat(trace_type).norm_Vm_CV(:,1),hists_stat(trace_type).norm_Vm_CV(:,2));
%         [hists_stat(trace_type).wilcoxon_p_Vm_CV,hists_stat(trace_type).wilcoxon_h_Vm_CV]= signrank(hists_stat(trace_type).Vm_CV(:,1),hists_stat(trace_type).Vm_CV(:,2));
%         clear norm_Vm_CV_noES norm_Vm_CV_ES diff_val
  end
%%   Plotting Vm median  
if print_flag==1;
  for trace_type=1:2;

Vm_median_Y=hists_stat(trace_type).Vm_median';
Vm_median_X(1,:)=ones(1,size(Vm_median_Y,2));
Vm_median_X(2,:)=2*ones(1,size(Vm_median_Y,2));
E =hists_stat(trace_type).Vm_median_std;

g1=figure;
hold on
line(Vm_median_X,Vm_median_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Vm_median_X(:,1), mean(Vm_median_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Median Vm [mV]', 'FontSize', 28,'fontname', 'arial');
        title([type_strng{trace_type},' activity n=' num2str(length(files_to_analyze)) ', p=' num2str(hists_stat(trace_type).wilcoxon_p_Vm_median)] ,'FontSize', 28,'fontname', 'arial');
  
        %% plotting Vm lower 5 percentile

        Vm_5prctile_Y=hists_stat(trace_type).Vm_5prctile';
Vm_5prctile_X(1,:)=ones(1,size(Vm_5prctile_Y,2));
Vm_5prctile_X(2,:)=2*ones(1,size(Vm_5prctile_Y,2));
E =hists_stat(trace_type).Vm_5prctile_std;
if hists_stat(trace_type).wilcoxon_p_Vm_5prctile >0.05 
    asterisk_sp='n.s.';
else if hists_stat(trace_type).wilcoxon_p_Vm_5prctile<0.05 && hists_stat(trace_type).wilcoxon_p_Vm_5prctile>0.01
    asterisk_sp='*';
    else if hists_stat(trace_type).wilcoxon_p_Vm_5prctile<0.01 && hists_stat(trace_type).wilcoxon_p_Vm_5prctile>0.001
            asterisk_sp='**';
    else if hists_stat(trace_type).wilcoxon_p_Vm_5prctile<0.001
             asterisk_sp='***';
        end
        end
    end
end
linex=[1;2];
% my=max(max(hists_stat(trace_type).Vm_5prctile))*1.1; 
liney=[-32 ;-32 ];

g2=figure;
hold on
line(Vm_5prctile_X,Vm_5prctile_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Vm_5prctile_X(:,1), mean(Vm_5prctile_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,-32,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim', [-85,-25], 'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('lower 5 percentile [mV]', 'FontSize', 28,'fontname', 'arial');
        title([type_strng{trace_type},' activity n=' num2str(length(files_to_analyze)) ', p=' num2str(hists_stat(trace_type).wilcoxon_p_Vm_5prctile)] ,'FontSize', 28,'fontname', 'arial');
  
 %%        %save figures  
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'; 
fn1=['Vm_median_', type_strng{trace_type}];
saveas(g1,['Vm_median_', type_strng{trace_type},'.fig']); 
print(g1,fn1,'-dpng','-r600','-opengl') 
fn2=['Vm_5prcentile_', type_strng{trace_type}];
saveas(g2,['Vm_5prcentile_', type_strng{trace_type},'.fig']); 
print(g2,fn2,'-dpng','-r600','-opengl') 
  end
end
 %% Plot power spectrum
% start_time = [0.4,6]; %[sec]
% duration = 2.5; %[sec] 
% x_value = 1;
% Y_abs = []; f = [];  interval = []; 
% 
% for t=1:length(start_time);
% DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = [];
% 
% start_sample = start_time(t).*sf{channel};
% if start_time(t)==0
%     start_sample = 1;
% end
% end_sample = start_sample+duration.*sf{channel}-1;
% interval(:,t) = start_sample:end_sample;
% spec_mat = raw_data{channel}(interval(:,t),:,x_value);
% DC= mean(spec_mat,1);
% l=size(spec_mat,1);
% l=l/2-1;
% DC=wextend('addrow', 'per', DC, l);
% spec_mat_noDC = spec_mat-DC;
% 
% [Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf{channel},0,0);
% end

%% Plot power spectrum without shaded error bar
% figure
%         xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
%         ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
%         set(gca,'xscale','log');
%         set(gca,'yscale','log');
%  for t=1:length(start_time);
%       hold on
%         plot(f(:,t),mean(Y_abs(:,:,t),2),'color', color_table(t,:))        
%         xlim([0 1000]); ylim([0.001 10])
% %         title('')
%         set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2)        
%       hold off
%  end
%  pause
%      end %temp file loop

 %% saving figures of power spectrum
%             if save_flag==1;
%                 cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'  
%                     saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_power_spectrum.fig']) 
%                     print(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_power_spectrum'],'-dpng','-r600','-opengl') 
%             end
%             end %temp file loop
%% visualizing raw data traces for ongoing and evoked NB- and NB+ together  
%  figure
%  hold on
%  plot(interval(:,1).*dt,raw_data{1}(interval(:,1),1:5,1),'k');
%   plot(interval(:,1).*dt,raw_data{1}(interval(:,2),1:5,1),'b');
%   hold off
%   pause
% 
%  figure
%  hold on
%  plot((stim2_X{2}(1,1):(stim2_X{2}(1,1)+sf{1}-1)).*dt,raw_data{1}(stim2_X{2}(1,1):stim2_X{2}(1,1)+sf{1}-1,1:3,2),'k');
%   plot((stim2_X{2}(1,1):(stim2_X{2}(1,1)+sf{1}-1)).*dt,raw_data{1}(stim2_X{2}(1,1):stim2_X{2}(1,1)+sf{1}-1,1:3,3),'b');
%   hold off
%   pause
%     end %temp

 
 %%
   cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'  
if save_flag==1;
   % filename='hists'; 
filename='hists'; 
save(filename, 'files_to_analyze', 'hists', 'hists_stat')
end