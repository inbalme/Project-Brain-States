%% Analyze Vm histograms and powerspectrum
% This file was created on 15/1/2017 based on
% "Analyze_NBES_Vm_hist_and_powerspec_10_20Hz_v2'
%This file is used for the analysis of files created with
%extract_NBES_Data_v3 or Extract_ChAT_Data_v3

%       protocol 5 (NBES+galvnao (5 x-values))
%       protocol 6 (ChAT intracellular current injections light+galvnao (15 x-values))
%       protocol 8 - ChAT intracellular: light only (3 x-values): 1) train
%       20 Hz, 10ms ondur, 50 pulses. 2) single pulse 5 sec ondur. 3)single pulse 10 ms
%       protocol 9 - IC: I-V curve
%       protocol 10 (NBES/ChAT+galvnao train+test (3 x-values))
%                 11 - NBES: 3 currents+galvano train+test (6 x-values)
%% for opening workspace saved 
clear all
close all
global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data Ch2_data stim2_X stim1_X 
 global exp_type
exp_type=2; %1-NBES, 2-ChAT

trace_type_input=1; %
analyze_time_before_train=0.1;
analyze_train_only_flag=0;
save_flag= 0;
print_flag=1;
norm_flag=0;
BP50HzLFP_flag=1; %removing 50Hz noise from LFP signal
BP50HzVm_flag=1; %removing 50Hz noise from Vm signal
BPLFP_flag=0; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=[1,200]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=0; %filtering LFP and Vm same as LFP was filtered in the multiclamp
bp_manual_Vm=[0,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
type_strng = {'ongoing', 'evoked'};
data_vec_all{1} = []; data_vec_all{2} = []; data_vec_residual_all{1} = []; data_vec_residual_all{2} = [];
%%
switch exp_type
    case 1
        files_to_analyze =[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75,82,84]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,44,46,48,50,52,56,58,62,72,75]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};
        path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm Histograms and Powerspec\Spontaneous';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end        

    case 2
        files_to_analyze =[74,76,77,80,82,84,87]; %
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light On', 'Light Off'};
        path_output= 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm Histograms and Powerspec\Spontaneous';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end        
end
        
    for fileind=1:length(files_to_analyze) ;
    close all
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
   %%
   Ch2_data= raw_data{3}./20; %dividing by the LFP gain           
    current_data=data_no_spikes{channel};    
    galvano_nstim = Param.facade(6);
    galvano_freq = Param.facade(7);

data_preprocessing   

 if ~isempty(current_data_filt)
     current_data=current_data_filt;
 end
 
 clear color_table
        whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;
        switch exp_type
            case 1 
                color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
            case 2
                color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
        end
 
%% Vm histogram - for spontaneous and evoked activity
 for trace_type=  trace_type_input   %1 or 3 for spont., 2 for evoked
     interval=[]; data_vec=[]; data_vec_residual=[]; data_vec_5traces=[]; data_vec_5traces_residual=[]; 
 intervals_to_analyze         
        %%
        for t=1:size(interval,2);
              data_vec(:,t)=reshape(current_data(interval(:,t),:,x_value(t)),numel(current_data(interval(:,t),:,x_value(t))),1); %take only data from x_value=1, and connect all traces into one long vector;
              data_vec_residual(:,t)=data_vec(:,t)-mean(data_vec(:,t));
              data_vec_5traces(:,t)=reshape(current_data(interval(:,t),1:5,x_value(t)),numel(current_data(interval(:,t),1:5,x_value(t))),1); %take only data from x_value=1, and connect all traces into one long vector;
              data_vec_5traces_residual(:,t)=data_vec_5traces(:,t)-mean(data_vec_5traces(:,1));
%               nbin=ceil(range(data_vec(:,1))); %set the range according to the interval before ES
            end
             nbin = 30;
             data_vec_median{trace_type}(fileind,:) = median(data_vec);
             data_vec_residual_median{trace_type}(fileind,:) = median(data_vec_residual);
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
            h1=bar(nbins(1,:),ncounts(1,:));
            h2=bar(nbins(2,:),ncounts(2,:));
            
%             obj = gmdistribution.fit(data_vec(:,1),2);
            set(h2,'FaceColor',[0 0 1],'EdgeColor','w')
            set(h1,'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
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
            h1=bar(nbins_residual(1,:),ncounts_residual(1,:));
            h2=bar(nbins_residual(2,:),ncounts_residual(2,:));
%             obj = gmdistribution.fit(data_vec(:,1),2);
            set(h2,'FaceColor',[0 0 1],'EdgeColor','w')
            set(h1,'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
            ylim1=get(gca,'ylim');
            median_residual(1)=line([data_vec_residual_median{trace_type}(fileind,1) data_vec_residual_median{trace_type}(fileind,1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
            median_residual(2)=line([data_vec_residual_median{trace_type}(fileind,2) data_vec_residual_median{trace_type}(fileind,2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
            hold off
%             pause
            
       clear data_no_spike_no_DC data_vec data_vec_residual
       %% saving figures
            if save_flag==1;
  cd(path_output)
                    saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_', type_strng{trace_type},'.fig']) 
                    print(Fig1,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_', type_strng{trace_type}],'-dpng','-r600','-opengl') 
                    saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_residuals_',type_strng{trace_type},'.fig']) 
                    print(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_residuals_', type_strng{trace_type}],'-dpng','-r600','-opengl') 
            end
 end
    end %temp for file loop
    end 
             
%% plotting a Vm histogram for all the cells 
% if print_flag==1;
% for trace_type=1:2;
%              data_vec_all_median{trace_type} = median(data_vec_all{trace_type});
%              data_vec_residual_all_median{trace_type} = median(data_vec_residual_all{trace_type});
%         Fig_all=figure;
%         hold on
% %             hist(data_vec_all(:,1),nbin) %problem because there are different number of traces from each cell - different weight for each cell.
% %             hist(data_vec_all(:,2),nbin)
%             [ncounts_all(1,:),nbins_all(1,:)]=hist(data_vec_all{trace_type}(:,1),nbin);
%             ncounts_all(1,:)=ncounts_all(1,:)./length(data_vec_all{trace_type}(:,1));
%             [ncounts_all(2,:),nbins_all(2,:)]=hist(data_vec_all{trace_type}(:,2),nbin);
%             ncounts_all(2,:)=ncounts_all(2,:)./length(data_vec_all{trace_type}(:,2));
%             h1=bar(nbins_all(1,:),ncounts_all(1,:));
%             h2=bar(nbins_all(2,:),ncounts_all(2,:));
%             set(h2,'FaceColor',[0 0 1],'EdgeColor','w')
%             set(h1,'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
%             ylim1=get(gca,'ylim');
%             median1=line([data_vec_all_median{trace_type}(1) data_vec_all_median{trace_type}(1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
%             median2=line([data_vec_all_median{trace_type}(2) data_vec_all_median{trace_type}(2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
%          hold off
%          
%   % plotting a mean subtracted (residuals) Vm histogram for all the cells             
%         Fig_all_residual=figure;
%         hold on
% %             hist(data_vec_residual_all(:,1),nbin) %problem because there are different number of traces from each cell - different weight for each cell.
% %             hist(data_vec_residual_all(:,2),nbin)
%             [ncounts_residual_all(1,:),nbins_residual_all(1,:)]=hist(data_vec_residual_all{trace_type}(:,1),nbin);
%             ncounts_residual_all(1,:)=ncounts_residual_all(1,:)./length(data_vec_residual_all{trace_type}(:,1));
%             [ncounts_residual_all(2,:),nbins_residual_all(2,:)]=hist(data_vec_residual_all{trace_type}(:,2),nbin);
%             ncounts_residual_all(2,:)=ncounts_residual_all(2,:)./length(data_vec_residual_all{trace_type}(:,2));
%             h1=bar(nbins_residual_all(1,:),ncounts_residual_all(1,:));
%            h2= bar(nbins_residual_all(2,:),ncounts_residual_all(2,:));
%             set(h2,'FaceColor',[0 0 1],'EdgeColor','w')
%             set(h1,'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
%             ylim1=get(gca,'ylim');
%             median1=line([data_vec_residual_all_median{trace_type}(1) data_vec_residual_all_median{trace_type}(1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
%             median2=line([data_vec_residual_all_median{trace_type}(2) data_vec_residual_all_median{trace_type}(2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
%          hold off
         
%           if save_flag==1;
%                 cd(path_output)
%                     saveas(Fig_all, ['Vm_histogram_all_', type_strng{trace_type},'.fig']) 
%                     print(Fig_all,['Vm_histogram_all_', type_strng{trace_type}],'-dpng','-r600','-opengl') 
%                     saveas(Fig_all_residual,['Vm_histogram_residuals_all_', type_strng{trace_type},'.fig']) 
%                     print(Fig_all_residual,['Vm_histogram_residuals_all_',type_strng{trace_type}],'-dpng','-r600','-opengl') 
%           end
% end
% end
%%   statistics       
          %stats for Vm median when x_value=1
          %it is meaningless to take the median or mean or std of the connected traces. the correct thing is to take the median/mean/std of each trace and then average.
  for trace_type=trace_type_input;

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
  for trace_type=trace_type_input;

Vm_median_Y=hists_stat(trace_type).Vm_median';
Vm_median_X(1,:)=ones(1,size(Vm_median_Y,2));
Vm_median_X(2,:)=2*ones(1,size(Vm_median_Y,2));
E =hists_stat(trace_type).Vm_median_std;
if hists_stat(trace_type).lillietest_h_Vm_median ==0 %data is normally dist. and we can use ttest
    p_median= hists_stat(trace_type).ttest_p_Vm_median;
else %data is not dist. noarmally and we need non-parametric test
     p_median=hists_stat(trace_type).wilcoxon_p_Vm_median;
end
        if p_median >0.05 
    asterisk_sp='n.s.';
else if p_median<0.05 && p_median>0.01
    asterisk_sp='*';
    else if p_median<0.01 && p_median>0.001
            asterisk_sp='**';
    else if p_median<0.001
             asterisk_sp='***';
        end
    end
end
        end
my=max(max(hists_stat(trace_type).Vm_median)); 
g1=figure;
hold on
line(Vm_median_X,Vm_median_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Vm_median_X(:,1), mean(Vm_median_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',20)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Median Vm [mV]', 'FontSize', 28,'fontname', 'arial');
        title([type_strng{trace_type},' activity n=' num2str(length(files_to_analyze)) ', p=' num2str(hists_stat(trace_type).wilcoxon_p_Vm_median)] ,'FontSize', 28,'fontname', 'arial');
  
        %% plotting Vm lower 5 percentile

        Vm_5prctile_Y=hists_stat(trace_type).Vm_5prctile';
Vm_5prctile_X(1,:)=ones(1,size(Vm_5prctile_Y,2));
Vm_5prctile_X(2,:)=2*ones(1,size(Vm_5prctile_Y,2));
E =hists_stat(trace_type).Vm_5prctile_std;
if hists_stat(trace_type).lillietest_h_Vm_5prctile ==0 %data is normally dist. and we can use ttest
    p_5prctile= hists_stat(trace_type).ttest_p_Vm_5prctile;
else %data is not dist. noarmally and we need non-parametric test
     p_5prctile=hists_stat(trace_type).wilcoxon_p_Vm_5prctile;
end
    
    if p_5prctile >0.05 
    asterisk_sp='n.s.';
else if p_5prctile<0.05 && p_5prctile>0.01
    asterisk_sp='*';
    else if p_5prctile<0.01 && p_5prctile>0.001
            asterisk_sp='**';
    else if p_5prctile<0.001
             asterisk_sp='***';
        end
        end
    end
end
linex=[1;2];
my=max(max(hists_stat(trace_type).Vm_5prctile)); 
liney=[-32 ;-32 ];

g2=figure;
hold on
line(Vm_5prctile_X,Vm_5prctile_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Vm_5prctile_X(:,1), mean(Vm_5prctile_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',20)
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks, 'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 'ylim', [-85,-25],
        ylabel('lower 5 percentile [mV]', 'FontSize', 28,'fontname', 'arial');
        title([type_strng{trace_type},' activity n=' num2str(length(files_to_analyze)) ', p=' num2str(hists_stat(trace_type).wilcoxon_p_Vm_5prctile)] ,'FontSize', 28,'fontname', 'arial');
  
 %%        %save figures  
 if save_flag==1
  cd(path_output)
fn1=['Vm_median_', type_strng{trace_type}];
saveas(g1,fn1,'fig'); 
print(g1,fn1,'-dpng','-r600','-opengl') 
fn2=['Vm_5prcentile_', type_strng{trace_type}];
saveas(g2,fn2,'fig'); 
print(g2,fn2,'-dpng','-r600','-opengl') 
 end
  end
end
 %% Plot power spectrum

Y_abs = []; f = [];  interval = []; 

for t=1:length(start_time);
DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = [];

spec_mat = raw_data{channel}(interval(:,t),:,x_value);
DC= mean(spec_mat,1);
l=size(spec_mat,1);
l=l/2-1;
DC=wextend('addrow', 'per', DC, l);
spec_mat_noDC = spec_mat-DC;

[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf{channel},0,0);
end

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
%                 cd(path_output)
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
   cd(path_output)
if save_flag==1;
   % filename='hists'; 
filename='hists'; 
save(filename, 'files_to_analyze', 'hists', 'hists_stat')
end