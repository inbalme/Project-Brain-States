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
print_flag=1;
data_vec_all = []; data_vec_residual_all = [];
files_to_analyze =[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,44,46,48,50,52,56,58,62,72,75]; 
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
                stim_temp=stim2_X{2}(:,1);
                stim2_X=[];
                stim2_X{2}=stim_temp; stim2_X{3}=stim_temp;
                
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
%  start_time = [0.4,6]; %[sec] %[0,5]
%          duration = 2.5; %[sec] 
%          xx=1;
%             for t=1:length(start_time);
%              start_sample(:,t) = ceil(start_time(t).*sf{1});
%                 if start_time(t)==0
%                     start_sample(:,t) = 1;
%                 end
%               end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
%               interval(:,t) = start_sample(:,t):end_sample(:,t);  
%               data_vec(:,t)=reshape(data_no_spikes{1}(interval(:,t),:,1),numel(data_no_spikes{1}(interval(:,t),:,1)),1); %take only data from x_value=1, and connect all traces into one long vector;
%               data_vec_residual(:,t)=data_vec(:,t)-mean(data_vec(:,1));
%               data_vec_5traces(:,t)=reshape(data_no_spikes{1}(interval(:,t),1:5,1),numel(data_no_spikes{1}(interval(:,t),1:5,1)),1); %take only data from x_value=1, and connect all traces into one long vector;
%               data_vec_5traces_residual(:,t)=data_vec_5traces(:,t)-mean(data_vec_5traces(:,1));
% %               nbin=ceil(range(data_vec(:,1))); %set the range according to the interval before ES
%             end
%              nbin = 30;
%              data_vec_median(fileind,:) = median(data_vec);
%              data_vec_residual_median(fileind,:) = median(data_vec_residual);
%              data_vec_5prctile(fileind,:)  = prctile(data_vec,5,1);
%              data_vec_mean(fileind,:) = mean(data_vec,1);
%              data_vec_std(fileind,:) = std(data_vec,1);
%              data_vec_CV(fileind,:) =data_vec_std(fileind,:)./ abs(data_vec_mean(fileind,:));
%             data_vec_all=[data_vec_all; data_vec_5traces];
%             data_vec_residual_all=[data_vec_residual_all; data_vec_5traces_residual];
%             
%         ongoing(fileind).Vm_median=data_vec_median(fileind,:);  
%         ongoing(fileind).Vm_5prctile=data_vec_5prctile(fileind,:);   
%         ongoing(fileind).Vm_mean=data_vec_mean(fileind,:);
%         ongoing(fileind).Vm_std=data_vec_std(fileind,:);
%         ongoing(fileind).Vm_CV=data_vec_CV(fileind,:);
 %% plotting Vm histogram for each cell:
%             Fig1=figure;
%             hold on
%             [ncounts(1,:),nbins(1,:)]=hist(data_vec(:,1),nbin);
%             ncounts(1,:)=ncounts(1,:)./length(data_vec(:,1));
%             [ncounts(2,:),nbins(2,:)]=hist(data_vec(:,2),nbin);
%             ncounts(2,:)=ncounts(2,:)./length(data_vec(:,2));
%             bar(nbins(1,:),ncounts(1,:));
%             bar(nbins(2,:),ncounts(2,:));
%             
%             obj = gmdistribution.fit(data_vec(:,1),2);
%             h= findobj(gca,'Type','patch');
%             set(h(1),'FaceColor',[0 0 1],'EdgeColor','w')
%             set(h(2),'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
%             ylim1=get(gca,'ylim');
%             median1(1)=line([data_vec_median(fileind,1) data_vec_median(fileind,1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
%             median1(2)=line([data_vec_median(fileind,2) data_vec_median(fileind,2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
%             hold off
%             pause
%% plotting mean subtracted-Vm histogram for each cell:
%             Fig2=figure;
%             hold on
% %             hist(data_vec_residual(:,1),nbin)
% %             hist(data_vec_residual(:,2),nbin)
%             [ncounts_residual(1,:),nbins_residual(1,:)]=hist(data_vec_residual(:,1),nbin);
%             ncounts_residual(1,:)=ncounts_residual(1,:)./length(data_vec_residual(:,1));
%             [ncounts_residual(2,:),nbins_residual(2,:)]=hist(data_vec_residual(:,2),nbin);
%             ncounts_residual(2,:)=ncounts_residual(2,:)./length(data_vec_residual(:,2));
%             bar(nbins_residual(1,:),ncounts_residual(1,:));
%             bar(nbins_residual(2,:),ncounts_residual(2,:));
% %             obj = gmdistribution.fit(data_vec(:,1),2);
%             h= findobj(gca,'Type','patch');
%             set(h(1),'FaceColor',[0 0 1],'EdgeColor','w')
%             set(h(2),'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
%             ylim1=get(gca,'ylim');
%             median_residual(1)=line([data_vec_residual_median(fileind,1) data_vec_residual_median(fileind,1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
%             median_residual(2)=line([data_vec_residual_median(fileind,2) data_vec_residual_median(fileind,2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
%             hold off
% %             pause
%             
%        clear data_no_spike_no_DC data_vec data_vec_residual
       %% saving figures
%             if save_flag==1;
%                 cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'  
%                     saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram.fig']) 
%                     print(Fig1,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram'],'-dpng','-r600','-opengl') 
%                     saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_residuals.fig']) 
%                     print(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_residuals'],'-dpng','-r600','-opengl') 
%             end
%     end %temp for file loop
%              data_vec_all_median = median(data_vec_all);
%              data_vec_residual_all_median = median(data_vec_residual_all);
%% plotting a Vm histogram for all the cells             
%         Fig_all=figure;
%         hold on
% %             hist(data_vec_all(:,1),nbin) %problem because there are different number of traces from each cell - different weight for each cell.
% %             hist(data_vec_all(:,2),nbin)
%             [ncounts_all(1,:),nbins_all(1,:)]=hist(data_vec_all(:,1),nbin);
%             ncounts_all(1,:)=ncounts_all(1,:)./length(data_vec_all(:,1));
%             [ncounts_all(2,:),nbins_all(2,:)]=hist(data_vec_all(:,2),nbin);
%             ncounts_all(2,:)=ncounts_all(2,:)./length(data_vec_all(:,2));
%             bar(nbins_all(1,:),ncounts_all(1,:));
%             bar(nbins_all(2,:),ncounts_all(2,:));
%             h= findobj(gca,'Type','patch');
%             set(h(1),'FaceColor',[0 0 1],'EdgeColor','w')
%             set(h(2),'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
%             ylim1=get(gca,'ylim');
%             median1=line([data_vec_all_median(1) data_vec_all_median(1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
%             median2=line([data_vec_all_median(2) data_vec_all_median(2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
%          hold off
%          
%   % plotting a standardized Vm histogram for all the cells             
%         Fig_all_residual=figure;
%         hold on
% %             hist(data_vec_residual_all(:,1),nbin) %problem because there are different number of traces from each cell - different weight for each cell.
% %             hist(data_vec_residual_all(:,2),nbin)
%             [ncounts_residual_all(1,:),nbins_residual_all(1,:)]=hist(data_vec_residual_all(:,1),nbin);
%             ncounts_residual_all(1,:)=ncounts_residual_all(1,:)./length(data_vec_residual_all(:,1));
%             [ncounts_residual_all(2,:),nbins_residual_all(2,:)]=hist(data_vec_residual_all(:,2),nbin);
%             ncounts_residual_all(2,:)=ncounts_residual_all(2,:)./length(data_vec_residual_all(:,2));
%             bar(nbins_residual_all(1,:),ncounts_residual_all(1,:));
%             bar(nbins_residual_all(2,:),ncounts_residual_all(2,:));
%             h= findobj(gca,'Type','patch');
%             set(h(1),'FaceColor',[0 0 1],'EdgeColor','w')
%             set(h(2),'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
%             ylim1=get(gca,'ylim');
%             median1=line([data_vec_residual_all_median(1) data_vec_residual_all_median(1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
%             median2=line([data_vec_residual_all_median(2) data_vec_residual_all_median(2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
%          hold off
%          
%           if save_flag==1;
%                 cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'  
%                     saveas(Fig_all,['Vm_histogram_all.fig']) 
%                     print(Fig_all,['Vm_histogram_all'],'-dpng','-r600','-opengl') 
%                     saveas(Fig_all_residual,['Vm_histogram_residuals_all.fig']) 
%                     print(Fig_all_residual,['Vm_histogram_residuals_all'],'-dpng','-r600','-opengl') 
%           end
%%   statistics       
%           %stats for Vm median when x_value=1
%         ongoing_stat.Vm_median=data_vec_median;
%         norm_Vm_median_noES= data_vec_median(:,1)./ data_vec_median(:,1);
%         norm_Vm_median_ES= data_vec_median(:,2)./ data_vec_median(:,1);
%         ongoing_stat.norm_Vm_median=[norm_Vm_median_noES, norm_Vm_median_ES];
%         ongoing_stat.Vm_median_m=mean(ongoing_stat.Vm_median,1);
%         ongoing_stat.Vm_median_std=std(ongoing_stat.Vm_median,0,1);
%         %testing for normal distribution
%         diff_val= data_vec_median(:,2)- data_vec_median(:,1);
%         [ongoing_stat.lillietest_h_Vm_median, ongoing_stat.lillietest_p_Vm_median] = lillietest(diff_val);
%         %paired ttest 
%         [ongoing_stat.ttest_h_Vm_median, ongoing_stat.ttest_p_Vm_median]= ttest(ongoing_stat.Vm_median(:,1),ongoing_stat.Vm_median(:,2));
%         [ongoing_stat.ttest_h_norm_Vm_median, ongoing_stat.ttest_p_norm_Vm_median]= ttest(ongoing_stat.norm_Vm_median(:,1),ongoing_stat.norm_Vm_median(:,2));
%         [ongoing_stat.wilcoxon_p_Vm_median, ongoing_stat.wilcoxon_h_Vm_median]= signrank(ongoing_stat.Vm_median(:,1),ongoing_stat.Vm_median(:,2));
%         clear norm_Vm_median_noES norm_Vm_median_ES diff_val
%         
%            %stats for Vm 5 percentile when x_value=1
%         ongoing_stat.Vm_5prctile=[data_vec_5prctile(:,:)];
%         norm_Vm_5prctile_noES= data_vec_5prctile(:,1)./ data_vec_5prctile(:,1);
%         norm_Vm_5prctile_ES= data_vec_5prctile(:,2)./ data_vec_5prctile(:,1);
%         ongoing_stat.norm_Vm_5prctile=[norm_Vm_5prctile_noES, norm_Vm_5prctile_ES];
%         ongoing_stat.Vm_5prctile_m=mean(ongoing_stat.Vm_5prctile,1);
%         ongoing_stat.Vm_5prctile_std=std(ongoing_stat.Vm_5prctile,0,1);
%         %testing for normal distribution
%         diff_val= data_vec_5prctile(:,2)- data_vec_5prctile(:,1);
%         [ongoing_stat.lillietest_h_Vm_5prctile, ongoing_stat.lillietest_p_Vm_5prctile] = lillietest(diff_val);
%         %paired ttest 
%         [ongoing_stat.ttest_h_Vm_5prctile, ongoing_stat.ttest_p_Vm_5prctile]= ttest(ongoing_stat.Vm_5prctile(:,1),ongoing_stat.Vm_5prctile(:,2));
%         [ongoing_stat.ttest_h_norm_Vm_5prctile, ongoing_stat.ttest_p_norm_Vm_5prctile]= ttest(ongoing_stat.norm_Vm_5prctile(:,1),ongoing_stat.norm_Vm_5prctile(:,2));
%         [ongoing_stat.wilcoxon_p_Vm_5prctile, ongoing_stat.wilcoxon_h_Vm_5prctile]= signrank(ongoing_stat.Vm_5prctile(:,1),ongoing_stat.Vm_5prctile(:,2));
%         clear norm_Vm_5prctile_noES norm_Vm_5prctile_ES diff_val
%         
%         %stats for Vm mean when x_value=1
%         ongoing_stat.Vm_mean=[data_vec_mean(:,:)];
%         norm_Vm_mean_noES= data_vec_mean(:,1)./ data_vec_mean(:,1);
%         norm_Vm_mean_ES= data_vec_mean(:,2)./ data_vec_mean(:,1);
%         ongoing_stat.norm_Vm_mean=[norm_Vm_mean_noES, norm_Vm_mean_ES];
%         ongoing_stat.Vm_mean_m=mean(ongoing_stat.Vm_mean,1);
%         ongoing_stat.Vm_mean_std=std(ongoing_stat.Vm_mean,0,1);
%         %testing for normal distribution
%         diff_val= data_vec_mean(:,2)- data_vec_mean(:,1);
%         [ongoing_stat.lillietest_h_Vm_mean, ongoing_stat.lillietest_p_Vm_mean] = lillietest(diff_val);
%         %paired ttest 
%         [ongoing_stat.ttest_h_Vm_mean, ongoing_stat.ttest_p_Vm_mean]= ttest(ongoing_stat.Vm_mean(:,1),ongoing_stat.Vm_mean(:,2));
%         [ongoing_stat.ttest_h_norm_Vm_mean, ongoing_stat.ttest_p_norm_Vm_mean]= ttest(ongoing_stat.norm_Vm_mean(:,1),ongoing_stat.norm_Vm_mean(:,2));
%         [ongoing_stat.wilcoxon_p_Vm_mean, ongoing_stat.wilcoxon_h_Vm_mean]= signrank(ongoing_stat.Vm_mean(:,1),ongoing_stat.Vm_mean(:,2));
%         clear norm_Vm_mean_noES norm_Vm_mean_ES diff_val
%         
%          %stats for Vm std when x_value=1
%         ongoing_stat.Vm_std=[data_vec_std(:,:)];
%         norm_Vm_std_noES= data_vec_std(:,1)./ data_vec_std(:,1);
%         norm_Vm_std_ES= data_vec_std(:,2)./ data_vec_std(:,1);
%         ongoing_stat.norm_Vm_std=[norm_Vm_std_noES, norm_Vm_std_ES];
%         ongoing_stat.Vm_std_m=std(ongoing_stat.Vm_std,1);
%         ongoing_stat.Vm_std_std=std(ongoing_stat.Vm_std,0,1);
%         %testing for normal distribution
%         diff_val= data_vec_std(:,2)- data_vec_std(:,1);
%         [ongoing_stat.lillietest_h_Vm_std, ongoing_stat.lillietest_p_Vm_std] = lillietest(diff_val);
%         %paired ttest 
%         [ongoing_stat.ttest_h_Vm_std, ongoing_stat.ttest_p_Vm_std]= ttest(ongoing_stat.Vm_std(:,1),ongoing_stat.Vm_std(:,2));
%         [ongoing_stat.ttest_h_norm_Vm_std, ongoing_stat.ttest_p_norm_Vm_std]= ttest(ongoing_stat.norm_Vm_std(:,1),ongoing_stat.norm_Vm_std(:,2));
%         [ongoing_stat.wilcoxon_p_Vm_std, ongoing_stat.wilcoxon_h_Vm_std]= signrank(ongoing_stat.Vm_std(:,1),ongoing_stat.Vm_std(:,2));
%         clear norm_Vm_std_noES norm_Vm_std_ES diff_val
%         
%          %stats for Vm CV when x_value=1
%         ongoing_stat.Vm_CV=data_vec_CV;
%         norm_Vm_CV_noES= data_vec_CV(:,1)./ data_vec_CV(:,1);
%         norm_Vm_CV_ES= data_vec_CV(:,2)./ data_vec_CV(:,1);
%         ongoing_stat.norm_Vm_CV=[norm_Vm_CV_noES, norm_Vm_CV_ES];
%         ongoing_stat.Vm_CV_m=mean(ongoing_stat.Vm_CV,1);
%         ongoing_stat.Vm_CV_std=std(ongoing_stat.Vm_CV,0,1);
%         %testing for normal distribution
%         diff_val= data_vec_CV(:,2)- data_vec_CV(:,1);
%         [ongoing_stat.lillietest_h_Vm_CV, ongoing_stat.lillietest_p_Vm_CV] = lillietest(diff_val);
%         %paired ttest 
%         [ongoing_stat.ttest_h_Vm_CV, ongoing_stat.ttest_p_Vm_CV]= ttest(ongoing_stat.Vm_CV(:,1),ongoing_stat.Vm_CV(:,2));
%         [ongoing_stat.ttest_h_norm_Vm_CV, ongoing_stat.ttest_p_norm_Vm_CV]= ttest(ongoing_stat.norm_Vm_CV(:,1),ongoing_stat.norm_Vm_CV(:,2));
%         [ongoing_stat.wilcoxon_p_Vm_CV, ongoing_stat.wilcoxon_h_Vm_CV]= signrank(ongoing_stat.Vm_CV(:,1),ongoing_stat.Vm_CV(:,2));
%         clear norm_Vm_CV_noES norm_Vm_CV_ES diff_val
%  
%%   Plotting Vm median     
% Vm_median_Y=ongoing_stat.Vm_median';
% Vm_median_X(1,:)=ones(1,size(Vm_median_Y,2));
% Vm_median_X(2,:)=2*ones(1,size(Vm_median_Y,2));
% E = ongoing_stat.Vm_median_std;
% 
% g1=figure;
% hold on
% line(Vm_median_X,Vm_median_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(Vm_median_X(:,1), mean(Vm_median_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
% %         y1limits = [0 1.1];
% %         y1ticks = [0,0.5,1];
%         set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Median Vm [mV]', 'FontSize', 28,'fontname', 'arial');
%         title(['Ongoing activity n=' num2str(length(files_to_analyze)) ', p=' num2str(ongoing_stat.wilcoxon_p_Vm_median)] ,'FontSize', 28,'fontname', 'arial');
%% plotting Vm lower 5 percentile
% Vm_5prctile_Y=ongoing_stat.Vm_5prctile';
% Vm_5prctile_X(1,:)=ones(1,size(Vm_5prctile_Y,2));
% Vm_5prctile_X(2,:)=2*ones(1,size(Vm_5prctile_Y,2));
% E = ongoing_stat.Vm_5prctile_std;
% 
% g2=figure;
% hold on
% line(Vm_5prctile_X,Vm_5prctile_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(Vm_5prctile_X(:,1), mean(Vm_5prctile_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
% 
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
% %         y1limits = [0 1.1];
% %         y1ticks = [0,0.5,1];
%         set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('lower 5 percentile [mV]', 'FontSize', 28,'fontname', 'arial');
%         title(['Ongoing activity n=' num2str(length(files_to_analyze)) ', p=' num2str(ongoing_stat.ttest_p_Vm_5prctile)] ,'FontSize', 28,'fontname', 'arial');

        %% plotting Vm mean
% Vm_mean_Y=ongoing_stat.Vm_mean';
% Vm_mean_X(1,:)=ones(1,size(Vm_mean_Y,2));
% Vm_mean_X(2,:)=2*ones(1,size(Vm_mean_Y,2));
% E = ongoing_stat.Vm_mean_std;
% 
% g3=figure;
% hold on
% line(Vm_mean_X,Vm_mean_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(Vm_mean_X(:,1), mean(Vm_mean_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
% 
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
% %         y1limits = [0 1.1];
% %         y1ticks = [0,0.5,1];
%         set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Vm mean [mV]', 'FontSize', 28,'fontname', 'arial');
%         title(['Ongoing activity n=' num2str(length(files_to_analyze)) ', p=' num2str(ongoing_stat.wilcoxon_p_Vm_mean)] ,'FontSize', 28,'fontname', 'arial');

                %% plotting Vm std
% Vm_std_Y=ongoing_stat.Vm_std';
% Vm_std_X(1,:)=ones(1,size(Vm_std_Y,2));
% Vm_std_X(2,:)=2*ones(1,size(Vm_std_Y,2));
% E = ongoing_stat.Vm_std_std;
% 
% g4=figure;
% hold on
% line(Vm_std_X,Vm_std_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(Vm_std_X(:,1), mean(Vm_std_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
% 
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
% %         y1limits = [0 1.1];
% %         y1ticks = [0,0.5,1];
%         set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Vm std [mV]', 'FontSize', 28,'fontname', 'arial');
%         title(['Ongoing activity n=' num2str(length(files_to_analyze)) ', p=' num2str(ongoing_stat.wilcoxon_p_Vm_std)] ,'FontSize', 28,'fontname', 'arial');
%%   Plotting Vm CV     
% Vm_CV_Y=ongoing_stat.Vm_CV';
% Vm_CV_X(1,:)=ones(1,size(Vm_CV_Y,2));
% Vm_CV_X(2,:)=2*ones(1,size(Vm_CV_Y,2));
% E = ongoing_stat.Vm_CV_std;
% 
% g5=figure;
% hold on
% line(Vm_CV_X,Vm_CV_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(Vm_CV_X(:,1), mean(Vm_CV_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% hold off
%         x1limits = [0.75 2.25];
%         x1ticks = [1,2];
% %         y1limits = [0 1.1];
% %         y1ticks = [0,0.5,1];
%         set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('CV Vm', 'FontSize', 28,'fontname', 'arial');
%         title(['Ongoing activity n=' num2str(length(files_to_analyze)) ', p=' num2str(ongoing_stat.wilcoxon_p_Vm_CV)] ,'FontSize', 28,'fontname', 'arial');

%%        %save figures  
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz'; 
% fn1='Ongoing_Vm_median';
% saveas(g1,'Ongoing_Vm_median.fig'); 
% print(g1,fn1,'-dpng','-r600','-opengl') 
% fn2='Ongoing_Vm_5prcentile';
% saveas(g2,'Ongoing_Vm_5prcentile.fig'); 
% print(g2,fn2,'-dpng','-r600','-opengl') 
% fn3='Ongoing_Vm_mean';
% saveas(g3,'Ongoing_Vm_mean.fig'); 
% print(g3,fn3,'-dpng','-r600','-opengl') 
% fn4='Ongoing_Vm_std';
% saveas(g4,'Ongoing_Vm_std.fig'); 
% print(g4,fn4,'-dpng','-r600','-opengl') 
% fn5='Ongoing_Vm_CV';
% saveas(g5,'Ongoing_Vm_CV.fig'); 
% print(g5,fn5,'-dpng','-r600','-opengl') 
 %% Plot power spectrum
start_time = [0.4,6]; %[sec]
duration = 2.5; %[sec] 
x_value = 1;
Y_abs = []; f = [];  interval = []; 

for t=1:length(start_time);
DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = [];

start_sample = start_time(t).*sf{channel};
if start_time(t)==0
    start_sample = 1;
end
end_sample = start_sample+duration.*sf{channel}-1;
interval(:,t) = start_sample:end_sample;
spec_mat = raw_data{channel}(interval(:,t),:,x_value);
DC= mean(spec_mat,1);
l=size(spec_mat,1);
l=l/2-1;
DC=wextend('addrow', 'per', DC, l);
spec_mat_noDC = spec_mat-DC;

[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf{channel},0,0);
end
%% Plot power spectrum with shaded error bar
% Fig2=figure;
%  for t=1:length(start_time);
%       hold on
%     fn_shadedErrorBar(f(:,t),mean(Y_abs(:,:,t),2),std(Y_abs(:,:,t),0,2),{'LineWidth',1,'color', color_table(t,:)},1);
%         xlim([0 1000]); ylim([-0.01 10])
%         title('')
%         set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2);
%         set(gca,'xscale','log');
%         set(gca,'yscale','log');
%         xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
%         ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
%        hold off
%  end
% pause
%% Plot power spectrum without shaded error bar
figure
        xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
        ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
        set(gca,'xscale','log');
        set(gca,'yscale','log');
 for t=1:length(start_time);
      hold on
        plot(f(:,t),mean(Y_abs(:,:,t),2),'color', color_table(t,:))        
        xlim([0 1000]); ylim([0.001 10])
%         title('')
        set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2)        
      hold off
 end
 pause
     end %temp file loop

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
% filename='ongoing'; 
% save(filename, 'files_to_analyze', 'ongoing', 'ongoing_stat')
