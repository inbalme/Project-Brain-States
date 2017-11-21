%% Analyze Vm histograms and powerspectrum
% This file was created on 15/1/2017 based on "Analyze_NBES_Vm_hist_and_powerspec_10_20Hz_v2'
%This file is used for the analysis of files created with extract_NBES_Data_v3 or Extract_ChAT_Data_v3
%when running: pay attention to change the exp_type and files_to_analyze.
%To analyze LFP, change current_data=Ch2_data. for powerspectrum this is in line 400.

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
exp_type=4; %1-NBES, 2-ChAT
data_type='Vm'; %'LFP', 'Vm'
trace_type_input=1; %
analyze_time_before_train=0;
analyze_train_only_flag=0;
analyze_hist_flag=1;
analyze_powerspec_flag=0;
save_flag= 1;
print_flag=0;
norm_flag=0;
no_DC_flag=1; %put 1 for exp_type=4
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
        files_to_analyze =[8,10,12,14,15,16,22,37,40,1,46,48,52,58,72,82,84]; %[46,52,58,62,72,75,84]; %[8,10,12,14,15,16,22,36,37,40,1,44,46,48,52,56,58,62,72,75,82,84]; %[46,52,58,62,72,75,84]; %one cell from each animal for LFP  
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB-', 'NB+'};
        path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm Histograms and Powerspec';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end        

    case 2
        files_to_analyze = 80; %[76,80,84,92,114]; %[74,76,77,80,82,84,87,90,92,112,114,115]; %[74,76,77,80,82,84,87]; %[76,80,84] one cell from each animal for LFP
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light Off', 'Light On'};
        path_output= 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm Histograms and Powerspec';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end   
     case 4
        files_to_analyze =  [118,120,121,122,126,127,128,129]; %119
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light Off', 'Light On'};
        path_output= 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm Histograms and Powerspec_awake';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end        
end

if analyze_hist_flag==1;
    for fileind=1:length(files_to_analyze) ;
    close all
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
   %%
if exp_type~=4
   Ch2_data= raw_data{3}./20; %dividing by the LFP gain        
end
    current_data=data_no_spikes{channel};   %Ch2_data; %data_no_spikes{channel};    
    galvano_nstim = Param.facade(6);
    galvano_freq = Param.facade(7);

    data_preprocessing       
   
 if isempty(current_data_filt)
     current_data_filt=current_data;
 end
 
%%
clear color_table
        whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;
        switch exp_type
            case 1 
                color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
            case 2
                color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [102, 172,255]./256; [30,75,14]/256; [112,172,90]/256];  
            case 4
                color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [102, 172,255]./256; [30,75,14]/256; [112,172,90]/256];  
        end
 
%% Vm histogram - for spontaneous and evoked activity
 for trace_type=  trace_type_input   %1 or 3 for spont., 2 for evoked
     interval=[]; data_vec=[]; data_vec_residual=[]; data_vec_5traces=[]; data_vec_5traces_residual=[]; data_mat=[]; data_mat_median=[];
 intervals_to_analyze         
        %%
        binsize=1; %[mV];
        for t=1:size(interval,2);
              data_vec(:,t)=reshape(current_data_filt(interval(:,t),:,x_value(t)),numel(current_data_filt(interval(:,t),:,x_value(t))),1); %take only data from x_value=1, and connect all traces into one long vector;
              data_mat{t}(:,:)=current_data_filt(interval(:,t),:,x_value(t));
              data_mat_median{t}(:,:)=median(data_mat{t}(:,:),1);
              data_vec_residual(:,t)=data_vec(:,t)-mean(data_vec(:,t));
              data_vec_5traces(:,t)=reshape(current_data_filt(interval(:,t),1:5,x_value(t)),numel(current_data_filt(interval(:,t),1:5,x_value(t))),1); %take only data from x_value=1, and connect all traces into one long vector;
              data_vec_5traces_residual(:,t)=data_vec_5traces(:,t)-mean(data_vec_5traces(:,1));
              nbin(t)=ceil(range(data_vec(:,t)))./binsize; %set the range according to the interval before ES
            end
%              nbin = 30;
             data_mat_median_m(fileind,:)=[mean(data_mat_median{1}(:,:)),mean(data_mat_median{2}(:,:))];
             data_mat_median_std(fileind,:)=[std(data_mat_median{1}(:,:),0,2),std(data_mat_median{2}(:,:),0,2)];
             data_vec_residual_median{trace_type}(fileind,:) = median(data_vec_residual);
             data_vec_5prctile{trace_type}(fileind,:)  = prctile(data_vec,5,1);
%              data_vec_mean{trace_type}(fileind,:) = mean(data_vec,1);
%              data_vec_std{trace_type}(fileind,:) = std(data_vec,1);
%              data_vec_CV{trace_type}(fileind,:) =data_vec_std{trace_type}(fileind,:)./ abs(data_vec_mean{trace_type}(fileind,:));
            data_vec_all{trace_type}=[data_vec_all{trace_type}; data_vec_5traces];
            data_vec_residual_all{trace_type}=[data_vec_residual_all{trace_type}; data_vec_5traces_residual];
  switch trace_type
      case 1
        ongoing(fileind).Vm_median=[data_mat_median{1}(:,:)',data_mat_median{2}(:,:)'];  
        ongoing(fileind).Vm_median_m=mean(ongoing(fileind).Vm_median);
        ongoing(fileind).Vm_median_std=std(ongoing(fileind).Vm_median,0,1);
        
        %testing for normal distribution
        diff_val= ongoing(fileind).Vm_median(:,2)- ongoing(fileind).Vm_median(:,1);
        [ongoing(fileind).lillietest_h_Vm_median,ongoing(fileind).lillietest_p_Vm_median] = lillietest(diff_val);
        %paired ttest 
        [ongoing(fileind).ttest_h_Vm_median,ongoing(fileind).ttest_p_Vm_median]= ttest(ongoing(fileind).Vm_median(:,1),ongoing(fileind).Vm_median(:,2));
        [ongoing(fileind).wilcoxon_p_Vm_median,ongoing(fileind).wilcoxon_h_Vm_median]= signrank(ongoing(fileind).Vm_median(:,1),ongoing(fileind).Vm_median(:,2));
        clear diff_val
        
        ongoing(fileind).Vm_5prctile=data_vec_5prctile{trace_type}(fileind,:);   

        hists.ongoing=ongoing;
      case 2
        evoked(fileind).Vm_median=data_mat_median_m(fileind,:);  
        evoked(fileind).Vm_5prctile=data_vec_5prctile{trace_type}(fileind,:);   

        hists.evoked=evoked;
  end
 %% plotting Vm histogram for each cell:
 if print_flag==1;
     ncounts=[]; nbins=[];
            Fig1=figure;
            hold on
            [ncounts{1}(1,:),nbins{1}(1,:)]=hist(data_vec(:,1),nbin(1));
            ncounts{1}(1,:)=ncounts{1}(1,:)./length(data_vec(:,1));
            [ncounts{2}(1,:),nbins{2}(1,:)]=hist(data_vec(:,2),nbins{1}(1,:)); %nbin(2)
            ncounts{2}(1,:)=ncounts{2}(1,:)./length(data_vec(:,2));
            h1=bar(nbins{1}(1,:),ncounts{1}(1,:));
            h2=bar(nbins{2}(1,:),ncounts{2}(1,:));
            
%             obj = gmdistribution.fit(data_vec(:,1),2);
            set(h2,'FaceColor',color_table(2,:),'EdgeColor','w','faceAlpha', 0.3)
            set(h1,'FaceColor',color_table(1,:),'EdgeColor','w','faceAlpha', 0.7)
            ylim1=get(gca,'ylim');
%             median1(1)=line([data_mat_median_m(fileind,1) data_mat_median_m(fileind,1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',color_table(1,:),'linewidth',1.5);
%             median1(2)=line([data_mat_median_m(fileind,2) data_mat_median_m(fileind,2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',color_table(2,:),'linewidth',1.5);
            hold off
            axis tight
            xlim1=get(gca,'xlim');
            xrim=diff(xlim).*0.1;
            xlim1(1)=xlim1(1)-xrim;
            xlim1(2)=xlim1(2)+xrim;
            set(gca,'xlim',xlim1,'FontSize', 20,'fontname', 'arial')
            xlabel('Vm [mV]', 'FontSize', 20,'fontname', 'arial');
            ylabel('Probability', 'FontSize', 20,'fontname', 'arial')
            [l,OBJH,OUTH,OUTM] = legend([h1 h2],legend_string, 'position',[0.8,0.9,0.2,0.1],'fontsize',10, 'box', 'off'); % returns a handle LEGH to the legend axes; a vector OBJH containing handles for the text, lines, and patches in the legend; a vector OUTH of handles to thelines and patches in the plot; and a cell array OUTM containingthe text in the legend.
lpatch=findobj(OBJH,'type','patch');
lpatch(1).Vertices(2:3,2)=lpatch(1).Vertices(1,2)+0.15;
lpatch(2).Vertices(2:3,2)=lpatch(2).Vertices(1,2)+0.15;
lpatch(1).Vertices(:,2)=lpatch(1).Vertices(:,2)+0.06;
lpatch(2).Vertices(:,2)=lpatch(2).Vertices(:,2)+0.08;
lpatch(1).FaceColor=color_table(1,:);
lpatch(1).FaceAlpha=0.7;
lpatch(2).FaceColor=color_table(2,:);
lpatch(2).FaceAlpha=0.3;
            %% plot the medians+ci on a separate figure
%             switch exp_type
%         case 1
%             fileind=12;
%             legend_string={'NB-','NB+'};
%             cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous'
%             load('Spontaneous activity.mat')
%             color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
%         case 2
%             fileind=4;
%             legend_string={'Light Off','Light On'};
%             cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous'
%             load('Spontaneous activity.mat')
%             color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
%     end

%% bootstrap CI for the median:  
%         for t=1:2;                 
%                     dataset=[];
% %                     dataset=data_vec(:,t)';
%                     dataset=data_mat{t}(:,:);
%                     k= size(dataset,2);     
%                     lower_bound=2.5;
%                     upper_bound=97.5;
%                     iterations=5000;
% %                      [prcntile1(t,:), prcntile2(t,:)]=fn_get_CI_w_bootstrap(dataset,length(dataset),iterations,2.5,97.5,'median');
%                      data_vec_med(1,t) = median(data_vec(:,t)); 
%         end
%             Fig2=figure;
%         hold on
%             e1=errorbar(1,data_vec_med(1,1),data_vec_med(1,1)-prcntile1(1,1),prcntile2(1,1)-data_vec_med(1,1),'color',color_table(3,:)); %
%             e2=errorbar(2,data_vec_med(1,2),data_vec_med(1,2)-prcntile1(2,1),prcntile2(2,1)-data_vec_med(1,2),'color',color_table(4,:)); % 
%             set(e1,'Marker','o','MarkerSize',5,'MarkerFaceColor',color_table(1,:),'MarkerEdgeColor',color_table(1,:),'linewidth',2)
%             set(e2,'Marker','o','MarkerSize',5,'MarkerFaceColor',color_table(2,:),'MarkerEdgeColor',color_table(2,:),'linewidth',2)
%             set(gca,'xtick',[1,2],'xticklabel',legend_string);
%             ylabel('Vm [mV]')
% 
%         hold off
%             
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
%             h1=bar(nbins_residual(1,:),ncounts_residual(1,:));
%             h2=bar(nbins_residual(2,:),ncounts_residual(2,:));
% %             obj = gmdistribution.fit(data_vec(:,1),2);
%             set(h2,'FaceColor',[0 0 1],'EdgeColor','w')
%             set(h1,'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
%             ylim1=get(gca,'ylim');
%             median_residual(1)=line([data_vec_residual_median{trace_type}(fileind,1) data_vec_residual_median{trace_type}(fileind,1)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0],'linewidth',1);
%             median_residual(2)=line([data_vec_residual_median{trace_type}(fileind,2) data_vec_residual_median{trace_type}(fileind,2)],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 1],'linewidth',1);
%             hold off
% %             pause
            
       clear data_no_spike_no_DC data_vec data_vec_residual
       %% saving figures
            if save_flag==1;
  cd(path_output)
                    saveas(Fig1,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_', type_strng{trace_type},'.fig']) 
                    print(Fig1,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_', type_strng{trace_type}],'-dpng','-r600','-opengl') 
%                      saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_median_ci_', type_strng{trace_type},'.fig']) 
%                     print(Fig2,['f' num2str(files_to_analyze(fileind)) '_median_ci_', type_strng{trace_type}],'-dpng','-r600','-opengl') 
%                     saveas(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_residuals_',type_strng{trace_type},'.fig']) 
%                     print(Fig2,['f' num2str(files_to_analyze(fileind)) '_Vm_histogram_residuals_', type_strng{trace_type}],'-dpng','-r600','-opengl') 
            end
 end
    end 
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

       hists_stat(trace_type).Vm_median=data_mat_median_m;
        norm_Vm_median_noES= data_mat_median_m(:,1)./ data_mat_median_m(:,1);
        norm_Vm_median_ES= data_mat_median_m(:,2)./ data_mat_median_m(:,1);
       hists_stat(trace_type).norm_Vm_median=[norm_Vm_median_noES, norm_Vm_median_ES];
       hists_stat(trace_type).Vm_median_m=mean(hists_stat(trace_type).Vm_median,1);
       hists_stat(trace_type).Vm_median_std=std(hists_stat(trace_type).Vm_median,0,1);
        %testing for normal distribution
        diff_val= data_mat_median_m(:,2)- data_mat_median_m(:,1);
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
        ylabel('Median Vm (mV)', 'FontSize', 28,'fontname', 'arial');
        title([type_strng{trace_type},' activity n=' num2str(length(files_to_analyze)) ', p=' num2str(hists_stat(trace_type).wilcoxon_p_Vm_median)] ,'FontSize', 28,'fontname', 'arial');
  
        %% plotting Vm lower 5 percentile

        Vm_5prctile_Y=hists_stat(trace_type).Vm_5prctile';
Vm_5prctile_X(1,:)=ones(1,size(Vm_5prctile_Y,2));
Vm_5prctile_X(2,:)=2*ones(1,size(Vm_5prctile_Y,2));
E =hists_stat(trace_type).Vm_5prctile_std;
% if hists_stat(trace_type).lillietest_h_Vm_5prctile ==0 %data is normally dist. and we can use ttest
%     p_5prctile= hists_stat(trace_type).ttest_p_Vm_5prctile;
% else %data is not dist. noarmally and we need non-parametric test
%      p_5prctile=hists_stat(trace_type).wilcoxon_p_Vm_5prctile;
% end
p_5prctile=hists_stat(trace_type).wilcoxon_p_Vm_5prctile;

    if p_5prctile >0.05 
    asterisk_sp='n.s.';
    a_fontsize=13;
else if p_5prctile<0.05 && p_5prctile>0.01
    asterisk_sp='*';
    a_fontsize=17;
    else if p_5prctile<0.01 && p_5prctile>0.001
            asterisk_sp='**';
            a_fontsize=17;
    else if p_5prctile<0.001
             asterisk_sp='***';
             a_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(hists_stat(trace_type).Vm_5prctile))+0.1*abs(max(max(hists_stat(trace_type).Vm_5prctile))); 
liney=[-32 ;-32 ];

g2=figure;
hold on
line(Vm_5prctile_X,Vm_5prctile_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Vm_5prctile_X(:,1), mean(Vm_5prctile_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','middle','fontsize',a_fontsize)
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = get(gca,'ylim');
        y1limits(2) = my;
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks, 'ylim', y1limits, 'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 'ylim', [-85,-25],
        ylabel('5%-ile (mV)', 'FontSize', 28,'fontname', 'arial');
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

 cd(path_output)
if save_flag==1;
   % filename='hists'; 
filename='hists'; 
save(filename, 'files_to_analyze', 'hists', 'hists_stat')
end  
end
 %% Plot power spectrum
if analyze_powerspec_flag==1
     for fileind=1:length(files_to_analyze) ;
    close all
    channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
   %%
   Ch2_data= raw_data{3}./20; %dividing by the LFP gain      
if strcmp(data_type,'LFP')
    current_data=Ch2_data;
else
    current_data= data_no_spikes{channel};  % Ch2_data;  data_no_spikes{channel};
end
    galvano_nstim = Param.facade(6);
    galvano_freq = Param.facade(7);

data_preprocessing   

 if isempty(current_data_filt)
     current_data_filt=current_data;
 end
 
 clear color_table
        whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;
        switch exp_type
            case 1 
                color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
            case 2
                color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
        end
 trace_type = trace_type_input;   %1 or 3 for spont., 2 for evoked
     interval=[]; data_vec=[]; data_vec_residual=[]; data_vec_5traces=[]; data_vec_5traces_residual=[]; 
 intervals_to_analyze 
        %%
Y_abs = []; f = [];  Y_abs_sum=[]; Y_abs_norm=[]; Y_abs_mean=[]; Y_abs_log=[];

for t=1:length(start_time);
DC=[]; spec_mat_noDC=[]; spec_mat = []; 
spec_mat = current_data_filt(interval(:,t),:,x_value(t));
DC= mean(spec_mat,1);
spec_mat_noDC=bsxfun(@minus,spec_mat,DC);
% l=size(spec_mat,1);
% l=l/2-1;
% DC=wextend('addrow', 'per', DC, l);
% spec_mat_noDC = spec_mat-DC;

[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf{channel},0,0);
Y_abs_sem(:,t)=std(Y_abs(:,:,t),0,2)./sqrt(size(Y_abs(:,:,t),2));
Y_abs_sum(t,:)=sum(Y_abs(:,:,t),1);
Y_abs_norm(:,:,t)=bsxfun(@rdivide,Y_abs(:,:,t),Y_abs_sum(t,:));
Y_abs_log(:,:,t)=log(Y_abs(:,:,t));

powerspec.fileind(fileind).fname=files_to_analyze(fileind);
powerspec.fileind(fileind).Y_abs=Y_abs;
powerspec.fileind(fileind).Y_abs_sum=Y_abs_sum;
powerspec.fileind(fileind).Y_abs_norm=Y_abs_norm;
powerspec.fileind(fileind).Y_abs_log=Y_abs_norm;
end
fround=round(f(:,1));
delta=[1,4];
theta=[4,8];
alpha=[8,12];
beta=[12,25];
gamma=[30,50];
lowfreq=[1,10];
a=fround>=delta(1) & fround<=delta(2);
b=fround>=theta(1) & fround<=theta(2);
c=fround>=alpha(1) & fround<=alpha(2);
d=fround>=beta(1) & fround<=beta(2);
e=fround>=gamma(1) & fround<=gamma(2);
e2=fround>=lowfreq(1) & fround<=lowfreq(2);
% Y_abs_mean=mean(Y_abs_norm,2);
% Y_abs_mean=mean(Y_abs_log,2);
Y_abs_mean=mean(Y_abs,2);
Y_abs_delta(1:2)=sum(Y_abs_mean(a,1,:));
Y_abs_theta(1:2)=sum(Y_abs_mean(b,1,:));
Y_abs_alpha(1:2)=sum(Y_abs_mean(c,1,:));
Y_abs_beta(1:2)=sum(Y_abs_mean(d,1,:));
Y_abs_gamma(1:2)=sum(Y_abs_mean(e,1,:));
Y_abs_lowfreq(1:2)=sum(Y_abs_mean(e2,1,:));
Y_delta_all(fileind,:)=Y_abs_delta;
Y_theta_all(fileind,:)=Y_abs_theta;
Y_alpha_all(fileind,:)=Y_abs_alpha;
Y_beta_all(fileind,:)=Y_abs_beta;
Y_gamma_all(fileind,:)=Y_abs_gamma;
Y_lowfreq_all(fileind,:)=Y_abs_lowfreq;
Y_diff(fileind,:)=[diff(Y_abs_lowfreq), diff(Y_abs_gamma)];
Y_sum(fileind,:)=mean(Y_abs_sum,2);
     end %end of files loop
powerspec.delta_range=delta;
powerspec.theta_range=theta;
powerspec.alpha_range=alpha;
powerspec.beta_range=beta;
powerspec.gamma_range=gamma;
powerspec.lowfreq_range=lowfreq;
powerspec.delta=Y_delta_all;
powerspec.theta=Y_theta_all;
powerspec.alpha=Y_alpha_all;
powerspec.beta=Y_beta_all;
powerspec.gamma=Y_gamma_all;
powerspec.lowfreq=Y_lowfreq_all;
powerspec.diff=Y_diff;
powerspec.Y_sum=Y_sum;
powerspec.BP50HzLFP_flag=1; %removing 50Hz noise from LFP signal
powerspec.BP50HzVm_flag=1; %removing 50Hz noise from Vm signal
powerspec.BPLFP_flag=0; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
powerspec.BPVm_flag=0; 

[powerspec.delta_ttest_h,powerspec.delta_ttest_p]=ttest(powerspec.delta(:,1),powerspec.delta(:,2));
[powerspec.theta_ttest_h,powerspec.theta_ttest_p]=ttest(powerspec.theta(:,1),powerspec.theta(:,2));
[powerspec.alpha_ttest_h,powerspec.alpha_ttest_p]=ttest(powerspec.alpha(:,1),powerspec.alpha(:,2));
[powerspec.beta_ttest_h,powerspec.beta_ttest_p]=ttest(powerspec.beta(:,1),powerspec.beta(:,2));
[powerspec.gamma_ttest_h,powerspec.gamma_ttest_p]=ttest(powerspec.gamma(:,1),powerspec.gamma(:,2));
[powerspec.lowfreq_ttest_h,powerspec.lowfreq_ttest_p]=ttest(powerspec.lowfreq(:,1),powerspec.lowfreq(:,2));
[powerspec.diff_ttest_h,powerspec.diff_ttest_p]=ttest(powerspec.diff(:,1),powerspec.diff(:,2));
[powerspec.Y_sum_ttest_h,powerspec.Y_sum_ttest_p]=ttest(powerspec.Y_sum(:,1),powerspec.Y_sum(:,2));

% %normalization to pre-NB instead of to the total power
% powerspec.change_delta=[(powerspec.delta(:,2)-powerspec.delta(:,1))./abs(powerspec.delta(:,1))].*100; %percent change
% powerspec.change_theta=[(powerspec.theta(:,2)-powerspec.theta(:,1))./abs(powerspec.theta(:,1))].*100; %percent change
% powerspec.change_alpha=[(powerspec.alpha(:,2)-powerspec.alpha(:,1))./abs(powerspec.alpha(:,1))].*100; %percent change
% powerspec.change_beta=[(powerspec.beta(:,2)-powerspec.beta(:,1))./abs(powerspec.beta(:,1))].*100; %percent change
% powerspec.change_gamma=[(powerspec.gamma(:,2)-powerspec.gamma(:,1))./abs(powerspec.gamma(:,1))].*100; %percent change
% powerspec.change_lowfreq=[(powerspec.lowfreq(:,2)-powerspec.lowfreq(:,1))./abs(powerspec.lowfreq(:,1))].*100; %percent change

if length(files_to_analyze)>4;
[powerspec.delta_lillietest_h,powerspec.delta_lillietest_p]=lillietest(powerspec.delta(:,2)-powerspec.delta(:,1));
[powerspec.delta_wilcoxon_p,powerspec.delta_wilcoxon_h]=signrank(powerspec.delta(:,1),powerspec.delta(:,2));
% [powerspec.change_delta_wilcoxon_p,powerspec.change_delta_wilcoxon_h]=signrank(powerspec.change_delta);

[powerspec.theta_lillietest_h,powerspec.theta_lillietest_p]=lillietest(powerspec.theta(:,2)-powerspec.theta(:,1));
[powerspec.theta_wilcoxon_p,powerspec.theta_wilcoxon_h]=signrank(powerspec.theta(:,1),powerspec.theta(:,2));
% [powerspec.change_theta_wilcoxon_p,powerspec.change_theta_wilcoxon_h]=signrank(powerspec.change_theta);

[powerspec.alpha_lillietest_h,powerspec.alpha_lillietest_p]=lillietest(powerspec.alpha(:,2)-powerspec.alpha(:,1));
[powerspec.alpha_wilcoxon_p,powerspec.alpha_wilcoxon_h]=signrank(powerspec.alpha(:,1),powerspec.alpha(:,2));
% [powerspec.change_alpha_wilcoxon_p,powerspec.change_alpha_wilcoxon_h]=signrank(powerspec.change_alpha);

[powerspec.beta_lillietest_h,powerspec.beta_lillietest_p]=lillietest(powerspec.beta(:,2)-powerspec.beta(:,1));
[powerspec.beta_wilcoxon_p,powerspec.beta_wilcoxon_h]=signrank(powerspec.beta(:,1),powerspec.beta(:,2));
% [powerspec.change_beta_wilcoxon_p,powerspec.change_beta_wilcoxon_h]=signrank(powerspec.change_beta);

[powerspec.gamma_lillietest_h,powerspec.gamma_lillietest_p]=lillietest(powerspec.gamma(:,2)-powerspec.gamma(:,1));
[powerspec.gamma_wilcoxon_p,powerspec.gamma_wilcoxon_h]=signrank(powerspec.gamma(:,1),powerspec.gamma(:,2));
% [powerspec.change_gamma_wilcoxon_p,powerspec.change_gamma_wilcoxon_h]=signrank(powerspec.change_gamma);

[powerspec.lowfreq_lillietest_h,powerspec.lowfreq_lillietest_p]=lillietest(powerspec.lowfreq(:,2)-powerspec.lowfreq(:,1));
[powerspec.lowfreq_wilcoxon_p,powerspec.lowfreq_wilcoxon_h]=signrank(powerspec.lowfreq(:,1),powerspec.lowfreq(:,2));
% [powerspec.change_lowfreq_wilcoxon_p,powerspec.change_lowfreq_wilcoxon_h]=signrank(powerspec.change_lowfreq);

[powerspec.diff_lillietest_h,powerspec.diff_lillietest_p]=lillietest(powerspec.diff(:,2)-powerspec.diff(:,1));
[powerspec.diff_wilcoxon_p,powerspec.diff_wilcoxon_h]=signrank(powerspec.diff(:,1),powerspec.diff(:,2));

[powerspec.Y_sum_lillietest_h,powerspec.Y_sum_lillietest_p]=lillietest(powerspec.Y_sum(:,2)-powerspec.Y_sum(:,1));
[powerspec.Y_sum_wilcoxon_p,powerspec.Y_sum_wilcoxon_h]=signrank(powerspec.Y_sum(:,1),powerspec.Y_sum(:,2));
end
%% plots
%% Powerspec
j1=figure;
for t=1:length(start_time)
      hold on
   d1(t)= fn_shadedErrorBar(f(f(:,t)<49,t),mean(Y_abs(f(:,t)<49,:,t),2),Y_abs_sem(f(:,t)<49,t),{'LineWidth',2,'color', color_table(t,:)});  
end
    hold off

 xlim([1 100]); 
 if strcmp(data_type,'LFP')
     ylim([0 0.01])
     y1tick=[10e-7,10e-6, 10e-5, 10e-4, 10e-3, 10e-2];
 else
   ylim([0 1])
     y1tick=[10e-6, 10e-5, 10e-4, 10e-3, 10e-2,10e-1];  
 end
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2, 'xtick',[1, 10,100],'ytick', y1tick,'yminortick','on') 
         xlabel('Frequency (Hz)','fontsize',20, 'fontname', 'arial')
        ylabel('PSD (mV^2/Hz)','fontsize',20, 'fontname', 'arial'); %Power spectral density
l=legend([d1(1).mainLine d1(2).mainLine ],legend_string,'Location','northeast', 'box', 'off');
l.LineWidth=1.5;
l.FontSize=12;

            if save_flag==1;
                cd(path_output)
                 if strcmp(data_type,'LFP') 
                     filename=['f' num2str(files_to_analyze(fileind)) '_LFP_power_spectrum'];
                 else
                     filename=['f' num2str(files_to_analyze(fileind)) '_Vm_power_spectrum'];
                 end
                    saveas(j1,filename,'fig') 
                    print(j1,filename,'-dpng','-r600','-opengl') 
            end
%%
% Delta Power
tmp_Y= powerspec.delta';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(tmp_Y');
if isfield(powerspec,'delta_wilcoxon_p')
    p_val=powerspec.delta_wilcoxon_p;
else
    p_val=powerspec.delta_ttest_p;
end
if p_val >0.05 
    asterisk_sp='n.s.';
    a_fontsize=13;
else if p_val<0.05 && p_val>0.01
    asterisk_sp='*';
    a_fontsize=17;
    else if p_val<0.01 && p_val>0.001
            asterisk_sp='**';
            a_fontsize=17;
    else if p_val<0.001
             asterisk_sp='***';
             a_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(powerspec.delta))*1.1; 
liney=[my ;my ];

g1=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
        y1limits(2)=max(max(powerspec.delta))*1.2; %my;
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Normalized PSD [A.U.]', 'FontSize', 28,'fontname', 'arial');
        ylabel('PSD [mV^2/Hz]','fontsize',28, 'fontname', 'arial');
        title(['Delta Power,  p=' num2str(p_val)] ,'FontSize', 20,'fontname', 'arial');  
        %% Theta Power
tmp_Y= powerspec.theta';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(tmp_Y');
if isfield(powerspec,'theta_wilcoxon_p')
    p_val=powerspec.theta_wilcoxon_p;
else
    p_val=powerspec.theta_ttest_p;
end
if p_val >0.05 
    asterisk_sp='n.s.';
    a_fontsize=13;
else if p_val<0.05 && p_val>0.01
    asterisk_sp='*';
    a_fontsize=17;
    else if p_val<0.01 && p_val>0.001
            asterisk_sp='**';
            a_fontsize=17;
    else if p_val<0.001
             asterisk_sp='***';
             a_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(powerspec.theta))*1.1; 
liney=[my ;my ];

g2=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
        y1limits(2)=max(max(powerspec.theta))*1.2; %my;
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Normalized PSD [A.U.]', 'FontSize', 28,'fontname', 'arial');
        ylabel('PSD [mV^2/Hz]','fontsize',20, 'fontname', 'arial');
        title(['Theta Power,  p=' num2str(p_val)] ,'FontSize', 20,'fontname', 'arial');  
        
%% Alpha Power
tmp_Y= powerspec.alpha';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(tmp_Y');
if isfield(powerspec,'alpha_wilcoxon_p')
    p_val=powerspec.alpha_wilcoxon_p;
else
    p_val=powerspec.alpha_ttest_p;
end
if p_val >0.05 
    asterisk_sp='n.s.';
    a_fontsize=13;
else if p_val<0.05 && p_val>0.01
    asterisk_sp='*';a_fontsize=17;
    else if p_val<0.01 && p_val>0.001
            asterisk_sp='**';a_fontsize=17;
    else if p_val<0.001
             asterisk_sp='***';a_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(powerspec.alpha))*1.1; 
liney=[my ;my ];

g3=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
        y1limits(2)=max(max(powerspec.alpha))*1.2; %my;
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Normalized PSD [A.U.]', 'FontSize', 28,'fontname', 'arial');
        ylabel('PSD [mV^2/Hz]','fontsize',20, 'fontname', 'arial');
        title(['Alpha Power,  p=' num2str(p_val)] ,'FontSize', 20,'fontname', 'arial');   
        %% Beta Power
tmp_Y= powerspec.beta';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(tmp_Y');
if isfield(powerspec,'beta_wilcoxon_p')
    p_val=powerspec.beta_wilcoxon_p;
else
    p_val=powerspec.beta_ttest_p;
end
if p_val >0.05 
    asterisk_sp='n.s.';
    a_fontsize=13;
else if p_val<0.05 && p_val>0.01
    asterisk_sp='*';
    a_fontsize=17;
    else if p_val<0.01 && p_val>0.001
            asterisk_sp='**';
            a_fontsize=17;
    else if p_val<0.001
             asterisk_sp='***';
             a_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(powerspec.beta))*1.1; 
liney=[my ;my ];

g4=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
        y1limits(2)=max(max(powerspec.beta))*1.2; %my;
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Normalized PSD [A.U.]', 'FontSize', 28,'fontname', 'arial');
        ylabel('PSD [mV^2/Hz]','fontsize',20, 'fontname', 'arial');
        title(['Beta Power,  p=' num2str(p_val)] ,'FontSize', 20,'fontname', 'arial');  
   %% Gamma Power
tmp_Y= powerspec.gamma';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(tmp_Y');
if isfield(powerspec,'gamma_wilcoxon_p')
    p_val=powerspec.gamma_wilcoxon_p;
else
    p_val=powerspec.gamma_ttest_p;
end
if p_val >0.05 
    asterisk_sp='n.s.';
    a_fontsize=13;
else if p_val<0.05 && p_val>0.01
    asterisk_sp='*';
    a_fontsize=17;
    else if p_val<0.01 && p_val>0.001
            asterisk_sp='**';
            a_fontsize=17;
    else if p_val<0.001
             asterisk_sp='***';
             a_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(powerspec.gamma))*1.1; 
liney=[my ;my ];

g5=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
        y1limits(2)=max(max(powerspec.gamma))*1.2; %my;
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Normalized PSD [A.U.]', 'FontSize', 28,'fontname', 'arial');
        ylabel('PSD [mV^2/Hz]','fontsize',20, 'fontname', 'arial');
        title(['Gamma Power,  p=' num2str(p_val)] ,'FontSize', 20,'fontname', 'arial');   
       
% Sum of Power
tmp_Y= powerspec.Y_sum';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(tmp_Y');
if isfield(powerspec,'Y_sum_wilcoxon_p')
    p_val=powerspec.Y_sum_wilcoxon_p;
else
    p_val=powerspec.Y_sum_ttest_p;
end
if p_val >0.05 
    asterisk_sp='n.s.'; a_fontsize=13;
else if p_val<0.05 && p_val>0.01
    asterisk_sp='*';a_fontsize=17;
    else if p_val<0.01 && p_val>0.001
            asterisk_sp='**';a_fontsize=17;
    else if p_val<0.001
             asterisk_sp='***';a_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(powerspec.Y_sum))*1.1; 
liney=[my ;my ];

g6=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
        y1limits(2)=max(max(powerspec.Y_sum))*1.2; %my;
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('PSD [mV^2/Hz]', 'FontSize', 28,'fontname', 'arial');
        title(['Total Power,  p=' num2str(p_val)] ,'FontSize', 20,'fontname', 'arial');  
        %% Low-frequency [1-10] Power
tmp_Y= powerspec.lowfreq';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(tmp_Y');
if isfield(powerspec,'lowfreq_wilcoxon_p')
    p_val=powerspec.lowfreq_wilcoxon_p;
else
    p_val=powerspec.lowfreq_ttest_p;
end
if p_val >0.05 
    asterisk_sp='n.s.';
    a_fontsize=13;
else if p_val<0.05 && p_val>0.01
    asterisk_sp='*';
    a_fontsize=17;
    else if p_val<0.01 && p_val>0.001
            asterisk_sp='**';
            a_fontsize=17;
    else if p_val<0.001
             asterisk_sp='***';
             a_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(powerspec.lowfreq))*1.1; 
liney=[my ;my ];

g7=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
        y1limits(2)=max(max(powerspec.lowfreq))*1.2; %my;
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
%         ylabel('Normalized PSD [A.U.]', 'FontSize', 28,'fontname', 'arial');
        ylabel('PSD [mV^2/Hz]','fontsize',20, 'fontname', 'arial');
        title(['Low-frequency Power,  p=' num2str(p_val)] ,'FontSize', 20,'fontname', 'arial');   
  %% Diff Power
tmp_Y= powerspec.diff';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(tmp_Y');
if isfield(powerspec,'diff_wilcoxon_p')
    p_val=powerspec.diff_wilcoxon_p;
else
    p_val=powerspec.diff_ttest_p;
end
if p_val >0.05 
    asterisk_sp='n.s.';
    a_fontsize=13;
else if p_val<0.05 && p_val>0.01
    asterisk_sp='*';
    a_fontsize=17;
    else if p_val<0.01 && p_val>0.001
            asterisk_sp='**';
            a_fontsize=17;
    else if p_val<0.001
             asterisk_sp='***';
             a_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(powerspec.diff))*1.1; 
liney=[my ;my ];

g8=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
        y1limits(2)=max(max(powerspec.diff))*1.2; %my;
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Normalized PSD [A.U.]', 'FontSize', 28,'fontname', 'arial');
        title(['Diff Power,  p=' num2str(p_val)] ,'FontSize', 20,'fontname', 'arial');        
%% Save figures and data
 if save_flag==1
        cd(path_output)
        print(g1,'Delta_Power','-dpng','-r600','-opengl')
        saveas(g1,'Delta_Power','fig') 
        print(g2,'Theta_Power','-dpng','-r600','-opengl')
        saveas(g2,'Theta_Power','fig') 
        print(g3,'Alpha_Power','-dpng','-r600','-opengl')
        saveas(g3,'Alpha_Power','fig') 
        print(g4,'Beta_Power','-dpng','-r600','-opengl')
        saveas(g4,'Beta_Power','fig') 
        print(g5,'Gamma_Power','-dpng','-r600','-opengl')
        saveas(g5,'Gamma_Power','fig') 
        print(g6,'Total_Power','-dpng','-r600','-opengl')
        saveas(g6,'Total_Power','fig') 
        save('Powerspec', 'files_to_analyze', 'powerspec')
        print(g7,'Lowfreq_Power','-dpng','-r600','-opengl')
        saveas(g7,'Lowfreq_Power','fig') 
    end
end
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
  