%% for opening workspace saved 
clear all
global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data Ch2_data stim2_X stim1_X 
 global exp_type
exp_type=1; %1-NBES, 2-ChAT
trace_type_input=1; %
analyze_time_before_train=0;
analyze_train_only_flag=1;
save_flag= 1;
print_flag=0;
norm_flag=0;
no_DC_flag=0;
manual_validation_flag=0;
plot_eventsoverlay_flag=0;
clamp_flag=[]; %[]; %3; %clamp_flag=1 for hyperpolarization traces, clamp_flag=2 for depolarization traces and clamp_flag=3 for no current traces (only clamp to resting Vm)
BP50HzLFP_flag=0; %removing 50Hz noise from LFP signal
BP50HzVm_flag=1; %removing 50Hz noise from Vm signal
BPLFP_flag=0; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=[1,200]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=1; %filtering Vm
bp_manual_Vm=[0,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
%%
switch exp_type
    case 1
        files_to_analyze =[8,10,12,14,15,16,22,37,40,1,46,48,52,58,72,82,84]; %[8,10,12,14,15,16,22,37,40,1,44,46,48,52,56,58,62,72,75,82,84];  %[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75,82,84]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB-', 'NB+'};  y_ax_label={'Vm'}; y_ax_units={'mV'};  
        path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end        

    case 2
        files_to_analyze =[74,76,77,80,82,84,87,90,92,115]; %[74,76,77,80,82,84,87]; %,84,87];112,114,
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light Off', 'Light On'};  y_ax_label={'Vm'}; y_ax_units={'mV'};  
        path_output= 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end   
    case 3 
            files_to_analyze =[31,38,42,51,69,71,74]; %[31,38,42,51,61,64,67,69,71,74,77];61,64,
            cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
            load NBES_Files_v2
            legend_string={'NB-', 'NB+'};    y_ax_label={'Im'}; y_ax_units={'pA'};    
            path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis_VC\Evoked';   
            a = exist(path_output,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
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
    if isempty(data_no_spikes)
            current_data=raw_data{channel};
            if clamp_flag==1
                current_data=(-1).*raw_data{channel};
            end
            data_used='raw_data';
        else
        current_data=data_no_spikes{channel}; %raw_data{channel}; 
        if clamp_flag==1
            current_data=(-1).*data_no_spikes{channel}; %raw_data{channel};
        end
        data_used='data_no_spikes';
        end
data_preprocessing 

 if isempty(current_data_filt)
     current_data_filt=current_data;
 end
 
current_data_smooth=sgolayfilt(current_data_filt, 1,81); 

clear color_table    
    whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;
    switch exp_type
            case 1 
                color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
            case 2
                color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
           case 3 
                color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];        
        end
 %%    
   for trace_type=trace_type_input %1 for spont., 2 for evoked
         interval=[];  interval_temp=[];     
    intervals_to_analyze 
            epoch_length=0.5; %[sec]
                 epochs=round(duration/epoch_length);    
                  finalAmp_Thres =1.5;
                  if exp_type==3
                       finalAmp_Thres = 4;
                  end
                                                  
        for t=1:2;
            event_start_vec{t}=[];  event_onVal_vec{t} = [];  event_amplitude_vec{t} =  [];    event_ampPos_vec{t} =  []; 
            event_ampVal_vec{t} = [];     event_halfWidth_vec{t} =  [];   event_halfWidthS_vec{t} =  []; event_peak90prcntTime_vec{t}=[];
            event_halfWidthE_vec{t} =  [];       
          for trace=1:size(current_data,2)    
              event_start{t,trace}=[];                      event_onVal{t,trace} = [];
            event_amplitude{t,trace} =  [];             event_ampPos{t,trace} =  [];             event_ampVal{t,trace} = [];
            event_halfWidth{t,trace} =  [];            event_halfWidthS{t,trace} =  [];             event_halfWidthE{t,trace} =  []; event_peak90prcntTime{t,trace}=[];
            
              for epoch=1:epochs
                  starting=[]; amplitude=[]; ampPos=[]; halfWidth=[]; halfWidthS=[]; halfWidthE=[];
                  peak_start_int = interval(1+sf{1}*epoch_length*(epoch-1),t);
                  peak_end_int = interval(sf{1}*epoch_length*epoch,t);
        %         voltages_input = current_data(peak_start_int:peak_end_int,trace,x_value(t));
                   voltages_input = current_data_smooth(peak_start_int:peak_end_int,trace,x_value(t));

                   doPlot = 0;
                    I_temp=0;
                    [voltages, tmp_starting,tmp_amplitude, tmp_ampPos, halfWidth,tmp_halfWidthS, tmp_halfWidthE] = fn_EventDetector_v2(voltages_input, dt, finalAmp_Thres, doPlot,I_temp);
            %         title(['t=', num2str(t), ' trace ',num2str(trace)]); 
            %         set(gca,'ylim',[-60 -30])
            %          pause

                %correcting to the real locations in the trace
                     starting=(tmp_starting+peak_start_int-1)'; 
                     ampPos=(tmp_ampPos+peak_start_int-1)'; 
                     halfWidthS=(tmp_halfWidthS+peak_start_int-1)'; 
                     halfWidthE=(tmp_halfWidthE+peak_start_int-1)'; 
                     amplitude=tmp_amplitude'; 

                     %if an event ends after a second event has started -> take halfWidthE to be the starting point of
                     %the next event:
                     tmp=find(starting(2:end)-halfWidthE(1:end-1)<0);
                     if ~isempty(tmp)
                          halfWidthE(tmp)=starting(tmp+1);
                     end         
                    halfWidth=(halfWidthE-halfWidthS).*dt.*1000; %in msec
                    event_start{t,trace} =  [event_start{t,trace};  starting];            
                    event_onVal{t,trace} = [event_onVal{t,trace}; current_data(starting,trace,x_value(t))];
        %             event_amplitude{t,trace} =  [event_amplitude{t,trace}; amplitude]; 
                    event_ampPos{t,trace} =  [event_ampPos{t,trace}; ampPos]; 
                    event_ampVal{t,trace} = [event_ampVal{t,trace}; current_data(ampPos,trace,x_value(t))];
                    event_halfWidth{t,trace} =  [event_halfWidth{t,trace}; halfWidth]; 
                    event_halfWidthS{t,trace} =  [event_halfWidthS{t,trace}; halfWidthS]; 
                    event_halfWidthE{t,trace} =  [event_halfWidthE{t,trace}; halfWidthE]; 
                    event_amplitude{t,trace} =  [event_amplitude{t,trace}; current_data(ampPos,trace,x_value(t))-current_data(starting,trace,x_value(t))]; 
              end        
         for j=1:length(event_start{t,trace})
             peak_int=event_start{t,trace}(j): event_ampPos{t,trace}(j);            
             peak90prcntLoc(:,j)=find(current_data(peak_int,trace,x_value(t))-event_onVal{t,trace}(j)>=0.9.*event_amplitude{t,trace}(j),1);
             peak90prcntTime(:,j)=peak90prcntLoc(:,j).*dt.*1000; %[time from event onset to 90% of peak, in mS]
         end
            event_peak90prcntTime{t,trace}=peak90prcntTime;
            %% visualize the detected events on the trace         
% %         if manual_validation_flag==1
%             ymax=max([current_data([start_time(1).*sf{1}:(start_time(1)+duration).*sf{1}-1],:,x_value(1)),current_data([start_time(2).*sf{1}:(start_time(2)+duration).*sf{1}-1],:,x_value(2))]);
%             ymin=min([current_data([start_time(1).*sf{1}:(start_time(1)+duration).*sf{1}-1],:,x_value(1)),current_data([start_time(2).*sf{1}:(start_time(2)+duration).*sf{1}-1],:,x_value(2))]);
%             ydiff=ymax-ymin;
%             yrange=ceil(1.2*max(ydiff));
%             figure(1); clf
%             hold on
%                      h1=plot(start_time(t):dt:(start_time(t)+duration)-dt, current_data([start_time(t).*sf{1}:(start_time(t)+duration).*sf{1}-1],trace,x_value(t)),'k','linewidth',1.2);
%                      h2=scatter(event_start{t,trace}*dt,current_data(event_start{t,trace},trace,x_value(t)),40,'c','fill'); %mark event onset
%                      h3=scatter(event_ampPos{t,trace}*dt,current_data(event_ampPos{t,trace},trace,x_value(t)),40,'b','fill'); %mark event peak
% %                      h4=scatter(event_halfWidthS{t,trace}*dt,current_data(event_halfWidthS{t,trace},trace,x_value(t)),'c','fill'); %mark event half-width start
% %                      h5=scatter(event_halfWidthE{t,trace}*dt,current_data(event_halfWidthE{t,trace},trace,x_value(t)),'g','fill'); %mark event  half-width end
%                     set(gca,'xlim',[start_time(t), start_time(t)+duration]);    
% ydata=get(gca,'ylim');
% set(gca,'ylim',[ydata(1), ydata(1)+yrange])
% %        plotting scale bar
%                 horiz_vert=1;        lengthh=0.5;     textit=[num2str(lengthh.*1000), ' mS'];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
%                         [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
%                 horiz_vert=0;        lengthh=5;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12; perc1=[]; perc2=[];
%                         [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
% 
%                 set(gca, 'visible', 'off') ;
% %                    set(gca,'ylim',[-70 -10]);
%                      hold off
% %                 pause         
%     if save_flag==1;
%         cd(path_output)
%          if t==1;
%              filename1=['detection example - ongoing NB-, trace ', num2str(trace)];
%             print(1,filename1,'-dpng','-r600','-opengl')
%             saveas(1,filename1,'fig') 
%              else if t==2
%                   filename1=['detection example - ongoing NB+, trace ', num2str(trace)];
%                 print(1,filename1,'-dpng','-r600','-opengl')
%                 saveas(1,filename1,'fig') 
%                  end
%              end
%     end     
%%              
            %make one vector for all traces
            event_start_vec{t}=[event_start_vec{t}(:); event_start{t,trace}(:)];
            event_onVal_vec{t} = [event_onVal_vec{t}(:); event_onVal{t,trace}(:)];
            event_amplitude_vec{t} =  [event_amplitude_vec{t}(:); event_amplitude{t,trace}(:)]; 
            event_ampPos_vec{t} =  [event_ampPos_vec{t}(:); event_ampPos{t,trace}(:)]; 
            event_ampVal_vec{t} = [event_ampVal_vec{t}(:); event_ampVal{t,trace}(:)];
            event_halfWidth_vec{t} =  [event_halfWidth_vec{t}(:); event_halfWidth{t,trace}(:)]; 
            event_halfWidthS_vec{t} =  [event_halfWidthS_vec{t}(:); event_halfWidthS{t,trace}(:)]; 
            event_halfWidthE_vec{t} =  [event_halfWidthE_vec{t}(:); event_halfWidthE{t,trace}(:)];       
            event_peak90prcntTime_vec{t}=[event_peak90prcntTime_vec{t}(:); event_peak90prcntTime{t,trace}(:)];
          end
           event_count(1,t)=length(event_start_vec{t}(:));
           event_freq(1,t) = event_count(1,t)/(duration*size(current_data,2));
        end
   end

   %% plotting overlayed PSPs aligned to event peak on x-axis and onset value on the y-axis
% 
if plot_eventsoverlay_flag==1;
     for t=1:2
         row_counter=1;
           for tr=1:size(current_data,2)      
               for event=1:length(event_ampPos{t,tr})
                   Events_forplot{t}(row_counter,:)=current_data(event_ampPos{t,tr}(event)-150:event_ampPos{t,tr}(event)+249,tr,x_value(t)); % mean(current_data(event_ampPos{t,tr}(event)-200:event_ampPos{t,tr}(event)+299,tr));
%                    Events_forplot{t}(row_counter,:)=Events_forplot{t}(row_counter,:)-Events_forplot{t}(row_counter,1);
                   Events_forplot{t}(row_counter,:)=Events_forplot{t}(row_counter,:)-current_data(event_ampPos{t,tr}(event),tr,x_value(t))+event_amplitude{t,tr}(event); %-Events_forplot{t}(row_counter,1); 
                   %                    Events_forplot{t}(row_counter,:)=current_data(event_start{t,tr}(event)-50:event_start{t,tr}(event)+449,tr)-event_onVal{t,tr}(event); % mean(current_data(event_ampPos{t,tr}(event)-200:event_ampPos{t,tr}(event)+299,tr));
                    row_counter=row_counter+1;
               end
           end
       eventnum(1,t)=row_counter-1;
     end
     if min(eventnum)>40;
         eventnum(:)=40;
     end
   if exp_type==1
       trace_color=[[255 102 102]./256,0.4];
   else if exp_type==2;
        trace_color=[[102 102 255]./256,0.3];   
       end
   end
   fg1=figure;
        hold on
              l1=plot([dt:dt:(size(Events_forplot{1}(:,:),2)).*dt],Events_forplot{1}(1:min(eventnum),:)','color',[0 0 0 0.3],'linewidth',1);
              l2=plot([dt:dt:(size(Events_forplot{1}(:,:),2)).*dt],Events_forplot{2}(1:min(eventnum),:)','color',trace_color,'linewidth',1);
              l3=plot([dt:dt:(size(Events_forplot{1}(:,:),2)).*dt],mean(Events_forplot{1}',2),'color',color_table(1,:),'linewidth',2);
              l4=plot([dt:dt:(size(Events_forplot{1}(:,:),2)).*dt],mean(Events_forplot{2}',2),'color',color_table(2,:),'linewidth',2);
%               ylim1=get(gca,'ylim');
              ylim1(1)=min(min([Events_forplot{1}(1:min(eventnum),:);Events_forplot{2}(1:min(eventnum),:)]));   
              ylim1(2)=max(max([Events_forplot{1}(1:min(eventnum),:);Events_forplot{2}(1:min(eventnum),:)]));    

              peak=line([150*dt 150*dt],[ylim1(1) ylim1(2)],'linestyle','-.','color',[0 0 0 0.3],'linewidth',1);
              hold on
              onset=line([dt 399*dt],[0 0],'linestyle','-.','color',[0 0 0 0.3],'linewidth',1);
          hold off
% %plotting scale bar
horiz_vert=1;        lengthh=0.01;     textit=[num2str(lengthh.*1000), ' mS'];     c=[0,0,0];  fonsizes=12; perc1=-0.04; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=2;     textit=[num2str(lengthh),y_ax_units{1}];     c=[0,0,0];  fonsizes=12;perc1=0.07; perc2=0.05;
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
% set(gca, 'visible', 'off') ;
if save_flag==1;
    cd(path_output)
    filename1=['f' num2str(files_to_analyze(fileind)) '_events overlay'];    
        saveas(fg1,filename1,'fig') 
        print(fg1,filename1,'-dpng','-r600','-opengl') 
end
end
%% Calculating Vm mean and std - not used anymore. I am using the same code in the script "Analyze_response_parameters_10_20Hz_v2"
interval=[]; Vm_spont=[]; VmM=[]; VmSTD=[];
%taking 
switch exp_type             
    case 1
         start_time = [0.4,5.6]; %[sec] %[0,5]
         duration = 1; %[sec]                
    case 2
        start_time=[0.4, stim1_X{x_value(1)}(1,1).*dt+0.4];
        duration =1;                                                
end
for t=1:2
         start_sample(:,t) = ceil(start_time(t).*sf{1}); 
         end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
         interval(:,t) = start_sample(:,t):end_sample(:,t);    
% from intervals to analyze
                Vm_spont(:,:,t) = current_data(interval(:,t),:,x_value(t));
                VmM(:,t) =mean(mean(Vm_spont(:,:,t),2));
                VmSTD(:,t) = mean(std(Vm_spont(:,:,t),0,2)); %mean std across traces (trial-to-trial)
end               
%%
        %organizing the data in structure:  
            event_ongoing(fileind).cells=files_to_analyze;
            event_ongoing(fileind).analysis_mfile='Analyze_by_single_trials_ongoing_v2.m';
            event_ongoing(fileind).fname = fname;
            event_ongoing(fileind).x_value=x_value;
            event_ongoing(fileind).used_data=data_used;
            event_ongoing(fileind).clamp_flag=clamp_flag;
            event_ongoing(fileind).trace_type_input = trace_type_input;
            event_ongoing(fileind).BP50HzLFP=BP50HzLFP_flag; %removing 50Hz noise from LFP signal
            event_ongoing(fileind).BP50HzVm=BP50HzVm_flag; %removing 50Hz noise from Vm signal
            event_ongoing(fileind).BPLFP_flag=BPLFP_flag; %filtering LFP 
            event_ongoing(fileind).BPLFP=bp_filt_LFP; %the BP frequency filter for  LFP, used if BPLFP_flag=1
            event_ongoing(fileind).BPVm_flag=BPVm_flag; %filtering Vm 
            event_ongoing(fileind).BPVm=bp_filt_Vm; %the BP frequency filter for Vm, used if BPVm_flag=1
            event_ongoing(fileind).finalAmp_Thres=finalAmp_Thres;
            event_ongoing(fileind).start_sample=start_sample;
            event_ongoing(fileind).end_sample=end_sample;
            event_ongoing(fileind).start_time=start_time;
            event_ongoing(fileind).duration=duration;
            event_ongoing(fileind).sf=sf{1};
            event_ongoing(fileind).VmM=VmM;
            event_ongoing(fileind).VmSTD=VmSTD;
            event_ongoing(fileind).count = event_count;
            event_ongoing(fileind).freq = event_freq;
            event_ongoing(fileind).start = event_start_vec;
            event_ongoing(fileind).onVal = event_onVal_vec;
            event_ongoing(fileind).ampPos = event_ampPos_vec;           
            event_ongoing(fileind).ampVal = event_ampVal_vec; 
            event_ongoing(fileind).amplitude = event_amplitude_vec; 
            event_ongoing(fileind).halfWidth = event_halfWidth_vec; 
            event_ongoing(fileind).halfWidthS = event_halfWidthS_vec; 
            event_ongoing(fileind).halfWidthE = event_halfWidthE_vec;          
            event_ongoing(fileind).peak90prcntTime = event_peak90prcntTime_vec;
             
            event_count_cells(fileind,:)=event_count;
            event_freq_cells(fileind,:)=event_freq;
            event_VmM(fileind,:)=VmM;
            event_VmSTD(fileind,:)=VmSTD;
            event_onVal_m(fileind,:)=cellfun(@mean,event_onVal_vec);
            event_onVal_std(fileind,:)=cellfun(@std,event_onVal_vec);
            event_ampVal_m(fileind,:)=cellfun(@mean,event_ampVal_vec);
            event_ampVal_std(fileind,:)=cellfun(@std,event_ampVal_vec);
            event_amplitude_m(fileind,:)=cellfun(@mean,event_amplitude_vec);
            event_halfWidth_m(fileind,:)=cellfun(@mean,event_halfWidth_vec);
            event_halfWidth_std(fileind,:)=cellfun(@std,event_halfWidth_vec);
            event_peak90prcntTime_m(fileind,:)=cellfun(@mean,event_peak90prcntTime_vec); 

            %% cell stats:
%  onVal_vec_tmp= event_onVal_vec{1};   
% onVal_mat=cell2mat(onVal_vec_tmp); 
% amplitude_mat=cell2mat(event_amplitude); 
% ampVal_mat=cell2mat(event_ampVal); 
% halfWidth_mat=cell2mat(event_halfWidth); 
% tmp=[];
%         
%     %onset
% 
%         event_evoked_stat.fileind(fileind).onset=tmp';
%         norm_peak_onset_noES(:,1)= tmp(1,:)./tmp(1,:);
%         norm_peak_onset_ES(:,1)= tmp(2,:)./tmp(1,:);
% %         event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_onset = [norm_peak_onset_noES, norm_peak_onset_ES];
%         event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset,1);
%         event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset,0,1);  
%         event_evoked_stat.stim_num(stim_num).onset_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_m;
%         event_evoked_stat.stim_num(stim_num).onset_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_std;
%         %testing for normal distribution
%         diff_onset= tmp(2,:)- tmp(1,:);
%         [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_onset, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_onset] = lillietest(diff_onset);
%         %unpaired ttest 
% %         [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_norm_onset, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_norm_onset]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_onset(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_onset(:,2));
%         [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_onset, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_onset]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).onset(:,2));
% 
%         clear norm_peak_onset_noES norm_peak_onset_ES diff_onset tmp        
    end
 %% population statistics
 
    %frequency of spontaneous events
        tmp = event_freq_cells;  
        event_ongoing_stat.freq = event_freq_cells;
        norm_event_freq_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_freq_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_freq = [norm_event_freq_noES, norm_event_freq_ES];
        event_ongoing_stat.freq_m= nanmean(tmp,1);
        event_ongoing_stat.freq_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_freq= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_freq, event_ongoing_stat.lillietest_p_freq] = lillietest(diff_event_freq);
        %paired ttest
        [event_ongoing_stat.ttest_h_norm_freq, event_ongoing_stat.ttest_p_norm_freq]= ttest(event_ongoing_stat.norm_freq(:,1),event_ongoing_stat.norm_freq(:,2));
        [event_ongoing_stat.ttest_h_freq, event_ongoing_stat.ttest_p_freq]= ttest(event_ongoing_stat.freq(:,1),event_ongoing_stat.freq(:,2));
        [event_ongoing_stat.wilcoxon_p_freq, event_ongoing_stat.wilcoxon_h_freq]= signrank(event_ongoing_stat.freq(:,1),event_ongoing_stat.freq(:,2));

        clear norm_event_freq_noES norm_event_freq_ES diff_event_freq tmp
 
        %VmM
        tmp = event_VmM;  
        event_ongoing_stat.VmM = event_VmM;
        norm_event_VmM_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_VmM_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_VmM = [norm_event_VmM_noES, norm_event_VmM_ES];
        event_ongoing_stat.VmM_m= nanmean(tmp,1);
        event_ongoing_stat.VmM_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_VmM= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_VmM, event_ongoing_stat.lillietest_p_VmM] = lillietest(diff_event_VmM);
        %paired ttest
        [event_ongoing_stat.ttest_h_norm_VmM, event_ongoing_stat.ttest_p_norm_VmM]= ttest(event_ongoing_stat.norm_VmM(:,1),event_ongoing_stat.norm_VmM(:,2));
        [event_ongoing_stat.ttest_h_VmM, event_ongoing_stat.ttest_p_VmM]= ttest(event_ongoing_stat.VmM(:,1),event_ongoing_stat.VmM(:,2));
        [event_ongoing_stat.wilcoxon_p_VmM, event_ongoing_stat.wilcoxon_h_VmM]= signrank(event_ongoing_stat.VmM(:,1),event_ongoing_stat.VmM(:,2));

        clear norm_event_VmM_noES norm_event_VmM_ES diff_event_VmM tmp

        %VmSTD
        tmp = event_VmSTD;  
        event_ongoing_stat.VmSTD = event_VmSTD;
        norm_event_VmSTD_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_VmSTD_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_VmSTD = [norm_event_VmSTD_noES, norm_event_VmSTD_ES];
        event_ongoing_stat.VmSTD_m= nanmean(tmp,1);
        event_ongoing_stat.VmSTD_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_VmSTD= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_VmSTD, event_ongoing_stat.lillietest_p_VmSTD] = lillietest(diff_event_VmSTD);
        %paired ttest
        [event_ongoing_stat.ttest_h_norm_VmSTD, event_ongoing_stat.ttest_p_norm_VmSTD]= ttest(event_ongoing_stat.norm_VmSTD(:,1),event_ongoing_stat.norm_VmSTD(:,2));
        [event_ongoing_stat.ttest_h_VmSTD, event_ongoing_stat.ttest_p_VmSTD]= ttest(event_ongoing_stat.VmSTD(:,1),event_ongoing_stat.VmSTD(:,2));
        [event_ongoing_stat.wilcoxon_p_VmSTD, event_ongoing_stat.wilcoxon_h_VmSTD]= signrank(event_ongoing_stat.VmSTD(:,1),event_ongoing_stat.VmSTD(:,2));

        clear norm_event_VmSTD_noES norm_event_VmSTD_ES diff_event_VmSTD tmp
        
 %mean onset values events
        tmp = event_onVal_m;  
        event_ongoing_stat.onVal_m = event_onVal_m;
        norm_event_onVal_m_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_onVal_m_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_onVal_m = [norm_event_onVal_m_noES, norm_event_onVal_m_ES];
        event_ongoing_stat.onVal_m_m= nanmean(tmp,1);
        event_ongoing_stat.onVal_m_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_onVal_m= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_onVal_m, event_ongoing_stat.lillietest_p_onVal_m] = lillietest(diff_event_onVal_m);
        %paired ttest
        [event_ongoing_stat.ttest_h_norm_onVal_m, event_ongoing_stat.ttest_p_norm_onVal_m]= ttest(event_ongoing_stat.norm_onVal_m(:,1),event_ongoing_stat.norm_onVal_m(:,2));
        [event_ongoing_stat.ttest_h_onVal_m, event_ongoing_stat.ttest_p_onVal_m]= ttest(event_ongoing_stat.onVal_m(:,1),event_ongoing_stat.onVal_m(:,2));
        [event_ongoing_stat.wilcoxon_p_onVal_m,event_ongoing_stat.wilcoxon_h_onVal_m]= signrank(event_ongoing_stat.onVal_m(:,1),event_ongoing_stat.onVal_m(:,2));

        clear norm_event_onVal_m_noES norm_event_onVal_m_ES diff_event_onVal_m tmp

        %std of onset values events
        tmp = event_onVal_std;  
        event_ongoing_stat.onVal_std = event_onVal_std;
        norm_event_onVal_std_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_onVal_std_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_onVal_std = [norm_event_onVal_std_noES, norm_event_onVal_std_ES];
        event_ongoing_stat.onVal_std_m= nanmean(tmp,1);
        event_ongoing_stat.onVal_std_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_onVal_std= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_onVal_std, event_ongoing_stat.lillietest_p_onVal_std] = lillietest(diff_event_onVal_std);
        %paired ttest
        [event_ongoing_stat.ttest_h_norm_onVal_std, event_ongoing_stat.ttest_p_norm_onVal_std]= ttest(event_ongoing_stat.norm_onVal_std(:,1),event_ongoing_stat.norm_onVal_std(:,2));
        [event_ongoing_stat.ttest_h_onVal_std, event_ongoing_stat.ttest_p_onVal_std]= ttest(event_ongoing_stat.onVal_std(:,1),event_ongoing_stat.onVal_std(:,2));
        [event_ongoing_stat.wilcoxon_p_onVal_std, event_ongoing_stat.wilcoxon_h_onVal_std]= signrank(event_ongoing_stat.onVal_std(:,1),event_ongoing_stat.onVal_std(:,2));

        clear norm_event_onVal_std_noES norm_event_onVal_std_ES diff_event_onVal_std tmp
        
        %mean peak values of spontaneous events
        tmp = event_ampVal_m;  
        event_ongoing_stat.ampVal_m = event_ampVal_m;
        norm_event_ampVal_m_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_ampVal_m_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_ampVal_m = [norm_event_ampVal_m_noES, norm_event_ampVal_m_ES];
        event_ongoing_stat.ampVal_m_m= nanmean(tmp,1);
        event_ongoing_stat.ampVal_m_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_ampVal_m= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_ampVal_m, event_ongoing_stat.lillietest_p_ampVal_m] = lillietest(diff_event_ampVal_m);
        %paired ttest
        [event_ongoing_stat.ttest_h_norm_ampVal_m, event_ongoing_stat.ttest_p_norm_ampVal_m]= ttest(event_ongoing_stat.norm_ampVal_m(:,1),event_ongoing_stat.norm_ampVal_m(:,2));
        [event_ongoing_stat.ttest_h_ampVal_m, event_ongoing_stat.ttest_p_ampVal_m]= ttest(event_ongoing_stat.ampVal_m(:,1),event_ongoing_stat.ampVal_m(:,2));
        [event_ongoing_stat.wilcoxon_p_ampVal_m,event_ongoing_stat.wilcoxon_h_ampVal_m]= signrank(event_ongoing_stat.ampVal_m(:,1),event_ongoing_stat.ampVal_m(:,2));

        clear norm_event_ampVal_m_noES norm_event_ampVal_m_ES diff_event_ampVal_m tmp
 
         %std of peak values of spontaneous events
        tmp = event_ampVal_std;  
        event_ongoing_stat.ampVal_std = event_ampVal_std;
        norm_event_ampVal_std_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_ampVal_std_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_ampVal_std = [norm_event_ampVal_std_noES, norm_event_ampVal_std_ES];
        event_ongoing_stat.ampVal_std_m= nanmean(tmp,1);
        event_ongoing_stat.ampVal_std_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_ampVal_std= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_ampVal_std, event_ongoing_stat.lillietest_p_ampVal_std] = lillietest(diff_event_ampVal_std);
        %paired ttest
        [event_ongoing_stat.ttest_h_norm_ampVal_std, event_ongoing_stat.ttest_p_norm_ampVal_std]= ttest(event_ongoing_stat.norm_ampVal_std(:,1),event_ongoing_stat.norm_ampVal_std(:,2));
        [event_ongoing_stat.ttest_h_ampVal_std, event_ongoing_stat.ttest_p_ampVal_std]= ttest(event_ongoing_stat.ampVal_std(:,1),event_ongoing_stat.ampVal_std(:,2));
        [event_ongoing_stat.wilcoxon_p_ampVal_std, event_ongoing_stat.wilcoxon_h_ampVal_std]= signrank(event_ongoing_stat.ampVal_std(:,1),event_ongoing_stat.ampVal_std(:,2));

        clear norm_event_ampVal_std_noES norm_event_ampVal_std_ES diff_event_ampVal_std tmp
 
        %mean amplitude of spontaneous events
        tmp = event_amplitude_m;  
        event_ongoing_stat.amplitude_m = event_amplitude_m;
        norm_event_amplitude_m_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_amplitude_m_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_amplitude_m = [norm_event_amplitude_m_noES, norm_event_amplitude_m_ES];
        event_ongoing_stat.amplitude_m_m= nanmean(tmp,1);
        event_ongoing_stat.amplitude_m_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_amplitude_m= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_amplitude_m, event_ongoing_stat.lillietest_p_amplitude_m] = lillietest(diff_event_amplitude_m);
        %paired ttest
        [event_ongoing_stat.ttest_h_norm_amplitude_m, event_ongoing_stat.ttest_p_norm_amplitude_m]= ttest(event_ongoing_stat.norm_amplitude_m(:,1),event_ongoing_stat.norm_amplitude_m(:,2));
        [event_ongoing_stat.ttest_h_amplitude_m, event_ongoing_stat.ttest_p_amplitude_m]= ttest(event_ongoing_stat.amplitude_m(:,1),event_ongoing_stat.amplitude_m(:,2));
        [event_ongoing_stat.wilcoxon_p_amplitude_m,event_ongoing_stat.wilcoxon_h_amplitude_m]= signrank(event_ongoing_stat.amplitude_m(:,1),event_ongoing_stat.amplitude_m(:,2));

        clear norm_event_amplitude_m_noES norm_event_amplitude_m_ES diff_event_amplitude_m tmp
        
        %half width of spontaneous events
        tmp = event_halfWidth_m;  
        event_ongoing_stat.halfWidth_m = event_halfWidth_m;
        norm_event_halfWidth_m_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_halfWidth_m_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_halfWidth_m = [norm_event_halfWidth_m_noES, norm_event_halfWidth_m_ES];
        event_ongoing_stat.halfWidth_m_m= nanmean(tmp,1);
        event_ongoing_stat.halfWidth_m_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_halfWidth_m= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_halfWidth_m, event_ongoing_stat.lillietest_p_halfWidth_m] = lillietest(diff_event_halfWidth_m);
        %paired ttest 
        [event_ongoing_stat.ttest_h_norm_halfWidth_m, event_ongoing_stat.ttest_p_norm_halfWidth_m]= ttest(event_ongoing_stat.norm_halfWidth_m(:,1),event_ongoing_stat.norm_halfWidth_m(:,2));
        [event_ongoing_stat.ttest_h_halfWidth_m, event_ongoing_stat.ttest_p_halfWidth_m]= ttest(event_ongoing_stat.halfWidth_m(:,1),event_ongoing_stat.halfWidth_m(:,2));
         [event_ongoing_stat.wilcoxon_p_halfWidth_m, event_ongoing_stat.wilcoxon_h_halfWidth_m]= signrank(event_ongoing_stat.halfWidth_m(:,1),event_ongoing_stat.halfWidth_m(:,2));

        clear norm_event_halfWidth_m_noES norm_event_halfWidth_m_ES diff_event_halfWidth_m tmp
         
        %half width of spontaneous events
        tmp = event_peak90prcntTime_m;  
        event_ongoing_stat.peak90prcntTime_m = event_peak90prcntTime_m;
        norm_event_peak90prcntTime_m_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_peak90prcntTime_m_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_peak90prcntTime_m = [norm_event_peak90prcntTime_m_noES, norm_event_peak90prcntTime_m_ES];
        event_ongoing_stat.peak90prcntTime_m_m= nanmean(tmp,1);
        event_ongoing_stat.peak90prcntTime_m_std= nanstd(tmp,0,1);
        %testing for normal distribution
        diff_event_peak90prcntTime_m= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_peak90prcntTime_m, event_ongoing_stat.lillietest_p_peak90prcntTime_m] = lillietest(diff_event_peak90prcntTime_m);
        %paired ttest 
        [event_ongoing_stat.ttest_h_norm_peak90prcntTime_m, event_ongoing_stat.ttest_p_norm_peak90prcntTime_m]= ttest(event_ongoing_stat.norm_peak90prcntTime_m(:,1),event_ongoing_stat.norm_peak90prcntTime_m(:,2));
        [event_ongoing_stat.ttest_h_peak90prcntTime_m, event_ongoing_stat.ttest_p_peak90prcntTime_m]= ttest(event_ongoing_stat.peak90prcntTime_m(:,1),event_ongoing_stat.peak90prcntTime_m(:,2));
         [event_ongoing_stat.wilcoxon_p_peak90prcntTime_m, event_ongoing_stat.wilcoxon_h_peak90prcntTime_m]= signrank(event_ongoing_stat.peak90prcntTime_m(:,1),event_ongoing_stat.peak90prcntTime_m(:,2));

        clear norm_event_peak90prcntTime_m_noES norm_event_peak90prcntTime_m_ES diff_event_peak90prcntTime_m tmp

 %% plots
tmp_Y= event_ongoing_stat.freq';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = event_ongoing_stat.freq_std;
if event_ongoing_stat.wilcoxon_p_freq >0.05 
    asterisk_sp='n.s.';
else if event_ongoing_stat.wilcoxon_p_freq<0.05 && event_ongoing_stat.wilcoxon_p_freq>0.01
    asterisk_sp='*';
    else if event_ongoing_stat.wilcoxon_p_freq<0.01 && event_ongoing_stat.wilcoxon_p_freq>0.001
            asterisk_sp='**';
    else if event_ongoing_stat.wilcoxon_p_freq<0.001
             asterisk_sp='***';
        end
        end
    end
end
linex=[1;2];
my=max(max(event_ongoing_stat.freq))*1.1; 
liney=[my ;my ];
g1=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off
y1limits=get(gca,'ylim');
        x1limits = [0.75 2.25];  x1ticks = [1,2];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 'ylim',[-1,20]
        ylabel('Frequency [Hz]', 'FontSize', 28,'fontname', 'arial');
        title(['Spontaneous event Frequency, p=' num2str(event_ongoing_stat.wilcoxon_p_freq)] ,'FontSize', 20,'fontname', 'arial');   
%%
tmp_Y= event_ongoing_stat.amplitude_m';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = event_ongoing_stat.amplitude_m_std;
if event_ongoing_stat.wilcoxon_p_amplitude_m >0.05 
    asterisk_sp='n.s.';
else if event_ongoing_stat.wilcoxon_p_amplitude_m<0.05 && event_ongoing_stat.wilcoxon_p_amplitude_m>0.01
    asterisk_sp='*';
    else if event_ongoing_stat.wilcoxon_p_amplitude_m<0.01 && event_ongoing_stat.wilcoxon_p_amplitude_m>0.001
            asterisk_sp='**';
    else if event_ongoing_stat.wilcoxon_p_amplitude_m<0.001
             asterisk_sp='***';
        end
        end
    end
end
linex=[1;2];
my=max(max(event_ongoing_stat.amplitude_m))*1.1-2; 
liney=[my ;my ];

g2=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel(['Amplitude [', y_ax_units{1},']'], 'FontSize', 28,'fontname', 'arial');
        title(['Spontaneous event Amplitude,  p=' num2str(event_ongoing_stat.wilcoxon_p_amplitude_m)] ,'FontSize', 20,'fontname', 'arial');   
%%
tmp_Y= event_ongoing_stat.halfWidth_m';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = event_ongoing_stat.halfWidth_m_std;
if event_ongoing_stat.wilcoxon_p_halfWidth_m >0.05 
    asterisk_sp='n.s.';
else if event_ongoing_stat.wilcoxon_p_halfWidth_m<0.05 && event_ongoing_stat.wilcoxon_p_halfWidth_m>0.01
    asterisk_sp='*';
    else if event_ongoing_stat.wilcoxon_p_halfWidth_m<0.01 && event_ongoing_stat.wilcoxon_p_halfWidth_m>0.001
            asterisk_sp='**';
    else if event_ongoing_stat.wilcoxon_p_halfWidth_m<0.001
             asterisk_sp='***';
        end
        end
    end
end
linex=[1;2];
my=max(max(event_ongoing_stat.halfWidth_m))+2; 
liney=[my ;my ];
g3=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)

hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Half-width [mS]', 'FontSize', 28,'fontname', 'arial');
        title(['Spontaneous event Half-width,  p=' num2str(event_ongoing_stat.wilcoxon_p_halfWidth_m)] ,'FontSize', 20,'fontname', 'arial');   
%%
tmp_Y= event_ongoing_stat.onVal_m';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = event_ongoing_stat.onVal_m_std;
g4=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel(['Event Onset value [', y_ax_units{1},']'], 'FontSize', 28,'fontname', 'arial');
        title(['Spontaneous event onset value,  p=' num2str(event_ongoing_stat.wilcoxon_p_onVal_m)] ,'FontSize', 20,'fontname', 'arial');
        
%% peak90prcntTime_m
tmp_Y= event_ongoing_stat.peak90prcntTime_m';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = event_ongoing_stat.peak90prcntTime_m_std;
if event_ongoing_stat.wilcoxon_p_peak90prcntTime_m >0.05 
    asterisk_sp='n.s.';
else if event_ongoing_stat.wilcoxon_p_peak90prcntTime_m<0.05 && event_ongoing_stat.wilcoxon_p_peak90prcntTime_m>0.01
    asterisk_sp='*';
    else if event_ongoing_stat.wilcoxon_p_peak90prcntTime_m<0.01 && event_ongoing_stat.wilcoxon_p_peak90prcntTime_m>0.001
            asterisk_sp='**';
    else if event_ongoing_stat.wilcoxon_p_peak90prcntTime_m<0.001
             asterisk_sp='***';
        end
        end
    end
end
my=max(max(event_ongoing_stat.peak90prcntTime_m))-1; 
g5=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Event Rise-Time [mS]', 'FontSize', 28,'fontname', 'arial');
        title(['Spontaneous event rise-time to 90% peak,  p=' num2str(event_ongoing_stat.wilcoxon_p_peak90prcntTime_m)] ,'FontSize', 20,'fontname', 'arial');

        %% VmM
tmp_Y= event_ongoing_stat.VmM';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = event_ongoing_stat.VmM_std;
if event_ongoing_stat.wilcoxon_p_VmM >0.05 
    asterisk_sp='n.s.';
else if event_ongoing_stat.wilcoxon_p_VmM<0.05 && event_ongoing_stat.wilcoxon_p_VmM>0.01
    asterisk_sp='*';
    else if event_ongoing_stat.wilcoxon_p_VmM<0.01 && event_ongoing_stat.wilcoxon_p_VmM>0.001
            asterisk_sp='**';
    else if event_ongoing_stat.wilcoxon_p_VmM<0.001
             asterisk_sp='***';
        end
        end
    end
end
my=max(max(event_ongoing_stat.VmM))+3; %+5; 
g6=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        switch exp_type
            case 1
                y1limits = [-80 -30];
            case 2
                y1limits = [-70 -45];
        end
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks, 'ylim', y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Vm mean [mV]', 'FontSize', 28,'fontname', 'arial');
        title(['VmM,  p=' num2str(event_ongoing_stat.wilcoxon_p_VmM)] ,'FontSize', 20,'fontname', 'arial');

        %% VmSTD
tmp_Y= event_ongoing_stat.VmSTD';
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = event_ongoing_stat.VmSTD_std;
if event_ongoing_stat.wilcoxon_p_VmSTD >0.05 
    asterisk_sp='n.s.';
else if event_ongoing_stat.wilcoxon_p_VmSTD<0.05 && event_ongoing_stat.wilcoxon_p_VmSTD>0.01
    asterisk_sp='*';
    else if event_ongoing_stat.wilcoxon_p_VmSTD<0.01 && event_ongoing_stat.wilcoxon_p_VmSTD>0.001
            asterisk_sp='**';
    else if event_ongoing_stat.wilcoxon_p_VmSTD<0.001
             asterisk_sp='***';
        end
        end
    end
end
my=max(max(event_ongoing_stat.VmSTD))+0.5; 
g7=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 11];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim', y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Vm std [mV]', 'FontSize', 28,'fontname', 'arial');
        title(['VmSTD,  p=' num2str(event_ongoing_stat.wilcoxon_p_VmSTD)] ,'FontSize', 20,'fontname', 'arial');

        %% save figures
if save_flag==1
cd(path_output)
print(g1,'Spontaneous event frequency','-dpng','-r600','-opengl')
saveas(g1,'Spontaneous event frequency','fig') 
print(g2,'Spontaneous event amplitude','-dpng','-r600','-opengl')
saveas(g2,'Spontaneous event amplitude','fig') 
print(g3,'Spontaneous event half-width','-dpng','-r600','-opengl')
saveas(g3,'Spontaneous event half-width','fig') 
print(g4,'Spontaneous event onset value','-dpng','-r600','-opengl')
saveas(g4,'Spontaneous event onset value','fig') 
print(g5,'Spontaneous event rise-time','-dpng','-r600','-opengl')
saveas(g5,'Spontaneous event rise-time','fig') 
print(g6,'Spontaneous VmM','-dpng','-r600','-opengl')
saveas(g6,'Spontaneous VmM','fig') 
print(g7,'Spontaneous VmSTD','-dpng','-r600','-opengl')
saveas(g7,'Spontaneous VmSTD','fig') 
filename='Spontaneous activity'; 
save(filename, 'files_to_analyze', 'event_ongoing', 'event_ongoing_stat')
end
