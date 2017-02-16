%% for opening workspace saved 
clear all
global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data Ch2_data stim2_X stim1_X 
 global exp_type
exp_type=1; %1-NBES, 2-ChAT, 3-NBES VC
trace_type_input=2; %
analyze_time_before_train=0.1;
analyze_train_only_flag=0;
save_flag= 1;
print_flag=0;
paired_plot_flag=0; 
norm_flag=0;
clamp_flag=3; %[]; %3; %clamp_flag=1 for hyperpolarization traces, clamp_flag=2 for depolarization traces and clamp_flag=3 for no current traces (only clamp to resting Vm)
BP50HzLFP_flag=0; %removing 50Hz noise from LFP signal
BP50HzVm_flag=1; %removing 50Hz noise from Vm signal
BPLFP_flag=0; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=[1,200]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=1; %filtering LFP and Vm same as LFP was filtered in the multiclamp
bp_manual_Vm=[0,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)

%%
switch exp_type
    case 1
        files_to_analyze =[8,10,12,14,15,16,22,36,37,40,1,44,46,48,52,56,58,62,72,75,82,84];  %[16,22,36,37,40,44,46,48,52,56,58,62,72,75,82,84]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};  y_ax_label={'Vm'}; y_ax_units={'mV'};   
        path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Evoked';
        if paired_plot_flag==1; 
            path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Evoked paired plots';
        end
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end        

    case 2
        files_to_analyze =[74,76,77,80,82,84,87]; %
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light On', 'Light Off'};  y_ax_label={'Vm'}; y_ax_units={'mV'};   
        path_output= 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Evoked';
         if paired_plot_flag==1; 
             path_output= 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Evoked paired plots';
         end
        a = exist(path_output,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end   
    case 3 
        files_to_analyze =[31,38,42,51,69,71,74]; %[31,38,42,51,61,64,67,69,71,74,77];
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};    y_ax_label={'Im'}; y_ax_units={'pA'};    
        path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis_VC\Evoked';   
        if paired_plot_flag==1; 
            path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis_VC\Evoked paired plots';   
        end
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
   for trace_type= trace_type_input      
          interval=[];   interval_temp=[];     
     intervals_to_analyze 
            finalAmp_Thres = -1 ; %negative values are interpreted in terms of std. positive values are fixed amplitude in mV       
%% detect events
% for i=1:10;
% voltages_input = current_data(i/dt:(i+1)/dt,1,2);
        for t=1:2;
          for trace=1:size(current_data,2)             
              for stim_num=1:11; %1:galvano_nstim;
                  count=0;
                  event_on_nonspecific{t,stim_num,trace}=[]; event_amplitude_nonspecific{t,stim_num,trace}=[]; event_ampPos_nonspecific{t,stim_num,trace}=[];
                  event_halfWidth_nonspecific{t,stim_num,trace}=[]; event_halfWidthS_nonspecific{t,stim_num,trace}=[]; event_halfWidthE_nonspecific{t,stim_num,trace}=[];
                  event_start_nonspecific{t,stim_num,trace}=[]; event_ampDel_nonspecific{t,stim_num,trace}=[];
                  stim_ISI =1/galvano_freq.*sf{1};
                  if stim_num>size(stim2{t},2); %if there is no test stimulus in this cell
                        event_on{t,stim_num,trace} = nan;
                        event_onVal{t,stim_num,trace} = nan;
                        event_start{t,stim_num,trace} = nan;
                        event_amplitude{t,stim_num,trace} = nan;
                        event_ampPos{t,stim_num,trace} = nan;
                        event_ampDel{t,stim_num,trace} = nan;
                        event_ampVal{t,stim_num,trace} = nan;
                        event_halfWidth{t,stim_num,trace} = nan;
                        event_halfWidthS{t,stim_num,trace} = nan;
                        event_halfWidthE{t,stim_num,trace} = nan;
                        failures{t,stim_num,trace} = nan;
                        event_peak10prcntTime{t,stim_num,trace}=nan;
                        event_peak90prcntTime{t,stim_num,trace}=nan;
                        event_peak10to90Time{t,stim_num,trace}=nan;
                        event_nonspecific_count{t,stim_num,trace}=nan;
                                  continue
                 else
                      peak_start_int = round(stim2{t}(1,stim_num)+0.003*sf{1});
                      peak_end_int = round(stim2{t}(1,stim_num)+2*stim_ISI-0.003*sf{1});    
                  end
%                   peak_start_int = interval(1+stim_ISI*(stim_num-1)+0.003*sf{1},t); %not suitable for the test stim
%                   peak_end_int = interval(stim_ISI*(stim_num+1),t);
%                   peak_start_int = interval(1,t);
%                   peak_end_int = interval(end,t);
        voltages_input = current_data_filt(peak_start_int:peak_end_int,trace,x_value(t));
%         finalAmp_Thres = 1 ;
        doPlot = 0;
        I_temp=0;
        [voltages, starting, amplitude, ampPos, halfWidth,halfWidthS, halfWidthE] = fn_EventDetector_v2(voltages_input, dt, finalAmp_Thres, doPlot,I_temp);
%  figure(1); clf
%         title(['t=', num2str(t), ' trace ',num2str(trace),' stim. ',num2str(stim_num)]);  
%             hold on
%                      h1=plot([1:size(voltages_input,1)].*dt, voltages_input,'k');
%                      h2=scatter(starting*dt,voltages_input(starting),'r','fill'); %mark event onset
%                      h3=scatter(ampPos*dt,voltages_input(ampPos),'b','fill'); %mark event peak
% %                 set(gca,'ylim',[-50 -20]); set(gca,'xlim',[0 1.1]);  
%           hold off
% %         plot whisker stim
%       if stim_num+2>size(stim2{t},2); %if there is no test stimulus in this cell
%       else
%         ylim_data=[get(gca,'ylim')]';
%         patch_xdata=[stim2{t}(:,stim_num+1:stim_num+2); flipud(stim2{t}(:,stim_num+1:stim_num+2))];
%         patch_xdata = (round(patch_xdata)-peak_start_int).*dt;
%         yex=wextend('1D','sym',ylim_data,1);
%         l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
%         temp_y=wextend('ac','sym',yex,l);
%         patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
%         patch_cdata=ones(size(patch_xdata));
%         p=patch(patch_xdata,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
%         set(gca,'linewidth',1.2)
%       end
%           hold off           
%          pause
%               end   %temporary
%           end     %temporary
%         end     %temporary. real END is in line 273
%  cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz\single trials'
%  if t==1;
% print(1,['detection example - evoked NB-, trace ', num2str(trace)],'-dpng','-r600','-opengl')
% % saveas(1,'detection example - evoked NB-.fig') 
%  else if t==2
%   print(1,['detection example - evoked NB+, trace ', num2str(trace)],'-dpng','-r600','-opengl')
% % saveas(1,'detection example - evoked NB+.fig') 
%      end
%  end
    %correcting to the real locations in the trace
         starting=starting+peak_start_int-1; 
         ampPos=ampPos+peak_start_int-1; 
         halfWidthS=halfWidthS+peak_start_int-1; 
         halfWidthE=halfWidthE+peak_start_int-1; 
    % taking only the events during stim_num    
%     next_stim_pos =stim2{t}(1,1)+stim_ISI*(stim_num);
%     if stim_num==11;
        next_stim_pos =stim2{t}(1,stim_num)+stim_ISI;
%     end
         ampPos(find(starting>=next_stim_pos))=[];   
         amplitude(find(starting>=next_stim_pos))=[];
         halfWidth(find(starting>=next_stim_pos))=[];
         halfWidthS(find(starting>=next_stim_pos))=[];   
         halfWidthE(find(starting>=next_stim_pos))=[];   
         starting(find(starting>=next_stim_pos))=[];   
        %if en event ends after the next stim had started - take halfWidthE to be the beginning of the next
        %stim.
%          tmp=find(halfWidthE>=next_stim_pos);
%           if ~isempty(tmp)
%              halfWidthE(tmp)=next_stim_pos;
%              halfWidth(tmp)=halfWidthE(tmp)-halfWidthS(tmp)*dt;
%              [amplitude(tmp),ampPos(tmp)]=max(current_data(starting(tmp): halfWidthE(tmp),trace,x_value(t))); 
%       %need to recalculate the halfwidth according to the new amplitude. I leave it for now...
%         end

         %if an event ends after a second event has started -> take halfWidthE to be the starting point of
         %the next event:
         tmp=find(starting(2:end)-halfWidthE(1:end-1)<0);
         if ~isempty(tmp)
              halfWidthE(tmp)=starting(tmp+1);
         end  
       halfWidth=(halfWidthE-halfWidthS).*dt.*1000; %in msec 
       if stim_num>size(stim2_X{x_value(t)},2) %will happen if in this cell there was no "test" stimulus.
           starting=[];
       end
        if isempty(starting);
            failures{t,stim_num,trace}=1;
%             if stim_num>galvano_nstim %if this is a cell with no test stimulus then this is not a failure
%                failures{t,stim_num,trace}=nan;
%             end
            %if this was a failure then all response parameters will be NaNs:
            event_on{t,stim_num,trace} = nan;
            event_onVal{t,stim_num,trace} = nan;
            event_start{t,stim_num,trace} = nan;
            event_amplitude{t,stim_num,trace} = nan;
            event_ampPos{t,stim_num,trace} = nan;
            event_ampDel{t,stim_num,trace} = nan;
            event_ampVal{t,stim_num,trace} = nan;
            event_halfWidth{t,stim_num,trace} = nan;
            event_halfWidthS{t,stim_num,trace} = nan;
            event_halfWidthE{t,stim_num,trace} = nan;
            event_peak10prcntTime{t,stim_num,trace}=nan;
            event_peak90prcntTime{t,stim_num,trace}=nan;
            event_peak10to90Time{t,stim_num,trace}=nan;
        elseif starting(1)>(peak_start_int+0.03*sf{1}) %if the delay to the event was too large (more than 30ms) then this is not a response but a non-specific event
             failures{t,stim_num,trace}=1;
            event_start{t,stim_num,trace} = nan;
            event_on{t,stim_num,trace} = nan;
            event_onVal{t,stim_num,trace} = nan;
            event_amplitude{t,stim_num,trace} = nan;
            event_ampPos{t,stim_num,trace} = nan;
            event_ampDel{t,stim_num,trace} = nan;
            event_ampVal{t,stim_num,trace} = nan;
            event_halfWidth{t,stim_num,trace} = nan;
            event_halfWidthS{t,stim_num,trace} = nan;
            event_halfWidthE{t,stim_num,trace} = nan;
            event_peak90prcntTime{t,stim_num,trace}=nan;
            event_peak10prcntTime{t,stim_num,trace}=nan;
            event_peak10to90Time{t,stim_num,trace}=nan;
            event_start_nonspecific{t,stim_num,trace}= starting(1);
            event_on_nonspecific{t,stim_num,trace}= (starting(1)-stim2{t}(1,stim_num))*dt*1000; %in msec
            event_amplitude_nonspecific{t,stim_num,trace}= amplitude(1);
            event_ampPos_nonspecific{t,stim_num,trace}= ampPos(1);
            event_ampDel_nonspecific{t,stim_num,trace}= (ampPos(1)-stim2{t}(1,stim_num))*dt*1000; %in msec
            event_halfWidth_nonspecific{t,stim_num,trace} =  halfWidth(1); 
            event_halfWidthS_nonspecific{t,stim_num,trace} =  halfWidthS(1); 
            event_halfWidthE_nonspecific{t,stim_num,trace} =  halfWidthE(1); 
            count = count+1;
%need to take as response amplitude the values at ampPos minus at starting,
%in the Raw Data and not in the filtered data. and then test it.
        else %in this case it is a real response
            failures{t,stim_num,trace}=0;
            event_start{t,stim_num,trace} =  starting(1);
            event_on{t,stim_num,trace} =  (starting(1)-stim2{t}(1,stim_num))*dt*1000; %in msec
            event_onVal{t,stim_num,trace} = current_data(starting(1),trace,x_value(t));
            event_ampPos{t,stim_num,trace} =  ampPos(1); 
            event_ampDel{t,stim_num,trace} =  (ampPos(1)-stim2{t}(1,stim_num))*dt*1000;  %in msec
            event_ampVal{t,stim_num,trace} = current_data(ampPos(1),trace,x_value(t));
            event_halfWidth{t,stim_num,trace} =  halfWidth(1); 
            event_halfWidthS{t,stim_num,trace} =  halfWidthS(1); 
            event_halfWidthE{t,stim_num,trace} =  halfWidthE(1); 
            event_amplitude{t,stim_num,trace} = current_data(ampPos(1),trace,x_value(t))-current_data(starting(1),trace,x_value(t)); 
           
%             event_amplitude{t,stim_num,trace} =  amplitude(1); 
% finding the delay to 10%peak and 90%peak and the rise-time (time from 10%to90% peak)
             peak_int=event_start{t,stim_num,trace}: event_ampPos{t,stim_num,trace};            
             event_peak90prcntLoc{t,stim_num,trace}=find((current_data(peak_int,trace,x_value(t))-event_onVal{t,stim_num,trace})>=0.9.*event_amplitude{t,stim_num,trace},1);
             event_peak10prcntLoc{t,stim_num,trace}=find((current_data(peak_int,trace,x_value(t))-event_onVal{t,stim_num,trace})>=0.1.*event_amplitude{t,stim_num,trace},1);
             event_peak90prcntTime{t,stim_num,trace}=event_peak90prcntLoc{t,stim_num,trace}.*dt.*1000+3; %[time from stimulus onset (hence added 3ms) to 90% of peak, in mS]
             event_peak10prcntTime{t,stim_num,trace}=event_peak10prcntLoc{t,stim_num,trace}.*dt.*1000+3; %[time from stimulus onset (hence added 3ms) to 10% of peak, in mS]
             event_peak10to90Time{t,stim_num,trace}=event_peak90prcntTime{t,stim_num,trace}-event_peak10prcntTime{t,stim_num,trace};

              if event_amplitude{t,stim_num,trace}<=0
                event_amplitude{t,stim_num,trace}=nan;
              end
            
        end
        if length(starting)>1
            event_start_nonspecific{t,stim_num,trace}=[event_start_nonspecific{t,stim_num,trace}, starting(2:end)];
            event_on_nonspecific{t,stim_num,trace}=[event_on_nonspecific{t,stim_num,trace}, (starting(2:end)-stim2{t}(1,stim_num)).*dt*1000]; %in msec
            event_amplitude_nonspecific{t,stim_num,trace}= [event_amplitude_nonspecific{t,stim_num,trace}, amplitude(2:end)];
            event_ampPos_nonspecific{t,stim_num,trace}= [event_ampPos_nonspecific{t,stim_num,trace}, ampPos(2:end)]; 
            event_ampDel_nonspecific{t,stim_num,trace}= [event_ampDel_nonspecific{t,stim_num,trace}, (ampPos(2:end)-stim2{t}(1,stim_num)).*dt*1000]; %in msec
            event_halfWidth_nonspecific{t,stim_num,trace} =  [event_halfWidth_nonspecific{t,stim_num,trace}, halfWidth(2:end)]; 
            event_halfWidthS_nonspecific{t,stim_num,trace} =  [event_halfWidthS_nonspecific{t,stim_num,trace}, halfWidthS(2:end)]; 
            event_halfWidthE_nonspecific{t,stim_num,trace} =  [event_halfWidthE_nonspecific{t,stim_num,trace}, halfWidthE(2:end)]; 
            count = count+length(starting(2:end));
        end
        event_nonspecific_count{t,stim_num,trace}=count;
%    figure(2); %stops whenever there is nan... need to fix
%             hold on
%                      h1=plot([1:size(current_data,1)].*dt, current_data(:,trace,x_value(t)),'k');
%                      h2=scatter(event_start{t,stim_num,trace}(:)*dt,current_data(event_start{t,stim_num,trace}(:),trace,x_value(t)),'r','fill'); %mark event onset
%                      h3=scatter(event_ampPos{t,stim_num,trace}(:)*dt,current_data(event_ampPos{t,stim_num,trace}(:),trace,x_value(t)),'b','fill'); %mark event peak
%                      hold off
% pause
              end
           end
        end
   end

        %organizing the data in structure:     
            event_evoked(fileind).cells=files_to_analyze;
            event_evoked(fileind).x_value=x_value;
            event_evoked(fileind).used_data=data_used;
            event_evoked(fileind).clamp_flag=clamp_flag;
            event_evoked(fileind).analysis_mfile='Analyze_by_single_trials_evoked_v2.m';
            event_evoked(fileind).fname = fname;
            event_evoked(fileind).win_duration = duration;
            event_evoked(fileind).trace_type_input = trace_type_input;
            event_evoked(fileind).analyze_time_before_train = analyze_time_before_train;
            event_evoked(fileind).analyze_train_only_flag = analyze_train_only_flag;
            event_evoked(fileind).BP50HzLFP=BP50HzLFP_flag; %removing 50Hz noise from LFP signal
            event_evoked(fileind).BP50HzVm=BP50HzVm_flag; %removing 50Hz noise from Vm signal
            event_evoked(fileind).BPLFP_flag=BPLFP_flag; %filtering LFP 
            event_evoked(fileind).BPLFP=bp_filt_LFP; %the BP frequency filter for  LFP, used if BPLFP_flag=1
            event_evoked(fileind).BPVm_flag=BPVm_flag; %filtering Vm 
            event_evoked(fileind).BPVm=bp_filt_Vm; %the BP frequency filter for Vm, used if BPVm_flag=1
            event_evoked(fileind).finalAmp_Thres=finalAmp_Thres;
            event_evoked(fileind).start_sample=start_sample;
            event_evoked(fileind).end_sample=end_sample;
            event_evoked(fileind).sf=sf{1};
            event_evoked(fileind).failures = failures;
            event_evoked(fileind).start = event_start;
            event_evoked(fileind).onset = event_on;
            event_evoked(fileind).onVal = event_onVal;
            event_evoked(fileind).ampPos = event_ampPos; 
            event_evoked(fileind).ampDel = event_ampDel; 
            event_evoked(fileind).ampVal = event_ampVal; 
            event_evoked(fileind).amplitude = event_amplitude; 
            event_evoked(fileind).halfWidth = event_halfWidth; 
            event_evoked(fileind).halfWidthS = event_halfWidthS; 
            event_evoked(fileind).halfWidthE = event_halfWidthE; 
            event_evoked(fileind).peak90prcntTime=event_peak90prcntTime; 
            event_evoked(fileind).peak10prcntTime=event_peak10prcntTime; 
            event_evoked(fileind).peak10to90Time=event_peak10to90Time; 
            event_evoked(fileind).nonspecific_count =event_nonspecific_count;
            event_evoked(fileind).onset_nonspecific = event_on_nonspecific;
            event_evoked(fileind).amplitude_nonspecific = event_amplitude_nonspecific;
            event_evoked(fileind).ampPos_nonspecific = event_ampPos_nonspecific;
            event_evoked(fileind).ampDel_nonspecific = event_ampDel_nonspecific;
            event_evoked(fileind).halfWidth_nonspecific = event_halfWidth_nonspecific; 
            event_evoked(fileind).halfWidthS_nonspecific = event_halfWidthS_nonspecific; 
            event_evoked(fileind).halfWidthE_nonspecific = event_halfWidthE_nonspecific; 
            
            %% cell stats: within cell there is no logic to use normalized values
clear onset_mat onVal_mat   amplitude_mat    ampDel_mat  ampVal_mat  halfWidth_mat peak90prcntTime_mat peak10prcntTime_mat peak10to90Time_mat...
    failures_mat nonspecific_count_mat  amplitude_mat_M adapt_amp adapt_amp_M
onset_mat=cell2mat(event_on); 
onVal_mat=cell2mat(event_onVal); 
amplitude_mat=cell2mat(event_amplitude); 
ampDel_mat=cell2mat(event_ampDel); 
ampVal_mat=cell2mat(event_ampVal); 
halfWidth_mat=cell2mat(event_halfWidth); 
peak90prcntTime_mat=cell2mat(event_peak90prcntTime); 
peak10prcntTime_mat=cell2mat(event_peak10prcntTime); 
peak10to90Time_mat=cell2mat(event_peak10to90Time); 
failures_mat=cell2mat(failures); 
nonspecific_count_mat = cell2mat(event_nonspecific_count);
amplitude_mat_M=nanmean(amplitude_mat,3);
adapt_amp=nanmean(amplitude_mat(:,8:10,:),2)./amplitude_mat(:,1,:);
% adapt_amp_M=nanmean(amplitude_mat_M(:,8:10),2)./nanmean(amplitude_mat_M(:,1),2);
adapt_amp_M=nanmean(adapt_amp,3);
event_evoked(fileind).adapt_amp(:,:)=permute(adapt_amp(1:2,1,:),[1 3 2])'; 

tmp=[];

    for stim_num=1:11;      
%adaptaion amplitude ratio 
         tmp(:,:) = adapt_amp(1:2,1,:);  
        event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp=tmp';       
        event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp,0,1);  
        event_evoked_stat.stim_num(stim_num).adapt_amp_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp_m;
        event_evoked_stat.stim_num(stim_num).adapt_amp_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp_std;
        %testing for normal distribution
        if (nnz(diff(tmp))-sum(isnan(diff(tmp))))>4; %if there are at least 5 values that are non-zero and not NaNs
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_adapt_amp, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_adapt_amp] = lillietest(diff(tmp));        
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_adapt_amp, event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_adapt_amp]= signrank(event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp(:,2));
        end
        %paired ttest 
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_adapt_amp, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_adapt_amp]= ttest(event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).adapt_amp(:,2));

        clear  tmp
     
     %failures
        tmp(:,:) = failures_mat(1:2,stim_num,:);  
        event_evoked_stat.fileind(fileind).stim_num(stim_num).failures=tmp';       
        event_evoked_stat.fileind(fileind).stim_num(stim_num).failures_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).failures,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).failures_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).failures,0,1);  
        event_evoked_stat.stim_num(stim_num).failures_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).failures_m;
        event_evoked_stat.stim_num(stim_num).failures_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).failures_std;
        %testing for normal distribution
        if (nnz(diff(tmp))-sum(isnan(diff(tmp))))>4; %if there are at least 5 values that are non-zero and not NaNs
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_failures, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_failures] = lillietest(diff(tmp));        
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_failures, event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_failures]= signrank(event_evoked_stat.fileind(fileind).stim_num(stim_num).failures(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).failures(:,2));
        end
        %paired ttest 
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_failures, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_failures]= ttest(event_evoked_stat.fileind(fileind).stim_num(stim_num).failures(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).failures(:,2));      

        clear  tmp
        
     %nonspecific events
        tmp(:,:) = nonspecific_count_mat(1:2,stim_num,:);  
        event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count=tmp';       
        event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count,0,1);  
        event_evoked_stat.stim_num(stim_num).nonspecific_count_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count_m;
        event_evoked_stat.stim_num(stim_num).nonspecific_count_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count_std;
        %testing for normal distribution
        if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_nonspecific_count, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_nonspecific_count] = lillietest(diff(tmp));        
        end
        %paired ttest 
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_nonspecific_count, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_nonspecific_count]= ttest(event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count(:,2));
       if sum(~isnan(event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count(:,1)))>4
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_nonspecific_count, event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_nonspecific_count]= signrank(event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count(:,2));
       end
        clear  tmp
        
    %onset
        tmp(:,:) = onset_mat(1:2,stim_num,:);          
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onset=tmp';       
       event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset,0,1);  
        event_evoked_stat.stim_num(stim_num).onset_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_m;
        event_evoked_stat.stim_num(stim_num).onset_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_std;
        %testing for normal distribution
        if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_onset, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_onset] = lillietest(diff(tmp));        
           [event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_onset, event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_onset]= signrank(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).onset(:,2));
        end
        %paired ttest 
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_onset, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_onset]= ttest(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).onset(:,2));

        clear  tmp
        
        %onset value
        tmp(:,:) = onVal_mat(1:2,stim_num,:);  
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal=tmp';
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal,0,1);  
        event_evoked_stat.stim_num(stim_num).onVal_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal_m;
        event_evoked_stat.stim_num(stim_num).onVal_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal_std;
        %testing for normal distribution
        if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_onVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_onVal] = lillietest(diff(tmp));        
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_onVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_onVal]= signrank(event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal(:,2));
        end
        %paired ttest 
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_onVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_onVal]= ttest(event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal(:,2));

        clear  tmp
        
    %amplitude
        tmp(:,:) = amplitude_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude=tmp';
        event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude,0,1);  
        event_evoked_stat.stim_num(stim_num).amplitude_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude_m;
        event_evoked_stat.stim_num(stim_num).amplitude_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude_std;
        %testing for normal distribution
        if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_amplitude, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_amplitude] = lillietest(diff(tmp));        
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_amplitude, event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_amplitude]= signrank(event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude(:,2));
       end
        %paired ttest 
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_amplitude, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_amplitude]= ttest(event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude(:,2));

        clear  tmp

    %peak value
        tmp(:,:) = ampVal_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal=tmp';
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal,0,1);  
        event_evoked_stat.stim_num(stim_num).ampVal_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal_m;
        event_evoked_stat.stim_num(stim_num).ampVal_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal_std;
        %testing for normal distribution
        if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_ampVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_ampVal] = lillietest(diff(tmp));        
           [event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_ampVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_ampVal]= signrank(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal(:,2));
        end
        %paired ttest 
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_ampVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_ampVal]= ttest(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal(:,2));

        clear  tmp
        
    %amplitude latency
        tmp(:,:) = ampDel_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel=tmp';
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel,0,1);  
        event_evoked_stat.stim_num(stim_num).ampDel_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel_m;
        event_evoked_stat.stim_num(stim_num).ampDel_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel_std;
        %testing for normal distribution
        if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_ampDel, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_ampDel] = lillietest(diff(tmp));        
           [event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_ampDel, event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_ampDel]= signrank(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel(:,2));
        end
        %paired ttest 
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_ampDel, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_ampDel]= ttest(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel(:,2));

        clear  tmp
 
    %half width
        tmp(:,:) = halfWidth_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth=tmp';
       event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth,0,1);  
        event_evoked_stat.stim_num(stim_num).halfWidth_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth_m;
        event_evoked_stat.stim_num(stim_num).halfWidth_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth_std;
        %testing for normal distribution
        if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
            [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_halfWidth, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_halfWidth] = lillietest(diff(tmp));        
          [event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_halfWidth, event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_halfWidth]= signrank(event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth(:,2));
        end
        %paired ttest 
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_halfWidth, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_halfWidth]= ttest(event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth(:,2));

        clear  tmp
        
    %delay to 90%peak
        tmp(:,:) = peak90prcntTime_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).peak90prcntTime=tmp';
        event_evoked_stat.fileind(fileind).stim_num(stim_num).peak90prcntTime_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).peak90prcntTime,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).peak90prcntTime_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).peak90prcntTime,0,1);  
        event_evoked_stat.stim_num(stim_num).peak90prcntTime_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).peak90prcntTime_m;
        event_evoked_stat.stim_num(stim_num).peak90prcntTime_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).peak90prcntTime_std;

        clear tmp
        
    %delay to 10%peak
        tmp(:,:) = peak10prcntTime_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10prcntTime=tmp';
        event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10prcntTime_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10prcntTime,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10prcntTime_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10prcntTime,0,1);  
        event_evoked_stat.stim_num(stim_num).peak10prcntTime_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10prcntTime_m;
        event_evoked_stat.stim_num(stim_num).peak10prcntTime_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10prcntTime_std;

        clear tmp
        
    %rise-time from 10%-90% peak
        tmp(:,:) = peak10to90Time_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10to90Time=tmp';
        event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10to90Time_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10to90Time,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10to90Time_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10to90Time,0,1);  
        event_evoked_stat.stim_num(stim_num).peak10to90Time_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10to90Time_m;
        event_evoked_stat.stim_num(stim_num).peak10to90Time_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).peak10to90Time_std;

%        if any(stim_num==[1,5,11])
%                 cell_stat_amplitude.wilk_h{fileind,stim_num}=event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_amplitude;         
%                 cell_stat_amplitude.wilk_p{fileind,stim_num}=event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_amplitude;   
% %               cell_stat_amplitude.wilk_h{fileind,stim_num}=event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_amplitude;
% %                 cell_stat_amplitude.wilk_p{fileind,stim_num}=event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_amplitude;
%                 cell_stat_amplitude.diff{fileind,stim_num}=diff(event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude_m,1,2); %first derivative along the 2 dimension
%                if cell_stat_amplitude.wilk_h{fileind,stim_num}==1
%                     if cell_stat_amplitude.diff{fileind,stim_num}>0
%                          cell_stat_amplitude.effect{fileind,stim_num}=1; %'increase'
%                     else
%                         cell_stat_amplitude.effect{fileind,stim_num}=-1; %'decrease'
%                     end
%                else 
%                    cell_stat_amplitude.effect{fileind,stim_num}=0; %'same'
%                end
%                 cell_stat_onset.wilk_h{fileind,stim_num}=event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_h_onset;
%                 cell_stat_onset.wilk_p{fileind,stim_num}=event_evoked_stat.fileind(fileind).stim_num(stim_num).wilcoxon_p_onset;   
% %               cell_stat_amplitude.wilk_h{fileind,stim_num}=event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_h_onset;
% %                 cell_stat_amplitude.wilk_p{fileind,stim_num}=event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest_p_onset;
%                 cell_stat_onset.diff{fileind,stim_num}=diff(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_m,1,2); %first derivative along the 2 dimension
%                if cell_stat_onset.wilk_h{fileind,stim_num}==1
%                     if cell_stat_onset.diff{fileind,stim_num}>0
%                          cell_stat_onset.effect{fileind,stim_num}=1; %'increase'
%                     else
%                         cell_stat_onset.effect{fileind,stim_num}=-1; %'decrease'
%                     end
%                else 
%                    cell_stat_onset.effect{fileind,stim_num}=0; %'same'
%                end               
%        end
    end
    end
 %% population statistics
 clear tmp
 
     h_failures_m=[];     h_onset_m=[];     h_onset_std=[];     h_amplitude_m=[];
     h_amplitude_std=[];     h_onVal_m=[];     h_ampVal_m=[];     h_ampDel_m=[];     h_ampDel_std=[];
     h_halfWidth_m=[];    h_halfWidth_std=[];     h_adapt_amp=[];     h_nonspecific_count_m=[];
    
 for stim_num=1:11;

    %average number of failures per stim
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).failures_m;           
        event_evoked_stat.stim_num(stim_num).change_failures_m = (tmp(:,2)-tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).failures_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).failures_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_failures_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_failures_m,1);
        event_evoked_stat.stim_num(stim_num).change_failures_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_failures_m,0,1);
       change_failures_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_failures_m;
        %testing for normal distribution
        diff_failures_m= tmp(:,2)- tmp(:,1);
    if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
        [event_evoked_stat.stim_num(stim_num).lillietest_h_failures_m, event_evoked_stat.stim_num(stim_num).lillietest_p_failures_m] = lillietest(diff_failures_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_failures_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_failures_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_failures_m);
        %paired ttest 
         [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_failures_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_failures_m]= signrank(event_evoked_stat.stim_num(stim_num).change_failures_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_failures_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_failures_m]= signrank(event_evoked_stat.stim_num(stim_num).failures_m(:,1),event_evoked_stat.stim_num(stim_num).failures_m(:,2));
    end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_failures_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_failures_m]= ttest(event_evoked_stat.stim_num(stim_num).change_failures_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_failures_m, event_evoked_stat.stim_num(stim_num).ttest_p_failures_m]= ttest(event_evoked_stat.stim_num(stim_num).failures_m(:,1),event_evoked_stat.stim_num(stim_num).failures_m(:,2));
      
        clear  diff_failures_m tmp
        
        %average number of nonspecific events per stim
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).nonspecific_count_m;           
        event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m = (tmp(:,2)-tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).nonspecific_count_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).nonspecific_count_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m,1);
        event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m,0,1);
       change_nonspecific_count_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m;
        %testing for normal distribution
        diff_nonspecific_count_m= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
        [event_evoked_stat.stim_num(stim_num).lillietest_h_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).lillietest_p_nonspecific_count_m] = lillietest(diff_nonspecific_count_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_nonspecific_count_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m);
        %paired ttest 
         [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_nonspecific_count_m]= signrank(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m);
         [event_evoked_stat.stim_num(stim_num).wilcoxon_p_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_nonspecific_count_m]= signrank(event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,1),event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_nonspecific_count_m]= ttest(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).ttest_p_nonspecific_count_m]= ttest(event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,1),event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,2));
       
        clear  diff_nonspecific_count_m tmp
        
    % onset   
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).onset_m;           
        event_evoked_stat.stim_num(stim_num).change_onset_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).onset_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).onset_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_onset_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_onset_m,1);
        event_evoked_stat.stim_num(stim_num).change_onset_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_onset_m,0,1);
        change_onset_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_onset_m;
        %testing for normal distribution
        diff_onset_m= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs
        [event_evoked_stat.stim_num(stim_num).lillietest_h_onset_m, event_evoked_stat.stim_num(stim_num).lillietest_p_onset_m] = lillietest(diff_onset_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_onset_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_onset_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_onset_m);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_onset_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onset_m]= signrank(event_evoked_stat.stim_num(stim_num).change_onset_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_onset_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_onset_m]= signrank(event_evoked_stat.stim_num(stim_num).onset_m(:,1),event_evoked_stat.stim_num(stim_num).onset_m(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_onset_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_onset_m]= ttest(event_evoked_stat.stim_num(stim_num).change_onset_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_onset_m, event_evoked_stat.stim_num(stim_num).ttest_p_onset_m]= ttest(event_evoked_stat.stim_num(stim_num).onset_m(:,1),event_evoked_stat.stim_num(stim_num).onset_m(:,2));
        
        clear  diff_onset_m tmp
 
    % onset std
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).onset_std;           
        event_evoked_stat.stim_num(stim_num).change_onset_std = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).onset_std_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).onset_std_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_onset_std_m= nanmean(event_evoked_stat.stim_num(stim_num).change_onset_std,1);
        event_evoked_stat.stim_num(stim_num).change_onset_std_std= nanstd(event_evoked_stat.stim_num(stim_num).change_onset_std,0,1);
        change_onset_std_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_onset_std;
        %testing for normal distribution
        diff_onset_std= tmp(:,2)- tmp(:,1);
  if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs        
        [event_evoked_stat.stim_num(stim_num).lillietest_h_onset_std, event_evoked_stat.stim_num(stim_num).lillietest_p_onset_std] = lillietest(diff_onset_std);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_onset_std, event_evoked_stat.stim_num(stim_num).lillietest_p_change_onset_std] = lillietest(event_evoked_stat.stim_num(stim_num).change_onset_std);
        %paired ttest 
         [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_onset_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onset_std]= signrank(event_evoked_stat.stim_num(stim_num).change_onset_std);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_onset_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_onset_std]= signrank(event_evoked_stat.stim_num(stim_num).onset_std(:,1),event_evoked_stat.stim_num(stim_num).onset_std(:,2));
  end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_onset_std, event_evoked_stat.stim_num(stim_num).ttest_p_change_onset_std]= ttest(event_evoked_stat.stim_num(stim_num).change_onset_std);
        [event_evoked_stat.stim_num(stim_num).ttest_h_onset_std, event_evoked_stat.stim_num(stim_num).ttest_p_onset_std]= ttest(event_evoked_stat.stim_num(stim_num).onset_std(:,1),event_evoked_stat.stim_num(stim_num).onset_std(:,2));
       
        clear  diff_onset_std tmp
        
        % onset value
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).onVal_m;           
        event_evoked_stat.stim_num(stim_num).change_onVal_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).onVal_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).onVal_m_std= nanstd(tmp,0,1);
        event_evoked_stat.stim_num(stim_num).change_onVal_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_onVal_m,1);
        event_evoked_stat.stim_num(stim_num).change_onVal_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_onVal_m,0,1);
        change_onVal_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_onVal_m;
        %testing for normal distribution
        diff_onVal_m= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs        
        [event_evoked_stat.stim_num(stim_num).lillietest_h_onVal_m, event_evoked_stat.stim_num(stim_num).lillietest_p_onVal_m] = lillietest(diff_onVal_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_onVal_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_onVal_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_onVal_m);
        %paired ttest 
         [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_onVal_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onVal_m]= signrank(event_evoked_stat.stim_num(stim_num).change_onVal_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_onVal_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_onVal_m]= signrank(event_evoked_stat.stim_num(stim_num).onVal_m(:,1),event_evoked_stat.stim_num(stim_num).onVal_m(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_onVal_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_onVal_m]= ttest(event_evoked_stat.stim_num(stim_num).change_onVal_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_onVal_m, event_evoked_stat.stim_num(stim_num).ttest_p_onVal_m]= ttest(event_evoked_stat.stim_num(stim_num).onVal_m(:,1),event_evoked_stat.stim_num(stim_num).onVal_m(:,2));
       
        clear  diff_onVal_m tmp
        
     % amplitude mean
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).amplitude_m;           
        event_evoked_stat.stim_num(stim_num).change_amplitude_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).amplitude_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).amplitude_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_amplitude_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_amplitude_m,1);
        event_evoked_stat.stim_num(stim_num).change_amplitude_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_amplitude_m,0,1);
       change_amplitude_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_amplitude_m;
        %testing for normal distribution
        diff_amplitude_m= tmp(:,2)- tmp(:,1);
        if (nnz(diff_amplitude_m)-sum(sum(isnan(diff_amplitude_m))))>4;
            [event_evoked_stat.stim_num(stim_num).lillietest_h_amplitude_m, event_evoked_stat.stim_num(stim_num).lillietest_p_amplitude_m] = lillietest(diff_amplitude_m);
            [event_evoked_stat.stim_num(stim_num).lillietest_h_change_amplitude_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_amplitude_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_amplitude_m);
        %paired ttest
            [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_amplitude_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_amplitude_m]= signrank(event_evoked_stat.stim_num(stim_num).change_amplitude_m);
            [event_evoked_stat.stim_num(stim_num).wilcoxon_p_amplitude_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_amplitude_m]= signrank(event_evoked_stat.stim_num(stim_num).amplitude_m(:,1),event_evoked_stat.stim_num(stim_num).amplitude_m(:,2));
        end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_amplitude_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_amplitude_m]= ttest(event_evoked_stat.stim_num(stim_num).change_amplitude_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_amplitude_m, event_evoked_stat.stim_num(stim_num).ttest_p_amplitude_m]= ttest(event_evoked_stat.stim_num(stim_num).amplitude_m(:,1),event_evoked_stat.stim_num(stim_num).amplitude_m(:,2));

        clear  diff_amplitude_m tmp
        
    % amplitude std
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).amplitude_std;           
        event_evoked_stat.stim_num(stim_num).change_amplitude_std = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).amplitude_std_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).amplitude_std_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_amplitude_std_m= nanmean(event_evoked_stat.stim_num(stim_num).change_amplitude_std,1);
        event_evoked_stat.stim_num(stim_num).change_amplitude_std_std= nanstd(event_evoked_stat.stim_num(stim_num).change_amplitude_std,0,1);
        change_amplitude_std_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_amplitude_std;
        %testing for normal distribution
        diff_amplitude_std= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs        
        [event_evoked_stat.stim_num(stim_num).lillietest_h_amplitude_std, event_evoked_stat.stim_num(stim_num).lillietest_p_amplitude_std] = lillietest(diff_amplitude_std);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_amplitude_std, event_evoked_stat.stim_num(stim_num).lillietest_p_change_amplitude_std] = lillietest(event_evoked_stat.stim_num(stim_num).change_amplitude_std);
        %paired ttest
         [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_amplitude_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_amplitude_std]= signrank(event_evoked_stat.stim_num(stim_num).change_amplitude_std);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_amplitude_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_amplitude_std]= signrank(event_evoked_stat.stim_num(stim_num).amplitude_std(:,1),event_evoked_stat.stim_num(stim_num).amplitude_std(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_amplitude_std, event_evoked_stat.stim_num(stim_num).ttest_p_change_amplitude_std]= ttest(event_evoked_stat.stim_num(stim_num).change_amplitude_std);
        [event_evoked_stat.stim_num(stim_num).ttest_h_amplitude_std, event_evoked_stat.stim_num(stim_num).ttest_p_amplitude_std]= ttest(event_evoked_stat.stim_num(stim_num).amplitude_std(:,1),event_evoked_stat.stim_num(stim_num).amplitude_std(:,2));
       
        clear  diff_amplitude_std tmp
        
        % peak value
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).ampVal_m;           
        event_evoked_stat.stim_num(stim_num).change_ampVal_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).ampVal_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).ampVal_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_ampVal_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_ampVal_m,1);
        event_evoked_stat.stim_num(stim_num).change_ampVal_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_ampVal_m,0,1);
        change_ampVal_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_ampVal_m;
        %testing for normal distribution
        diff_ampVal_m= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs        
        [event_evoked_stat.stim_num(stim_num).lillietest_h_ampVal_m, event_evoked_stat.stim_num(stim_num).lillietest_p_ampVal_m] = lillietest(diff_ampVal_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampVal_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_ampVal_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_ampVal_m);
        %paired ttest 
          [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_ampVal_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampVal_m]= signrank(event_evoked_stat.stim_num(stim_num).change_ampVal_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_ampVal_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampVal_m]= signrank(event_evoked_stat.stim_num(stim_num).ampVal_m(:,1),event_evoked_stat.stim_num(stim_num).ampVal_m(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_ampVal_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_ampVal_m]= ttest(event_evoked_stat.stim_num(stim_num).change_ampVal_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_ampVal_m, event_evoked_stat.stim_num(stim_num).ttest_p_ampVal_m]= ttest(event_evoked_stat.stim_num(stim_num).ampVal_m(:,1),event_evoked_stat.stim_num(stim_num).ampVal_m(:,2));
      
        clear  diff_ampVal_m tmp
        
    % amplitude latency mean
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).ampDel_m;           
        event_evoked_stat.stim_num(stim_num).change_ampDel_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).ampDel_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).ampDel_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_ampDel_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_ampDel_m,1);
        event_evoked_stat.stim_num(stim_num).change_ampDel_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_ampDel_m,0,1);
        change_ampDel_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_ampDel_m;
        %testing for normal distribution
        diff_ampDel_m= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs        
        [event_evoked_stat.stim_num(stim_num).lillietest_h_ampDel_m, event_evoked_stat.stim_num(stim_num).lillietest_p_ampDel_m] = lillietest(diff_ampDel_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampDel_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_ampDel_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_ampDel_m);
        %paired ttest
          [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_ampDel_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampDel_m]= signrank(event_evoked_stat.stim_num(stim_num).change_ampDel_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_ampDel_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampDel_m]= signrank(event_evoked_stat.stim_num(stim_num).ampDel_m(:,1),event_evoked_stat.stim_num(stim_num).ampDel_m(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_ampDel_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_ampDel_m]= ttest(event_evoked_stat.stim_num(stim_num).change_ampDel_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_ampDel_m, event_evoked_stat.stim_num(stim_num).ttest_p_ampDel_m]= ttest(event_evoked_stat.stim_num(stim_num).ampDel_m(:,1),event_evoked_stat.stim_num(stim_num).ampDel_m(:,2));
      
        clear  diff_ampDel_m tmp
        
     % amplitude latency std
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).ampDel_std;           
        event_evoked_stat.stim_num(stim_num).change_ampDel_std = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).ampDel_std_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).ampDel_std_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_ampDel_std_m= nanmean(event_evoked_stat.stim_num(stim_num).change_ampDel_std,1);
        event_evoked_stat.stim_num(stim_num).change_ampDel_std_std= nanstd(event_evoked_stat.stim_num(stim_num).change_ampDel_std,0,1);
        change_ampDel_std_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_ampDel_std;
        %testing for normal distribution
        diff_ampDel_std= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs        
        [event_evoked_stat.stim_num(stim_num).lillietest_h_ampDel_std, event_evoked_stat.stim_num(stim_num).lillietest_p_ampDel_std] = lillietest(diff_ampDel_std);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampDel_std, event_evoked_stat.stim_num(stim_num).lillietest_p_change_ampDel_std] = lillietest(event_evoked_stat.stim_num(stim_num).change_ampDel_std);
        %paired ttest
         [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_ampDel_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampDel_std]= signrank(event_evoked_stat.stim_num(stim_num).change_ampDel_std);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_ampDel_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampDel_std]= signrank(event_evoked_stat.stim_num(stim_num).ampDel_std(:,1),event_evoked_stat.stim_num(stim_num).ampDel_std(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_ampDel_std, event_evoked_stat.stim_num(stim_num).ttest_p_change_ampDel_std]= ttest(event_evoked_stat.stim_num(stim_num).change_ampDel_std);
        [event_evoked_stat.stim_num(stim_num).ttest_h_ampDel_std, event_evoked_stat.stim_num(stim_num).ttest_p_ampDel_std]= ttest(event_evoked_stat.stim_num(stim_num).ampDel_std(:,1),event_evoked_stat.stim_num(stim_num).ampDel_std(:,2));
       
        clear  diff_ampDel_std tmp
        
    % half width mean
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).halfWidth_m;           
        event_evoked_stat.stim_num(stim_num).change_halfWidth_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).halfWidth_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).halfWidth_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_halfWidth_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_halfWidth_m,1);
        event_evoked_stat.stim_num(stim_num).change_halfWidth_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_halfWidth_m,0,1);
        change_halfWidth_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_halfWidth_m;
        %testing for normal distribution
        diff_halfWidth_m= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs        
        [event_evoked_stat.stim_num(stim_num).lillietest_h_halfWidth_m, event_evoked_stat.stim_num(stim_num).lillietest_p_halfWidth_m] = lillietest(diff_halfWidth_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_halfWidth_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_halfWidth_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_halfWidth_m);
%         paired ttest 
         [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_halfWidth_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_halfWidth_m]= signrank(event_evoked_stat.stim_num(stim_num).change_halfWidth_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_halfWidth_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_halfWidth_m]= signrank(event_evoked_stat.stim_num(stim_num).halfWidth_m(:,1),event_evoked_stat.stim_num(stim_num).halfWidth_m(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_halfWidth_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_halfWidth_m]= ttest(event_evoked_stat.stim_num(stim_num).change_halfWidth_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_halfWidth_m, event_evoked_stat.stim_num(stim_num).ttest_p_halfWidth_m]= ttest(event_evoked_stat.stim_num(stim_num).halfWidth_m(:,1),event_evoked_stat.stim_num(stim_num).halfWidth_m(:,2));
       
        clear  diff_halfWidth_m tmp
        
    % half width std
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).halfWidth_std;           
        event_evoked_stat.stim_num(stim_num).change_halfWidth_std = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)) ;
        event_evoked_stat.stim_num(stim_num).halfWidth_std_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).halfWidth_std_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_halfWidth_std_m= nanmean(event_evoked_stat.stim_num(stim_num).change_halfWidth_std,1);
        event_evoked_stat.stim_num(stim_num).change_halfWidth_std_std= nanstd(event_evoked_stat.stim_num(stim_num).change_halfWidth_std,0,1);
        change_halfWidth_std_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_halfWidth_std;
        %testing for normal distribution
        diff_halfWidth_std= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs        
        [event_evoked_stat.stim_num(stim_num).lillietest_h_halfWidth_std, event_evoked_stat.stim_num(stim_num).lillietest_p_halfWidth_std] = lillietest(diff_halfWidth_std);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_halfWidth_std, event_evoked_stat.stim_num(stim_num).lillietest_p_change_halfWidth_std] = lillietest(event_evoked_stat.stim_num(stim_num).change_halfWidth_std);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_halfWidth_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_halfWidth_std]= signrank(event_evoked_stat.stim_num(stim_num).change_halfWidth_std);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_halfWidth_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_halfWidth_std]= signrank(event_evoked_stat.stim_num(stim_num).halfWidth_std(:,1),event_evoked_stat.stim_num(stim_num).halfWidth_std(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_halfWidth_std, event_evoked_stat.stim_num(stim_num).ttest_p_change_halfWidth_std]= ttest(event_evoked_stat.stim_num(stim_num).change_halfWidth_std);
        [event_evoked_stat.stim_num(stim_num).ttest_h_halfWidth_std, event_evoked_stat.stim_num(stim_num).ttest_p_halfWidth_std]= ttest(event_evoked_stat.stim_num(stim_num).halfWidth_std(:,1),event_evoked_stat.stim_num(stim_num).halfWidth_std(:,2));
            
        clear  diff_halfWidth_std tmp

%delay from stimulus onset to 10% peak
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).peak10prcntTime_m;
        event_evoked_stat.stim_num(stim_num).peak10prcntTime_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).peak10prcntTime_m_std= nanstd(tmp,0,1);

%delay from stimulus onset to 90% peak
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).peak90prcntTime_m;
        event_evoked_stat.stim_num(stim_num).peak90prcntTime_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).peak90prcntTime_m_std= nanstd(tmp,0,1);
        
%rise-time from 10% to 90% peak
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).peak10to90Time_m;
        event_evoked_stat.stim_num(stim_num).peak10to90Time_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).peak10to90Time_m_std= nanstd(tmp,0,1);
        
        %adaptation amplitude ratio
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).adapt_amp_m;           
        event_evoked_stat.stim_num(stim_num).change_adapt_amp = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).adapt_amp_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).adapt_amp_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_adapt_amp_m= nanmean(event_evoked_stat.stim_num(stim_num).change_adapt_amp,1);
        event_evoked_stat.stim_num(stim_num).change_adapt_amp_std= nanstd(event_evoked_stat.stim_num(stim_num).change_adapt_amp,0,1);
        change_adapt_amp_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_adapt_amp;
        %testing for normal distribution
        diff_adapt_amp= tmp(:,2)- tmp(:,1);
      if (nnz(diff(tmp))-sum(sum(isnan(diff(tmp)))))>4; %if there are at least 5 values that are non-zero and not NaNs        
        [event_evoked_stat.stim_num(stim_num).lillietest_h_adapt_amp, event_evoked_stat.stim_num(stim_num).lillietest_p_adapt_amp] = lillietest(diff_adapt_amp);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_adapt_amp, event_evoked_stat.stim_num(stim_num).lillietest_p_change_adapt_amp] = lillietest(event_evoked_stat.stim_num(stim_num).change_adapt_amp);
        %paired ttest 
         [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_adapt_amp, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_adapt_amp]= signrank(event_evoked_stat.stim_num(stim_num).change_adapt_amp);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_adapt_amp, event_evoked_stat.stim_num(stim_num).wilcoxon_h_adapt_amp]= signrank(event_evoked_stat.stim_num(stim_num).adapt_amp_m(:,1),event_evoked_stat.stim_num(stim_num).adapt_amp_m(:,2));
      end
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_adapt_amp, event_evoked_stat.stim_num(stim_num).ttest_p_change_adapt_amp]= ttest(event_evoked_stat.stim_num(stim_num).change_adapt_amp);
        [event_evoked_stat.stim_num(stim_num).ttest_h_adapt_amp, event_evoked_stat.stim_num(stim_num).ttest_p_adapt_amp]= ttest(event_evoked_stat.stim_num(stim_num).adapt_amp_m(:,1),event_evoked_stat.stim_num(stim_num).adapt_amp_m(:,2));
       
        clear  diff_adapt_amp tmp
%taking the p-values for non-normalized data:        
%         if event_evoked_stat.stim_num(stim_num).lillietest_h_failures_m==0      
%             h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_failures_m; 
%   else h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_failures_m; 
%         end
%      if event_evoked_stat.stim_num(stim_num).lillietest_h_onset_m==0      
%             h_onset_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_onset_m; 
%   else h_onset_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_onset_m; 
%      end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_onset_std==0      
%             h_onset_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_onset_std; 
%   else h_onset_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_onset_std; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_amplitude_m==0      
%             h_amplitude_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_amplitude_m; 
%   else h_amplitude_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_amplitude_m; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_amplitude_std==0      
%             h_amplitude_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_amplitude_std; 
%   else h_amplitude_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_amplitude_std; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_onVal_m==0      
%             h_onVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_onVal_m; 
%   else h_onVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_onVal_m; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_ampVal_m==0      
%             h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_ampVal_m; 
%   else h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampVal_m; 
%    end
%  if event_evoked_stat.stim_num(stim_num).lillietest_h_ampDel_m==0      
%             h_ampDel_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_ampDel_m; 
%   else h_ampDel_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampDel_m; 
%  end
%     if event_evoked_stat.stim_num(stim_num).lillietest_h_ampDel_std==0      
%             h_ampDel_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_ampDel_std; 
%   else h_ampDel_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampDel_std; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_halfWidth_m==0      
%             h_halfWidth_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_halfWidth_m; 
%   else h_halfWidth_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_halfWidth_m; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_halfWidth_std==0      
%             h_halfWidth_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_halfWidth_std; 
%   else h_halfWidth_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_halfWidth_std; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_adapt_amp==0      
%             h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_adapt_amp; 
%   else h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_adapt_amp; 
%    end
%  if event_evoked_stat.stim_num(stim_num).lillietest_h_nonspecific_count_m==0      
%             h_nonspecific_count_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_nonspecific_count_m; 
%   else h_nonspecific_count_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_nonspecific_count_m; 
%    end

%taking the p-values for normalized data:    
if isfield(event_evoked_stat,'lillietest_h_change_failures_m');
        if event_evoked_stat.stim_num(stim_num).lillietest_h_change_failures_m==0      
            h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_failures_m; 
  else h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_failures_m; 
        end
end
     if event_evoked_stat.stim_num(stim_num).lillietest_h_change_onset_m==0      
            h_onset_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_onset_m; 
  else h_onset_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onset_m; 
     end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_change_onset_std==0      
            h_onset_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_onset_std; 
  else h_onset_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onset_std; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_change_amplitude_m==0      
            h_amplitude_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_amplitude_m; 
  else h_amplitude_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_amplitude_m; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_change_amplitude_std==0      
            h_amplitude_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_amplitude_std; 
  else h_amplitude_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_amplitude_std; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_change_onVal_m==0      
            h_onVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_onVal_m; 
  else h_onVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onVal_m; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampVal_m==0      
            h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_ampVal_m; 
  else h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampVal_m; 
   end
 if event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampDel_m==0      
            h_ampDel_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_ampDel_m; 
  else h_ampDel_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampDel_m; 
 end
    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampDel_std==0      
            h_ampDel_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_ampDel_std; 
  else h_ampDel_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampDel_std; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_change_halfWidth_m==0      
            h_halfWidth_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_halfWidth_m; 
  else h_halfWidth_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_halfWidth_m; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_change_halfWidth_std==0      
            h_halfWidth_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_halfWidth_std; 
  else h_halfWidth_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_halfWidth_std; 
   end
   if isfield(event_evoked_stat,'lillietest_h_change_adapt_amp')
       if event_evoked_stat.stim_num(stim_num).lillietest_h_change_adapt_amp==0      
                h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_adapt_amp; 
        else h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_adapt_amp; 
       end
   end
 if event_evoked_stat.stim_num(stim_num).lillietest_h_change_nonspecific_count_m==0      
            h_nonspecific_count_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_nonspecific_count_m; 
  else h_nonspecific_count_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_nonspecific_count_m; 
   end
    end
    %getting the x-axis (stimulus number) and y-axis for the asterisks of significance
     s_failures_m=find(h_failures_m==1);
     s_onset_m=find(h_onset_m==1);
     s_onset_std=find(h_onset_std==1); 
     s_amplitude_m=find(h_amplitude_m==1);
     s_amplitude_std=find(h_amplitude_std==1);
     s_onVal_m=find(h_onVal_m==1);
     s_ampVal_m=find(h_ampVal_m==1);
     s_ampDel_m=find(h_ampDel_m==1);
     s_ampDel_std=find(h_ampDel_std==1);
     s_halfWidth_m=find(h_halfWidth_m==1);
     s_halfWidth_std=find(h_halfWidth_std==1);
     s_adapt_amp=find(h_adapt_amp==1);
     s_nonspecific_count_m=find(h_nonspecific_count_m==1);
 %% Plot parameters along the train stim - version 3:bars+error bars of the mean values
 close all
     
        barwidth1=0.3;
    for stim_num=1:11      
     failures_m_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).failures_m(:,1);
     failures_m_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).failures_m(:,2);
     nonspecific_count_m_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,1);
     nonspecific_count_m_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,2);
     onset_m_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).onset_m(:,1);
     onset_m_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).onset_m(:,2);
     onset_std_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).onset_std(:,1);
     onset_std_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).onset_std(:,2);
     amplitude_m_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).amplitude_m(:,1);
     amplitude_m_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).amplitude_m(:,2);
     amplitude_std_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).amplitude_std(:,1);
     amplitude_std_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).amplitude_std(:,2);
     ampVal_m_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).ampVal_m(:,1);
     ampVal_m_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).ampVal_m(:,2);
     onVal_m_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).onVal_m(:,1);
     onVal_m_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).onVal_m(:,2);
     ampDel_m_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).ampDel_m(:,1);
     ampDel_m_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).ampDel_m(:,2);
     ampDel_std_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).ampDel_std(:,1);
     ampDel_std_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).ampDel_std(:,2);
     halfWidth_m_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).halfWidth_m(:,1);
     halfWidth_m_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).halfWidth_m(:,2);
     halfWidth_std_mat{1}(:,stim_num)= event_evoked_stat.stim_num(stim_num).halfWidth_std(:,1);
     halfWidth_std_mat{2}(:,stim_num)= event_evoked_stat.stim_num(stim_num).halfWidth_std(:,2);
     peak10prcntTime_m_mat{1}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak10prcntTime_m(:,1);
     peak10prcntTime_m_mat{2}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak10prcntTime_m(:,2);
     peak10prcntTime_std_mat{1}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak10prcntTime_std(:,1);
     peak10prcntTime_std_mat{2}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak10prcntTime_std(:,2);
     peak90prcntTime_m_mat{1}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak90prcntTime_m(:,1);
     peak90prcntTime_m_mat{2}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak90prcntTime_m(:,2);
     peak90prcntTime_std_mat{1}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak90prcntTime_std(:,1);
     peak90prcntTime_std_mat{2}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak90prcntTime_std(:,2);
     peak10to90Time_m_mat{1}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak10to90Time_m(:,1);
     peak10to90Time_m_mat{2}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak10to90Time_m(:,2);   
     peak10to90Time_std_mat{1}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak10to90Time_std(:,1);
     peak10to90Time_std_mat{2}(:,stim_num)=event_evoked_stat.stim_num(stim_num).peak10to90Time_std(:,2);   
    end     
%% repeated measures ANOVA for amplitude over the train stim.
[stat_failures]=fn_sensory_train_response_2_way_rmANOVA(failures_m_mat{1}(:,:), failures_m_mat{2}(:,:));
[stat_nonspecific_count]=fn_sensory_train_response_2_way_rmANOVA(nonspecific_count_m_mat{1}(:,:), nonspecific_count_m_mat{2}(:,:));
[stat_onset]=fn_sensory_train_response_2_way_rmANOVA(onset_m_mat{1}(:,:), onset_m_mat{2}(:,:));
[stat_onset_std]=fn_sensory_train_response_2_way_rmANOVA(onset_std_mat{1}(:,:), onset_std_mat{2}(:,:));
[stat_amp]=fn_sensory_train_response_2_way_rmANOVA(amplitude_m_mat{1}(:,:), amplitude_m_mat{2}(:,:));
[stat_amp_std]=fn_sensory_train_response_2_way_rmANOVA(amplitude_std_mat{1}(:,:), amplitude_std_mat{2}(:,:));
[stat_ampVal]=fn_sensory_train_response_2_way_rmANOVA(ampVal_m_mat{1}(:,:), ampVal_m_mat{2}(:,:));
[stat_onVal]=fn_sensory_train_response_2_way_rmANOVA(onVal_m_mat{1}(:,:), onVal_m_mat{2}(:,:));
[stat_ampDel]=fn_sensory_train_response_2_way_rmANOVA(ampDel_m_mat{1}(:,:), ampDel_m_mat{2}(:,:));
[stat_ampDel_std]=fn_sensory_train_response_2_way_rmANOVA(ampDel_std_mat{1}(:,:), ampDel_std_mat{2}(:,:));
[stat_halfWidth]=fn_sensory_train_response_2_way_rmANOVA(halfWidth_m_mat{1}(:,:), halfWidth_m_mat{2}(:,:));
[stat_halfWidth_std]=fn_sensory_train_response_2_way_rmANOVA(halfWidth_std_mat{1}(:,:), halfWidth_std_mat{2}(:,:));
[stat_peak10prcntTime]=fn_sensory_train_response_2_way_rmANOVA(peak10prcntTime_m_mat{1}(:,:), peak10prcntTime_m_mat{2}(:,:));
[stat_peak90prcntTime]=fn_sensory_train_response_2_way_rmANOVA(peak90prcntTime_m_mat{1}(:,:), peak90prcntTime_m_mat{2}(:,:));
[stat_peak10to90Time]=fn_sensory_train_response_2_way_rmANOVA(peak10to90Time_m_mat{1}(:,:), peak10to90Time_m_mat{2}(:,:));
[stat_peak10prcntTime_std]=fn_sensory_train_response_2_way_rmANOVA(peak10prcntTime_std_mat{1}(:,:), peak10prcntTime_std_mat{2}(:,:));
[stat_peak90prcntTime_std]=fn_sensory_train_response_2_way_rmANOVA(peak90prcntTime_std_mat{1}(:,:), peak90prcntTime_std_mat{2}(:,:));
[stat_peak10to90Time_std]=fn_sensory_train_response_2_way_rmANOVA(peak10to90Time_std_mat{1}(:,:), peak10to90Time_std_mat{2}(:,:));
%% Plots with p-values from repeated-measures ANOVA
color_table_lines = rand(length(files_to_analyze),3);
%% failures
      h1=figure; 
      stimz=size(failures_m_mat{2}(:,:),2);
for stim=1:stimz;
        if stat_failures.p_stim{stim,1}>0.05 
            asterisk_failures{stim,1}='n.s.';
        else if stat_failures.p_stim{stim,1}<0.05 && stat_failures.p_stim{stim,1}>0.01
            asterisk_failures{stim,1}='*';
            else if stat_failures.p_stim{stim,1}<0.01 && stat_failures.p_stim{stim,1}>0.001
                    asterisk_failures{stim,1}='**';
            else if stat_failures.p_stim{stim,1}<0.001
                     asterisk_failures{stim,1}='***';
                end
                end
            end
        end
end

hold on      
        errbar_h=errorbar([1:size(failures_m_mat{1}(:,:),2)]-0.3,nanmean(failures_m_mat{1}(:,:),1),zeros(1,size(failures_m_mat{1}(:,:),2)), nanstd(failures_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(failures_m_mat{2}(:,:),2)],nanmean(failures_m_mat{2}(:,:),1),zeros(1,size(failures_m_mat{2}(:,:),2)), nanstd(failures_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(failures_m_mat{1}(:,:),2)]-0.3, nanmean(failures_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(failures_m_mat{2}(:,:),2)], nanmean(failures_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))

        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-0.05;
        for stim=1:size(failures_m_mat{2}(:,:),2);
            if strcmp('n.s.',asterisk_failures{stim,1})==0          
             text(stim,my,asterisk_failures{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(failures_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Failures rate' ,'FontSize', 16);    
        title(['Mean Rate of Failures,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);      
    %% nonspecific_count
      h2=figure; 
for stim=1:size(nonspecific_count_m_mat{2}(:,:),2);
        if stat_nonspecific_count.p_stim{stim,1}>0.05 
            asterisk_nonspecific_count{stim,1}='n.s.';
        else if stat_nonspecific_count.p_stim{stim,1}<0.05 && stat_nonspecific_count.p_stim{stim,1}>0.01
            asterisk_nonspecific_count{stim,1}='*';
            else if stat_nonspecific_count.p_stim{stim,1}<0.01 && stat_nonspecific_count.p_stim{stim,1}>0.001
                    asterisk_nonspecific_count{stim,1}='**';
            else if stat_nonspecific_count.p_stim{stim,1}<0.001
                     asterisk_nonspecific_count{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(nonspecific_count_m_mat{1}(:,:),2)]-0.3,nanmean(nonspecific_count_m_mat{1}(:,:),1),zeros(1,size(nonspecific_count_m_mat{1}(:,:),2)), nanstd(nonspecific_count_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(nonspecific_count_m_mat{2}(:,:),2)],nanmean(nonspecific_count_m_mat{2}(:,:),1),zeros(1,size(nonspecific_count_m_mat{2}(:,:),2)), nanstd(nonspecific_count_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(nonspecific_count_m_mat{1}(:,:),2)]-0.3, nanmean(nonspecific_count_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(nonspecific_count_m_mat{2}(:,:),2)], nanmean(nonspecific_count_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-0.05;
        for stim=1:size(nonspecific_count_m_mat{2}(:,:),2);
            if strcmp('n.s.',asterisk_nonspecific_count{stim,1})==0          
             text(stim,my,asterisk_nonspecific_count{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end 
        hold off
        set(gca,'xlim',[0,size(nonspecific_count_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Non-specific response rate' ,'FontSize', 16);    
        title(['Mean Rate of Non-specific responses,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
         %% onset_m
      h3=figure; 
      stimz=size(onset_m_mat{2}(:,:),2);
for stim=1:stimz
        if stat_onset.p_stim{stim,1}>0.05 
            asterisk_onset{stim,1}='n.s.';
        else if stat_onset.p_stim{stim,1}<0.05 && stat_onset.p_stim{stim,1}>0.01
            asterisk_onset{stim,1}='*';
            else if stat_onset.p_stim{stim,1}<0.01 && stat_onset.p_stim{stim,1}>0.001
                    asterisk_onset{stim,1}='**';
            else if stat_onset.p_stim{stim,1}<0.001
                     asterisk_onset{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(onset_m_mat{1}(:,:),2)]-0.3,nanmean(onset_m_mat{1}(:,:),1),zeros(1,size(onset_m_mat{1}(:,:),2)), nanstd(onset_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(onset_m_mat{2}(:,:),2)],nanmean(onset_m_mat{2}(:,:),1),zeros(1,size(onset_m_mat{2}(:,:),2)), nanstd(onset_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(onset_m_mat{1}(:,:),2)]-0.3, nanmean(onset_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(onset_m_mat{2}(:,:),2)], nanmean(onset_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
%for line pair-plots:
if paired_plot_flag==1; 
          for stim=1:stimz
                F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).onset_m;
                F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
                F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
           end
                for cell=1:length(files_to_analyze)
                    for stim=1:stimz
                        if isnan(F1_Y{stim}(cell,:))
                        else
                             plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                           pause
                        end
                    end
                end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-1;
        for stim=1:stimz
            if strcmp('n.s.',asterisk_onset{stim,1})==0          
             text(stim,my,asterisk_onset{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(onset_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Onset Latency [mS]' ,'FontSize', 16);    
        title(['Mean Onset Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        %% onset std
        F1_Y=[]; F1_X=[];
      h4=figure; 
      stimz=size(onset_std_mat{2}(:,:),2);
for stim=1:stimz;
        if stat_onset_std.p_stim{stim,1}>0.05 
            asterisk_onset_std{stim,1}='n.s.';
        else if stat_onset_std.p_stim{stim,1}<0.05 && stat_onset_std.p_stim{stim,1}>0.01
            asterisk_onset_std{stim,1}='*';
            else if stat_onset_std.p_stim{stim,1}<0.01 && stat_onset_std.p_stim{stim,1}>0.001
                    asterisk_onset_std{stim,1}='**';
            else if stat_onset_std.p_stim{stim,1}<0.001
                     asterisk_onset_std{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(onset_std_mat{1}(:,:),2)]-0.3,nanmean(onset_std_mat{1}(:,:),1),zeros(1,size(onset_std_mat{1}(:,:),2)), nanstd(onset_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(onset_std_mat{2}(:,:),2)],nanmean(onset_std_mat{2}(:,:),1),zeros(1,size(onset_std_mat{2}(:,:),2)), nanstd(onset_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(onset_std_mat{1}(:,:),2)]-0.3, nanmean(onset_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(onset_std_mat{2}(:,:),2)], nanmean(onset_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
   %for line pair-plots:
if paired_plot_flag==1; 
   for stim=1:stimz
        F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).onset_std;
        F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
        F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
   end
        for cell=1:length(files_to_analyze)
            for stim=1:stimz
                 if isnan(F1_Y{stim}(cell,:))
                 else
                   plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                   pause
                 end
            end
        end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-1;
        for stim=1:stimz;
            if strcmp('n.s.',asterisk_onset_std{stim,1})==0          
             text(stim,my,asterisk_onset_std{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(onset_std_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Onset Latency STD [mS]' ,'FontSize', 16);    
        title(['Mean Response Jitter ,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        %%  amplitude 
      h5=figure; 
      stimz=size(amplitude_m_mat{2}(:,:),2);
for stim=1:stimz
        if stat_amp.p_stim{stim,1}>0.05 
            asterisk_amp{stim,1}='n.s.';
        else if stat_amp.p_stim{stim,1}<0.05 && stat_amp.p_stim{stim,1}>0.01
            asterisk_amp{stim,1}='*';
            else if stat_amp.p_stim{stim,1}<0.01 && stat_amp.p_stim{stim,1}>0.001
                    asterisk_amp{stim,1}='**';
            else if stat_amp.p_stim{stim,1}<0.001
                     asterisk_amp{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(amplitude_m_mat{1}(:,:),2)]-0.3,nanmean(amplitude_m_mat{1}(:,:),1),zeros(1,size(amplitude_m_mat{1}(:,:),2)), nanstd(amplitude_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(amplitude_m_mat{2}(:,:),2)],nanmean(amplitude_m_mat{2}(:,:),1),zeros(1,size(amplitude_m_mat{2}(:,:),2)), nanstd(amplitude_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(amplitude_m_mat{1}(:,:),2)]-0.3, nanmean(amplitude_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(amplitude_m_mat{2}(:,:),2)], nanmean(amplitude_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
%for line pair-plots:
if paired_plot_flag==1; 
   for stim=1:stimz
        F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).amplitude_m;
        F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
        F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
   end
        for cell=1:length(files_to_analyze)
            for stim=1:stimz
                 if isnan(F1_Y{stim}(cell,:))
                    else
                        plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                         pause
                 end
            end
        end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-2;
        for stim=1:stimz
            if strcmp('n.s.',asterisk_amp{stim,1})==0          
             text(stim,my,asterisk_amp{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(amplitude_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
       ylabel(['Amplitude [', y_ax_units{1},']'],'FontSize', 16);    
        title(['Mean Response Amplitude, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);     
        %%
        data_vec=[];
        nbin=20;
        for stim=1:stimz
           data_vec=[data_vec,event_evoked_stat.stim_num(stim).change_amplitude_m'];
        end
         Fig1=figure;
            hold on
            [ncounts,nbins]=hist(data_vec,nbin);
%             ncounts(2,:)=ncounts(2,:)./length(data_vec(:,2));
            b1=bar(nbins,ncounts);
            
%             obj = gmdistribution.fit(data_vec(:,1),2);
%             set(h2,'FaceColor',[0 0 1],'EdgeColor','w')
%             set(h1,'FaceColor',[0 0 0],'EdgeColor','w','faceAlpha', 0.3)
            ylim1=get(gca,'ylim');
       hold off
      %%  amplitude std
      h6=figure; 
      stimz=size(amplitude_std_mat{2}(:,:),2);
for stim=1:stimz
        if stat_amp_std.p_stim{stim,1}>0.05 
            asterisk_amp_std{stim,1}='n.s.';
        else if stat_amp_std.p_stim{stim,1}<0.05 && stat_amp_std.p_stim{stim,1}>0.01
            asterisk_amp_std{stim,1}='*';
            else if stat_amp_std.p_stim{stim,1}<0.01 && stat_amp_std.p_stim{stim,1}>0.001
                    asterisk_amp_std{stim,1}='**';
            else if stat_amp_std.p_stim{stim,1}<0.001
                     asterisk_amp_std{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(amplitude_std_mat{1}(:,:),2)]-0.3,nanmean(amplitude_std_mat{1}(:,:),1),zeros(1,size(amplitude_std_mat{1}(:,:),2)), nanstd(amplitude_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(amplitude_std_mat{2}(:,:),2)],nanmean(amplitude_std_mat{2}(:,:),1),zeros(1,size(amplitude_std_mat{2}(:,:),2)), nanstd(amplitude_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(amplitude_std_mat{1}(:,:),2)]-0.3, nanmean(amplitude_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(amplitude_std_mat{2}(:,:),2)], nanmean(amplitude_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
    %for line pair-plots:
if paired_plot_flag==1; 
   for stim=1:stimz
        F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).amplitude_std;
        F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
        F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
   end
        for cell=1:length(files_to_analyze)
            for stim=1:stimz
                 if isnan(F1_Y{stim}(cell,:))
                    else
                        plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                         pause
                 end
            end
        end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-3;
        for stim=1:size(amplitude_std_mat{2}(:,:),2);
            if strcmp('n.s.',asterisk_amp_std{stim,1})==0          
             text(stim,my,asterisk_amp_std{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(amplitude_std_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel(['Amplitude STD [', y_ax_units{1},']'] ,'FontSize', 16);    
        title(['Mean Response Amplitude Jitter,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
         %%  amplitude peak value
         stimz=size(ampVal_m_mat{2}(:,:),2);
      h7=figure; 
for stim=1:stimz
        if stat_ampVal.p_stim{stim,1}>0.05 
            asterisk_ampVal{stim,1}='n.s.';
        else if stat_ampVal.p_stim{stim,1}<0.05 && stat_ampVal.p_stim{stim,1}>0.01
            asterisk_ampVal{stim,1}='*';
            else if stat_ampVal.p_stim{stim,1}<0.01 && stat_ampVal.p_stim{stim,1}>0.001
                    asterisk_ampVal{stim,1}='**';
            else if stat_ampVal.p_stim{stim,1}<0.001
                     asterisk_ampVal{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(ampVal_m_mat{1}(:,:),2)]-0.3,nanmean(ampVal_m_mat{1}(:,:),1),nanstd(ampVal_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k',zeros(1,size(ampVal_m_mat{1}(:,:),2)), 
        errbar_h=errorbar([1:size(ampVal_m_mat{2}(:,:),2)],nanmean(ampVal_m_mat{2}(:,:),1),nanstd(ampVal_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k',zeros(1,size(ampVal_m_mat{2}(:,:),2)), 
        bar([1:size(ampVal_m_mat{1}(:,:),2)]-0.3, nanmean(ampVal_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(ampVal_m_mat{2}(:,:),2)], nanmean(ampVal_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
    %for line pair-plots:
if paired_plot_flag==1; 
    for stim=1:stimz
        F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).ampVal_m;
        F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
        F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
   end
        for cell=1:length(files_to_analyze)
            for stim=1:stimz
                 if isnan(F1_Y{stim}(cell,:))
                    else
                        plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                         pause
                 end
            end
        end
end
        axis tight
        ylim_data=[get(gca,'ylim')]';
 if ylim_data(2)<=0
    ylim_data=[ylim_data(1)-5 ylim_data(1)+30]; %get(gca,'ylim'); %[0 8.5];
    my=min(ylim_data)+2;
else
    my=max(ylim_data)-5;
    ylim_data(1)=0;
end
        for stim=1:size(ampVal_m_mat{2}(:,:),2);
            if strcmp('n.s.',asterisk_ampVal{stim,1})==0          
             text(stim,my,asterisk_ampVal{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(ampVal_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11],'ylim',ylim_data)
        xlabel('Stim. serial number' ,'FontSize', 16);
       ylabel(['Peak Value [', y_ax_units{1},']'] ,'FontSize', 16);    
        title(['Mean Response Peak Value,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);      
                %%  Response Onset Value
                stimz=size(onVal_m_mat{2}(:,:),2);
      h8=figure; 
for stim=1:stimz
        if stat_onVal.p_stim{stim,1}>0.05 
            asterisk_onVal{stim,1}='n.s.';
        else if stat_onVal.p_stim{stim,1}<0.05 && stat_onVal.p_stim{stim,1}>0.01
            asterisk_onVal{stim,1}='*';
            else if stat_onVal.p_stim{stim,1}<0.01 && stat_onVal.p_stim{stim,1}>0.001
                    asterisk_onVal{stim,1}='**';
            else if stat_onVal.p_stim{stim,1}<0.001
                     asterisk_onVal{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(onVal_m_mat{1}(:,:),2)]-0.3,nanmean(onVal_m_mat{1}(:,:),1),nanstd(onVal_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k',zeros(1,size(onVal_m_mat{1}(:,:),2)), 
        errbar_h=errorbar([1:size(onVal_m_mat{2}(:,:),2)],nanmean(onVal_m_mat{2}(:,:),1),nanstd(onVal_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k',zeros(1,size(onVal_m_mat{2}(:,:),2)), 
        bar([1:size(onVal_m_mat{1}(:,:),2)]-0.3, nanmean(onVal_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(onVal_m_mat{2}(:,:),2)], nanmean(onVal_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
       
    %for line pair-plots:
if paired_plot_flag==1;     
   for stim=1:stimz
        F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).onVal_m;
        F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
        F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
   end
        for cell=1:length(files_to_analyze)
            for stim=1:stimz
                 if isnan(F1_Y{stim}(cell,:))
                    else
                        plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                         pause
                 end
            end
        end
end
        axis tight
        ylim_data=[get(gca,'ylim')]';
 if ylim_data(2)<=0
    ylim_data=[ylim_data(1)-5 ylim_data(1)+25]; %get(gca,'ylim'); %[0 8.5];
    my=min(ylim_data)+2;
else
    my=max(ylim_data)-5;
    ylim_data(1)=0;
end
        for stim=1:size(onVal_m_mat{2}(:,:),2);
            if strcmp('n.s.',asterisk_amp{stim,1})==0          
             text(stim,my,asterisk_amp{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(onVal_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11],'ylim',ylim_data)
        xlabel('Stim. serial number' ,'FontSize', 16);
       ylabel(['Onset Value [', y_ax_units{1},']'],'FontSize', 16);    
        title(['Mean Response Onset Value,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
                %%  Response Peak Latency
                stimz=size(ampDel_m_mat{2}(:,:),2);
      h9=figure; 
for stim=1:stimz
        if stat_ampDel.p_stim{stim,1}>0.05 
            asterisk_ampDel{stim,1}='n.s.';
        else if stat_ampDel.p_stim{stim,1}<0.05 && stat_ampDel.p_stim{stim,1}>0.01
            asterisk_ampDel{stim,1}='*';
            else if stat_ampDel.p_stim{stim,1}<0.01 && stat_ampDel.p_stim{stim,1}>0.001
                    asterisk_ampDel{stim,1}='**';
            else if stat_ampDel.p_stim{stim,1}<0.001
                     asterisk_ampDel{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(ampDel_m_mat{1}(:,:),2)]-0.3,nanmean(ampDel_m_mat{1}(:,:),1),zeros(1,size(ampDel_m_mat{1}(:,:),2)), nanstd(ampDel_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(ampDel_m_mat{2}(:,:),2)],nanmean(ampDel_m_mat{2}(:,:),1),zeros(1,size(ampDel_m_mat{2}(:,:),2)), nanstd(ampDel_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(ampDel_m_mat{1}(:,:),2)]-0.3, nanmean(ampDel_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(ampDel_m_mat{2}(:,:),2)], nanmean(ampDel_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
    %for line pair-plots:
if paired_plot_flag==1;        
   for stim=1:stimz
        F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).ampDel_m;
        F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
        F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
   end
        for cell=1:length(files_to_analyze)
            for stim=1:stimz
                 if isnan(F1_Y{stim}(cell,:))
                    else
                        plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                         pause
                 end
            end
        end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-5;
        for stim=1:size(ampDel_m_mat{2}(:,:),2);
            if strcmp('n.s.',asterisk_ampDel{stim,1})==0          
             text(stim,my,asterisk_ampDel{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(ampDel_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
       ylabel('Mean Peak Latency [mS]' ,'FontSize', 16);    
        title(['Mean Response Peak Latency, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);     
                %%  Response Peak Latency STD
                stimz=size(ampDel_std_mat{2}(:,:),2);
      h10=figure; 
for stim=1:stimz
        if stat_ampDel_std.p_stim{stim,1}>0.05 
            asterisk_ampDel_std{stim,1}='n.s.';
        else if stat_ampDel_std.p_stim{stim,1}<0.05 && stat_ampDel_std.p_stim{stim,1}>0.01
            asterisk_ampDel_std{stim,1}='*';
            else if stat_ampDel_std.p_stim{stim,1}<0.01 && stat_ampDel_std.p_stim{stim,1}>0.001
                    asterisk_ampDel_std{stim,1}='**';
            else if stat_ampDel_std.p_stim{stim,1}<0.001
                     asterisk_ampDel_std{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(ampDel_std_mat{1}(:,:),2)]-0.3,nanmean(ampDel_std_mat{1}(:,:),1),zeros(1,size(ampDel_std_mat{1}(:,:),2)), nanstd(ampDel_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(ampDel_std_mat{2}(:,:),2)],nanmean(ampDel_std_mat{2}(:,:),1),zeros(1,size(ampDel_std_mat{2}(:,:),2)), nanstd(ampDel_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(ampDel_std_mat{1}(:,:),2)]-0.3, nanmean(ampDel_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(ampDel_std_mat{2}(:,:),2)], nanmean(ampDel_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
    %for line pair-plots:
if paired_plot_flag==1;       
   for stim=1:stimz
        F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).ampDel_std;
        F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
        F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
   end
        for cell=1:length(files_to_analyze)
            for stim=1:stimz
                 if isnan(F1_Y{stim}(cell,:))
                    else
                        plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                         pause
                 end
            end
        end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-5;
        for stim=1:size(ampDel_std_mat{2}(:,:),2);
            if strcmp('n.s.',asterisk_ampDel_std{stim,1})==0          
             text(stim,my,asterisk_ampDel_std{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(ampDel_std_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
       ylabel('Peak Latency STD [mS]' ,'FontSize', 16);    
        title(['Mean Response Peak Latency Jitter, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);     
                        %%  Response Half-Width
                        stimz=size(halfWidth_m_mat{2}(:,:),2);
      h11=figure; 
for stim=1:stimz
        if stat_halfWidth.p_stim{stim,1}>0.05 
            asterisk_halfWidth{stim,1}='n.s.';
        else if stat_halfWidth.p_stim{stim,1}<0.05 && stat_halfWidth.p_stim{stim,1}>0.01
            asterisk_halfWidth{stim,1}='*';
            else if stat_halfWidth.p_stim{stim,1}<0.01 && stat_halfWidth.p_stim{stim,1}>0.001
                    asterisk_halfWidth{stim,1}='**';
            else if stat_halfWidth.p_stim{stim,1}<0.001
                     asterisk_halfWidth{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(halfWidth_m_mat{1}(:,:),2)]-0.3,nanmean(halfWidth_m_mat{1}(:,:),1),zeros(1,size(halfWidth_m_mat{1}(:,:),2)), nanstd(halfWidth_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(halfWidth_m_mat{2}(:,:),2)],nanmean(halfWidth_m_mat{2}(:,:),1),zeros(1,size(halfWidth_m_mat{2}(:,:),2)), nanstd(halfWidth_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(halfWidth_m_mat{1}(:,:),2)]-0.3, nanmean(halfWidth_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(halfWidth_m_mat{2}(:,:),2)], nanmean(halfWidth_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
    %for line pair-plots:
if paired_plot_flag==1;
   for stim=1:stimz
        F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).halfWidth_m;
        F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
        F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
   end
        for cell=1:length(files_to_analyze)
            for stim=1:stimz
                 if isnan(F1_Y{stim}(cell,:))
                    else
                        plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                         pause
                 end
            end
        end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-5;
        for stim=1:size(halfWidth_m_mat{2}(:,:),2);
            if strcmp('n.s.',asterisk_halfWidth{stim,1})==0          
             text(stim,my,asterisk_halfWidth{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(halfWidth_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
       ylabel('Half-Width [mS]' ,'FontSize', 16);    
        title(['Mean Response Half-Width, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
                                %%  Response Half-Width STD
    stimz=size(halfWidth_std_mat{2}(:,:),2);
      h12=figure; 
for stim=1:stimz
        if stat_halfWidth_std.p_stim{stim,1}>0.05 
            asterisk_halfWidth_std{stim,1}='n.s.';
        else if stat_halfWidth_std.p_stim{stim,1}<0.05 && stat_halfWidth_std.p_stim{stim,1}>0.01
            asterisk_halfWidth_std{stim,1}='*';
            else if stat_halfWidth_std.p_stim{stim,1}<0.01 && stat_halfWidth_std.p_stim{stim,1}>0.001
                    asterisk_halfWidth_std{stim,1}='**';
            else if stat_halfWidth_std.p_stim{stim,1}<0.001
                     asterisk_halfWidth_std{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(halfWidth_std_mat{1}(:,:),2)]-0.3,nanmean(halfWidth_std_mat{1}(:,:),1),zeros(1,size(halfWidth_std_mat{1}(:,:),2)), nanstd(halfWidth_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(halfWidth_std_mat{2}(:,:),2)],nanmean(halfWidth_std_mat{2}(:,:),1),zeros(1,size(halfWidth_std_mat{2}(:,:),2)), nanstd(halfWidth_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(halfWidth_std_mat{1}(:,:),2)]-0.3, nanmean(halfWidth_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(halfWidth_std_mat{2}(:,:),2)], nanmean(halfWidth_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
    %for line pair-plots:
if paired_plot_flag==1;
   for stim=1:stimz
        F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).halfWidth_std;
        F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
        F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
   end
        for cell=1:length(files_to_analyze)
            for stim=1:stimz
                 if isnan(F1_Y{stim}(cell,:))
                    else
                        plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                         pause
                 end
            end
        end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-5;
        for stim=1:size(halfWidth_std_mat{2}(:,:),2);
            if strcmp('n.s.',asterisk_halfWidth_std{stim,1})==0          
             text(stim,my,asterisk_halfWidth_std{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(halfWidth_std_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
       ylabel('Half-Width STD [mS]' ,'FontSize', 16);    
        title(['Mean Response Half-Width STD, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);   
        
                %% latency to 10% peak (mean)
      h13=figure; 
      stimz=size(peak10prcntTime_m_mat{2}(:,:),2);
for stim=1:stimz
        if stat_peak10prcntTime.p_stim{stim,1}>0.05 
            asterisk_peak10prcntTime{stim,1}='n.s.';
        else if stat_peak10prcntTime.p_stim{stim,1}<0.05 && stat_peak10prcntTime.p_stim{stim,1}>0.01
            asterisk_peak10prcntTime{stim,1}='*';
            else if stat_peak10prcntTime.p_stim{stim,1}<0.01 && stat_peak10prcntTime.p_stim{stim,1}>0.001
                    asterisk_peak10prcntTime{stim,1}='**';
            else if stat_peak10prcntTime.p_stim{stim,1}<0.001
                     asterisk_peak10prcntTime{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(peak10prcntTime_m_mat{1}(:,:),2)]-0.3,nanmean(peak10prcntTime_m_mat{1}(:,:),1),zeros(1,size(peak10prcntTime_m_mat{1}(:,:),2)), nanstd(peak10prcntTime_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak10prcntTime_m_mat{2}(:,:),2)],nanmean(peak10prcntTime_m_mat{2}(:,:),1),zeros(1,size(peak10prcntTime_m_mat{2}(:,:),2)), nanstd(peak10prcntTime_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak10prcntTime_m_mat{1}(:,:),2)]-0.3, nanmean(peak10prcntTime_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak10prcntTime_m_mat{2}(:,:),2)], nanmean(peak10prcntTime_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
%for line pair-plots:
if paired_plot_flag==1; 
          for stim=1:stimz
                F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).peak10prcntTime_m;
                F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
                F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
           end
                for cell=1:length(files_to_analyze)
                    for stim=1:stimz
                        if isnan(F1_Y{stim}(cell,:))
                        else
                             plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                           pause
                        end
                    end
                end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-1;
        for stim=1:stimz
            if strcmp('n.s.',asterisk_peak10prcntTime{stim,1})==0          
             text(stim,my,asterisk_peak10prcntTime{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(peak10prcntTime_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Peak 10% Latency [mS]' ,'FontSize', 16);    
        title(['Mean Latency to 10% of peak,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        %% jitter to 10% peak (latency std)
      h14=figure; 
      stimz=size(peak10prcntTime_std_mat{2}(:,:),2);
for stim=1:stimz
        if stat_peak10prcntTime_std.p_stim{stim,1}>0.05 
            asterisk_peak10prcntTime_std{stim,1}='n.s.';
        else if stat_peak10prcntTime_std.p_stim{stim,1}<0.05 && stat_peak10prcntTime_std.p_stim{stim,1}>0.01
            asterisk_peak10prcntTime_std{stim,1}='*';
            else if stat_peak10prcntTime_std.p_stim{stim,1}<0.01 && stat_peak10prcntTime_std.p_stim{stim,1}>0.001
                    asterisk_peak10prcntTime_std{stim,1}='**';
            else if stat_peak10prcntTime_std.p_stim{stim,1}<0.001
                     asterisk_peak10prcntTime_std{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(peak10prcntTime_std_mat{1}(:,:),2)]-0.3,nanmean(peak10prcntTime_std_mat{1}(:,:),1),zeros(1,size(peak10prcntTime_std_mat{1}(:,:),2)), nanstd(peak10prcntTime_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak10prcntTime_std_mat{2}(:,:),2)],nanmean(peak10prcntTime_std_mat{2}(:,:),1),zeros(1,size(peak10prcntTime_std_mat{2}(:,:),2)), nanstd(peak10prcntTime_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak10prcntTime_std_mat{1}(:,:),2)]-0.3, nanmean(peak10prcntTime_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak10prcntTime_std_mat{2}(:,:),2)], nanmean(peak10prcntTime_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
%for line pair-plots:
if paired_plot_flag==1; 
          for stim=1:stimz
                F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).peak10prcntTime_m;
                F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
                F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
           end
                for cell=1:length(files_to_analyze)
                    for stim=1:stimz
                        if isnan(F1_Y{stim}(cell,:))
                        else
                             plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                           pause
                        end
                    end
                end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-1;
        for stim=1:stimz
            if strcmp('n.s.',asterisk_peak10prcntTime_std{stim,1})==0          
             text(stim,my,asterisk_peak10prcntTime_std{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(peak10prcntTime_std_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Peak 10% Jitter [mS]' ,'FontSize', 16);    
        title(['Mean Latency STD to 10% of peak,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
                        %% latency to 90% peak (mean)
      h15=figure; 
      stimz=size(peak90prcntTime_m_mat{2}(:,:),2);
for stim=1:stimz
        if stat_peak90prcntTime.p_stim{stim,1}>0.05 
            asterisk_peak90prcntTime{stim,1}='n.s.';
        else if stat_peak90prcntTime.p_stim{stim,1}<0.05 && stat_peak90prcntTime.p_stim{stim,1}>0.01
            asterisk_peak90prcntTime{stim,1}='*';
            else if stat_peak90prcntTime.p_stim{stim,1}<0.01 && stat_peak90prcntTime.p_stim{stim,1}>0.001
                    asterisk_peak90prcntTime{stim,1}='**';
            else if stat_peak90prcntTime.p_stim{stim,1}<0.001
                     asterisk_peak90prcntTime{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(peak90prcntTime_m_mat{1}(:,:),2)]-0.3,nanmean(peak90prcntTime_m_mat{1}(:,:),1),zeros(1,size(peak90prcntTime_m_mat{1}(:,:),2)), nanstd(peak90prcntTime_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak90prcntTime_m_mat{2}(:,:),2)],nanmean(peak90prcntTime_m_mat{2}(:,:),1),zeros(1,size(peak90prcntTime_m_mat{2}(:,:),2)), nanstd(peak90prcntTime_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak90prcntTime_m_mat{1}(:,:),2)]-0.3, nanmean(peak90prcntTime_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak90prcntTime_m_mat{2}(:,:),2)], nanmean(peak90prcntTime_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
%for line pair-plots:
if paired_plot_flag==1; 
          for stim=1:stimz
                F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).peak90prcntTime_m;
                F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
                F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
           end
                for cell=1:length(files_to_analyze)
                    for stim=1:stimz
                        if isnan(F1_Y{stim}(cell,:))
                        else
                             plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',90,'markerfacecolor','k')
%                           pause
                        end
                    end
                end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-1;
        for stim=1:stimz
            if strcmp('n.s.',asterisk_peak90prcntTime{stim,1})==0          
             text(stim,my,asterisk_peak90prcntTime{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(peak90prcntTime_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Peak 90% Latency [mS]' ,'FontSize', 16);    
        title(['Mean Latency to 90% of peak,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        %% jitter to 90% peak (latency std)
      h16=figure; 
      stimz=size(peak90prcntTime_std_mat{2}(:,:),2);
for stim=1:stimz
        if stat_peak90prcntTime_std.p_stim{stim,1}>0.05 
            asterisk_peak90prcntTime_std{stim,1}='n.s.';
        else if stat_peak90prcntTime_std.p_stim{stim,1}<0.05 && stat_peak90prcntTime_std.p_stim{stim,1}>0.01
            asterisk_peak90prcntTime_std{stim,1}='*';
            else if stat_peak90prcntTime_std.p_stim{stim,1}<0.01 && stat_peak90prcntTime_std.p_stim{stim,1}>0.001
                    asterisk_peak90prcntTime_std{stim,1}='**';
            else if stat_peak90prcntTime_std.p_stim{stim,1}<0.001
                     asterisk_peak90prcntTime_std{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(peak90prcntTime_std_mat{1}(:,:),2)]-0.3,nanmean(peak90prcntTime_std_mat{1}(:,:),1),zeros(1,size(peak90prcntTime_std_mat{1}(:,:),2)), nanstd(peak90prcntTime_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak90prcntTime_std_mat{2}(:,:),2)],nanmean(peak90prcntTime_std_mat{2}(:,:),1),zeros(1,size(peak90prcntTime_std_mat{2}(:,:),2)), nanstd(peak90prcntTime_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak90prcntTime_std_mat{1}(:,:),2)]-0.3, nanmean(peak90prcntTime_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak90prcntTime_std_mat{2}(:,:),2)], nanmean(peak90prcntTime_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
%for line pair-plots:
if paired_plot_flag==1; 
          for stim=1:stimz
                F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).peak90prcntTime_m;
                F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
                F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
           end
                for cell=1:length(files_to_analyze)
                    for stim=1:stimz
                        if isnan(F1_Y{stim}(cell,:))
                        else
                             plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                           pause
                        end
                    end
                end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-1;
        for stim=1:stimz
            if strcmp('n.s.',asterisk_peak90prcntTime_std{stim,1})==0          
             text(stim,my,asterisk_peak90prcntTime_std{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(peak90prcntTime_std_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Peak 90% Jitter [mS]' ,'FontSize', 16);    
        title(['Mean Latency STD to 90% of peak,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
                        %% rise-time from 10% to 90% peak (mean)
      h17=figure; 
      stimz=size(peak10to90Time_m_mat{2}(:,:),2);
for stim=1:stimz
        if stat_peak10to90Time.p_stim{stim,1}>0.05 
            asterisk_peak10to90Time{stim,1}='n.s.';
        else if stat_peak10to90Time.p_stim{stim,1}<0.05 && stat_peak10to90Time.p_stim{stim,1}>0.01
            asterisk_peak10to90Time{stim,1}='*';
            else if stat_peak10to90Time.p_stim{stim,1}<0.01 && stat_peak10to90Time.p_stim{stim,1}>0.001
                    asterisk_peak10to90Time{stim,1}='**';
            else if stat_peak10to90Time.p_stim{stim,1}<0.001
                     asterisk_peak10to90Time{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(peak10to90Time_m_mat{1}(:,:),2)]-0.3,nanmean(peak10to90Time_m_mat{1}(:,:),1),zeros(1,size(peak10to90Time_m_mat{1}(:,:),2)), nanstd(peak10to90Time_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak10to90Time_m_mat{2}(:,:),2)],nanmean(peak10to90Time_m_mat{2}(:,:),1),zeros(1,size(peak10to90Time_m_mat{2}(:,:),2)), nanstd(peak10to90Time_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak10to90Time_m_mat{1}(:,:),2)]-0.3, nanmean(peak10to90Time_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak10to90Time_m_mat{2}(:,:),2)], nanmean(peak10to90Time_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
%for line pair-plots:
if paired_plot_flag==1; 
          for stim=1:stimz
                F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).peak10to90Time_m;
                F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
                F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
           end
                for cell=1:length(files_to_analyze)
                    for stim=1:stimz
                        if isnan(F1_Y{stim}(cell,:))
                        else
                             plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                           pause
                        end
                    end
                end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-1;
        for stim=1:stimz
            if strcmp('n.s.',asterisk_peak10to90Time{stim,1})==0          
             text(stim,my,asterisk_peak10to90Time{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(peak10to90Time_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Rise-Time [mS]' ,'FontSize', 16);    
        title(['Rise-time from 10% to 90% peak, n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        %% jitter in rise-time from 10% to 90% peak (std)
      h18=figure; 
      stimz=size(peak10to90Time_std_mat{2}(:,:),2);
for stim=1:stimz
        if stat_peak10to90Time_std.p_stim{stim,1}>0.05 
            asterisk_peak10to90Time_std{stim,1}='n.s.';
        else if stat_peak10to90Time_std.p_stim{stim,1}<0.05 && stat_peak10to90Time_std.p_stim{stim,1}>0.01
            asterisk_peak10to90Time_std{stim,1}='*';
            else if stat_peak10to90Time_std.p_stim{stim,1}<0.01 && stat_peak10to90Time_std.p_stim{stim,1}>0.001
                    asterisk_peak10to90Time_std{stim,1}='**';
            else if stat_peak10to90Time_std.p_stim{stim,1}<0.001
                     asterisk_peak10to90Time_std{stim,1}='***';
                end
                end
            end
        end
end
   
hold on      
        errbar_h=errorbar([1:size(peak10to90Time_std_mat{1}(:,:),2)]-0.3,nanmean(peak10to90Time_std_mat{1}(:,:),1),zeros(1,size(peak10to90Time_std_mat{1}(:,:),2)), nanstd(peak10to90Time_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        errbar_h=errorbar([1:size(peak10to90Time_std_mat{2}(:,:),2)],nanmean(peak10to90Time_std_mat{2}(:,:),1),zeros(1,size(peak10to90Time_std_mat{2}(:,:),2)), nanstd(peak10to90Time_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        bar([1:size(peak10to90Time_std_mat{1}(:,:),2)]-0.3, nanmean(peak10to90Time_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(peak10to90Time_std_mat{2}(:,:),2)], nanmean(peak10to90Time_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
%for line pair-plots:
if paired_plot_flag==1; 
          for stim=1:stimz
                F1_Y{stim}(:,:)=event_evoked_stat.stim_num(stim).peak10to90Time_m;
                F1_X{stim}(:,1)=(stim-0.3)*ones(size(F1_Y,1),1);
                F1_X{stim}(:,2)=stim*ones(size(F1_Y,1),1);
           end
                for cell=1:length(files_to_analyze)
                    for stim=1:stimz
                        if isnan(F1_Y{stim}(cell,:))
                        else
                             plot(F1_X{stim}(1,:),F1_Y{stim}(cell,:),'color',color_table_lines(cell,:),'linewidth',1.5); %,'markersize',10,'markerfacecolor','k')
%                           pause
                        end
                    end
                end
end
        ylim_data=[get(gca,'ylim')]';
        my=max(ylim_data)-1;
        for stim=1:stimz
            if strcmp('n.s.',asterisk_peak10to90Time_std{stim,1})==0          
             text(stim,my,asterisk_peak10to90Time_std{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
            end
        end
        hold off
        set(gca,'xlim',[0,size(peak10to90Time_std_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Rise-Time STD [mS]' ,'FontSize', 16);    
        title(['Rise-Time STD 10%-90% of peak,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%% rmANOVA with matlab function -amplitude. 
% peak_amp_noNB(:,:)= amplitude_m_mat{1}(:,:);
% peak_amp_NB(:,:)= amplitude_m_mat{2}(:,:);
% peak_amp_all(:,:)=[peak_amp_noNB;peak_amp_NB];
%  S_1(:,1)=1:length(files_to_analyze); 
%  
% % clear ta within rm_amp mauchly_amp ranovatbl_amp eps_amp eps NBbyTime_amp TimebyNB_amp
% 
% ta_vector_names={'peak_amp_noNB(:,1:10)','peak_amp_NB(:,1:10)'};
% ta=table(S_1,peak_amp_noNB(:,1),peak_amp_noNB(:,2),peak_amp_noNB(:,3),peak_amp_noNB(:,4),peak_amp_noNB(:,5),peak_amp_noNB(:,6),...
%     peak_amp_noNB(:,7),peak_amp_noNB(:,8),peak_amp_noNB(:,9),peak_amp_noNB(:,10),peak_amp_noNB(:,11),...
%     peak_amp_NB(:,1),peak_amp_NB(:,2),peak_amp_NB(:,3),peak_amp_NB(:,4),peak_amp_NB(:,5),peak_amp_NB(:,6),...
%     peak_amp_NB(:,7),peak_amp_NB(:,8),peak_amp_NB(:,9),peak_amp_NB(:,10),peak_amp_NB(:,11),...
%     'variablenames', {'cells','Y1','Y2','Y3','Y4','Y5','Y6','Y7','Y8','Y9','Y10','Y11','Y12','Y13','Y14','Y15','Y16','Y17','Y18','Y19','Y20','Y21','Y22'}); %,,'Y21','Y22'
% factorNames = {'NB','Time'};
% within = table({'N';'N';'N';'N';'N';'N';'N';'N';'N';'N';'N';'Y';'Y';'Y';'Y';'Y';'Y';'Y';'Y';'Y';'Y';'Y'},{'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11'},'VariableNames',factorNames); %the levels of each factor in each measurment witgin subject
% % fit the repeated measures model
% rm_amp = fitrm(ta,'Y1-Y22~1','WithinDesign',within);
% %test for sphericity
% mauchly_amp=mauchly(rm_amp);
% % run my repeated measures anova here
% [ranovatbl_amp] = ranova(rm_amp,'withinmodel','NB*Time');
% %multiple comparisons:
% NBbyTime_amp = multcompare(rm_amp,'NB','By','Time');
% TimebyNB_amp = multcompare(rm_amp,'Time','By','NB');
% eps_amp = epsilon(rm_amp);
% eps=table2cell(eps_amp(1,2));
% if eps{1}<0.75
%     amp_rmANOVA_NB_effect_p= table2cell(ranovatbl_amp(3,6));
%     amp_rmANOVA_Time_effect_p= table2cell(ranovatbl_amp(5,6));
%     amp_rmANOVA_interaction_effect_p= table2cell(ranovatbl_amp(7,6));
% else
%      amp_rmANOVA_NB_effect_p= table2cell(ranovatbl_amp(3,7));
%      amp_rmANOVA_Time_effect_p= table2cell(ranovatbl_amp(5,7));
%      amp_rmANOVA_interaction_effect_p= table2cell(ranovatbl_amp(7,7));
% end
% table_column1=table2cell(NBbyTime_amp(:,1));
% 
% for stim=1:size(peak_amp_noNB,2) 
%     row= find(strcmp(num2str(stim),table_column1),1);
%     amp_p_stim{1}(stim,:)=table2cell(NBbyTime_amp(row,6));
% end
%     
% stat_amp.table=ta;
% stat_amp.table_data_vecs=ta_vector_names;
% stat_amp.within_design=within;
% stat_amp.rm=rm_amp;
% stat_amp.mauchly=mauchly_amp;
% stat_amp.eps=eps_amp;
% stat_amp.ANOVA=ranovatbl_amp;
% stat_amp.multcomp_NBbyTime=NBbyTime_amp;
% stat_amp.multcomp_TimebyNB=TimebyNB_amp;
% stat_amp.rmANOVA_NB_effect_p= amp_rmANOVA_NB_effect_p;
% stat_amp.rmANOVA_Time_effect_p= amp_rmANOVA_Time_effect_p;
% stat_amp.rmANOVA_interaction_effect_p= amp_rmANOVA_interaction_effect_p;
% stat_amp.p_stim=amp_p_stim{1};
%% Plots with p-values from t-tests
% h1=figure;        
%         hold on
%         errbar_h=errorbar([1:size(failures_m_mat{1}(:,:),2)]-0.3,nanmean(failures_m_mat{1}(:,:),1),zeros(1,size(failures_m_mat{1}(:,:),2)), nanstd(failures_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         errbar_h=errorbar([1:size(failures_m_mat{2}(:,:),2)],nanmean(failures_m_mat{2}(:,:),1),zeros(1,size(failures_m_mat{2}(:,:),2)), nanstd(failures_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         bar([1:size(failures_m_mat{1}(:,:),2)]-0.3, nanmean(failures_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(failures_m_mat{2}(:,:),2)], nanmean(failures_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%         p_failures_m=max(ylim_data).*ones(size(s_failures_m));
%         plot(s_failures_m-0.15,p_failures_m,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Failures rate' ,'FontSize', 16);    
%         title(['Mean Rate of Failures,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%         h2=figure;        
%         hold on
%         errbar_h=errorbar([1:size(nonspecific_count_m_mat{1}(:,:),2)]-0.3,nanmean(nonspecific_count_m_mat{1}(:,:),1),zeros(1,size(nonspecific_count_m_mat{1}(:,:),2)), nanstd(nonspecific_count_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         errbar_h=errorbar([1:size(nonspecific_count_m_mat{2}(:,:),2)],nanmean(nonspecific_count_m_mat{2}(:,:),1),zeros(1,size(nonspecific_count_m_mat{2}(:,:),2)), nanstd(nonspecific_count_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         bar([1:size(nonspecific_count_m_mat{1}(:,:),2)]-0.3, nanmean(nonspecific_count_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(nonspecific_count_m_mat{2}(:,:),2)], nanmean(nonspecific_count_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%         p_nonspecific_count_m=max(ylim_data).*ones(size(s_nonspecific_count_m));
%         plot(s_nonspecific_count_m-0.15,p_nonspecific_count_m,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Non-specific response rate' ,'FontSize', 16);    
%         title(['Mean Rate of Non-specific responses,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
% 
%         h3=figure;        
%         hold on
%         errbar_h=errorbar([1:size(onset_m_mat{1}(:,:),2)]-0.3,nanmean(onset_m_mat{1}(:,:),1),zeros(1,size(onset_m_mat{1}(:,:),2)), nanstd(onset_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         errbar_h=errorbar([1:size(onset_m_mat{2}(:,:),2)],nanmean(onset_m_mat{2}(:,:),1),zeros(1,size(onset_m_mat{2}(:,:),2)), nanstd(onset_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         bar([1:size(onset_m_mat{1}(:,:),2)]-0.3, nanmean(onset_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(onset_m_mat{2}(:,:),2)], nanmean(onset_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%         p_onset_m=max(ylim_data).*ones(size(s_onset_m));
%         plot(s_onset_m-0.15,p_onset_m,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Onset Latency [mS]' ,'FontSize', 16);    
%         title(['Mean Onset Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
% 
%         h4=figure;   
%         hold on
%         errbar_h=errorbar([1:size(onset_std_mat{1}(:,:),2)]-0.3,nanmean(onset_std_mat{1}(:,:),1),zeros(1,size(onset_std_mat{1}(:,:),2)), nanstd(onset_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         errbar_h=errorbar([1:size(onset_std_mat{2}(:,:),2)],nanmean(onset_std_mat{2}(:,:),1),zeros(1,size(onset_std_mat{2}(:,:),2)), nanstd(onset_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         bar([1:size(onset_std_mat{1}(:,:),2)]-0.3, nanmean(onset_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(onset_std_mat{2}(:,:),2)], nanmean(onset_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%         p_onset_std=max(ylim_data).*ones(size(s_onset_std));
%         plot(s_onset_std-0.15,p_onset_std,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Onset Latency STD [mS]' ,'FontSize', 16);    
%         title(['Mean response jitter,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  

% 
%       h5=figure; 
% for stim=1:size(amplitude_m_mat{2}(:,:),2);
%         if stat_amp.p_stim{stim,1}>0.05 
%             asterisk_amp{stim,1}='n.s.';
%         else if stat_amp.p_stim{stim,1}<0.05 && stat_amp.p_stim{stim,1}>0.01
%             asterisk_amp{stim,1}='*';
%             else if stat_amp.p_stim{stim,1}<0.01 && stat_amp.p_stim{stim,1}>0.001
%                     asterisk_amp{stim,1}='**';
%             else if stat_amp.p_stim{stim,1}<0.001
%                      asterisk_amp{stim,1}='***';
%                 end
%                 end
%             end
%         end
% end
%    
% hold on      
%         errbar_h=errorbar([1:size(amplitude_m_mat{1}(:,:),2)]-0.3,nanmean(amplitude_m_mat{1}(:,:),1),zeros(1,size(amplitude_m_mat{1}(:,:),2)), nanstd(amplitude_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         errbar_h=errorbar([1:size(amplitude_m_mat{2}(:,:),2)],nanmean(amplitude_m_mat{2}(:,:),1),zeros(1,size(amplitude_m_mat{2}(:,:),2)), nanstd(amplitude_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         bar([1:size(amplitude_m_mat{1}(:,:),2)]-0.3, nanmean(amplitude_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(amplitude_m_mat{2}(:,:),2)], nanmean(amplitude_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))       
%         ylim_data=[get(gca,'ylim')]';
%         my=max(ylim_data)-1;
%         for stim=1:size(amplitude_m_mat{2}(:,:),2);
%             if strcmp('n.s.',asterisk_amp{stim,1})==0          
%              text(stim,my,asterisk_amp{stim,1},'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)      
%             end
%         end
% %         p_amplitude_m=max(ylim_data).*ones(size(s_amplitude_m));
% %         plot(s_amplitude_m-0.15,p_amplitude_m,'k*')
% %         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b') %change line to zero
%         hold off
%         set(gca,'xlim',[0,size(amplitude_m_mat{2}(:,:),2)+0.5], 'xtick',[1:1:11],'xticklabel',[1:1:11])
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Mean Response Amplitude [mV]' ,'FontSize', 16);    
%         title(['Mean Response Amplitude - local baseline,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%          h6=figure;        
%         hold on
%         errbar_h=errorbar([1:size(amplitude_std_mat{1}(:,:),2)]-0.3,nanmean(amplitude_std_mat{1}(:,:),1),zeros(1,size(amplitude_std_mat{1}(:,:),2)), nanstd(amplitude_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         errbar_h=errorbar([1:size(amplitude_std_mat{2}(:,:),2)],nanmean(amplitude_std_mat{2}(:,:),1),zeros(1,size(amplitude_std_mat{2}(:,:),2)), nanstd(amplitude_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         bar([1:size(amplitude_std_mat{1}(:,:),2)]-0.3, nanmean(amplitude_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(amplitude_std_mat{2}(:,:),2)], nanmean(amplitude_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%         p_amplitude_std=max(ylim_data).*ones(size(s_amplitude_std));
%         plot(s_amplitude_std-0.15,p_amplitude_std,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Amplitude STD [mV]' ,'FontSize', 16);    
%         title(['Mean Amplitude STD (trial-to-trial variability) ,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%  
%         h7=figure;        
%         hold on
% %         set(gca,'Ydir','reverse'); %added because the membrane pot. values are negative. more changes: error bars are symmetric and caps not removed, and asterisks yvalue is min instead of max
%         errbar_h=errorbar([1:size(ampVal_m_mat{1}(:,:),2)]-0.3,nanmean(ampVal_m_mat{1}(:,:),1), nanstd(ampVal_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         errbar_h=errorbar([1:size(ampVal_m_mat{2}(:,:),2)],nanmean(ampVal_m_mat{2}(:,:),1), nanstd(ampVal_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         bar([1:size(ampVal_m_mat{1}(:,:),2)]-0.3, nanmean(ampVal_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(ampVal_m_mat{2}(:,:),2)], nanmean(ampVal_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%          set(gca,'ylim',[min(ylim_data),min(ylim_data)+30])
%         p_ampVal_m=min(ylim_data).*ones(size(s_ampVal_m))+5;
%         plot(s_ampVal_m-0.15,p_ampVal_m,'k*')  
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Response Peak Value [mV]' ,'FontSize', 16);    
%         title(['Mean Response Peak Value ,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%                
%         h8=figure;        
%         hold on
% %         set(gca,'Ydir','reverse');
%         errbar_h=errorbar([1:size(onVal_m_mat{1}(:,:),2)]-0.3,nanmean(onVal_m_mat{1}(:,:),1), nanstd(onVal_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         errbar_h=errorbar([1:size(onVal_m_mat{2}(:,:),2)],nanmean(onVal_m_mat{2}(:,:),1), nanstd(onVal_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         bar([1:size(onVal_m_mat{1}(:,:),2)]-0.3, nanmean(onVal_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(onVal_m_mat{2}(:,:),2)], nanmean(onVal_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%      set(gca,'ylim',[min(ylim_data),min(ylim_data)+30])
%         p_onVal_m=min(ylim_data).*ones(size(s_onVal_m))+5; 
%         plot(s_onVal_m-0.15,p_onVal_m,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Onset Vm [mV]' ,'FontSize', 16);    
%         title(['Mean Onset Vm value,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%         h9=figure;        
%         hold on
%         errbar_h=errorbar([1:size(ampDel_m_mat{1}(:,:),2)]-0.3,nanmean(ampDel_m_mat{1}(:,:),1),zeros(1,size(ampDel_m_mat{1}(:,:),2)), nanstd(ampDel_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         errbar_h=errorbar([1:size(ampDel_m_mat{2}(:,:),2)],nanmean(ampDel_m_mat{2}(:,:),1),zeros(1,size(ampDel_m_mat{2}(:,:),2)), nanstd(ampDel_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         bar([1:size(ampDel_m_mat{1}(:,:),2)]-0.3, nanmean(ampDel_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(ampDel_m_mat{2}(:,:),2)], nanmean(ampDel_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%         p_ampDel_m=max(ylim_data).*ones(size(s_ampDel_m));
%         plot(s_ampDel_m-0.15,p_ampDel_m,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Peak Latency [mS]' ,'FontSize', 16);    
%         title(['Mean Peak Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%     
%          h10=figure;        
%         hold on
%         errbar_h=errorbar([1:size(ampDel_std_mat{1}(:,:),2)]-0.3,nanmean(ampDel_std_mat{1}(:,:),1),zeros(1,size(ampDel_std_mat{1}(:,:),2)), nanstd(ampDel_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         errbar_h=errorbar([1:size(ampDel_std_mat{2}(:,:),2)],nanmean(ampDel_std_mat{2}(:,:),1),zeros(1,size(ampDel_std_mat{2}(:,:),2)), nanstd(ampDel_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         bar([1:size(ampDel_std_mat{1}(:,:),2)]-0.3, nanmean(ampDel_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(ampDel_std_mat{2}(:,:),2)], nanmean(ampDel_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%         p_ampDel_std=max(ylim_data).*ones(size(s_ampDel_std));
%         plot(s_ampDel_std-0.15,p_ampDel_std,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Peak Latency STD [mS]' ,'FontSize', 16);    
%         title(['Mean Peak Latency Jitter,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%          h11=figure;        
%         hold on
%         errbar_h=errorbar([1:size(halfWidth_m_mat{1}(:,:),2)]-0.3,nanmean(halfWidth_m_mat{1}(:,:),1),zeros(1,size(halfWidth_m_mat{1}(:,:),2)), nanstd(halfWidth_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         errbar_h=errorbar([1:size(halfWidth_m_mat{2}(:,:),2)],nanmean(halfWidth_m_mat{2}(:,:),1),zeros(1,size(halfWidth_m_mat{2}(:,:),2)), nanstd(halfWidth_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         bar([1:size(halfWidth_m_mat{1}(:,:),2)]-0.3, nanmean(halfWidth_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(halfWidth_m_mat{2}(:,:),2)], nanmean(halfWidth_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%         p_halfWidth_m=max(ylim_data).*ones(size(s_halfWidth_m));
%         plot(s_halfWidth_m-0.15,p_halfWidth_m,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Half-Width [mS]' ,'FontSize', 16);    
%         title(['Mean Half-Width,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%          h12=figure;        
%         hold on
%         errbar_h=errorbar([1:size(halfWidth_std_mat{1}(:,:),2)]-0.3,nanmean(halfWidth_std_mat{1}(:,:),1),zeros(1,size(halfWidth_std_mat{1}(:,:),2)), nanstd(halfWidth_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         errbar_h=errorbar([1:size(halfWidth_std_mat{2}(:,:),2)],nanmean(halfWidth_std_mat{2}(:,:),1),zeros(1,size(halfWidth_std_mat{2}(:,:),2)), nanstd(halfWidth_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
% %         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
%         bar([1:size(halfWidth_std_mat{1}(:,:),2)]-0.3, nanmean(halfWidth_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
%         bar([1:size(halfWidth_std_mat{2}(:,:),2)], nanmean(halfWidth_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
%         ylim_data=[get(gca,'ylim')]';
%         p_halfWidth_std=max(ylim_data).*ones(size(s_halfWidth_std));
%         plot(s_halfWidth_std-0.15,p_halfWidth_std,'k*')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Half-Width STD [mS]' ,'FontSize', 16);    
%         title(['Mean Half-Width STD,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);        
%% Plot parameters along the train stim - version 1: line+error bars of percent change - preparations
% for stim_num=1:11;
%        if event_evoked_stat.stim_num(stim_num).lillietest_h_change_failures_m==0      
%             h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_failures_m; 
%   else h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_failures_m; 
%         end
%      if event_evoked_stat.stim_num(stim_num).lillietest_h_change_onset_m==0      
%             h_onset_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_onset_m; 
%   else h_onset_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onset_m; 
%      end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_onset_std==0      
%             h_onset_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_onset_std; 
%   else h_onset_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onset_std; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_amplitude_m==0      
%             h_amplitude_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_amplitude_m; 
%   else h_amplitude_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_amplitude_m; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_amplitude_std==0      
%             h_amplitude_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_amplitude_std; 
%   else h_amplitude_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_amplitude_std; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_onVal_m==0      
%             h_onVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_onVal_m; 
%   else h_onVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onVal_m; 
%    end
% %    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_onVal_std==0      
% %             h_onVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_onVal_std; 
% %   else h_onVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onVal_std; 
% %    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampVal_m==0      
%             h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_ampVal_m; 
%   else h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampVal_m; 
%    end
% %    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampVal_std==0      
% %             h_ampVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_ampVal_std; 
% %   else h_ampVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampVal_std; 
% %    end
%  if event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampDel_m==0      
%             h_ampDel_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_ampDel_m; 
%   else h_ampDel_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampDel_m; 
%  end
%     if event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampDel_std==0      
%             h_ampDel_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_ampDel_std; 
%   else h_ampDel_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampDel_std; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_halfWidth_m==0      
%             h_halfWidth_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_halfWidth_m; 
%   else h_halfWidth_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_halfWidth_m; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_halfWidth_std==0      
%             h_halfWidth_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_halfWidth_std; 
%   else h_halfWidth_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_halfWidth_std; 
%    end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_adapt_amp==0      
%             h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_adapt_amp; 
%   else h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_adapt_amp; 
%    end
%  if event_evoked_stat.stim_num(stim_num).lillietest_h_change_nonspecific_count_m==0      
%             h_nonspecific_count_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_nonspecific_count_m; 
%   else h_nonspecific_count_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_nonspecific_count_m; 
%    end
%     end
%     %getting the x-axis (stimulus number) and y-axis for the asterisks of significance
%      s_failures_m=find(h_failures_m==1);
%      s_onset_m=find(h_onset_m==1);
%      s_onset_std=find(h_onset_std==1); 
%      s_amplitude_m=find(h_amplitude_m==1);
%      s_amplitude_std=find(h_amplitude_std==1);
%      s_onVal_m=find(h_onVal_m==1);
% %      s_onVal_std=find(h_onVal_std==1);
%      s_ampVal_m=find(h_ampVal_m==1);
% %      s_ampVal_std=find(h_ampVal_std==1);
%      s_ampDel_m=find(h_ampDel_m==1);
%      s_ampDel_std=find(h_ampDel_std==1);
%      s_halfWidth_m=find(h_halfWidth_m==1);
%      s_halfWidth_std=find(h_halfWidth_std==1);
%      s_adapt_amp=find(h_adapt_amp==1);
%      s_nonspecific_count_m=find(h_nonspecific_count_m==1);
%%    Plot parameters along the train stim - version 1: line+error bars of percent change 
% close all
%     if print_flag==1;
% h1=figure;        
%         hold on
%         errorbar([1:size(change_failures_m_mat,2)], nanmean(change_failures_m_mat,1),nanstd(change_failures_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_failures_m=max(ylim_data).*ones(size(s_failures_m));
%         plot(s_failures_m,p_failures_m,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Failures difference[%]' ,'FontSize', 16);    
%         title(['Mean Rate of Failures,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%         h2=figure;        
%         hold on
%         errorbar([1:size(change_nonspecific_count_m_mat,2)], nanmean(change_nonspecific_count_m_mat,1),nanstd(change_nonspecific_count_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_nonspecific_count_m=max(ylim_data).*ones(size(s_nonspecific_count_m));
%         plot(s_nonspecific_count_m,p_nonspecific_count_m,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Non-specific responses difference[%]' ,'FontSize', 16);    
%         title(['Mean Rate of Non-specific responses,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
% 
%         h3=figure;        
%         hold on
%         errorbar([1:size(change_onset_m_mat,2)], nanmean(change_onset_m_mat,1),nanstd(change_onset_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_onset_m=max(ylim_data).*ones(size(s_onset_m));
%         plot(s_onset_m,p_onset_m,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Onset Latency [%change]' ,'FontSize', 16);    
%         title(['Mean Onset Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
% 
%         h4=figure;        
%         hold on
%         errorbar([1:size(change_onset_std_mat,2)], nanmean(change_onset_std_mat,1),nanstd(change_onset_std_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_onset_std=max(ylim_data).*ones(size(s_onset_std));
%         plot(s_onset_std,p_onset_std,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Onset Latency STD [%change]' ,'FontSize', 16);    
%         title(['Mean response jitter,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
% 
% h5=figure;        
%         hold on
%         errorbar([1:size(change_amplitude_m_mat,2)], nanmean(change_amplitude_m_mat,1),nanstd(change_amplitude_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_amplitude_m=max(ylim_data).*ones(size(s_amplitude_m));
%         plot(s_amplitude_m,p_amplitude_m,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Response Amplitude [%change]' ,'FontSize', 16);    
%         title(['Mean Response Amplitude,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%  h6=figure;        
%         hold on
%         errorbar([1:size(change_amplitude_std_mat,2)], nanmean(change_amplitude_std_mat,1),nanstd(change_amplitude_std_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_amplitude_std=max(ylim_data).*ones(size(s_amplitude_std));
%         plot(s_amplitude_std,p_amplitude_std,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Amplitude STD [%change]' ,'FontSize', 16);    
%         title(['Mean Amplitude STD (trial-to-trial variability) ,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%  
%         h7=figure;        
%         hold on
%         errorbar([1:size(change_ampVal_m_mat,2)], nanmean(change_ampVal_m_mat,1),nanstd(change_ampVal_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_ampVal_m=max(ylim_data).*ones(size(s_ampVal_m));
%         plot(s_ampVal_m,p_ampVal_m,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Response Peak Value [%change]' ,'FontSize', 16);    
%         title(['Mean Response Peak Value ,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%                
%         h8=figure;        
%         hold on
%         errorbar([1:size(change_onVal_m_mat,2)], nanmean(change_onVal_m_mat,1),nanstd(change_onVal_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_onVal_m=max(ylim_data).*ones(size(s_onVal_m));
%         plot(s_onVal_m,p_onVal_m,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Onset Vm [%change]' ,'FontSize', 16);    
%         title(['Mean Onset Vm value,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%         h9=figure;        
%         hold on
%         errorbar([1:size(change_ampDel_m_mat,2)], nanmean(change_ampDel_m_mat,1),nanstd(change_ampDel_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_ampDel_m=max(ylim_data).*ones(size(s_ampDel_m));
%         plot(s_ampDel_m,p_ampDel_m,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Peak Latency [%change]' ,'FontSize', 16);    
%         title(['Mean Peak Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%     
%          h10=figure;        
%         hold on
%         errorbar([1:size(change_ampDel_std_mat,2)], nanmean(change_ampDel_std_mat,1),nanstd(change_ampDel_std_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_ampDel_std=max(ylim_data).*ones(size(s_ampDel_std));
%         plot(s_ampDel_std,p_ampDel_std,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Peak Latency STD [%change]' ,'FontSize', 16);    
%         title(['Mean Peak Latency Jitter,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%          h11=figure;        
%         hold on
%         errorbar([1:size(change_halfWidth_m_mat,2)], nanmean(change_halfWidth_m_mat,1),nanstd(change_halfWidth_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_halfWidth_m=max(ylim_data).*ones(size(s_halfWidth_m));
%         plot(s_halfWidth_m,p_halfWidth_m,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Half-Width [%change]' ,'FontSize', 16);    
%         title(['Mean Half-Width,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%          h12=figure;        
%         hold on
%         errorbar([1:size(change_halfWidth_std_mat,2)], nanmean(change_halfWidth_std_mat,1),nanstd(change_halfWidth_std_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         ylim_data=[get(gca,'ylim')]';
%         p_halfWidth_std=max(ylim_data).*ones(size(s_halfWidth_std));
%         plot(s_halfWidth_std,p_halfWidth_std,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Half-Width STD [%change]' ,'FontSize', 16);    
%         title(['Mean Half-Width STD,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
%         
%     end 
 %% plots
%  close all
for stim_num=1;
tmp_Y=event_evoked_stat.stim_num(stim_num).adapt_amp_m';
nanmat=~isnan(tmp_Y);
nancount=sum(nanmat(1,:));

tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E =event_evoked_stat.stim_num(stim_num).adapt_amp_m_std;
if event_evoked_stat.stim_num(stim_num).wilcoxon_p_adapt_amp >0.05 
    asterisk_sp='n.s.';
else if event_evoked_stat.stim_num(stim_num).wilcoxon_p_adapt_amp<0.05 && event_evoked_stat.stim_num(stim_num).wilcoxon_p_adapt_amp>0.01
    asterisk_sp='*';
    else if event_evoked_stat.stim_num(stim_num).wilcoxon_p_adapt_amp<0.01 && event_evoked_stat.stim_num(stim_num).wilcoxon_p_adapt_amp>0.001
            asterisk_sp='**';
    else if event_evoked_stat.stim_num(stim_num).wilcoxon_p_adapt_amp<0.001
             asterisk_sp='***';
        end
        end
    end
end
my=max(max(event_evoked_stat.stim_num(stim_num).adapt_amp_m))*1; 
g1=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
% boxplot(tmp_Y')
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Adaptation Amplitude Ratio', 'FontSize', 28,'fontname', 'arial');
        title(['Adaptation Amplitude Ratio, n=' num2str(nancount) ', p=' num2str(event_evoked_stat.stim_num(stim_num).ttest_p_adapt_amp)] ,'FontSize', 20,'fontname', 'arial');   

%         pause 
% close all
end
        %% save figures
if save_flag==1
cd(path_output)
saveas(h1,'Evoked Failures','fig') 
print(h1,'Evoked Failures','-dpng','-r600','-opengl')
print(h2,'Evoked Nonspecific Responses','-dpng','-r600','-opengl')
saveas(h2,'Evoked Nonspecific Responses','fig') 
print(h3,'Evoked Onset Latency','-dpng','-r600','-opengl')
saveas(h3,'Evoked Onset Latency','fig') 
print(h4,'Evoked Onset Latency STD','-dpng','-r600','-opengl')
saveas(h4,'Evoked Onset Latency STD','fig') 
print(h5,'Evoked Amplitude','-dpng','-r600','-opengl')
saveas(h5,'Evoked Amplitude','fig') 
print(h6,'Evoked Amplitude STD','-dpng','-r600','-opengl')
saveas(h6,'Evoked Amplitude STD','fig') 
print(h7,'Evoked Peak Value','-dpng','-r600','-opengl')
saveas(h7,'Evoked Peak Value','fig') 
print(h8,'Evoked Onset Value','-dpng','-r600','-opengl')
saveas(h8,'Evoked Onset Value','fig') 
print(h9,'Evoked Peak Latency','-dpng','-r600','-opengl')
saveas(h9,'Evoked Peak Latency','fig') 
print(h10,'Evoked Peak Latency STD','-dpng','-r600','-opengl')
saveas(h10,'Evoked Peak Latency STD','fig') 
print(h11,'Evoked Half-Width','-dpng','-r600','-opengl')
saveas(h11,'Evoked Half-Width','fig') 
print(h12,'Evoked Half-Width STD','-dpng','-r600','-opengl')
saveas(h12,'Evoked Half-Width STD','fig') 
saveas(h13,'Peak 10% Latency','fig') 
print(h13,'Peak 10% Latency','-dpng','-r600','-opengl')
saveas(h14,'Peak 10% Latency STD','fig') 
print(h14,'Peak 10% Latency STD','-dpng','-r600','-opengl')
saveas(h15,'Peak 90% Latency','fig') 
print(h15,'Peak 90% Latency','-dpng','-r600','-opengl')
saveas(h16,'Peak 90% Latency STD','fig') 
print(h16,'Peak 90% Latency STD','-dpng','-r600','-opengl')
saveas(h17,'Rise Time 10-90% Peak','fig') 
print(h17,'Rise Time 10-90% Peak','-dpng','-r600','-opengl')
saveas(h18,'Rise Time 10-90% Peak STD','fig') 
print(h18,'Rise Time 10-90% Peak STD','-dpng','-r600','-opengl')


print(g1,'Evoked Adaptation Amplitude Ratio_paired plot','-dpng','-r600','-opengl')
saveas(g1,'Evoked Adaptation Amplitude Ratio_paired_plot','fig') 

filename='Evoked activity event detection'; 
save(filename, 'files_to_analyze', 'event_evoked', 'event_evoked_stat','stat_failures','stat_nonspecific_count','stat_onset',...
    'stat_onset_std','stat_amp','stat_amp_std','stat_ampVal','stat_onVal','stat_ampDel','stat_ampDel_std','stat_halfWidth','stat_halfWidth_std',...
    'stat_peak10prcntTime', 'stat_peak10prcntTime_std', 'stat_peak90prcntTime', 'stat_peak90prcntTime_std','stat_peak10to90Time', 'stat_peak10to90Time_std')
end
