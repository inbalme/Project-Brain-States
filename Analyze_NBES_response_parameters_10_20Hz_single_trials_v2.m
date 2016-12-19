%% for opening workspace saved 
clear all
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
peaks_for_xls=[]; peak_for_xls_mean=[]; 
save_flag= 0;
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
%                 stim_temp=stim2_X{2}(:,1);
%                 stim2_X=[];
%                 stim2_X{2}=stim_temp; stim2_X{3}=stim_temp;
%                 
%% low-pass filtering below 300Hz  
lp=300;
                for xx=1:3
    for trace= 1:size(data_no_spikes{channel},2)    
            jj=data_no_spikes{channel}(:,trace,xx);
            data_no_spikes{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jj,dt,-1,lp,0,0); 
    end
                end 
              %% 
   for trace_type= 2;%1:2; %1:2;  %1 for spont., 2 for evoked
     interval=[]; 
        x_value=[2:3]; %2:3; %[1,1];
                                          
                 start_sample =stim2_X{x_value(1)}(1,1);  %start with sensory stim
        %          duration = 1; %[sec]
                galvano_nstim = size(stim2_X{2},2); %Param.facade(6);
                galvano_freq = Param.facade(7);
                duration = galvano_nstim./galvano_freq+0.1;
%                 end_sample = start_sample+duration.*sf{1}-1;
                end_sample = stim2_X{x_value(1)}(1,end)+0.1.*sf{1}-1;
                interval(:,1) = round(start_sample:end_sample);
                interval(:,2) = interval(:,1);
                        finalAmp_Thres = -1 ; %negative values are interpreted in terms of std. positive values are fixed amplitude in mV       
%% detect events
% for i=1:10;
% voltages_input = data_no_spikes{channel}(i/dt:(i+1)/dt,1,2);
        for t=1:2;
          for trace=1:size(data_no_spikes{1},2)             
              for stim_num=1:11; %1:galvano_nstim;
                  count=0;
                  event_on_nonspecific{t,stim_num,trace}=[]; event_amplitude_nonspecific{t,stim_num,trace}=[]; event_ampPos_nonspecific{t,stim_num,trace}=[];
                  event_halfWidth_nonspecific{t,stim_num,trace}=[]; event_halfWidthS_nonspecific{t,stim_num,trace}=[]; event_halfWidthE_nonspecific{t,stim_num,trace}=[];
                  event_start_nonspecific{t,stim_num,trace}=[]; event_ampDel_nonspecific{t,stim_num,trace}=[];
                  stim_ISI =1/galvano_freq.*sf{1};
                  if stim_num<=galvano_nstim;
                      peak_start_int = stim2_X{2}(1,stim_num)+0.003*sf{1};
                      peak_end_int = stim2_X{2}(1,stim_num)+2*stim_ISI-0.003*sf{1};    
                  end
%                   peak_start_int = interval(1+stim_ISI*(stim_num-1)+0.003*sf{1},t); %not suitable for the test stim
%                   peak_end_int = interval(stim_ISI*(stim_num+1),t);
%                   peak_start_int = interval(1,t);
%                   peak_end_int = interval(end,t);
        voltages_input = data_no_spikes{1}(peak_start_int:peak_end_int,trace,x_value(t));
%         finalAmp_Thres = 1 ;
        doPlot = 0;
        I_temp=0;
        [voltages, starting, amplitude, ampPos, halfWidth,halfWidthS, halfWidthE] = EventDetector_v2(voltages_input, dt, finalAmp_Thres, doPlot,I_temp);
        title(['t=', num2str(t), ' trace ',num2str(trace),' stim. ',num2str(stim_num)]);  
%  figure(1); clf
%             hold on
%                      h1=plot([1:size(voltages_input,1)].*dt, voltages_input,'k');
%                      h2=scatter(starting*dt,voltages_input(starting),'r','fill'); %mark event onset
%                      h3=scatter(ampPos*dt,voltages_input(ampPos),'b','fill'); %mark event peak
% %                 set(gca,'ylim',[-50 -20]); set(gca,'xlim',[0 1.1]);  
%           hold off
%         plot whisker stim
%         ylim_data=[get(gca,'ylim')]';
%         patch_xdata=[stim2_X{x_value(2)}(:,1:10); flipud(stim2_X{x_value(2)}(:,1:10))];
%         patch_xdata = (patch_xdata-interval(1,t))*dt;
%         yex=wextend('1D','sym',ylim_data,1);
%         l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
%         temp_y=wextend('ac','sym',yex,l);
%         patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
%         patch_cdata=ones(size(patch_xdata));
%         p=patch(patch_xdata,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
%         set(gca,'linewidth',1.2)
%           hold off           
%          pause
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
    next_stim_pos = stim2_X{x_value(t)}(1,1)+stim_ISI*(stim_num);
    if galvano_nstim==11;
        next_stim_pos = stim2_X{x_value(t)}(1,stim_num)+stim_ISI;
    end
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
%              [amplitude(tmp),ampPos(tmp)]=max(data_no_spikes{1}(starting(tmp): halfWidthE(tmp),trace,x_value(t))); 
%       %need to recalculate the halfwidth according to the new amplitude. I leave it for now...
%         end

         %if an event ends after a second event has started -> take halfWidthE to be the starting point of
         %the next event:
         tmp=find(starting(2:end)-halfWidthE(1:end-1)<0);
         if ~isempty(tmp)
              halfWidthE(tmp)=starting(tmp+1);
         end  
       halfWidth=(halfWidthE-halfWidthS).*dt.*1000; %in msec 
       if stim_num>size(stim2_X{x_value(t)},2)
           starting=[];
       end
        if isempty(starting);
            failures{t,stim_num,trace}=1;
            if stim_num>galvano_nstim
               failures{t,stim_num,trace}=nan;
            end
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
        elseif starting(1)>(peak_start_int+0.03*sf{1})
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
            event_start_nonspecific{t,stim_num,trace}= starting(1);
            event_on_nonspecific{t,stim_num,trace}= (starting(1)-stim2_X{x_value(t)}(1,stim_num))*dt*1000; %in msec
            event_amplitude_nonspecific{t,stim_num,trace}= amplitude(1);
            event_ampPos_nonspecific{t,stim_num,trace}= ampPos(1);
            event_ampDel_nonspecific{t,stim_num,trace}= (ampPos(1)-stim2_X{x_value(t)}(1,stim_num))*dt*1000; %in msec
            event_halfWidth_nonspecific{t,stim_num,trace} =  halfWidth(1); 
            event_halfWidthS_nonspecific{t,stim_num,trace} =  halfWidthS(1); 
            event_halfWidthE_nonspecific{t,stim_num,trace} =  halfWidthE(1); 
            count = count+1;

        else
            failures{t,stim_num,trace}=0;
            event_start{t,stim_num,trace} =  starting(1);
            event_on{t,stim_num,trace} =  (starting(1)-stim2_X{x_value(t)}(1,stim_num))*dt*1000; %in msec
            event_onVal{t,stim_num,trace} = data_no_spikes{1}(starting(1),trace,x_value(t));
            event_amplitude{t,stim_num,trace} =  amplitude(1); 
            event_ampPos{t,stim_num,trace} =  ampPos(1); 
            event_ampDel{t,stim_num,trace} =  (ampPos(1)-stim2_X{x_value(t)}(1,stim_num))*dt*1000;  %in msec
            event_ampVal{t,stim_num,trace} = data_no_spikes{1}(ampPos(1),trace,x_value(t));
            event_halfWidth{t,stim_num,trace} =  halfWidth(1); 
            event_halfWidthS{t,stim_num,trace} =  halfWidthS(1); 
            event_halfWidthE{t,stim_num,trace} =  halfWidthE(1); 
        end
        if length(starting)>1
            event_start_nonspecific{t,stim_num,trace}=[event_start_nonspecific{t,stim_num,trace}, starting(2:end)];
            event_on_nonspecific{t,stim_num,trace}=[event_on_nonspecific{t,stim_num,trace}, (starting(2:end)-stim2_X{x_value(t)}(1,stim_num)).*dt*1000]; %in msec
            event_amplitude_nonspecific{t,stim_num,trace}= [event_amplitude_nonspecific{t,stim_num,trace}, amplitude(2:end)];
            event_ampPos_nonspecific{t,stim_num,trace}= [event_ampPos_nonspecific{t,stim_num,trace}, ampPos(2:end)]; 
            event_ampDel_nonspecific{t,stim_num,trace}= [event_ampDel_nonspecific{t,stim_num,trace}, (ampPos(2:end)-stim2_X{x_value(t)}(1,stim_num)).*dt*1000]; %in msec
            event_halfWidth_nonspecific{t,stim_num,trace} =  [event_halfWidth_nonspecific{t,stim_num,trace}, halfWidth(2:end)]; 
            event_halfWidthS_nonspecific{t,stim_num,trace} =  [event_halfWidthS_nonspecific{t,stim_num,trace}, halfWidthS(2:end)]; 
            event_halfWidthE_nonspecific{t,stim_num,trace} =  [event_halfWidthE_nonspecific{t,stim_num,trace}, halfWidthE(2:end)]; 
            count = count+length(starting(2:end));
        end
        event_nonspecific_count{t,stim_num,trace}=count;
%    figure(2); %stops whenever there is nan... need to fix
%             hold on
%                      h1=plot([1:size(data_no_spikes{channel},1)].*dt, data_no_spikes{channel}(:,trace,x_value(t)),'k');
%                      h2=scatter(event_start{t,stim_num,trace}(:)*dt,data_no_spikes{channel}(event_start{t,stim_num,trace}(:),trace,x_value(t)),'r','fill'); %mark event onset
%                      h3=scatter(event_ampPos{t,stim_num,trace}(:)*dt,data_no_spikes{channel}(event_ampPos{t,stim_num,trace}(:),trace,x_value(t)),'b','fill'); %mark event peak
%                      hold off
% pause
              end
           end
        end
   end
        %organizing the data in structure:     
            event_evoked(fileind).cells=files_to_analyze;
            event_evoked(fileind).analysis_mfile='Analyze_NBES_response_parameters_10_20Hz_single_trials_v2.m';
            event_evoked(fileind).lp=lp;
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
            event_evoked(fileind).nonspecific_count =event_nonspecific_count;
            event_evoked(fileind).onset_nonspecific = event_on_nonspecific;
            event_evoked(fileind).amplitude_nonspecific = event_amplitude_nonspecific;
            event_evoked(fileind).ampPos_nonspecific = event_ampPos_nonspecific;
            event_evoked(fileind).ampDel_nonspecific = event_ampDel_nonspecific;
            event_evoked(fileind).halfWidth_nonspecific = event_halfWidth_nonspecific; 
            event_evoked(fileind).halfWidthS_nonspecific = event_halfWidthS_nonspecific; 
            event_evoked(fileind).halfWidthE_nonspecific = event_halfWidthE_nonspecific; 
            
            %% cell stats:
onset_mat=cell2mat(event_on); 
onVal_mat=cell2mat(event_onVal); 
amplitude_mat=cell2mat(event_amplitude); 
ampDel_mat=cell2mat(event_ampDel); 
ampVal_mat=cell2mat(event_ampVal); 
halfWidth_mat=cell2mat(event_halfWidth); 
failures_mat=cell2mat(failures); 
nonspecific_count_mat = cell2mat(event_nonspecific_count);
amplitude_mat_M=nanmean(amplitude_mat,3);
adapt_amp=nanmean(amplitude_mat_M(:,8:10),2)./nanmean(amplitude_mat_M(:,1),2);
event_evoked(fileind).adapt_amp=adapt_amp';

tmp=[];

    for stim_num=1:11;      
%adaptaion amplitude ratio
 event_evoked_stat.stim_num(stim_num).adapt_amp(fileind,:)=event_evoked(fileind).adapt_amp;
     %failures
        tmp(:,:) = failures_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).failures=tmp';               
        event_evoked_stat.fileind(fileind).stim_num(stim_num).failures_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).failures,1);
        event_evoked_stat.stim_num(stim_num).failures_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).failures_m;
        clear tmp
        
     %nonspecific events
        tmp(:,:) = nonspecific_count_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count=tmp';               
        event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count,1);
        event_evoked_stat.stim_num(stim_num).nonspecific_count_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).nonspecific_count_m;
        clear tmp
        
    %onset
        tmp(:,:) = onset_mat(1:2,stim_num,:);  
%         if sum(~isnan(tmp(1,:)))<5 | sum(~isnan(tmp(2,:)))<5
             if sum(~isnan(tmp(1,(~isnan(tmp(2,:))==1))))<5
            continue
        end
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onset=tmp';       
        norm_peak_onset_noES(:,1)= tmp(1,:)./tmp(1,:);
        norm_peak_onset_ES(:,1)= tmp(2,:)./tmp(1,:);
%         event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_onset = [norm_peak_onset_noES, norm_peak_onset_ES];
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset,0,1);  
        event_evoked_stat.stim_num(stim_num).onset_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_m;
        event_evoked_stat.stim_num(stim_num).onset_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onset_std;
        %testing for normal distribution
        diff_onset= tmp(2,:)- tmp(1,:);
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_onset, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_onset] = lillietest(diff_onset);
        %unpaired ttest 
%         [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_norm_onset, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_norm_onset]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_onset(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_onset(:,2));
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_onset, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_onset]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).onset(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).onset(:,2));

        clear norm_peak_onset_noES norm_peak_onset_ES diff_onset tmp
        
        %onset value
        tmp(:,:) = onVal_mat(1:2,stim_num,:);  
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal=tmp';
        norm_peak_onVal_noES(:,1)= tmp(1,:)./tmp(1,:);
        norm_peak_onVal_ES(:,1)= tmp(2,:)./tmp(1,:);
%         event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_onVal = [norm_peak_onVal_noES, norm_peak_onVal_ES];
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal,0,1);  
        event_evoked_stat.stim_num(stim_num).onVal_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal_m;
        event_evoked_stat.stim_num(stim_num).onVal_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal_std;
        %testing for normal distribution
        diff_onVal= tmp(2,:)- tmp(1,:);
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_onVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_onVal] = lillietest(diff_onVal);
        %unpaired ttest 
%         [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_norm_onVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_norm_onVal]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_onVal(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_onVal(:,2));
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_onVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_onVal]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).onVal(:,2));

        clear norm_peak_onVal_noES norm_peak_onVal_ES diff_onVal tmp
        
    %amplitude
        tmp(:,:) = amplitude_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude=tmp';
        norm_peak_amplitude_noES(:,1)= tmp(1,:)./tmp(1,:);
        norm_peak_amplitude_ES(:,1)= tmp(2,:)./tmp(1,:);
%         event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_amplitude = [norm_peak_amplitude_noES, norm_peak_amplitude_ES];
        event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude,0,1);   
        event_evoked_stat.stim_num(stim_num).amplitude_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude_m;
        event_evoked_stat.stim_num(stim_num).amplitude_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude_std;
        %testing for normal distribution
        diff_amplitude= tmp(2,:)- tmp(1,:);
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_amplitude, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_amplitude] = lillietest(diff_amplitude);
        %paired ttest2 
%         [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_norm_amplitude, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_norm_amplitude]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_amplitude(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_amplitude(:,2));
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_amplitude, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_amplitude]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).amplitude(:,2));

        clear norm_peak_amplitude_noES norm_peak_amplitude_ES diff_amplitude tmp

    %peak value
        tmp(:,:) = ampVal_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal=tmp';
        norm_peak_ampVal_noES(:,1)= tmp(1,:)./tmp(1,:);
        norm_peak_ampVal_ES(:,1)= tmp(2,:)./tmp(1,:);
%         event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_ampVal = [norm_peak_ampVal_noES, norm_peak_ampVal_ES];
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal,0,1);   
        event_evoked_stat.stim_num(stim_num).ampVal_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal_m;
        event_evoked_stat.stim_num(stim_num).ampVal_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal_std;
        %testing for normal distribution
        diff_ampVal= tmp(2,:)- tmp(1,:);
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_ampVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_ampVal] = lillietest(diff_ampVal);
        %paired ttest2 
%         [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_norm_ampVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_norm_ampVal]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_ampVal(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_ampVal(:,2));
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_ampVal, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_ampVal]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).ampVal(:,2));

        clear norm_peak_ampVal_noES norm_peak_ampVal_ES diff_ampVal tmp
        
    %amplitude latency
        tmp(:,:) = ampDel_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel=tmp';
        norm_peak_ampDel_noES(:,1)= tmp(1,:)./tmp(1,:);
        norm_peak_ampDel_ES(:,1)= tmp(2,:)./tmp(1,:);
%         event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_ampDel = [norm_peak_ampDel_noES, norm_peak_ampDel_ES];
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel,0,1);    
        event_evoked_stat.stim_num(stim_num).ampDel_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel_m;
        event_evoked_stat.stim_num(stim_num).ampDel_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel_std;
        %testing for normal distribution
        diff_ampDel= tmp(2,:)- tmp(1,:);
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_ampDel, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_ampDel] = lillietest(diff_ampDel);
        %unpaired ttest 
%         [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_norm_ampDel, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_norm_ampDel]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_ampDel(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_ampDel(:,2));
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_ampDel, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_ampDel]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).ampDel(:,2));

        clear norm_peak_ampDel_noES norm_peak_ampDel_ES diff_ampDel tmp
 
    %half width
        tmp(:,:) = halfWidth_mat(1:2,stim_num,:);        
        event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth=tmp';
        norm_peak_halfWidth_noES(:,1)= tmp(1,:)./tmp(1,:);
        norm_peak_halfWidth_ES(:,1)= tmp(2,:)./tmp(1,:);
%         event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_halfWidth = [norm_peak_halfWidth_noES, norm_peak_halfWidth_ES];
        event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth_m=nanmean(event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth,1);
        event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth_std=nanstd(event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth,0,1);     
        event_evoked_stat.stim_num(stim_num).halfWidth_m(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth_m;
        event_evoked_stat.stim_num(stim_num).halfWidth_std(fileind,:)= event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth_std;
        %testing for normal distribution
        diff_halfWidth= tmp(2,:)- tmp(1,:);
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_h_halfWidth, event_evoked_stat.fileind(fileind).stim_num(stim_num).lillietest_p_halfWidth] = lillietest(diff_halfWidth);
        %unpaired ttest2 
%         [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_norm_halfWidth, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_norm_halfWidth]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_halfWidth(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).norm_halfWidth(:,2));
        [event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_h_halfWidth, event_evoked_stat.fileind(fileind).stim_num(stim_num).ttest2_p_halfWidth]= ttest2(event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth(:,1),event_evoked_stat.fileind(fileind).stim_num(stim_num).halfWidth(:,2));

        clear norm_peak_halfWidth_noES norm_peak_halfWidth_ES diff_halfWidth tmp
    end
    if galvano_nstim==11;
            stim11(fileind,1)=1;
    else stim11(fileind,1)=0;
    end
    end
 %% population statistics
 clear tmp
 for stim_num=1:11;
      if stim_num==11;           
            event_evoked_stat.stim_num(stim_num).nonspecific_count_m(stim11==0,:)=nan; 
            event_evoked_stat.stim_num(stim_num).onset_std(stim11==0,:)=nan;
            event_evoked_stat.stim_num(stim_num).onset_m(stim11==0,:)=nan;
            event_evoked_stat.stim_num(stim_num).onset_std(stim11==0,:)=nan;
             event_evoked_stat.stim_num(stim_num).failures_m(stim11==0,:)=nan;            
             event_evoked_stat.stim_num(stim_num).onVal_m(stim11==0,:)=nan;
             event_evoked_stat.stim_num(stim_num).amplitude_m(stim11==0,:)=nan;
             event_evoked_stat.stim_num(stim_num).amplitude_std(stim11==0,:)=nan;
             event_evoked_stat.stim_num(stim_num).ampVal_m(stim11==0,:)=nan;
             event_evoked_stat.stim_num(stim_num).ampDel_m(stim11==0,:)=nan;
             event_evoked_stat.stim_num(stim_num).ampDel_std(stim11==0,:)=nan;
             event_evoked_stat.stim_num(stim_num).halfWidth_m(stim11==0,:)=nan;
             event_evoked_stat.stim_num(stim_num).halfWidth_std(stim11==0,:)=nan;
        end
    %average number of failures per stim
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).failures_m;           
        event_evoked_stat.stim_num(stim_num).change_failures_m = (tmp(:,2)-tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).failures_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).failures_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_failures_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_failures_m,1);
        event_evoked_stat.stim_num(stim_num).change_failures_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_failures_m,0,1);
       change_failures_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_failures_m;
        %testing for normal distribution
        diff_failures_m= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_failures_m, event_evoked_stat.stim_num(stim_num).lillietest_p_failures_m] = lillietest(diff_failures_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_failures_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_failures_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_failures_m);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_failures_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_failures_m]= ttest(event_evoked_stat.stim_num(stim_num).change_failures_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_failures_m, event_evoked_stat.stim_num(stim_num).ttest_p_failures_m]= ttest(event_evoked_stat.stim_num(stim_num).failures_m(:,1),event_evoked_stat.stim_num(stim_num).failures_m(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_failures_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_failures_m]= signrank(event_evoked_stat.stim_num(stim_num).change_failures_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_failures_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_failures_m]= signrank(event_evoked_stat.stim_num(stim_num).failures_m(:,1),event_evoked_stat.stim_num(stim_num).failures_m(:,2));

        clear  diff_failures_m tmp
        
        %average number of nonspecific events per stim
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).nonspecific_count_m;           
        event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m = (tmp(:,2)-tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).nonspecific_count_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).nonspecific_count_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m,1);
        event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m,0,1);
       change_nonspecific_count_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m;
        %testing for normal distribution
        diff_nonspecific_count_m= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).lillietest_p_nonspecific_count_m] = lillietest(diff_nonspecific_count_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_nonspecific_count_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_nonspecific_count_m]= ttest(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).ttest_p_nonspecific_count_m]= ttest(event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,1),event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_nonspecific_count_m]= signrank(event_evoked_stat.stim_num(stim_num).change_nonspecific_count_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_nonspecific_count_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_nonspecific_count_m]= signrank(event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,1),event_evoked_stat.stim_num(stim_num).nonspecific_count_m(:,2));

        clear  diff_nonspecific_count_m tmp
        
    % onset   
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).onset_m;           
        event_evoked_stat.stim_num(stim_num).change_onset_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).onset_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).onset_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_onset_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_onset_m,1);
        event_evoked_stat.stim_num(stim_num).change_onset_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_onset_m,0,1);
        change_onset_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_onset_m;
        %testing for normal distribution
        diff_onset_m= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_onset_m, event_evoked_stat.stim_num(stim_num).lillietest_p_onset_m] = lillietest(diff_onset_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_onset_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_onset_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_onset_m);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_onset_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_onset_m]= ttest(event_evoked_stat.stim_num(stim_num).change_onset_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_onset_m, event_evoked_stat.stim_num(stim_num).ttest_p_onset_m]= ttest(event_evoked_stat.stim_num(stim_num).onset_m(:,1),event_evoked_stat.stim_num(stim_num).onset_m(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_onset_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onset_m]= signrank(event_evoked_stat.stim_num(stim_num).change_onset_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_onset_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_onset_m]= signrank(event_evoked_stat.stim_num(stim_num).onset_m(:,1),event_evoked_stat.stim_num(stim_num).onset_m(:,2));

        clear  diff_onset_m tmp
 
    % onset std
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).onset_std;           
        event_evoked_stat.stim_num(stim_num).change_onset_std = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).onset_std_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).onset_std_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_onset_std_m= nanmean(event_evoked_stat.stim_num(stim_num).change_onset_std,1);
        event_evoked_stat.stim_num(stim_num).change_onset_std_std= nanstd(event_evoked_stat.stim_num(stim_num).change_onset_std,0,1);
        change_onset_std_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_onset_std;
        %testing for normal distribution
        diff_onset_std= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_onset_std, event_evoked_stat.stim_num(stim_num).lillietest_p_onset_std] = lillietest(diff_onset_std);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_onset_std, event_evoked_stat.stim_num(stim_num).lillietest_p_change_onset_std] = lillietest(event_evoked_stat.stim_num(stim_num).change_onset_std);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_onset_std, event_evoked_stat.stim_num(stim_num).ttest_p_change_onset_std]= ttest(event_evoked_stat.stim_num(stim_num).change_onset_std);
        [event_evoked_stat.stim_num(stim_num).ttest_h_onset_std, event_evoked_stat.stim_num(stim_num).ttest_p_onset_std]= ttest(event_evoked_stat.stim_num(stim_num).onset_std(:,1),event_evoked_stat.stim_num(stim_num).onset_std(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_onset_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onset_std]= signrank(event_evoked_stat.stim_num(stim_num).change_onset_std);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_onset_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_onset_std]= signrank(event_evoked_stat.stim_num(stim_num).onset_std(:,1),event_evoked_stat.stim_num(stim_num).onset_std(:,2));

        clear  diff_onset_std tmp
        
        % onset value
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).onVal_m;           
        event_evoked_stat.stim_num(stim_num).change_onVal_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).onVal_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).onVal_m_std= nanstd(tmp,0,1);
        event_evoked_stat.stim_num(stim_num).change_onVal_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_onVal_m,1);
        event_evoked_stat.stim_num(stim_num).change_onVal_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_onVal_m,0,1);
        change_onVal_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_onVal_m;
        %testing for normal distribution
        diff_onVal_m= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_onVal_m, event_evoked_stat.stim_num(stim_num).lillietest_p_onVal_m] = lillietest(diff_onVal_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_onVal_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_onVal_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_onVal_m);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_onVal_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_onVal_m]= ttest(event_evoked_stat.stim_num(stim_num).change_onVal_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_onVal_m, event_evoked_stat.stim_num(stim_num).ttest_p_onVal_m]= ttest(event_evoked_stat.stim_num(stim_num).onVal_m(:,1),event_evoked_stat.stim_num(stim_num).onVal_m(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_onVal_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onVal_m]= signrank(event_evoked_stat.stim_num(stim_num).change_onVal_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_onVal_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_onVal_m]= signrank(event_evoked_stat.stim_num(stim_num).onVal_m(:,1),event_evoked_stat.stim_num(stim_num).onVal_m(:,2));

        clear  diff_onVal_m tmp
        
     % amplitude mean
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).amplitude_m;           
        event_evoked_stat.stim_num(stim_num).change_amplitude_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).amplitude_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).amplitude_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_amplitude_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_amplitude_m,1);
        event_evoked_stat.stim_num(stim_num).change_amplitude_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_amplitude_m,0,1);
       change_amplitude_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_amplitude_m;
        %testing for normal distribution
        diff_amplitude_m= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_amplitude_m, event_evoked_stat.stim_num(stim_num).lillietest_p_amplitude_m] = lillietest(diff_amplitude_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_amplitude_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_amplitude_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_amplitude_m);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_amplitude_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_amplitude_m]= ttest(event_evoked_stat.stim_num(stim_num).change_amplitude_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_amplitude_m, event_evoked_stat.stim_num(stim_num).ttest_p_amplitude_m]= ttest(event_evoked_stat.stim_num(stim_num).amplitude_m(:,1),event_evoked_stat.stim_num(stim_num).amplitude_m(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_amplitude_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_amplitude_m]= signrank(event_evoked_stat.stim_num(stim_num).change_amplitude_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_amplitude_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_amplitude_m]= signrank(event_evoked_stat.stim_num(stim_num).amplitude_m(:,1),event_evoked_stat.stim_num(stim_num).amplitude_m(:,2));

        clear  diff_amplitude_m tmp
        
    % amplitude std
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).amplitude_std;           
        event_evoked_stat.stim_num(stim_num).change_amplitude_std = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).amplitude_std_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).amplitude_std_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_amplitude_std_m= nanmean(event_evoked_stat.stim_num(stim_num).change_amplitude_std,1);
        event_evoked_stat.stim_num(stim_num).change_amplitude_std_std= nanstd(event_evoked_stat.stim_num(stim_num).change_amplitude_std,0,1);
        change_amplitude_std_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_amplitude_std;
        %testing for normal distribution
        diff_amplitude_std= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_amplitude_std, event_evoked_stat.stim_num(stim_num).lillietest_p_amplitude_std] = lillietest(diff_amplitude_std);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_amplitude_std, event_evoked_stat.stim_num(stim_num).lillietest_p_change_amplitude_std] = lillietest(event_evoked_stat.stim_num(stim_num).change_amplitude_std);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_amplitude_std, event_evoked_stat.stim_num(stim_num).ttest_p_change_amplitude_std]= ttest(event_evoked_stat.stim_num(stim_num).change_amplitude_std);
        [event_evoked_stat.stim_num(stim_num).ttest_h_amplitude_std, event_evoked_stat.stim_num(stim_num).ttest_p_amplitude_std]= ttest(event_evoked_stat.stim_num(stim_num).amplitude_std(:,1),event_evoked_stat.stim_num(stim_num).amplitude_std(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_amplitude_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_amplitude_std]= signrank(event_evoked_stat.stim_num(stim_num).change_amplitude_std);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_amplitude_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_amplitude_std]= signrank(event_evoked_stat.stim_num(stim_num).amplitude_std(:,1),event_evoked_stat.stim_num(stim_num).amplitude_std(:,2));

        clear  diff_amplitude_std tmp
        
        % peak value
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).ampVal_m;           
        event_evoked_stat.stim_num(stim_num).change_ampVal_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).ampVal_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).ampVal_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_ampVal_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_ampVal_m,1);
        event_evoked_stat.stim_num(stim_num).change_ampVal_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_ampVal_m,0,1);
        change_ampVal_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_ampVal_m;
        %testing for normal distribution
        diff_ampVal_m= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_ampVal_m, event_evoked_stat.stim_num(stim_num).lillietest_p_ampVal_m] = lillietest(diff_ampVal_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampVal_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_ampVal_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_ampVal_m);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_ampVal_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_ampVal_m]= ttest(event_evoked_stat.stim_num(stim_num).change_ampVal_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_ampVal_m, event_evoked_stat.stim_num(stim_num).ttest_p_ampVal_m]= ttest(event_evoked_stat.stim_num(stim_num).ampVal_m(:,1),event_evoked_stat.stim_num(stim_num).ampVal_m(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_ampVal_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampVal_m]= signrank(event_evoked_stat.stim_num(stim_num).change_ampVal_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_ampVal_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampVal_m]= signrank(event_evoked_stat.stim_num(stim_num).ampVal_m(:,1),event_evoked_stat.stim_num(stim_num).ampVal_m(:,2));

        clear  diff_onset_m tmp
        
    % amplitude latency mean
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).ampDel_m;           
        event_evoked_stat.stim_num(stim_num).change_ampDel_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).ampDel_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).ampDel_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_ampDel_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_ampDel_m,1);
        event_evoked_stat.stim_num(stim_num).change_ampDel_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_ampDel_m,0,1);
        change_ampDel_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_ampDel_m;
        %testing for normal distribution
        diff_ampDel_m= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_ampDel_m, event_evoked_stat.stim_num(stim_num).lillietest_p_ampDel_m] = lillietest(diff_ampDel_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampDel_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_ampDel_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_ampDel_m);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_ampDel_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_ampDel_m]= ttest(event_evoked_stat.stim_num(stim_num).change_ampDel_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_ampDel_m, event_evoked_stat.stim_num(stim_num).ttest_p_ampDel_m]= ttest(event_evoked_stat.stim_num(stim_num).ampDel_m(:,1),event_evoked_stat.stim_num(stim_num).ampDel_m(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_ampDel_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampDel_m]= signrank(event_evoked_stat.stim_num(stim_num).change_ampDel_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_ampDel_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampDel_m]= signrank(event_evoked_stat.stim_num(stim_num).ampDel_m(:,1),event_evoked_stat.stim_num(stim_num).ampDel_m(:,2));

        clear  diff_ampDel_m tmp
        
     % amplitude latency std
         tmp(:,:) = event_evoked_stat.stim_num(stim_num).ampDel_std;           
        event_evoked_stat.stim_num(stim_num).change_ampDel_std = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).ampDel_std_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).ampDel_std_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_ampDel_std_m= nanmean(event_evoked_stat.stim_num(stim_num).change_ampDel_std,1);
        event_evoked_stat.stim_num(stim_num).change_ampDel_std_std= nanstd(event_evoked_stat.stim_num(stim_num).change_ampDel_std,0,1);
        change_ampDel_std_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_ampDel_std;
        %testing for normal distribution
        diff_ampDel_std= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_ampDel_std, event_evoked_stat.stim_num(stim_num).lillietest_p_ampDel_std] = lillietest(diff_ampDel_std);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampDel_std, event_evoked_stat.stim_num(stim_num).lillietest_p_change_ampDel_std] = lillietest(event_evoked_stat.stim_num(stim_num).change_ampDel_std);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_ampDel_std, event_evoked_stat.stim_num(stim_num).ttest_p_change_ampDel_std]= ttest(event_evoked_stat.stim_num(stim_num).change_ampDel_std);
        [event_evoked_stat.stim_num(stim_num).ttest_h_ampDel_std, event_evoked_stat.stim_num(stim_num).ttest_p_ampDel_std]= ttest(event_evoked_stat.stim_num(stim_num).ampDel_std(:,1),event_evoked_stat.stim_num(stim_num).ampDel_std(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_ampDel_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampDel_std]= signrank(event_evoked_stat.stim_num(stim_num).change_ampDel_std);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_ampDel_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampDel_std]= signrank(event_evoked_stat.stim_num(stim_num).ampDel_std(:,1),event_evoked_stat.stim_num(stim_num).ampDel_std(:,2));

        clear  diff_ampDel_std tmp
        
    % half width mean
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).halfWidth_m;           
        event_evoked_stat.stim_num(stim_num).change_halfWidth_m = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).halfWidth_m_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).halfWidth_m_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_halfWidth_m_m= nanmean(event_evoked_stat.stim_num(stim_num).change_halfWidth_m,1);
        event_evoked_stat.stim_num(stim_num).change_halfWidth_m_std= nanstd(event_evoked_stat.stim_num(stim_num).change_halfWidth_m,0,1);
        change_halfWidth_m_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_halfWidth_m;
        %testing for normal distribution
        diff_halfWidth_m= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_halfWidth_m, event_evoked_stat.stim_num(stim_num).lillietest_p_halfWidth_m] = lillietest(diff_halfWidth_m);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_halfWidth_m, event_evoked_stat.stim_num(stim_num).lillietest_p_change_halfWidth_m] = lillietest(event_evoked_stat.stim_num(stim_num).change_halfWidth_m);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_halfWidth_m, event_evoked_stat.stim_num(stim_num).ttest_p_change_halfWidth_m]= ttest(event_evoked_stat.stim_num(stim_num).change_halfWidth_m);
        [event_evoked_stat.stim_num(stim_num).ttest_h_halfWidth_m, event_evoked_stat.stim_num(stim_num).ttest_p_halfWidth_m]= ttest(event_evoked_stat.stim_num(stim_num).halfWidth_m(:,1),event_evoked_stat.stim_num(stim_num).halfWidth_m(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_halfWidth_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_halfWidth_m]= signrank(event_evoked_stat.stim_num(stim_num).change_halfWidth_m);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_halfWidth_m, event_evoked_stat.stim_num(stim_num).wilcoxon_h_halfWidth_m]= signrank(event_evoked_stat.stim_num(stim_num).halfWidth_m(:,1),event_evoked_stat.stim_num(stim_num).halfWidth_m(:,2));

        clear  diff_halfWidth_m tmp
        
    % half width std
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).halfWidth_std;           
        event_evoked_stat.stim_num(stim_num).change_halfWidth_std = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).halfWidth_std_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).halfWidth_std_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_halfWidth_std_m= nanmean(event_evoked_stat.stim_num(stim_num).change_halfWidth_std,1);
        event_evoked_stat.stim_num(stim_num).change_halfWidth_std_std= nanstd(event_evoked_stat.stim_num(stim_num).change_halfWidth_std,0,1);
        change_halfWidth_std_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_halfWidth_std;
        %testing for normal distribution
        diff_halfWidth_std= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_halfWidth_std, event_evoked_stat.stim_num(stim_num).lillietest_p_halfWidth_std] = lillietest(diff_halfWidth_std);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_halfWidth_std, event_evoked_stat.stim_num(stim_num).lillietest_p_change_halfWidth_std] = lillietest(event_evoked_stat.stim_num(stim_num).change_halfWidth_std);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_halfWidth_std, event_evoked_stat.stim_num(stim_num).ttest_p_change_halfWidth_std]= ttest(event_evoked_stat.stim_num(stim_num).change_halfWidth_std);
        [event_evoked_stat.stim_num(stim_num).ttest_h_halfWidth_std, event_evoked_stat.stim_num(stim_num).ttest_p_halfWidth_std]= ttest(event_evoked_stat.stim_num(stim_num).halfWidth_std(:,1),event_evoked_stat.stim_num(stim_num).halfWidth_std(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_halfWidth_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_halfWidth_std]= signrank(event_evoked_stat.stim_num(stim_num).change_halfWidth_std);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_halfWidth_std, event_evoked_stat.stim_num(stim_num).wilcoxon_h_halfWidth_std]= signrank(event_evoked_stat.stim_num(stim_num).halfWidth_std(:,1),event_evoked_stat.stim_num(stim_num).halfWidth_std(:,2));

        clear  diff_halfWidth_std tmp
        
        %adaptation amplitude ratio
        tmp(:,:) = event_evoked_stat.stim_num(stim_num).adapt_amp;           
        event_evoked_stat.stim_num(stim_num).change_adapt_amp = (tmp(:,2)-tmp(:,1))./abs(tmp(:,1)).*100 ;
        event_evoked_stat.stim_num(stim_num).adapt_amp_m= nanmean(tmp,1);
        event_evoked_stat.stim_num(stim_num).adapt_amp_std= nanstd(tmp,0,1);
         event_evoked_stat.stim_num(stim_num).change_adapt_amp_m= nanmean(event_evoked_stat.stim_num(stim_num).change_adapt_amp,1);
        event_evoked_stat.stim_num(stim_num).change_adapt_amp_std= nanstd(event_evoked_stat.stim_num(stim_num).change_adapt_amp,0,1);
        change_adapt_amp_mat(:,stim_num)= event_evoked_stat.stim_num(stim_num).change_adapt_amp;
        %testing for normal distribution
        diff_adapt_amp= tmp(:,2)- tmp(:,1);
        [event_evoked_stat.stim_num(stim_num).lillietest_h_adapt_amp, event_evoked_stat.stim_num(stim_num).lillietest_p_adapt_amp] = lillietest(diff_adapt_amp);
         [event_evoked_stat.stim_num(stim_num).lillietest_h_change_adapt_amp, event_evoked_stat.stim_num(stim_num).lillietest_p_change_adapt_amp] = lillietest(event_evoked_stat.stim_num(stim_num).change_adapt_amp);
        %paired ttest 
        [event_evoked_stat.stim_num(stim_num).ttest_h_change_adapt_amp, event_evoked_stat.stim_num(stim_num).ttest_p_change_adapt_amp]= ttest(event_evoked_stat.stim_num(stim_num).change_adapt_amp);
        [event_evoked_stat.stim_num(stim_num).ttest_h_adapt_amp, event_evoked_stat.stim_num(stim_num).ttest_p_adapt_amp]= ttest(event_evoked_stat.stim_num(stim_num).adapt_amp(:,1),event_evoked_stat.stim_num(stim_num).adapt_amp(:,2));
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_change_adapt_amp, event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_adapt_amp]= signrank(event_evoked_stat.stim_num(stim_num).change_adapt_amp);
        [event_evoked_stat.stim_num(stim_num).wilcoxon_p_adapt_amp, event_evoked_stat.stim_num(stim_num).wilcoxon_h_adapt_amp]= signrank(event_evoked_stat.stim_num(stim_num).adapt_amp(:,1),event_evoked_stat.stim_num(stim_num).adapt_amp(:,2));

        clear  diff_adapt_amp tmp
        
        if event_evoked_stat.stim_num(stim_num).lillietest_h_failures_m==0      
            h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_failures_m; 
  else h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_failures_m; 
        end
     if event_evoked_stat.stim_num(stim_num).lillietest_h_onset_m==0      
            h_onset_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_onset_m; 
  else h_onset_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_onset_m; 
     end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_onset_std==0      
            h_onset_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_onset_std; 
  else h_onset_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_onset_std; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_amplitude_m==0      
            h_amplitude_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_amplitude_m; 
  else h_amplitude_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_amplitude_m; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_amplitude_std==0      
            h_amplitude_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_amplitude_std; 
  else h_amplitude_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_amplitude_std; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_onVal_m==0      
            h_onVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_onVal_m; 
  else h_onVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_onVal_m; 
   end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_onVal_std==0      
%             h_onVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_onVal_std; 
%   else h_onVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_onVal_std; 
%    end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_ampVal_m==0      
            h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_ampVal_m; 
  else h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampVal_m; 
   end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_ampVal_std==0      
%             h_ampVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_ampVal_std; 
%   else h_ampVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampVal_std; 
%    end
 if event_evoked_stat.stim_num(stim_num).lillietest_h_ampDel_m==0      
            h_ampDel_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_ampDel_m; 
  else h_ampDel_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampDel_m; 
 end
    if event_evoked_stat.stim_num(stim_num).lillietest_h_ampDel_std==0      
            h_ampDel_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_ampDel_std; 
  else h_ampDel_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_ampDel_std; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_halfWidth_m==0      
            h_halfWidth_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_halfWidth_m; 
  else h_halfWidth_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_halfWidth_m; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_halfWidth_std==0      
            h_halfWidth_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_halfWidth_std; 
  else h_halfWidth_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_halfWidth_std; 
   end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_adapt_amp==0      
            h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_adapt_amp; 
  else h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_adapt_amp; 
   end
 if event_evoked_stat.stim_num(stim_num).lillietest_h_nonspecific_count_m==0      
            h_nonspecific_count_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_nonspecific_count_m; 
  else h_nonspecific_count_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_nonspecific_count_m; 
   end
    end
    %getting the x-axis (stimulus number) and y-axis for the asterisks of significance
     s_failures_m=find(h_failures_m==1);
     s_onset_m=find(h_onset_m==1);
     s_onset_std=find(h_onset_std==1); 
     s_amplitude_m=find(h_amplitude_m==1);
     s_amplitude_std=find(h_amplitude_std==1);
     s_onVal_m=find(h_onVal_m==1);
%      s_onVal_std=find(h_onVal_std==1);
     s_ampVal_m=find(h_ampVal_m==1);
%      s_ampVal_std=find(h_ampVal_std==1);
     s_ampDel_m=find(h_ampDel_m==1);
     s_ampDel_std=find(h_ampDel_std==1);
     s_halfWidth_m=find(h_halfWidth_m==1);
     s_halfWidth_std=find(h_halfWidth_std==1);
     s_adapt_amp=find(h_adapt_amp==1);
     s_nonspecific_count_m=find(h_nonspecific_count_m==1);

 %% Plot parameters along the train stim - version 3:bars+error bars of the mean values
 close all
     clear color_table
    color_table=[0 0 0; [216 22 22]/256; [136 137 138]/256; [255 153 153]/256];
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
    end 
    
    capincrease=-0.2; %#units to increase length of errorbar cap
    asymmetric=1;
    barwidth1=0.3;
    
    h1=figure;        
        hold on
        errbar_h=errorbar([1:size(failures_m_mat{1}(:,:),2)]-0.3,nanmean(failures_m_mat{1}(:,:),1),zeros(1,size(failures_m_mat{1}(:,:),2)), nanstd(failures_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(failures_m_mat{2}(:,:),2)],nanmean(failures_m_mat{2}(:,:),1),zeros(1,size(failures_m_mat{2}(:,:),2)), nanstd(failures_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(failures_m_mat{1}(:,:),2)]-0.3, nanmean(failures_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(failures_m_mat{2}(:,:),2)], nanmean(failures_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_failures_m=max(ylim_data).*ones(size(s_failures_m));
        plot(s_failures_m-0.15,p_failures_m,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Failures rate' ,'FontSize', 16);    
        title(['Mean Rate of Failures,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
        h2=figure;        
        hold on
        errbar_h=errorbar([1:size(nonspecific_count_m_mat{1}(:,:),2)]-0.3,nanmean(nonspecific_count_m_mat{1}(:,:),1),zeros(1,size(nonspecific_count_m_mat{1}(:,:),2)), nanstd(nonspecific_count_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(nonspecific_count_m_mat{2}(:,:),2)],nanmean(nonspecific_count_m_mat{2}(:,:),1),zeros(1,size(nonspecific_count_m_mat{2}(:,:),2)), nanstd(nonspecific_count_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(nonspecific_count_m_mat{1}(:,:),2)]-0.3, nanmean(nonspecific_count_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(nonspecific_count_m_mat{2}(:,:),2)], nanmean(nonspecific_count_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_nonspecific_count_m=max(ylim_data).*ones(size(s_nonspecific_count_m));
        plot(s_nonspecific_count_m-0.15,p_nonspecific_count_m,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Non-specific response rate' ,'FontSize', 16);    
        title(['Mean Rate of Non-specific responses,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  

        h3=figure;        
        hold on
        errbar_h=errorbar([1:size(onset_m_mat{1}(:,:),2)]-0.3,nanmean(onset_m_mat{1}(:,:),1),zeros(1,size(onset_m_mat{1}(:,:),2)), nanstd(onset_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(onset_m_mat{2}(:,:),2)],nanmean(onset_m_mat{2}(:,:),1),zeros(1,size(onset_m_mat{2}(:,:),2)), nanstd(onset_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(onset_m_mat{1}(:,:),2)]-0.3, nanmean(onset_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(onset_m_mat{2}(:,:),2)], nanmean(onset_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_onset_m=max(ylim_data).*ones(size(s_onset_m));
        plot(s_onset_m-0.15,p_onset_m,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Onset Latency [mS]' ,'FontSize', 16);    
        title(['Mean Onset Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  

        h4=figure;        
        hold on
        errbar_h=errorbar([1:size(onset_std_mat{1}(:,:),2)]-0.3,nanmean(onset_std_mat{1}(:,:),1),zeros(1,size(onset_std_mat{1}(:,:),2)), nanstd(onset_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(onset_std_mat{2}(:,:),2)],nanmean(onset_std_mat{2}(:,:),1),zeros(1,size(onset_std_mat{2}(:,:),2)), nanstd(onset_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(onset_std_mat{1}(:,:),2)]-0.3, nanmean(onset_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(onset_std_mat{2}(:,:),2)], nanmean(onset_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_onset_std=max(ylim_data).*ones(size(s_onset_std));
        plot(s_onset_std-0.15,p_onset_std,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Onset Latency STD [mS]' ,'FontSize', 16);    
        title(['Mean response jitter,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
      h5=figure;      
        hold on      
        errbar_h=errorbar([1:size(amplitude_m_mat{1}(:,:),2)]-0.3,nanmean(amplitude_m_mat{1}(:,:),1),zeros(1,size(amplitude_m_mat{1}(:,:),2)), nanstd(amplitude_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(amplitude_m_mat{2}(:,:),2)],nanmean(amplitude_m_mat{2}(:,:),1),zeros(1,size(amplitude_m_mat{2}(:,:),2)), nanstd(amplitude_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(amplitude_m_mat{1}(:,:),2)]-0.3, nanmean(amplitude_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(amplitude_m_mat{2}(:,:),2)], nanmean(amplitude_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_amplitude_m=max(ylim_data).*ones(size(s_amplitude_m));
        plot(s_amplitude_m-0.15,p_amplitude_m,'k*')
%         line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b') %change line to zero
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Mean Response Amplitude [mV]' ,'FontSize', 16);    
        title(['Mean Response Amplitude - local baseline,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
         h6=figure;        
        hold on
        errbar_h=errorbar([1:size(amplitude_std_mat{1}(:,:),2)]-0.3,nanmean(amplitude_std_mat{1}(:,:),1),zeros(1,size(amplitude_std_mat{1}(:,:),2)), nanstd(amplitude_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(amplitude_std_mat{2}(:,:),2)],nanmean(amplitude_std_mat{2}(:,:),1),zeros(1,size(amplitude_std_mat{2}(:,:),2)), nanstd(amplitude_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(amplitude_std_mat{1}(:,:),2)]-0.3, nanmean(amplitude_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(amplitude_std_mat{2}(:,:),2)], nanmean(amplitude_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_amplitude_std=max(ylim_data).*ones(size(s_amplitude_std));
        plot(s_amplitude_std-0.15,p_amplitude_std,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Amplitude STD [mV]' ,'FontSize', 16);    
        title(['Mean Amplitude STD (trial-to-trial variability) ,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
 
        h7=figure;        
        hold on
        set(gca,'Ydir','reverse'); %added because the membrane pot. values are negative. more changes: error bars are symmetric and caps not removed, and asterisks yvalue is min instead of max
        errbar_h=errorbar([1:size(ampVal_m_mat{1}(:,:),2)]-0.3,nanmean(ampVal_m_mat{1}(:,:),1), nanstd(ampVal_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(ampVal_m_mat{2}(:,:),2)],nanmean(ampVal_m_mat{2}(:,:),1), nanstd(ampVal_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(ampVal_m_mat{1}(:,:),2)]-0.3, nanmean(ampVal_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(ampVal_m_mat{2}(:,:),2)], nanmean(ampVal_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_ampVal_m=min(ylim_data).*ones(size(s_ampVal_m));
        plot(s_ampVal_m-0.15,p_ampVal_m,'k*')      
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Response Peak Value [mV]' ,'FontSize', 16);    
        title(['Mean Response Peak Value ,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
               
        h8=figure;        
        hold on
        set(gca,'Ydir','reverse');
        errbar_h=errorbar([1:size(onVal_m_mat{1}(:,:),2)]-0.3,nanmean(onVal_m_mat{1}(:,:),1), nanstd(onVal_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(onVal_m_mat{2}(:,:),2)],nanmean(onVal_m_mat{2}(:,:),1), nanstd(onVal_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
%         fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(onVal_m_mat{1}(:,:),2)]-0.3, nanmean(onVal_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(onVal_m_mat{2}(:,:),2)], nanmean(onVal_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_onVal_m=min(ylim_data).*ones(size(s_onVal_m)); 
        plot(s_onVal_m-0.15,p_onVal_m,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Onset Vm [mV]' ,'FontSize', 16);    
        title(['Mean Onset Vm value,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
        h9=figure;        
        hold on
        errbar_h=errorbar([1:size(ampDel_m_mat{1}(:,:),2)]-0.3,nanmean(ampDel_m_mat{1}(:,:),1),zeros(1,size(ampDel_m_mat{1}(:,:),2)), nanstd(ampDel_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(ampDel_m_mat{2}(:,:),2)],nanmean(ampDel_m_mat{2}(:,:),1),zeros(1,size(ampDel_m_mat{2}(:,:),2)), nanstd(ampDel_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(ampDel_m_mat{1}(:,:),2)]-0.3, nanmean(ampDel_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(ampDel_m_mat{2}(:,:),2)], nanmean(ampDel_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_ampDel_m=max(ylim_data).*ones(size(s_ampDel_m));
        plot(s_ampDel_m-0.15,p_ampDel_m,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Peak Latency [mS]' ,'FontSize', 16);    
        title(['Mean Peak Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
    
         h10=figure;        
        hold on
        errbar_h=errorbar([1:size(ampDel_std_mat{1}(:,:),2)]-0.3,nanmean(ampDel_std_mat{1}(:,:),1),zeros(1,size(ampDel_std_mat{1}(:,:),2)), nanstd(ampDel_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(ampDel_std_mat{2}(:,:),2)],nanmean(ampDel_std_mat{2}(:,:),1),zeros(1,size(ampDel_std_mat{2}(:,:),2)), nanstd(ampDel_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(ampDel_std_mat{1}(:,:),2)]-0.3, nanmean(ampDel_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(ampDel_std_mat{2}(:,:),2)], nanmean(ampDel_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_ampDel_std=max(ylim_data).*ones(size(s_ampDel_std));
        plot(s_ampDel_std-0.15,p_ampDel_std,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Peak Latency STD [mS]' ,'FontSize', 16);    
        title(['Mean Peak Latency Jitter,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
         h11=figure;        
        hold on
        errbar_h=errorbar([1:size(halfWidth_m_mat{1}(:,:),2)]-0.3,nanmean(halfWidth_m_mat{1}(:,:),1),zeros(1,size(halfWidth_m_mat{1}(:,:),2)), nanstd(halfWidth_m_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(halfWidth_m_mat{2}(:,:),2)],nanmean(halfWidth_m_mat{2}(:,:),1),zeros(1,size(halfWidth_m_mat{2}(:,:),2)), nanstd(halfWidth_m_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(halfWidth_m_mat{1}(:,:),2)]-0.3, nanmean(halfWidth_m_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(halfWidth_m_mat{2}(:,:),2)], nanmean(halfWidth_m_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_halfWidth_m=max(ylim_data).*ones(size(s_halfWidth_m));
        plot(s_halfWidth_m-0.15,p_halfWidth_m,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Half-Width [mS]' ,'FontSize', 16);    
        title(['Mean Half-Width,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
         h12=figure;        
        hold on
        errbar_h=errorbar([1:size(halfWidth_std_mat{1}(:,:),2)]-0.3,nanmean(halfWidth_std_mat{1}(:,:),1),zeros(1,size(halfWidth_std_mat{1}(:,:),2)), nanstd(halfWidth_std_mat{1}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        errbar_h=errorbar([1:size(halfWidth_std_mat{2}(:,:),2)],nanmean(halfWidth_std_mat{2}(:,:),1),zeros(1,size(halfWidth_std_mat{2}(:,:),2)), nanstd(halfWidth_std_mat{2}(:,:),0,1),'.k', 'LineWidth',1.5,'marker','none'); %'markerfacecolor','k'
        fn_errorbar_capsize(errbar_h,capincrease,asymmetric)
        bar([1:size(halfWidth_std_mat{1}(:,:),2)]-0.3, nanmean(halfWidth_std_mat{1}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(1,:),'edgecolor', color_table(1,:))
        bar([1:size(halfWidth_std_mat{2}(:,:),2)], nanmean(halfWidth_std_mat{2}(:,:),1),'barwidth',barwidth1,'facecolor', color_table(2,:),'edgecolor', color_table(2,:))
        ylim_data=[get(gca,'ylim')]';
        p_halfWidth_std=max(ylim_data).*ones(size(s_halfWidth_std));
        plot(s_halfWidth_std-0.15,p_halfWidth_std,'k*')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Half-Width STD [mS]' ,'FontSize', 16);    
        title(['Mean Half-Width STD,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
%% Plot parameters along the train stim - version 1: line+error bars of percent change - preparations
for stim_num=1:11;
       if event_evoked_stat.stim_num(stim_num).lillietest_h_change_failures_m==0      
            h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_failures_m; 
  else h_failures_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_failures_m; 
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
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_onVal_std==0      
%             h_onVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_onVal_std; 
%   else h_onVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_onVal_std; 
%    end
   if event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampVal_m==0      
            h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_ampVal_m; 
  else h_ampVal_m(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampVal_m; 
   end
%    if event_evoked_stat.stim_num(stim_num).lillietest_h_change_ampVal_std==0      
%             h_ampVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_ampVal_std; 
%   else h_ampVal_std(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_ampVal_std; 
%    end
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
   if event_evoked_stat.stim_num(stim_num).lillietest_h_change_adapt_amp==0      
            h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).ttest_h_change_adapt_amp; 
  else h_adapt_amp(1,stim_num)=event_evoked_stat.stim_num(stim_num).wilcoxon_h_change_adapt_amp; 
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
%      s_onVal_std=find(h_onVal_std==1);
     s_ampVal_m=find(h_ampVal_m==1);
%      s_ampVal_std=find(h_ampVal_std==1);
     s_ampDel_m=find(h_ampDel_m==1);
     s_ampDel_std=find(h_ampDel_std==1);
     s_halfWidth_m=find(h_halfWidth_m==1);
     s_halfWidth_std=find(h_halfWidth_std==1);
     s_adapt_amp=find(h_adapt_amp==1);
     s_nonspecific_count_m=find(h_nonspecific_count_m==1);
%%    Plot parameters along the train stim - version 1: line+error bars of percent change 
close all
    if print_flag==1;
h1=figure;        
        hold on
        errorbar([1:size(change_failures_m_mat,2)], nanmean(change_failures_m_mat,1),nanstd(change_failures_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_failures_m=max(ylim_data).*ones(size(s_failures_m));
        plot(s_failures_m,p_failures_m,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Failures difference[%]' ,'FontSize', 16);    
        title(['Mean Rate of Failures,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
        h2=figure;        
        hold on
        errorbar([1:size(change_nonspecific_count_m_mat,2)], nanmean(change_nonspecific_count_m_mat,1),nanstd(change_nonspecific_count_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_nonspecific_count_m=max(ylim_data).*ones(size(s_nonspecific_count_m));
        plot(s_nonspecific_count_m,p_nonspecific_count_m,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Non-specific responses difference[%]' ,'FontSize', 16);    
        title(['Mean Rate of Non-specific responses,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  

        h3=figure;        
        hold on
        errorbar([1:size(change_onset_m_mat,2)], nanmean(change_onset_m_mat,1),nanstd(change_onset_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_onset_m=max(ylim_data).*ones(size(s_onset_m));
        plot(s_onset_m,p_onset_m,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Onset Latency [%change]' ,'FontSize', 16);    
        title(['Mean Onset Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  

        h4=figure;        
        hold on
        errorbar([1:size(change_onset_std_mat,2)], nanmean(change_onset_std_mat,1),nanstd(change_onset_std_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_onset_std=max(ylim_data).*ones(size(s_onset_std));
        plot(s_onset_std,p_onset_std,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Onset Latency STD [%change]' ,'FontSize', 16);    
        title(['Mean response jitter,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  

h5=figure;        
        hold on
        errorbar([1:size(change_amplitude_m_mat,2)], nanmean(change_amplitude_m_mat,1),nanstd(change_amplitude_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_amplitude_m=max(ylim_data).*ones(size(s_amplitude_m));
        plot(s_amplitude_m,p_amplitude_m,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Response Amplitude [%change]' ,'FontSize', 16);    
        title(['Mean Response Amplitude,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
 h6=figure;        
        hold on
        errorbar([1:size(change_amplitude_std_mat,2)], nanmean(change_amplitude_std_mat,1),nanstd(change_amplitude_std_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_amplitude_std=max(ylim_data).*ones(size(s_amplitude_std));
        plot(s_amplitude_std,p_amplitude_std,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Amplitude STD [%change]' ,'FontSize', 16);    
        title(['Mean Amplitude STD (trial-to-trial variability) ,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
 
        h7=figure;        
        hold on
        errorbar([1:size(change_ampVal_m_mat,2)], nanmean(change_ampVal_m_mat,1),nanstd(change_ampVal_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_ampVal_m=max(ylim_data).*ones(size(s_ampVal_m));
        plot(s_ampVal_m,p_ampVal_m,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Response Peak Value [%change]' ,'FontSize', 16);    
        title(['Mean Response Peak Value ,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
               
        h8=figure;        
        hold on
        errorbar([1:size(change_onVal_m_mat,2)], nanmean(change_onVal_m_mat,1),nanstd(change_onVal_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_onVal_m=max(ylim_data).*ones(size(s_onVal_m));
        plot(s_onVal_m,p_onVal_m,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Onset Vm [%change]' ,'FontSize', 16);    
        title(['Mean Onset Vm value,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
        h9=figure;        
        hold on
        errorbar([1:size(change_ampDel_m_mat,2)], nanmean(change_ampDel_m_mat,1),nanstd(change_ampDel_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_ampDel_m=max(ylim_data).*ones(size(s_ampDel_m));
        plot(s_ampDel_m,p_ampDel_m,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Peak Latency [%change]' ,'FontSize', 16);    
        title(['Mean Peak Latency,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
    
         h10=figure;        
        hold on
        errorbar([1:size(change_ampDel_std_mat,2)], nanmean(change_ampDel_std_mat,1),nanstd(change_ampDel_std_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_ampDel_std=max(ylim_data).*ones(size(s_ampDel_std));
        plot(s_ampDel_std,p_ampDel_std,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Peak Latency STD [%change]' ,'FontSize', 16);    
        title(['Mean Peak Latency Jitter,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
         h11=figure;        
        hold on
        errorbar([1:size(change_halfWidth_m_mat,2)], nanmean(change_halfWidth_m_mat,1),nanstd(change_halfWidth_m_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_halfWidth_m=max(ylim_data).*ones(size(s_halfWidth_m));
        plot(s_halfWidth_m,p_halfWidth_m,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Half-Width [%change]' ,'FontSize', 16);    
        title(['Mean Half-Width,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
         h12=figure;        
        hold on
        errorbar([1:size(change_halfWidth_std_mat,2)], nanmean(change_halfWidth_std_mat,1),nanstd(change_halfWidth_std_mat,0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
        ylim_data=[get(gca,'ylim')]';
        p_halfWidth_std=max(ylim_data).*ones(size(s_halfWidth_std));
        plot(s_halfWidth_std,p_halfWidth_std,'k*')
        line([0;12],[0;0],'linestyle','--','linewidth',2,'color','b')
        hold off
        xlabel('Stim. serial number' ,'FontSize', 16);
        ylabel('Half-Width STD [%change]' ,'FontSize', 16);    
        title(['Mean Half-Width STD,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);  
        
    end 
 %% plots
%  close all
for stim_num=1;
tmp_Y=event_evoked_stat.stim_num(stim_num).adapt_amp';
nanmat=~isnan(tmp_Y);
nancount=sum(nanmat(1,:));

tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E =event_evoked_stat.stim_num(stim_num).adapt_amp_std;
g1=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot(tmp_Y')
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Adaptation Amplitude Ratio]', 'FontSize', 28,'fontname', 'arial');
        title(['Adaptation Amplitude Ratio, n=' num2str(nancount) ', p=' num2str(event_evoked_stat.stim_num(stim_num).ttest_p_adapt_amp)] ,'FontSize', 20,'fontname', 'arial');   

%         pause 
% close all
end
        %% save figures
if save_flag==1
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz\single trials'
print(h1,'Evoked Failures','-dpng','-r600','-opengl')
saveas(h1,'Evoked Failures.fig') 
print(h2,'Evoked Nonspecific Responses','-dpng','-r600','-opengl')
saveas(h2,'Evoked Nonspecific Responses.fig') 
print(h3,'Evoked Onset Latency','-dpng','-r600','-opengl')
saveas(h3,'Evoked Onset Latency.fig') 
print(h4,'Evoked Onset Latency STD','-dpng','-r600','-opengl')
saveas(h4,'Evoked Onset Latency STD.fig') 
print(h5,'Evoked Amplitude','-dpng','-r600','-opengl')
saveas(h5,'Evoked Amplitude.fig') 
print(h6,'Evoked Amplitude STD','-dpng','-r600','-opengl')
saveas(h6,'Evoked Amplitude STD.fig') 
print(h7,'Evoked Peak Value','-dpng','-r600','-opengl')
saveas(h7,'Evoked Peak Value.fig') 
print(h8,'Evoked Onset Value','-dpng','-r600','-opengl')
saveas(h8,'Evoked Onset Value.fig') 
print(h9,'Evoked Peak Latency','-dpng','-r600','-opengl')
saveas(h9,'Evoked Peak Latency.fig') 
print(h10,'Evoked Peak Latency STD','-dpng','-r600','-opengl')
saveas(h10,'Evoked Peak Latency STD.fig') 
print(h11,'Evoked Half-Width','-dpng','-r600','-opengl')
saveas(h11,'Evoked Half-Width.fig') 
print(h12,'Evoked Half-Width STD','-dpng','-r600','-opengl')
saveas(h12,'Evoked Half-Width STD.fig') 

print(g1,'Evoked Adaptation Amplitude Ratio_paired plot','-dpng','-r600','-opengl')
saveas(g1,'Evoked Adaptation Amplitude Ratio_paired_plot.fig') 

filename='Evoked activity event detection'; 
save(filename, 'files_to_analyze', 'event_evoked', 'event_evoked_stat')
end
