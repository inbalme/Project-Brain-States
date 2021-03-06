%% for opening workspace saved 
close all
clear all
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
peaks_for_xls=[]; peak_for_xls_mean=[]; 
save_flag= 0;
print_flag=0;
data_vec_all = []; data_vec_residual_all = [];
files_to_analyze =[8,10,11,12,14,15,16,22,36,37,40,1,44,46,48,50,52,56,58,62,72,75,82,84]; %[8,10,11,12,14,15,16,22,36,37,40]; %[1,44,46,48,50,52,56,58,62,72,75]; 
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
   for trace_type= 1;%1:2; %1:2;  %1 for spont., 2 for evoked
     interval=[];       
            x_value=[1,1];    
            
                 start_time = [0.4,5.6]; %[sec] %[0,5]
                 duration = 2.5; %[sec] 
                 epoch_length=0.5; %[sec]
                 epochs=round(duration/epoch_length); 
                    for t=1:length(start_time);
                     start_sample(:,t) = ceil(start_time(t).*sf{1});
                        if start_time(t)==0
                            start_sample(:,t) = 1;
                        end
                      end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
                      interval(:,t) = start_sample(:,t):end_sample(:,t);
                      finalAmp_Thres = 2 ;
                    end                  
%% detect events
% for i=1:10;
% voltages_input = data_no_spikes{channel}(i/dt:(i+1)/dt,1,2);
        for t=1:2;
            event_start_vec{t}=[];  event_onVal_vec{t} = [];  event_amplitude_vec{t} =  [];    event_ampPos_vec{t} =  []; 
            event_ampVal_vec{t} = [];     event_halfWidth_vec{t} =  [];   event_halfWidthS_vec{t} =  []; 
            event_halfWidthE_vec{t} =  [];       
          for trace=1:size(data_no_spikes{1},2)    
              event_start{t,trace}=[];                      event_onVal{t,trace} = [];
            event_amplitude{t,trace} =  [];             event_ampPos{t,trace} =  [];             event_ampVal{t,trace} = [];
            event_halfWidth{t,trace} =  [];            event_halfWidthS{t,trace} =  [];             event_halfWidthE{t,trace} =  []; 
            
              for epoch=1:epochs
                  starting=[]; amplitude=[]; ampPos=[]; halfWidth=[]; halfWidthS=[]; halfWidthE=[];
                  peak_start_int = interval(1+sf{1}*epoch_length*(epoch-1),t);
                  peak_end_int = interval(sf{1}*epoch_length*epoch,t);
        voltages_input = data_no_spikes{1}(peak_start_int:peak_end_int,trace,x_value(t));
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
            event_onVal{t,trace} = [event_onVal{t,trace}; data_no_spikes{1}(starting,trace,x_value(t))];
            event_amplitude{t,trace} =  [event_amplitude{t,trace}; amplitude]; 
            event_ampPos{t,trace} =  [event_ampPos{t,trace}; ampPos]; 
            event_ampVal{t,trace} = [event_ampVal{t,trace}; data_no_spikes{1}(ampPos,trace,x_value(t))];
            event_halfWidth{t,trace} =  [event_halfWidth{t,trace}; halfWidth]; 
            event_halfWidthS{t,trace} =  [event_halfWidthS{t,trace}; halfWidthS]; 
            event_halfWidthE{t,trace} =  [event_halfWidthE{t,trace}; halfWidthE]; 
    end        
            %% visualize the detected events on the trace     
      if print_flag==1
            figure(1); clf
            hold on
                     h1=plot([1:size(data_no_spikes{channel},1)].*dt, data_no_spikes{channel}(:,trace,x_value(t)),'k');
                     h2=scatter(event_start{t,trace}*dt,data_no_spikes{channel}(event_start{t,trace},trace,x_value(t)),'r','fill'); %mark event onset
                     h3=scatter(event_ampPos{t,trace}*dt,data_no_spikes{channel}(event_ampPos{t,trace},trace,x_value(t)),'b','fill'); %mark event peak
                     h4=scatter(event_halfWidthS{t,trace}*dt,data_no_spikes{channel}(event_halfWidthS{t,trace},trace,x_value(t)),'c','fill'); %mark event half-width start
                     h5=scatter(event_halfWidthE{t,trace}*dt,data_no_spikes{channel}(event_halfWidthE{t,trace},trace,x_value(t)),'g','fill'); %mark event  half-width end

                   if t==1;  set(gca,'xlim',[0.4 2.9]); end
                   if t==2;  set(gca,'xlim',[5.6 8.1]); end
%                    set(gca,'ylim',[-50 -10]);
                     hold off
        pause
      end
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz\single trials'
%  if t==1;
% print(1,['detection example - ongoing NB-, trace ', num2str(trace)],'-dpng','-r600','-opengl')
% % saveas(1,'detection example - evoked NB-.fig') 
%  else if t==2
%   print(1,['detection example - ongoing NB+, trace ', num2str(trace)],'-dpng','-r600','-opengl')
% % saveas(1,'detection example - evoked NB+.fig') 
%      end
%  end
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
          end
           event_count(1,t)=length(event_start_vec{t}(:));
           event_freq(1,t) = event_count(1,t)/(duration*size(data_no_spikes{1},2));
        end
   end
        %organizing the data in structure:  
            event_ongoing(fileind).cells=files_to_analyze;
            event_ongoing(fileind).analysis_mfile='Analyze_NBES_ongoing_10_20Hz_single_trials_v2.m';
            event_ongoing(fileind).lp=lp;
            event_ongoing(fileind).finalAmp_Thres=finalAmp_Thres;
            event_ongoing(fileind).start_sample=start_sample;
            event_ongoing(fileind).end_sample=end_sample;
            event_ongoing(fileind).start_time=start_time;
            event_ongoing(fileind).duration=duration;
            event_ongoing(fileind).sf=sf{1};
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
            
            event_count_cells(fileind,:)=event_count;
            event_freq_cells(fileind,:)=event_freq;
            event_onVal_m(fileind,:)=cellfun(@mean,event_onVal_vec);
            event_onVal_std(fileind,:)=cellfun(@std,event_onVal_vec);
            event_ampVal_m(fileind,:)=cellfun(@mean,event_ampVal_vec);
            event_ampVal_std(fileind,:)=cellfun(@std,event_ampVal_vec);
            event_amplitude_m(fileind,:)=cellfun(@mean,event_amplitude_vec);
            event_halfWidth_m(fileind,:)=cellfun(@mean,event_halfWidth_vec);
            event_halfWidth_std(fileind,:)=cellfun(@std,event_halfWidth_vec);

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
        event_ongoing_stat.freq_m= mean(tmp,1);
        event_ongoing_stat.freq_std= std(tmp,0,1);
        %testing for normal distribution
        diff_event_freq= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_freq, event_ongoing_stat.lillietest_p_freq] = lillietest(diff_event_freq);
        %paired ttest
        [event_ongoing_stat.ttest_h_norm_freq, event_ongoing_stat.ttest_p_norm_freq]= ttest(event_ongoing_stat.norm_freq(:,1),event_ongoing_stat.norm_freq(:,2));
        [event_ongoing_stat.ttest_h_freq, event_ongoing_stat.ttest_p_freq]= ttest(event_ongoing_stat.freq(:,1),event_ongoing_stat.freq(:,2));
        [event_ongoing_stat.wilcoxon_p_freq, event_ongoing_stat.wilcoxon_h_freq]= signrank(event_ongoing_stat.freq(:,1),event_ongoing_stat.freq(:,2));

        clear norm_event_freq_noES norm_event_freq_ES diff_event_freq tmp
        
 %mean onset values events
        tmp = event_onVal_m;  
        event_ongoing_stat.onVal_m = event_onVal_m;
        norm_event_onVal_m_noES(:,1)= tmp(:,1)./tmp(:,1);
        norm_event_onVal_m_ES(:,1)= tmp(:,2)./tmp(:,1);
        event_ongoing_stat.norm_onVal_m = [norm_event_onVal_m_noES, norm_event_onVal_m_ES];
        event_ongoing_stat.onVal_m_m= mean(tmp,1);
        event_ongoing_stat.onVal_m_std= std(tmp,0,1);
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
        event_ongoing_stat.onVal_std_m= mean(tmp,1);
        event_ongoing_stat.onVal_std_std= std(tmp,0,1);
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
        event_ongoing_stat.ampVal_m_m= mean(tmp,1);
        event_ongoing_stat.ampVal_m_std= std(tmp,0,1);
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
        event_ongoing_stat.ampVal_std_m= mean(tmp,1);
        event_ongoing_stat.ampVal_std_std= std(tmp,0,1);
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
        event_ongoing_stat.amplitude_m_m= mean(tmp,1);
        event_ongoing_stat.amplitude_m_std= std(tmp,0,1);
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
        event_ongoing_stat.halfWidth_m_m= mean(tmp,1);
        event_ongoing_stat.halfWidth_m_std= std(tmp,0,1);
        %testing for normal distribution
        diff_event_halfWidth_m= tmp(:,2)- tmp(:,1);
        [event_ongoing_stat.lillietest_h_halfWidth_m, event_ongoing_stat.lillietest_p_halfWidth_m] = lillietest(diff_event_halfWidth_m);
        %paired ttest 
        [event_ongoing_stat.ttest_h_norm_halfWidth_m, event_ongoing_stat.ttest_p_norm_halfWidth_m]= ttest(event_ongoing_stat.norm_halfWidth_m(:,1),event_ongoing_stat.norm_halfWidth_m(:,2));
        [event_ongoing_stat.ttest_h_halfWidth_m, event_ongoing_stat.ttest_p_halfWidth_m]= ttest(event_ongoing_stat.halfWidth_m(:,1),event_ongoing_stat.halfWidth_m(:,2));
         [event_ongoing_stat.wilcoxon_p_halfWidth_m, event_ongoing_stat.wilcoxon_h_halfWidth_m]= signrank(event_ongoing_stat.halfWidth_m(:,1),event_ongoing_stat.halfWidth_m(:,2));

        clear norm_event_halfWidth_m_noES norm_event_halfWidth_m_ES diff_event_halfWidth_m tmp
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
errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off
y1limits=get(gca,'ylim');
        x1limits = [0.75 2.25];  x1ticks = [1,2];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'ylim',[-1,20],'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
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
my=max(max(event_ongoing_stat.amplitude_m))*1.1; 
liney=[my ;my ];

g2=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=get(gca,'ylim');
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Amplitude [mV]', 'FontSize', 28,'fontname', 'arial');
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
my=max(max(event_ongoing_stat.halfWidth_m))*1.1; 
liney=[my ;my ];
g3=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)

hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
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
errorbar(tmp_X(:,1), mean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Event Onset value [mV]', 'FontSize', 28,'fontname', 'arial');
        title(['Spontaneous event onset value,  p=' num2str(event_ongoing_stat.wilcoxon_p_onVal_m)] ,'FontSize', 20,'fontname', 'arial');   
        %% save figures
if save_flag==1
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Traces+std+mean+summary\10_20Hz\single trials\ongoing'
print(g1,'Spontaneous event frequency','-dpng','-r600','-opengl')
saveas(g1,'Spontaneous event frequency.fig') 
print(g2,'Spontaneous event amplitude','-dpng','-r600','-opengl')
saveas(g2,'Spontaneous event amplitude.fig') 
print(g3,'Spontaneous event half-width','-dpng','-r600','-opengl')
saveas(g3,'Spontaneous event half-width.fig') 
print(g4,'Spontaneous event onset value','-dpng','-r600','-opengl')
saveas(g4,'Spontaneous event onset value.fig') 
filename='Spontaneous activity'; 
save(filename, 'files_to_analyze', 'event_ongoing', 'event_ongoing_stat')
end
%% Plotting parameters along the stim train
% color_table_rand = rand(length(files_to_analyze),3);
%     if print_flag==1;
% h1=figure;
%         figure(h1);     
%         hold on
%         errorbar([1:size(peak_amp{3},2)], mean(peak_amp{3}./peak_amp{2},1),std(peak_amp{3}./peak_amp{2},0,1),'k-', 'LineWidth',1.5,'marker','o','markerfacecolor','k')
%         line([0;12],[1;1],'linestyle','--','linewidth',2,'color','b')
%         hold off
%         xlabel('Stim. serial number' ,'FontSize', 16);
%         ylabel('Mean Response Amplitude [mV]' ,'FontSize', 16);    
%         title(['Mean Response Amplitude - local baseline,n=' num2str(length(files_to_analyze))] ,'FontSize', 16);     
%         
