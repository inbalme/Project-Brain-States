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
%                 stim_temp=stim2_X{2}(:,1);
%                 stim2_X=[];
%                 stim2_X{2}=stim_temp; stim2_X{3}=stim_temp;
%                 
%% low-pass filtering below 300Hz                
                for xx=1:3
    for trace= 1:size(data_no_spikes{channel},2)    
            jj=data_no_spikes{channel}(:,trace,xx);
            data_no_spikes{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jj,dt,-1,300,0,0); 
    end
                end
                %% for spontaneous activity
%  start_time = [0.4,6]; %[sec] %[0,5]
%          duration = 2.5; %[sec] 
%          x_value=1;
%             for t=1:length(start_time);
%              start_sample(:,t) = ceil(start_time(t).*sf{1});
%                 if start_time(t)==0
%                     start_sample(:,t) = 1;
%                 end
%               end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
%               interval(:,t) = start_sample(:,t):end_sample(:,t);  
              %% 
   for trace_type= 2;%1:2; %1:2;  %1 for spont., 2 for evoked
     interval=[]; 
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
                      finalAmp_Thres = 1 ;

                    end 
            case 2
                 start_sample =stim2_X{x_value(1)}(1,1);  %start with sensory stim
        %          duration = 1; %[sec]
                galvano_nstim = Param.facade(6);
                galvano_freq = Param.facade(7);
                duration = galvano_nstim./galvano_freq+0.1;
                end_sample = start_sample+duration.*sf{1}-1;
                interval(:,1) = round(start_sample:end_sample);
                interval(:,2) = interval(:,1);
                        finalAmp_Thres = -1 ;

       end
%% detect events
% for i=1:10;
% voltages_input = data_no_spikes{channel}(i/dt:(i+1)/dt,1,2);
        for t=1:2;
          for trace=1:size(data_no_spikes{1},2)             
              for stim_num=1:galvano_nstim;
                  stim_ISI =1/galvano_freq.*sf{1};
                  peak_start_int = interval(1+stim_ISI*(stim_num-1)+0.003*sf{1},t);
                  peak_end_int = interval(stim_ISI*(stim_num+1),t);
        voltages_input = data_no_spikes{1}(peak_start_int:peak_end_int,trace,x_value(t));
%         finalAmp_Thres = 1 ;
        doPlot = 0;
        I_temp=0;
        [voltages, starting, amplitude, ampPos, halfWidth,halfWidthS, halfWidthE] = EventDetector_v2(voltages_input, dt, finalAmp_Thres, doPlot,I_temp);
        title(['t=', num2str(t), ' trace ',num2str(trace),' stim. ',num2str(stim_num)]);  
       
        %plot whisker stim
%         ylim_data=[get(gca,'ylim')]';
%         patch_xdata=[stim2_X{x_value(2)}; flipud(stim2_X{x_value(2)})];
%         patch_xdata = patch_xdata-interval(1,t);
%         yex=wextend('1D','sym',ylim_data,1);
%         l=ceil((size(patch_xdata,2)-1)/2);        %if the number of stim is odd
%         temp_y=wextend('ac','sym',yex,l);
%         patch_ydata=temp_y(1:size(patch_xdata,1),1:size(patch_xdata,2));
%         patch_cdata=ones(size(patch_xdata));
%         p=patch(patch_xdata,patch_ydata,patch_cdata,'faceColor','r','edgecolor','none','faceAlpha', 0.3);
%         set(gca,'linewidth',1.2)
%          pause
            end
          end
        end
   end
    end
        end
   end
