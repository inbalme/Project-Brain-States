
%%
close all
clear all
 global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data
 global exp_type
exp_type=4; %1-NBES, 2-ChAT, 4-ChAT awake
channel = 1;
save_flag= 0;
print_flag=1;

switch exp_type
    case 1
        files_to_analyze =[1,8,17,22,23,26,27,28,33,55,60]; %[10]
        % whisker stim onset: [23 6.5; 27 5; 28 5; 33 5.5; 55 5.5; 60 6];
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB-', 'NB+'};
        path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        spike_det_threshold = [1 6 nan; 8 6 nan; 10 5 nan; 17 5 nan; 22 6 nan; 23 2.5 nan; 26 3 nan; 27 1.5 nan; 28 1.5 nan; 33 2 nan; 55 5 nan; 60 2 nan]; %35 2 nan; %29 2 nan; 

    case 2
        files_to_analyze =125; %[74,77];
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light Off', 'Light On'};
        path_output= 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Raster+PSTH';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
         spike_det_threshold =[108 5 nan] ; %[74 10 nan; 77 10 nan];

case 4
        files_to_analyze =124; %[74,77];
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light Off', 'Light On'};
        path_output= 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Raster+PSTH';
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
         spike_det_threshold = [124 5 nan] ;%[108 5 nan] ; %[74 10 nan; 77 10 nan];
         end

    for fileind=1:length(files_to_analyze);
        spikes_stat=[];
    clearvars -except files spike_det_threshold fileind i spikes_stat  spikes channel save_flag print_flag files_to_analyze...
        legend_string path_output exp_type
  
    if channel==1||channel==2
        threshold = spike_det_threshold(spike_det_threshold(:,1)==files_to_analyze(fileind),2);
        cell_ch=1;
    else
        threshold = spike_det_threshold(spike_det_threshold(:,1)==files_to_analyze(fileind),3);
        cell_ch=2;
    end
    
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
    
    clear color_table    
    whisker_stim_color(1,:)=[255 153 153]/256; %[239 188 86]/256;
    switch exp_type
        case 1
            color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256]; 
        case 2
            color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
        case 4
            color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256]; 
    end
%% For cell attached

sf = Param.sf_Vm;
dt = 1/sf;
dt_galvano = 1/Param.sf_galvano;

        high_pass_freq = 300; %[Hz]
        low_pass_freq = 3000; %[Hz]
        BaselineTime = 1000; %[ms]
        robustness = 5;
        peaks_value = [];peaks_location = [];
        interpeak_time = 1; %[ms]
        PSTH=[]; raster=[];
        raster = zeros(size(raw_data{channel}));
        

 for x_value = 1:size(data.x_value,2) 

    data_HP(:,:,x_value) = fn_High_Pass (raw_data{channel}(:,:,x_value), sf, high_pass_freq);
    data_HP(:,:,x_value) = fn_Low_Pass (data_HP(:,:,x_value), sf, low_pass_freq);    
        if exp_type==1&&ES_flag(x_value)==1
            data_HP(29990:35010,:,x_value) =nan; %Ignoring the ES artifact
        end
 
for trace = 1:size(data_HP(:,:,x_value),2)
    data_HP_trace(:,:,x_value) = data_HP(:,trace,x_value);
    
%     threshold_peak(1,trace,x_value) = fn_Threshold_Peak (data_HP_trace(:,:,x_value), sf, BaselineTime, robustness);
    threshold_peak(1,trace,x_value)= threshold;
    
%find peaks that are above the threshold and with minimal distance between
%them specified by "threshold_peak" and "interpeak_time":
    [peaks_value{trace,x_value}(:,1),peaks_location{trace,x_value}(:,1)] = fn_Detect_Spike_HPF ...
        (data_HP_trace(:,:,x_value), sf, interpeak_time, threshold_peak(1,trace,x_value));
    %% making rasters
% raster is a 3D matrix where rasters (x and y axes) are stored. the matrices are sparse (all 0,
% except for 1 at spike times.

        for peak_ind = 1:length(peaks_location{trace,x_value}(:,1))
            raster(peaks_location{trace,x_value}(peak_ind,1),trace,x_value) = 1;
        end
 end
 end
%% Making PSTH from raster
bin_time = 10; %[ms]
bin_size = ceil(bin_time./1000.*sf);

 for x_value = 1:size(data.x_value,2) 
    A = raster(:,:,x_value)'; % convert from sparse to full
    % Plot a line on each spike location 
    [M, N] = size(A);
    [X,Y] = meshgrid(1:N,0:M-1);
    locations_X(1,:) = X(:);
    locations_X(2,:) = X(:);
    locations_Y(1,:) = [Y(:)+1].*A(:);
    locations_Y(2,:) = [Y(:)+1.5].*A(:);
    indxs = find(locations_Y(1,:) ~= 0);
    locations_x{x_value} = locations_X(:,indxs);
    locations_y{x_value} = locations_Y(:,indxs);
    locations_y_galvano_raster = ones(size(stim2_X{x_value})).*(max(locations_Y(2,:))+6);   
    locations_y_ES_raster = ones(size(stim1_X{channel})).*(max(locations_Y(2,:))+6);   

        clear locations_X locations_Y
        
    PSTH(:,:,x_value) = fn_PSTH(bin_time, bin_size, raster(:,:,x_value));
    locations_y_galvano_PSTH = ones(size( stim2_X{x_value})).*(max(PSTH(:,:,x_value))+4);  
 end
%% Making raster+PSTH plots - a figure for spont. activity and a figure for evoked activity
if print_flag==1;
    %spontaneous:
x_value = 1; 
switch exp_type
    case 1
            x1limits(1,:) = [0 3];
            x1limits(2,:) = [3.5 6.5];
    case 2
            x1limits(1,:) = [0.5 4.5];
            x1limits(2,:) = [4.5 8.5];            
    case 4
            x1limits(1,:) = [0 4];
            x1limits(2,:) = [4 8];
end
        
            x1ticks = [];
            y1limits = [0 size(A,1)+8];
            y1ticks = [];
            
for t=1:2;
        f(t)=figure;
        set(gcf,'color','w')
        position_raster = [0.15 , 0.55 , 0.75 , 0.3];
        position_PSTH = [0.15 0.2 0.75 0.3];

        %Raster plot axes:
            
    axes('position', position_raster);           
        hold on
            line(locations_x{x_value}(:,:).*dt,locations_y{x_value}(:,:),'LineWidth',1,'color', color_table(t,:))    %'Color','k')
            line(stim1_X{channel}.*dt,locations_y_ES_raster,'LineWidth',6,'Color','c')
        hold off
    set( gca, 'xlim', x1limits(t,:), 'ylim', y1limits,'xtick', x1ticks,'LineWidth',2,...
            'ticklength', [0.010 0.010],'fontname','arial','fontsize',20,'box','off'); %'fontweight', 'bold', 
    title(['file', num2str(files_to_analyze(fileind)), ' spontaneous activity'] ,'FontSize', 20,'fontname', 'arial');
    %         xlabel('Time [sec]' ,'FontSize', 12);
    ylabel('Trial no.', 'FontSize', 20,'fontname', 'arial');
%             x2limits = [0 size(A,2).*dt];
%             x2ticks = [0:size(A,2).*dt/10:size(A,2).*dt]; %[0:size(A,2).*dt/10:size(A,2).*dt]; x2limits;
            x2ticks(1,:) = [0:1:3];
            x2ticks(2,:) = [3.5:1:6.5];
            x2tickslabel = [0:1:3];
            y2limits = [0 max(PSTH(:,:,x_value))+6];
            y2ticks = [];

     axes('position', position_PSTH);        
        hold on
            bar(time_axis{channel}(bin_size:bin_size:end,x_value),PSTH(:,:,x_value), 'facecolor', color_table(t,:),'edgecolor', color_table(t,:)) %'k'                   
    set(gca, 'xlim', x1limits(t,:), 'ylim', y2limits,'xtick', x2ticks(t,:),'xticklabel',x2tickslabel,'LineWidth',2,...
            'ticklength', [0.010 0.010],'fontname','arial','fontsize',20,'box','off');      %'fontweight', 'bold', 
     xlabel('Time [sec]' ,'FontSize', 20,'fontname', 'arial');
    ylabel('#Spikes/sec', 'FontSize', 20,'fontname', 'arial');
    hold off
    end
%         pause
                %save figure  
                    if save_flag==1;
                        cd(path_output)
                        switch exp_type
                            case 1
                                filename=['file', num2str(files_to_analyze(fileind)), ' spontaneous activity NB-'];
                                saveas(f(1),filename,'fig'); 
                                print(f(1),filename,'-dpng','-r600','-opengl')    
                                filename=['file', num2str(files_to_analyze(fileind)), ' spontaneous activity NB+'];
                                saveas(f(2),filename,'fig'); 
                                print(f(2),filename,'-dpng','-r600','-opengl')  
                            case 2
                                filename=['file', num2str(files_to_analyze(fileind)), ' spontaneous activity Light Off'];
                                saveas(f(1),filename,'fig'); 
                                print(f(1),filename,'-dpng','-r600','-opengl')    
                                filename=['file', num2str(files_to_analyze(fileind)), ' spontaneous activity Light On'];
                                saveas(f(2),filename,'fig'); 
                                print(f(2),filename,'-dpng','-r600','-opengl')  
                            case 4
                                filename=['file', num2str(files_to_analyze(fileind)), ' spontaneous activity Light Off'];
                                saveas(f(1),filename,'fig'); 
                                print(f(1),filename,'-dpng','-r600','-opengl')    
                                filename=['file', num2str(files_to_analyze(fileind)), ' spontaneous activity Light On'];
                                saveas(f(2),filename,'fig'); 
                                print(f(2),filename,'-dpng','-r600','-opengl')
                        end
                    end
%     close all
%% sensory evoked:
for x_value = 2:3; %size(data.x_value,2)
    % x_value=3;
        f(x_value+1)=figure;
        set(gcf,'color','w')
        position_raster = [0.15 , 0.55 , 0.75 , 0.3];
        position_PSTH = [0.15 0.2 0.75 0.3];

        %Raster plot axes:
       
%             x1limits = [stim2_time{1,2}(1,1)-1, stim2_time{1,2}+1]; %[5.5,8]; %[4, 6]; %[5.3 6.8]; [6, 7.5];
            x1limits = [(stim1_X{1,1}(2,1)+1)./sf, stim2_time{1,2}(2,end)+1]; %[5.5,8]; %[4, 6]; %[5.3 6.8]; [6, 7.5];
            x1ticks = [];          
    
            y1limits = [0 size(A,1)+8];
            y1ticks = [];
    axes('position', position_raster);            
        hold on
            line(locations_x{x_value}(:,:).*dt,locations_y{x_value}(:,:),'LineWidth',1,'color', color_table(x_value-1,:))    %'Color','k')
            if x_value==3
                 line(stim1_X{channel}.*dt,locations_y_ES_raster,'LineWidth',6,'Color','c')
            end
            line(stim2_X{x_value}.*dt,locations_y_galvano_raster,'LineWidth',6,'Color','k')      
        hold off
    set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'LineWidth',2,...
            'ticklength', [0.010 0.010],'fontname','arial','fontsize',20,'box','off'); %'fontweight', 'bold', 
            title(['file', num2str(files_to_analyze(fileind)), ' whisker evoked activity'] ,'FontSize', 20,'fontname', 'arial');
    %         xlabel('Time [sec]' ,'FontSize', 12);
            ylabel('Trial no.', 'FontSize', 20,'fontname', 'arial');
            
%             x2limits = [3.5 size(A,2).*dt];
            x2limits = x1limits;
            x2ticks = [x2limits(1):0.5:x2limits(2)]; %[0:size(A,2).*dt/10:size(A,2).*dt]; [3.5:1:size(A,2).*dt];
            y2limits = [0 max(PSTH(:,:,x_value))+6];
            y2ticks = [];
            switch exp_type
                case 1
                    x2tickslabel = [0:0.5:x2limits(2)-x2limits(1)];
                case 2
                    x2tickslabel = [0:0.5:x2limits(2)-x2limits(1)];
            end
   axes('position', position_PSTH);
   hold on
             set(gca, 'xlim', x2limits, 'ylim', y2limits,'xtick', x2ticks,'xticklabel',x2tickslabel,'LineWidth',2,...
                'ticklength', [0.010 0.010],'fontname','arial','fontsize',20,'box','off');      %'fontweight', 'bold', 
        xlabel('Time [sec]' ,'FontSize', 20,'fontname', 'arial');
        ylabel('#Spikes/sec', 'FontSize', 20,'fontname', 'arial');
        
            bar(time_axis{channel}(bin_size:bin_size:end,x_value),PSTH(:,:,x_value),'facecolor', color_table(x_value-1,:),'edgecolor', color_table(x_value-1,:)) %'k'
            % line(locations_x_galvano.*dt_laser,locations_y_galvano_PSTH,'LineWidth',6,'Color','c')
        hold off
        end
%         pause
                %save figure  
                    if save_flag==1;
                        cd(path_output)
                         switch exp_type
                            case 1
                                filename=['file', num2str(files_to_analyze(fileind)), ' whisker evoked activity NB-'];
                                saveas(f(3),filename,'fig'); 
                                print(f(3),filename,'-dpng','-r600','-opengl')    
                                filename=['file', num2str(files_to_analyze(fileind)), ' whisker evoked activity NB+'];
                                saveas(f(4),filename,'fig'); 
                                print(f(4),filename,'-dpng','-r600','-opengl')  
                            case 2
                                filename=['file', num2str(files_to_analyze(fileind)), ' whisker evoked activityLight Off'];
                                saveas(f(3),filename,'fig'); 
                                print(f(3),filename,'-dpng','-r600','-opengl')    
                                filename=['file', num2str(files_to_analyze(fileind)), ' whisker evoked activity Light On'];
                                saveas(f(4),filename,'fig'); 
                                print(f(4),filename,'-dpng','-r600','-opengl')
                         end                    
                    end
%     close all
end
%% Making raster+PSTH plots - one x_value in each figure
% if print_flag==1;
%     for x_value = 1:3; %size(data.x_value,2)
%     % x_value=3;
%         f1=figure;
%         set(gcf,'color','w')
%         position_raster = [0.1 , 0.42 , 0.8 , 0.3];
%         position_PSTH = [0.1 0.1 0.8 0.3];
% 
%         %Raster plot axes:
%             x1limits = [0 size(A,2).*dt];
%             x1ticks = [];
%             y1limits = [0 size(A,1)+8];
%             y1ticks = [];
%             axes('position', position_raster);
%             set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,...
%             'ticklength', [0.010 0.010],'fontname', 'arial','box', 'off'); %'fontweight', 'bold', 
%             title(['file', num2str(Param.name)] ,'FontSize', 20,'fontname', 'arial');
%     %         xlabel('Time [sec]' ,'FontSize', 12);
%             ylabel('Trial no.', 'FontSize', 20,'fontname', 'arial');
%         hold on
%         line(locations_x{x_value}.*dt,locations_y{x_value},'LineWidth',1,'color', color_table(x_value,:))    %'Color','k')
%         line(stim1_X{channel}.*dt,locations_y_ES_raster,'LineWidth',6,'Color','c')
%         if isempty(stim2_X{x_value})
%         else
%          line(stim2_X{x_value}.*dt,locations_y_galvano_raster,'LineWidth',6,'Color','k')
%         end
% 
%         hold off
% 
%             x2limits = [0 size(A,2).*dt];
%             x2ticks = [0:size(A,2).*dt/10:size(A,2).*dt]; %[0:size(A,2).*dt/10:size(A,2).*dt]; x2limits;
%             y2limits = [0 max(PSTH(:,:,x_value))+6];
%             y2ticks = [];
% 
%            axes('position', position_PSTH);
%            set(gca, 'xlim', x2limits, 'ylim', y2limits,'xtick', x2ticks,'XMinorTick','on','fontsize',20,...
%                 'ticklength', [0.010 0.010],'fontname', 'myriad pro','box', 'off' );      %'fontweight', 'bold', 
%         % title(['file', num2str(Param.name), ' header', num2str(data.header(header).x_value(x_value).original_header), ' PSTH Spikes ', ...
%         %           'bin size [ms]=', num2str(bin_time)] ,'FontSize', 14);
%         xlabel('Time [sec]' ,'FontSize', 20,'fontname', 'arial');
%         ylabel('#Spikes/sec', 'FontSize', 20,'fontname', 'arial');
%         hold on
%             bar(time_axis{channel}(bin_size:bin_size:end,x_value),PSTH(:,:,x_value), 'k')
%             % line(locations_x_galvano.*dt_laser,locations_y_galvano_PSTH,'LineWidth',6,'Color','c')
%         hold off
%         pause
%                 %save figure  
%                     if save_flag==1;
%                         cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Raster+PSTH';
%                         filename='jhg';
%                         saveas(gcf,'mjb.fig'); 
%                         print(gcf,filename,'-dpng','-r600','-opengl') 
%                         saveas(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x2+3_mean-subt.fig']) 
%                       print(Fig9,['f' num2str(files_to_analyze(fileind)) '_var_x2+3_mean-subt'],'-dpng','-r600','-opengl')     
%                     end
%     end
%     close all
% end
%% Making raster+PSTH plots - 3 x_values in each figure
% color_raster=[0 0 0; 0 0 0; 0 0 1; 0 0 0; 0 0 1];
% color_bar = ['k';'k';'b';'k';'b'];
% for x_value = 1:3   %1:size(data.x_value,2)
% % x_value=3;
%     figure
%     set(gcf,'color','w')
%     
%     position_raster = [0.1 , 0.42 , 0.8 , 0.3];
%     position_PSTH = [0.1 0.1 0.8 0.3];
%         
%         x1limits = [0 size(A,2).*dt];
%         x1ticks = [];
%         y1limits = [0 size(A,1)+8];
%         y1ticks = [];
% 
%         axes('position', position_raster);
%         set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'arial','box', 'off'); %'fontweight', 'bold', 
%         title(['file', num2str(Param.name)] ,'FontSize', 20,'fontname', 'arial');
% %         xlabel('Time [sec]' ,'FontSize', 12);
%         ylabel('Trial no.', 'FontSize', 20,'fontname', 'arial');
% 
%     hold on
%     line(locations_x{x_value}.*dt,locations_y{x_value},'LineWidth',1,'color', color_raster(x_value,:))   
%     line(stim1_X{channel}.*dt,locations_y_ES_raster,'LineWidth',6,'Color','k')
%     if isempty(stim2_X{x_value})
%     else
%      line(stim2_X{x_value}.*dt,locations_y_galvano_raster,'LineWidth',6,'Color','k')
%     end
% 
%     hold off
%     
%         x2limits = [0 size(A,2).*dt];
%         x2ticks = [0:size(A,2).*dt/10:size(A,2).*dt]; %[0:size(A,2).*dt/10:size(A,2).*dt]; x2limits;
%         y2limits = [0 max(PSTH(:,:,x_value))+6];
%         y2ticks = [];
%    
%    axes('position', position_PSTH);
%    set(gca, 'xlim', x2limits, 'ylim', y2limits,'xtick', x2ticks,'XMinorTick','on','fontsize',20,'linewidth',1,...
%         'ticklength', [0.010 0.010],'fontname', 'myriad pro','box', 'off' );      %'fontweight', 'bold', 
% % title(['file', num2str(Param.name), ' header', num2str(data.header(header).x_value(x_value).original_header), ' PSTH Spikes ', ...
% %           'bin size [ms]=', num2str(bin_time)] ,'FontSize', 14);
% xlabel('Time [sec]' ,'FontSize', 20,'fontname', 'arial');
% ylabel('#Spikes/sec', 'FontSize', 20,'fontname', 'arial');
% hold on
% bar(time_axis{channel}(bin_size:bin_size:end,x_value),PSTH(:,:,x_value),color_bar(x_value))
% % line(locations_x_galvano.*dt_laser,locations_y_galvano_PSTH,'LineWidth',6,'Color','c')
% hold off
% end
 %% Calculating SNR  and Modulation Index
 %According to "Vinck_Cardin_2014_Arousal and locomotion make distinct contribution..."
 %the MI was calculated in the same manner in Cohen&Maunsel 2009. test significance for MI>0.
 % SNR = (MEANsignal - MEANbackground/(MEANsignal + MEANbackground)
 %MEANsignal = mean response following stimulus onset, calculated from PSTH.
 %MEANbackground = mean response prior to stimulus onset, calculated from PSTH.
 %MI = (R_ES-R_ctr)/(R_ES+R_ctr)
 %R_ES = sensory evoked response with ES
 %R_ctr = sensory evoked response without ES
  response_window = 500; %[ms]
 for x_value = 1:size(data.x_value,2) 
    if isempty(stim2_X{x_value})
       mean_signal{x_value}=[];
       mean_background{x_value}=[];
    else  

       stim_begin_sample = stim2_X{x_value}(1);
       response_end_sample = stim_begin_sample+response_window*sf/1000;
       background_begin_sample = stim_begin_sample-response_window*sf/1000;
       background_end_sample = stim_begin_sample;
       
       mean_signal{x_value}=mean(mean(raster(stim_begin_sample:response_end_sample,:,x_value),2));
       mean_background{x_value}=mean(mean(raster(background_begin_sample:background_end_sample,:,x_value),2));
       SNR{x_value}= (mean_signal{x_value}-mean_background{x_value})/(mean_signal{x_value}+mean_background{x_value});       
       
                % change in firing:
        nspikes_signal{x_value}=sum(mean(raster(stim_begin_sample:response_end_sample,:,x_value),2));
        nspikes_background{x_value}=sum(mean(raster(background_begin_sample:background_end_sample,:,x_value),2));
        nspikes_signal_net{x_value}=nspikes_signal{x_value}-nspikes_background{x_value};
        nspikes_signal_Hz{x_value}=nspikes_signal{x_value}.*(1000/response_window);
        nspikes_background_Hz{x_value}=nspikes_background{x_value}.*(1000/response_window);
        nspikes_signal_net_Hz{x_value}=nspikes_signal_net{x_value}.*(1000/response_window);

    end
 end
 MI = (mean_signal{3}-mean_signal{2})/(mean_signal{3}+mean_signal{2});
 

 %% Calculating SNR 
%  %According to Stuttgen&Schwarz 2010 'Integration of vibrotactile signals...'
%  % Glass's Delta = (MEANsignal - MEANbackground)/SDbackground
%  %MEANsignal = mean response 20ms following stimulus onset, calculated from PSTH.
%  %MEANbackground = mean response 20ms prior to stimulus onset, calculated from PSTH.
%  %SDbackground = standard deviation of the same period as MEANbackground.
%  for x_value = 1:size(data.x_value,2) 
%     if isempty(stim2_X{x_value})
%        mean_signal{x_value}=[];
%        mean_background{x_value}=[];
%        std_background{x_value}=[];
%     else  
% 
%        stim_begin_sample = stim2_X{x_value}(1);
%        response_window = 50; %[ms]
%        response_end_sample = stim_begin_sample+response_window*sf/1000;
%        background_begin_sample = stim_begin_sample-response_window*sf/1000;
%        background_end_sample = stim_begin_sample;
%        
%        mean_signal{x_value}=mean(mean(raster(stim_begin_sample:response_end_sample,:,x_value),2));
%        mean_background{x_value}=mean(mean(raster(background_begin_sample:background_end_sample,:,x_value),2));
%        std_background{x_value}=std(mean(raster(background_begin_sample:background_end_sample,:,x_value),2));
%        SNR{x_value}= (mean_signal{x_value}-mean_background{x_value})/std_background{x_value};
%     end
%  end
%% Calculating latency to first stim in the train + failures rate
 lat_response_window = 50; %[ms]
 train_flag = 0;
whisker_response_flag=0;
for x_value = 2:3; %size(data.x_value,2)
    if train_flag==1;
         stim_num = size(stim2_X{x_value},2);
    else
        stim_num = 1;
    end
            for stim_serial_num = 1:stim_num; %The number of stim in the train. 
               
                stim_begin_sample = stim2_X{x_value}(stim_serial_num,1);
                response_refractory_time = 4; %[ms] minimal time for whisker responses in the cortex
                response_refractory_sample = response_refractory_time*sf/1000;
                response_begin_sample = stim_begin_sample+response_refractory_sample;
                response_end_sample = stim_begin_sample+lat_response_window*sf/1000;
            
         for trace = 1:size(data_HP(:,:,x_value),2)
             if isempty(find(raster(response_begin_sample:response_end_sample,trace,x_value)==1,1,'first'))
                 response_latency{stim_serial_num,x_value}(trace,1) =nan;
             else
                     response_latency{stim_serial_num,x_value}(trace,1) =( find(raster(response_begin_sample:response_end_sample,trace,x_value)==1,1,'first')+response_refractory_sample)*dt*1000; %ms]
                    
             end
         end
          response_latency_mean{stim_serial_num,x_value} = nanmean(response_latency{stim_serial_num,x_value});
          response_latency_std{stim_serial_num,x_value} = nanstd(response_latency{stim_serial_num,x_value});
          success_rate{stim_serial_num,x_value} = 1-sum(isnan(response_latency{stim_serial_num, x_value}))/length(response_latency{stim_serial_num,x_value});
end
end
         if (success_rate{1,2}>0.5 || success_rate{1,3}>0.5)
            whisker_response_flag=1;
         end
        
%%

if fileind==1
    spikes = [];     
end
% load('extracellular_spikes_500ms')

% cell=find(spike_det_threshold(:,1)==fileind);   %this is the serial number of cells that are included in the extracellular analysis 
spikes(fileind).files_to_analyze = files_to_analyze;
spikes(fileind).Ch = cell_ch; %this is not the recording channel but preamp1 or preamp 2
spikes(fileind).res_time_win = response_window;
spikes(fileind).SNR = SNR;
spikes(fileind).modulation_ind = MI;
spikes(fileind).mean_signal = mean_signal; 
spikes(fileind).mean_background = mean_background;
spikes(fileind).nspikes_signal=nspikes_signal; %sum of spikes along the train averaged across trials
spikes(fileind).nspikes_background=nspikes_background; %sum of spikes along the background averaged across trials
spikes(fileind).nspikes_signal_net=nspikes_signal_net; %nspikes _signal - nspikes_background
spikes(fileind).nspikes_signal_Hz=nspikes_signal_Hz; %same but in Hz
spikes(fileind).nspikes_background_Hz=nspikes_background_Hz; %same but in Hz
spikes(fileind).nspikes_signal_net_Hz=nspikes_signal_net_Hz; %same but in Hz
spikes(fileind).raster = raster;
spikes(fileind).PSTH = PSTH;
spikes(fileind).PSTH_binsize = bin_size;
spikes(fileind).file = Param.name;
spikes(fileind).spike_threshold = threshold;
spikes(fileind).sf = sf;
spikes(fileind).lat_response_win = lat_response_window;
spikes(fileind).latency = response_latency;
spikes(fileind).latency_mean = response_latency_mean;
spikes(fileind).latency_std = response_latency_std;
spikes(fileind).success_rate = success_rate;
spikes(fileind).whisker_response_flag = whisker_response_flag;
% save('extracellular_spikes_500ms','spikes')

    end
%% plotting a trace
% x_value=2
%   [Fig4,h5]= fn_Plot_Trace_v2(data_HP(:,7,x_value), dt, dt, stim1_X, dt_galvano, stim2_X{x_value});
%    [Fig3,h3]= fn_Plot_Trace_v2(PSTH(:,1,x_value), dt*bin_size, dt*bin_size, stim1_X, dt_galvano, stim2_X{x_value});
%%
clear SNR nspikes_signal nspikes_background nspikes_signal_Hz nspikes_background_Hz
for cell_ind = 1:size(spikes,2)
    SNR(cell_ind,:) = cell2mat(spikes(cell_ind).SNR(1,2:3));
    MI(cell_ind,:) = spikes(cell_ind).modulation_ind;
    R_ctrl_NBES_temp(cell_ind,:) = cell2mat(spikes(cell_ind).nspikes_signal_net(1,2:3));
    nspikes_signal(cell_ind,:) = cell2mat(spikes(cell_ind).nspikes_signal(1,2:3));
    nspikes_background(cell_ind,:) = cell2mat(spikes(cell_ind).nspikes_background(1,2:3));
    R_ctrl_NBES_Hz_temp(cell_ind,:) = cell2mat(spikes(cell_ind).nspikes_signal_net_Hz(1,2:3));
    nspikes_signal_Hz(cell_ind,:) = cell2mat(spikes(cell_ind).nspikes_signal_Hz(1,2:3));
    nspikes_background_Hz(cell_ind,:) = cell2mat(spikes(cell_ind).nspikes_background_Hz(1,2:3));
    Latency_mean(cell_ind,:) = cell2mat(spikes(cell_ind).latency_mean(1,2:3));
    Latency_std(cell_ind,:) = cell2mat(spikes(cell_ind).latency_std(1,2:3));
    Success_rate(cell_ind,:) = cell2mat(spikes(cell_ind).success_rate(1,2:3));
    Whisker_response(cell_ind,:) = spikes(cell_ind).whisker_response_flag;
end
%% statistics
      
%SNR
        spikes_stat(stim_num).SNR=SNR;
        spikes_stat(stim_num).SNR_m=nanmean(spikes_stat(stim_num).SNR,1);
        spikes_stat(stim_num).SNR_std=nanstd(spikes_stat(stim_num).SNR,0,1);
        spikes_stat(stim_num).change_SNR=[(spikes_stat(stim_num).SNR(:,2)-spikes_stat(stim_num).SNR(:,1))./abs(spikes_stat(stim_num).SNR(:,1))].*100; %percent change
        spikes_stat(stim_num).change_SNR_m=nanmean(spikes_stat(stim_num).change_SNR,1);
        spikes_stat(stim_num).change_SNR_std=nanstd(spikes_stat(stim_num).change_SNR,0,1);
        change_SNR_mat(:,stim_num)= spikes_stat(stim_num).change_SNR;
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_SNR, spikes_stat(stim_num).lillietest_p_SNR] = lillietest(spikes_stat(stim_num).SNR(:,2)- spikes_stat(stim_num).SNR(:,1));
        [spikes_stat(stim_num).lillietest_h_change_SNR, spikes_stat(stim_num).lillietest_p_change_SNR] = lillietest(spikes_stat(stim_num).change_SNR);
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_SNR, spikes_stat(stim_num).ttest_p_SNR,spikes_stat(stim_num).ttest_stats_SNR]= ttest(spikes_stat(stim_num).SNR(:,1),spikes_stat(stim_num).SNR(:,2));   
        [spikes_stat(stim_num).wilcoxon_p_SNR, spikes_stat(stim_num).wilcoxon_h_SNR,spikes_stat(stim_num).wilcoxon_stats_SNR]= signrank(spikes_stat(stim_num).SNR(:,1),spikes_stat(stim_num).SNR(:,2));
        [spikes_stat(stim_num).ttest_h_change_SNR, spikes_stat(stim_num).ttest_p_change_SNR]= ttest(spikes_stat(stim_num).change_SNR);
        [spikes_stat(stim_num).wilcoxon_p_change_SNR, spikes_stat(stim_num).wilcoxon_h_change_SNR]= signrank(spikes_stat(stim_num).change_SNR);
        
        %Response Modulation
        spikes_stat(stim_num).res_modulation=R_ctrl_NBES_temp;
       spikes_stat(stim_num).res_modulation_m=nanmean(spikes_stat(stim_num).res_modulation,1);
        spikes_stat(stim_num).res_modulation_std=nanstd(spikes_stat(stim_num).res_modulation,0,1);
        spikes_stat(stim_num).change_res_modulation=[(spikes_stat(stim_num).res_modulation(:,2)-spikes_stat(stim_num).res_modulation(:,1))./abs(spikes_stat(stim_num).res_modulation(:,1))].*100; %percent change
        spikes_stat(stim_num).change_res_modulation_m=nanmean(spikes_stat(stim_num).change_res_modulation,1);
        spikes_stat(stim_num).change_res_modulation_std=nanstd(spikes_stat(stim_num).change_res_modulation,0,1);
        change_res_modulation_mat(:,stim_num)= spikes_stat(stim_num).change_res_modulation;
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_res_modulation, spikes_stat(stim_num).lillietest_p_res_modulation] = lillietest(spikes_stat(stim_num).res_modulation(:,2)- spikes_stat(stim_num).res_modulation(:,1));
        [spikes_stat(stim_num).lillietest_h_change_res_modulation, spikes_stat(stim_num).lillietest_p_change_res_modulation] = lillietest(spikes_stat(stim_num).change_res_modulation);
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_res_modulation, spikes_stat(stim_num).ttest_p_res_modulation,spikes_stat(stim_num).ttest_stats_res_modulation]= ttest(spikes_stat(stim_num).res_modulation(:,1),spikes_stat(stim_num).res_modulation(:,2));   
        [spikes_stat(stim_num).wilcoxon_p_res_modulation, spikes_stat(stim_num).wilcoxon_h_res_modulation,spikes_stat(stim_num).wilcoxon_stats_res_modulation]= signrank(spikes_stat(stim_num).res_modulation(:,1),spikes_stat(stim_num).res_modulation(:,2));
        [spikes_stat(stim_num).ttest_h_change_res_modulation, spikes_stat(stim_num).ttest_p_change_res_modulation]= ttest(spikes_stat(stim_num).change_res_modulation);
        [spikes_stat(stim_num).wilcoxon_p_change_res_modulation, spikes_stat(stim_num).wilcoxon_h_change_res_modulation]= signrank(spikes_stat(stim_num).change_res_modulation);

    %Response Modulation [Hz]
        spikes_stat(stim_num).res_modulation_Hz=R_ctrl_NBES_Hz_temp;
       spikes_stat(stim_num).res_modulation_Hz_m=nanmean(spikes_stat(stim_num).res_modulation_Hz,1);
        spikes_stat(stim_num).res_modulation_Hz_std=nanstd(spikes_stat(stim_num).res_modulation_Hz,0,1);
        spikes_stat(stim_num).change_res_modulation_Hz=[(spikes_stat(stim_num).res_modulation_Hz(:,2)-spikes_stat(stim_num).res_modulation_Hz(:,1))./abs(spikes_stat(stim_num).res_modulation_Hz(:,1))].*100; %percent change
        spikes_stat(stim_num).change_res_modulation_Hz_m=nanmean(spikes_stat(stim_num).change_res_modulation_Hz,1);
        spikes_stat(stim_num).change_res_modulation_Hz_std=nanstd(spikes_stat(stim_num).change_res_modulation_Hz,0,1);
        change_res_modulation_Hz_mat(:,stim_num)= spikes_stat(stim_num).change_res_modulation_Hz;
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_res_modulation_Hz, spikes_stat(stim_num).lillietest_p_res_modulation_Hz] = lillietest(spikes_stat(stim_num).res_modulation_Hz(:,2)- spikes_stat(stim_num).res_modulation_Hz(:,1));
        [spikes_stat(stim_num).lillietest_h_change_res_modulation_Hz, spikes_stat(stim_num).lillietest_p_change_res_modulation_Hz] = lillietest(spikes_stat(stim_num).change_res_modulation_Hz);
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_res_modulation_Hz, spikes_stat(stim_num).ttest_p_res_modulation_Hz,spikes_stat(stim_num).ttest_stats_res_modulation_Hz]= ttest(spikes_stat(stim_num).res_modulation_Hz(:,1),spikes_stat(stim_num).res_modulation_Hz(:,2));   
        [spikes_stat(stim_num).wilcoxon_p_res_modulation_Hz, spikes_stat(stim_num).wilcoxon_h_res_modulation_Hz,spikes_stat(stim_num).wilcoxon_stats_res_modulation_Hz]= signrank(spikes_stat(stim_num).res_modulation_Hz(:,1),spikes_stat(stim_num).res_modulation_Hz(:,2));
        [spikes_stat(stim_num).ttest_h_change_res_modulation_Hz, spikes_stat(stim_num).ttest_p_change_res_modulation_Hz]= ttest(spikes_stat(stim_num).change_res_modulation_Hz);
        [spikes_stat(stim_num).wilcoxon_p_change_res_modulation_Hz, spikes_stat(stim_num).wilcoxon_h_change_res_modulation_Hz]= signrank(spikes_stat(stim_num).change_res_modulation_Hz);

%mean response spike count over the specified time window
        spikes_stat(stim_num).nspikes_signal=nspikes_signal;
        spikes_stat(stim_num).nspikes_signal_m=nanmean(spikes_stat(stim_num).nspikes_signal,1);
        spikes_stat(stim_num).nspikes_signal_std=nanstd(spikes_stat(stim_num).nspikes_signal,0,1);
        spikes_stat(stim_num).change_nspikes_signal=[(spikes_stat(stim_num).nspikes_signal(:,2)-spikes_stat(stim_num).nspikes_signal(:,1))./abs(spikes_stat(stim_num).nspikes_signal(:,1))].*100; %percent change
        spikes_stat(stim_num).change_nspikes_signal_m=nanmean(spikes_stat(stim_num).change_nspikes_signal,1);
        spikes_stat(stim_num).change_nspikes_signal_std=nanstd(spikes_stat(stim_num).change_nspikes_signal,0,1);
        change_nspikes_signal_mat(:,stim_num)= spikes_stat(stim_num).change_nspikes_signal;
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_nspikes_signal, spikes_stat(stim_num).lillietest_p_nspikes_signal] = lillietest(spikes_stat(stim_num).nspikes_signal(:,2)- spikes_stat(stim_num).nspikes_signal(:,1));
        [spikes_stat(stim_num).lillietest_h_change_nspikes_signal, spikes_stat(stim_num).lillietest_p_change_nspikes_signal] = lillietest(spikes_stat(stim_num).change_nspikes_signal);
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_nspikes_signal, spikes_stat(stim_num).ttest_p_nspikes_signal,spikes_stat(stim_num).ttest_stats_nspikes_signal]= ttest(spikes_stat(stim_num).nspikes_signal(:,1),spikes_stat(stim_num).nspikes_signal(:,2));   
        [spikes_stat(stim_num).wilcoxon_p_nspikes_signal, spikes_stat(stim_num).wilcoxon_h_nspikes_signal,spikes_stat(stim_num).wilcoxon_stats_nspikes_signal]= signrank(spikes_stat(stim_num).nspikes_signal(:,1),spikes_stat(stim_num).nspikes_signal(:,2));
        [spikes_stat(stim_num).ttest_h_change_nspikes_signal, spikes_stat(stim_num).ttest_p_change_nspikes_signal]= ttest(spikes_stat(stim_num).change_nspikes_signal);
        [spikes_stat(stim_num).wilcoxon_p_change_nspikes_signal, spikes_stat(stim_num).wilcoxon_h_change_nspikes_signal]= signrank(spikes_stat(stim_num).change_nspikes_signal);
                
% mean background spike count over the specified time window
        spikes_stat(stim_num).nspikes_background=nspikes_background;
        spikes_stat(stim_num).nspikes_background_m=nanmean(spikes_stat(stim_num).nspikes_background,1);
        spikes_stat(stim_num).nspikes_background_std=nanstd(spikes_stat(stim_num).nspikes_background,0,1);
        spikes_stat(stim_num).change_nspikes_background=[(spikes_stat(stim_num).nspikes_background(:,2)-spikes_stat(stim_num).nspikes_background(:,1))./abs(spikes_stat(stim_num).nspikes_background(:,1))].*100; %percent change
        spikes_stat(stim_num).change_nspikes_background_m=nanmean(spikes_stat(stim_num).change_nspikes_background,1);
        spikes_stat(stim_num).change_nspikes_background_std=nanstd(spikes_stat(stim_num).change_nspikes_background,0,1);
        change_nspikes_background_mat(:,stim_num)= spikes_stat(stim_num).change_nspikes_background;
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_nspikes_background, spikes_stat(stim_num).lillietest_p_nspikes_background] = lillietest(spikes_stat(stim_num).nspikes_background(:,2)- spikes_stat(stim_num).nspikes_background(:,1));
        [spikes_stat(stim_num).lillietest_h_change_nspikes_background, spikes_stat(stim_num).lillietest_p_change_nspikes_background] = lillietest(spikes_stat(stim_num).change_nspikes_background);
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_nspikes_background, spikes_stat(stim_num).ttest_p_nspikes_background,spikes_stat(stim_num).ttest_stats_nspikes_background]= ttest(spikes_stat(stim_num).nspikes_background(:,1),spikes_stat(stim_num).nspikes_background(:,2));   
        [spikes_stat(stim_num).wilcoxon_p_nspikes_background, spikes_stat(stim_num).wilcoxon_h_nspikes_background,spikes_stat(stim_num).wilcoxon_stats_nspikes_background]= signrank(spikes_stat(stim_num).nspikes_background(:,1),spikes_stat(stim_num).nspikes_background(:,2));
        [spikes_stat(stim_num).ttest_h_change_nspikes_background, spikes_stat(stim_num).ttest_p_change_nspikes_background]= ttest(spikes_stat(stim_num).change_nspikes_background);
        [spikes_stat(stim_num).wilcoxon_p_change_nspikes_background, spikes_stat(stim_num).wilcoxon_h_change_nspikes_background]= signrank(spikes_stat(stim_num).change_nspikes_background);

 % mean background spike count over the specified time window [Hz]
        spikes_stat(stim_num).nspikes_background_Hz=nspikes_background_Hz;
        spikes_stat(stim_num).nspikes_background_Hz_m=nanmean(spikes_stat(stim_num).nspikes_background_Hz,1);
        spikes_stat(stim_num).nspikes_background_Hz_std=nanstd(spikes_stat(stim_num).nspikes_background_Hz,0,1);
        spikes_stat(stim_num).change_nspikes_background_Hz=[(spikes_stat(stim_num).nspikes_background_Hz(:,2)-spikes_stat(stim_num).nspikes_background_Hz(:,1))./abs(spikes_stat(stim_num).nspikes_background_Hz(:,1))].*100; %percent change
        spikes_stat(stim_num).change_nspikes_background_Hz_m=nanmean(spikes_stat(stim_num).change_nspikes_background_Hz,1);
        spikes_stat(stim_num).change_nspikes_background_Hz_std=nanstd(spikes_stat(stim_num).change_nspikes_background_Hz,0,1);
        change_nspikes_background_Hz_mat(:,stim_num)= spikes_stat(stim_num).change_nspikes_background_Hz;
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_nspikes_background_Hz, spikes_stat(stim_num).lillietest_p_nspikes_background_Hz] = lillietest(spikes_stat(stim_num).nspikes_background_Hz(:,2)- spikes_stat(stim_num).nspikes_background_Hz(:,1));
        [spikes_stat(stim_num).lillietest_h_change_nspikes_background_Hz, spikes_stat(stim_num).lillietest_p_change_nspikes_background_Hz] = lillietest(spikes_stat(stim_num).change_nspikes_background_Hz);
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_nspikes_background_Hz, spikes_stat(stim_num).ttest_p_nspikes_background_Hz,spikes_stat(stim_num).ttest_stats_nspikes_background_Hz]= ttest(spikes_stat(stim_num).nspikes_background_Hz(:,1),spikes_stat(stim_num).nspikes_background_Hz(:,2));   
        [spikes_stat(stim_num).wilcoxon_p_nspikes_background_Hz, spikes_stat(stim_num).wilcoxon_h_nspikes_background_Hz,spikes_stat(stim_num).wilcoxon_stats_nspikes_background_Hz]= signrank(spikes_stat(stim_num).nspikes_background_Hz(:,1),spikes_stat(stim_num).nspikes_background_Hz(:,2));
        [spikes_stat(stim_num).ttest_h_change_nspikes_background_Hz, spikes_stat(stim_num).ttest_p_change_nspikes_background_Hz]= ttest(spikes_stat(stim_num).change_nspikes_background_Hz);
        [spikes_stat(stim_num).wilcoxon_p_change_nspikes_background_Hz, spikes_stat(stim_num).wilcoxon_h_change_nspikes_background_Hz]= signrank(spikes_stat(stim_num).change_nspikes_background_Hz);
       
%Latency Mean (only cells which faithfully responded to the 1st stim in the train)        
        spikes_stat(stim_num).latency=Latency_mean(Whisker_response==1,:);        
        spikes_stat(stim_num).latency_m=nanmean(spikes_stat(stim_num).latency,1);
        spikes_stat(stim_num).latency_std=nanstd(spikes_stat(stim_num).latency,0,1);
        spikes_stat(stim_num).change_latency=[(spikes_stat(stim_num).latency(:,2)-spikes_stat(stim_num).latency(:,1))./abs(spikes_stat(stim_num).latency(:,1))].*100; %percent change
        spikes_stat(stim_num).change_latency_m=nanmean(spikes_stat(stim_num).change_latency,1);
        spikes_stat(stim_num).change_latency_std=nanstd(spikes_stat(stim_num).change_latency,0,1);
        change_latency_mat(:,stim_num)= spikes_stat(stim_num).change_latency;
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_latency, spikes_stat(stim_num).lillietest_p_latency] = lillietest(spikes_stat(stim_num).latency(:,2)- spikes_stat(stim_num).latency(:,1));
        [spikes_stat(stim_num).lillietest_h_change_latency, spikes_stat(stim_num).lillietest_p_change_latency] = lillietest(spikes_stat(stim_num).change_latency);
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_latency, spikes_stat(stim_num).ttest_p_latency,spikes_stat(stim_num).ttest_stats_latency]= ttest(spikes_stat(stim_num).latency(:,1),spikes_stat(stim_num).latency(:,2));   
        [spikes_stat(stim_num).wilcoxon_p_latency, spikes_stat(stim_num).wilcoxon_h_latency,spikes_stat(stim_num).wilcoxon_stats_latency]= signrank(spikes_stat(stim_num).latency(:,1),spikes_stat(stim_num).latency(:,2));
        [spikes_stat(stim_num).ttest_h_change_latency, spikes_stat(stim_num).ttest_p_change_latency]= ttest(spikes_stat(stim_num).change_latency);
        [spikes_stat(stim_num).wilcoxon_p_change_latency, spikes_stat(stim_num).wilcoxon_h_change_latency]= signrank(spikes_stat(stim_num).change_latency);

%Latency STD (jitter) (only cells which faithfully responded to the 1st stim in the train)       
        spikes_stat(stim_num).latency_std=Latency_std(Whisker_response==1,:);
        spikes_stat(stim_num).latency_std_m=nanmean(spikes_stat(stim_num).latency_std,1);
        spikes_stat(stim_num).latency_std_std=nanstd(spikes_stat(stim_num).latency_std,0,1);
        spikes_stat(stim_num).change_latency_std=[(spikes_stat(stim_num).latency_std(:,2)-spikes_stat(stim_num).latency_std(:,1))./abs(spikes_stat(stim_num).latency_std(:,1))].*100; %percent change
        spikes_stat(stim_num).change_latency_std_m=nanmean(spikes_stat(stim_num).change_latency_std,1);
        spikes_stat(stim_num).change_latency_std_std=nanstd(spikes_stat(stim_num).change_latency_std,0,1);
        change_latency_std_mat(:,stim_num)= spikes_stat(stim_num).change_latency_std;
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_latency_std, spikes_stat(stim_num).lillietest_p_latency_std] = lillietest(spikes_stat(stim_num).latency_std(:,2)- spikes_stat(stim_num).latency_std(:,1));
        [spikes_stat(stim_num).lillietest_h_change_latency_std, spikes_stat(stim_num).lillietest_p_change_latency_std] = lillietest(spikes_stat(stim_num).change_latency_std);
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_latency_std, spikes_stat(stim_num).ttest_p_latency_std,spikes_stat(stim_num).ttest_stats_latency_std]= ttest(spikes_stat(stim_num).latency_std(:,1),spikes_stat(stim_num).latency_std(:,2));   
        [spikes_stat(stim_num).wilcoxon_p_latency_std, spikes_stat(stim_num).wilcoxon_h_latency_std,spikes_stat(stim_num).wilcoxon_stats_latency_std]= signrank(spikes_stat(stim_num).latency_std(:,1),spikes_stat(stim_num).latency_std(:,2));
        [spikes_stat(stim_num).ttest_h_change_latency_std, spikes_stat(stim_num).ttest_p_change_latency_std]= ttest(spikes_stat(stim_num).change_latency_std);
        [spikes_stat(stim_num).wilcoxon_p_change_latency_std, spikes_stat(stim_num).wilcoxon_h_change_latency_std]= signrank(spikes_stat(stim_num).change_latency_std);
        
         spikes_stat(stim_num).whisker_response=Whisker_response;
         
         %Success Rate
        spikes_stat(stim_num).success_rate=Success_rate;
        spikes_stat(stim_num).success_rate_m=nanmean(spikes_stat(stim_num).success_rate,1);
        spikes_stat(stim_num).success_rate_std=nanstd(spikes_stat(stim_num).success_rate,0,1);
        spikes_stat(stim_num).change_success_rate=[(spikes_stat(stim_num).success_rate(:,2)-spikes_stat(stim_num).success_rate(:,1))./abs(spikes_stat(stim_num).success_rate(:,1))].*100; %percent change
        spikes_stat(stim_num).change_success_rate_m=nanmean(spikes_stat(stim_num).change_success_rate,1);
        spikes_stat(stim_num).change_success_rate_std=nanstd(spikes_stat(stim_num).change_success_rate,0,1);
        change_success_rate_mat(:,stim_num)= spikes_stat(stim_num).change_success_rate;
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_success_rate, spikes_stat(stim_num).lillietest_p_success_rate] = lillietest(spikes_stat(stim_num).success_rate(:,2)- spikes_stat(stim_num).success_rate(:,1));
        [spikes_stat(stim_num).lillietest_h_change_success_rate, spikes_stat(stim_num).lillietest_p_change_success_rate] = lillietest(spikes_stat(stim_num).change_success_rate);
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_success_rate, spikes_stat(stim_num).ttest_p_success_rate,spikes_stat(stim_num).ttest_stats_success_rate]= ttest(spikes_stat(stim_num).success_rate(:,1),spikes_stat(stim_num).success_rate(:,2));   
        [spikes_stat(stim_num).wilcoxon_p_success_rate, spikes_stat(stim_num).wilcoxon_h_success_rate,spikes_stat(stim_num).wilcoxon_stats_success_rate]= signrank(spikes_stat(stim_num).success_rate(:,1),spikes_stat(stim_num).success_rate(:,2));
        [spikes_stat(stim_num).ttest_h_change_success_rate, spikes_stat(stim_num).ttest_p_change_success_rate]= ttest(spikes_stat(stim_num).change_success_rate);
        [spikes_stat(stim_num).wilcoxon_p_change_success_rate, spikes_stat(stim_num).wilcoxon_h_change_success_rate]= signrank(spikes_stat(stim_num).change_success_rate);
 
        %Modulation Index
        spikes_stat(stim_num).MI=MI;
        spikes_stat(stim_num).MI_m=nanmean(spikes_stat(stim_num).MI,1);
        spikes_stat(stim_num).MI_std=nanstd(spikes_stat(stim_num).MI,0,1);       
        %testing for normal distribution       
        [spikes_stat(stim_num).lillietest_h_MI, spikes_stat(stim_num).lillietest_p_MI] = lillietest(spikes_stat(stim_num).MI(:,1));
        %paired ttest 
        [spikes_stat(stim_num).ttest_h_MI, spikes_stat(stim_num).ttest_p_MI,spikes_stat(stim_num).ttest_stats_MI]= ttest(spikes_stat(stim_num).MI(:,1));   
        [spikes_stat(stim_num).wilcoxon_p_MI, spikes_stat(stim_num).wilcoxon_h_MI,spikes_stat(stim_num).wilcoxon_stats_MI]= signrank(spikes_stat(stim_num).MI(:,1));

         cd(path_output)
        save('extracellular_spikes_500ms','spikes','spikes_stat')

%% SNR
SNR_Y= spikes_stat(stim_num).SNR';
SNR_X(1,:)=ones(1,size(SNR_Y,2));
SNR_X(2,:)=2*ones(1,size(SNR_Y,2));
E = std(SNR_Y,0,2);
linex=[1;2];
my=max(max(SNR_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_SNR>0.05 
    asterisk='n.s.';
    a_fontsize=10;
else if spikes_stat.wilcoxon_p_SNR<0.05 && spikes_stat.wilcoxon_p_SNR>0.01
    asterisk='*';
    a_fontsize=17;
    else if spikes_stat.wilcoxon_p_SNR<0.01 && spikes_stat.wilcoxon_p_SNR>0.001
            asterisk='**';
            a_fontsize=17;
    else if spikes_stat.wilcoxon_p_SNR<0.001
             asterisk='***';
             a_fontsize=17;
        end
        end
    end
end
figure
hold on
line(SNR_X,SNR_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(SNR_X(:,1), mean(SNR_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','fontsize',a_fontsize) %'verticalAlignment','bottom',
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 1.2];
        y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('SNR index', 'FontSize', 28,'fontname', 'arial');
        title(['SNR Index, n=', num2str(size(SNR_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_SNR)])
        
        %save figure  
cd(path_output)
filename='SNR';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 
%%  Response Modulation  
Res_mod_Y =  spikes_stat(stim_num).res_modulation';
Res_mod_X(1,:)=ones(1,size(Res_mod_Y,2));
Res_mod_X(2,:)=2*ones(1,size(Res_mod_Y,2));
E=std(Res_mod_Y,0,2);
linex=[1;2];
my=max(max(Res_mod_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_res_modulation>0.05 
    asterisk='n.s.';
    a_fontsize=10;
else if spikes_stat.wilcoxon_p_res_modulation<0.05 && spikes_stat.wilcoxon_p_res_modulation>0.01
    asterisk='*';
    a_fontsize=17;
    else if spikes_stat.wilcoxon_p_res_modulation<0.01 && spikes_stat.wilcoxon_p_res_modulation>0.001
            asterisk='**';
            a_fontsize=17;
    else if spikes_stat.wilcoxon_p_res_modulation<0.001
             asterisk='***';
             a_fontsize=17;
        end
        end
    end
end
figure
hold on
line(Res_mod_X,Res_mod_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Res_mod_X(:,1), mean(Res_mod_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [-5 30];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'xtick', x1ticks,'ylim', y1limits,'fontsize',28,'linewidth',1,...%'ylim', y1limits,
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Response [#spikes/train]', 'FontSize', 28,'fontname', 'arial');
        title(['Response Modulation, n=', num2str(size(Res_mod_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_res_modulation)])
%% save figure  
cd(path_output)
filename='Response_modulation';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 
%%  Response Modulation  [Hz]
Res_mod_Y =  spikes_stat(stim_num).res_modulation_Hz';
Res_mod_X(1,:)=ones(1,size(Res_mod_Y,2));
Res_mod_X(2,:)=2*ones(1,size(Res_mod_Y,2));
E=std(Res_mod_Y,0,2);
linex=[1;2];
my=max(max(Res_mod_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_res_modulation_Hz>0.05 
    asterisk='n.s.';
    a_fontsize=10;
else if spikes_stat.wilcoxon_p_res_modulation_Hz<0.05 && spikes_stat.wilcoxon_p_res_modulation_Hz>0.01
    asterisk='*';
    a_fontsize=17;
    else if spikes_stat.wilcoxon_p_res_modulation_Hz<0.01 && spikes_stat.wilcoxon_p_res_modulation_Hz>0.001
            asterisk='**';
            a_fontsize=17;
    else if spikes_stat.wilcoxon_p_res_modulation_Hz<0.001
             asterisk='***';
             a_fontsize=17;
        end
        end
    end
end
figure
hold on
line(Res_mod_X,Res_mod_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(Res_mod_X(:,1), mean(Res_mod_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [-5 60];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'xtick', x1ticks,'ylim', y1limits,'fontsize',28,'linewidth',1,...%'ylim', y1limits,
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Spikes (Hz)', 'FontSize', 28,'fontname', 'arial');
        title(['Response Modulation, n=', num2str(size(Res_mod_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_res_modulation_Hz)])
%% save figure  
cd(path_output)
filename='Response_modulation_Hz';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 
%%  Modulation Index 
MI_Y =  spikes_stat(stim_num).MI';
MI_X(1,:)=ones(1,size(MI_Y,2));
% MI_X(2,:)=2*ones(1,size(MI_Y,2));
E=std(MI_Y,0,2);
linex=[0.95;1.05];
my=max(max(MI_Y))*1.2; 
liney=[my;my];
if spikes_stat.wilcoxon_p_MI>0.05 
    asterisk='n.s.';
    a_fontsize=10;
    a_fontsize=17;
else if spikes_stat.wilcoxon_p_MI<0.05 && spikes_stat.wilcoxon_p_MI>0.01
    asterisk='*';
    a_fontsize=17;
    else if spikes_stat.wilcoxon_p_MI<0.01 && spikes_stat.wilcoxon_p_MI>0.001
            asterisk='**';
            a_fontsize=17;
    else if spikes_stat.wilcoxon_p_MI<0.001
             asterisk='***';
             a_fontsize=17;
        end
        end
    end
end
figure
hold on
% line(MI_X,MI_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
% errorbar(MI_X(:,1), mean(MI_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% boxplot( spikes_stat(stim_num).MI);
scatter(MI_X,MI_Y,50,'markerfacecolor',color_table(1,:),'markeredgecolor',color_table(1,:),'markerfacealpha',0.3)
scatter(MI_X(1),mean(MI_Y),100,'markeredgecolor',color_table(2,:),'linewidth',2); %'markerfacecolor',color_table(2,:),
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a_fontsize)
hold off
        x1limits = [0.9 1.1];
        x1ticks = [];
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'xtick', x1ticks,'fontsize',20,'linewidth',1,...%'ylim', y1limits,
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{''},'box', 'off'); %'fontweight', 'bold', 
        ylabel('Modulation Index', 'FontSize', 20,'fontname', 'arial');
%         title(['Modulation Index, n=', num2str(size(MI_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_MI)])
%% save figure  
cd(path_output)
filename='Modulation index';
saveas(gcf,'Modulation index.fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 
%% signal nspikes
nspikes_signal_Y= spikes_stat(stim_num).nspikes_signal';
nspikes_signal_X(1,:)=ones(1,size(nspikes_signal_Y,2));
nspikes_signal_X(2,:)=2*ones(1,size(nspikes_signal_Y,2));
E = std(nspikes_signal_Y,0,2);
linex=[1;2];
my=max(max(nspikes_signal_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_nspikes_signal>0.05 
    asterisk='n.s.';
    a_fontsize=10;
else if spikes_stat.wilcoxon_p_nspikes_signal<0.05 && spikes_stat.wilcoxon_p_nspikes_signal>0.01
    asterisk='*';
    a_fontsize=17;
    else if spikes_stat.wilcoxon_p_nspikes_signal<0.01 && spikes_stat.wilcoxon_p_nspikes_signal>0.001
            asterisk='**';
            a_fontsize=17;
    else if spikes_stat.wilcoxon_p_nspikes_signal<0.001
             asterisk='***';
             a_fontsize=17;
        end
        end
    end
end
figure
hold on
line(nspikes_signal_X,nspikes_signal_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(nspikes_signal_X(:,1), mean(nspikes_signal_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','fontsize',a_fontsize) %'verticalAlignment','bottom',
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [-5 35];
        y1ticks = [-5,0,10,20,30];
        set( gca, 'xlim', x1limits, 'xtick', x1ticks,'ylim', y1limits,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('#Spikes/Train', 'FontSize', 28,'fontname', 'arial');
        title(['Signal #Spikes, n=', num2str(size(nspikes_signal_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_nspikes_signal)])
        
        %save figure  
cd(path_output)
filename='Signal Spikes';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 
%% background nspikes
nspikes_background_Y= spikes_stat(stim_num).nspikes_background';
nspikes_background_X(1,:)=ones(1,size(nspikes_background_Y,2));
nspikes_background_X(2,:)=2*ones(1,size(nspikes_background_Y,2));
E = std(nspikes_background_Y,0,2);
linex=[1;2];
my=max(max(nspikes_background_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_nspikes_background>0.05 
    asterisk='n.s.';
    a_fontsize=10;
else if spikes_stat.wilcoxon_p_nspikes_background<0.05 && spikes_stat.wilcoxon_p_nspikes_background>0.01
    asterisk='*';
    a_fontsize=17;
    else if spikes_stat.wilcoxon_p_nspikes_background<0.01 && spikes_stat.wilcoxon_p_nspikes_background>0.001
            asterisk='**';
            a_fontsize=17;
    else if spikes_stat.wilcoxon_p_nspikes_background<0.001
             asterisk='***';
             a_fontsize=17;
        end
        end
    end
end
figure
hold on
line(nspikes_background_X,nspikes_background_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(nspikes_background_X(:,1), mean(nspikes_background_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','fontsize',a_fontsize) %'verticalAlignment','bottom',
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [-5 20];
        y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('#Spikes/Train', 'FontSize', 28,'fontname', 'arial');
        title(['Background #Spikes, n=', num2str(size(nspikes_background_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_nspikes_background)])
        
        %% save figure  
cd(path_output)
filename='Background Spikes';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 
%% background nspikes (Hz)
nspikes_background_Y= spikes_stat(stim_num).nspikes_background_Hz';
nspikes_background_X(1,:)=ones(1,size(nspikes_background_Y,2));
nspikes_background_X(2,:)=2*ones(1,size(nspikes_background_Y,2));
E = std(nspikes_background_Y,0,2);
linex=[1;2];
my=max(max(nspikes_background_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_nspikes_background_Hz>0.05 
    asterisk='n.s.';
    a_fontsize=10;
else if spikes_stat.wilcoxon_p_nspikes_background_Hz<0.05 && spikes_stat.wilcoxon_p_nspikes_background_Hz>0.01
    asterisk='*';
    a_fontsize=17;
    else if spikes_stat.wilcoxon_p_nspikes_background_Hz<0.01 && spikes_stat.wilcoxon_p_nspikes_background_Hz>0.001
            asterisk='**';
            a_fontsize=17;
    else if spikes_stat.wilcoxon_p_nspikes_background_Hz<0.001
             asterisk='***';
             a_fontsize=17;
        end
        end
    end
end
figure
hold on
line(nspikes_background_X,nspikes_background_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(nspikes_background_X(:,1), mean(nspikes_background_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','fontsize',a_fontsize) %'verticalAlignment','bottom',
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [-5 40];
        y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Spikes (Hz)', 'FontSize', 28,'fontname', 'arial');
        title(['Background #Spikes, n=', num2str(size(nspikes_background_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_nspikes_background_Hz)])
        
        %% save figure  
cd(path_output)
filename='Background Spikes_Hz';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 
%% Latency (only for cells which faithfully responded to the 1st stim in the train)
latency_Y= spikes_stat(stim_num).latency';
latency_X(1,:)=ones(1,size(latency_Y,2));
latency_X(2,:)=2*ones(1,size(latency_Y,2));
E = std(latency_Y,0,2);
linex=[1;2];
my=max(max(latency_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_latency>0.05 
    asterisk='n.s.';
else if spikes_stat.wilcoxon_p_latency<0.05 && spikes_stat.wilcoxon_p_latency>0.01
    asterisk='*';
    else if spikes_stat.wilcoxon_p_latency<0.01 && spikes_stat.wilcoxon_p_latency>0.001
            asterisk='**';
    else if spikes_stat.wilcoxon_p_latency<0.001
             asterisk='***';
        end
        end
    end
end
figure
hold on
line(latency_X,latency_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(latency_X(:,1), mean(latency_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','fontsize',17) %'verticalAlignment','bottom',
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 1.1];
        y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Latency [mS]', 'FontSize', 28,'fontname', 'arial');
        title(['Response Latency, n=', num2str(size(latency_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_latency)])
        
        %save figure  
cd(path_output)
filename='Latency';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 
%% Latency STD (Jitter)
latency_std_Y= spikes_stat(stim_num).latency_std';
latency_std_X(1,:)=ones(1,size(latency_std_Y,2));
latency_std_X(2,:)=2*ones(1,size(latency_std_Y,2));
E = std(latency_std_Y,0,2);
linex=[1;2];
my=max(max(latency_std_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_latency_std>0.05 
    asterisk='n.s.';
else if spikes_stat.wilcoxon_p_latency_std<0.05 && spikes_stat.wilcoxon_p_latency_std>0.01
    asterisk='*';
    else if spikes_stat.wilcoxon_p_latency_std<0.01 && spikes_stat.wilcoxon_p_latency_std>0.001
            asterisk='**';
    else if spikes_stat.wilcoxon_p_latency_std<0.001
             asterisk='***';
        end
        end
    end
end
figure
hold on
line(latency_std_X,latency_std_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(latency_std_X(:,1), mean(latency_std_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [0 1.1];
        y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Latency STD [mS]', 'FontSize', 28,'fontname', 'arial');
        title(['Response Jitter, n=', num2str(size(latency_std_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_latency_std)])
        
        %save figure  
cd (path_output)
filename='Latency STD';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 
%% Success rate
success_rate_Y= spikes_stat(stim_num).success_rate';
success_rate_X(1,:)=ones(1,size(success_rate_Y,2));
success_rate_X(2,:)=2*ones(1,size(success_rate_Y,2));
E = std(success_rate_Y,0,2);
linex=[1;2];
my=max(max(success_rate_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_success_rate>0.05 
    asterisk='n.s.';
else if spikes_stat.wilcoxon_p_success_rate<0.05 && spikes_stat.wilcoxon_p_success_rate>0.01
    asterisk='*';
    else if spikes_stat.wilcoxon_p_success_rate<0.01 && spikes_stat.wilcoxon_p_success_rate>0.001
            asterisk='**';
    else if spikes_stat.wilcoxon_p_success_rate<0.001
             asterisk='***';
        end
        end
    end
end
figure
hold on
line(success_rate_X,success_rate_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(success_rate_X(:,1), mean(success_rate_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',17)
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [-0.1 1.2];
        y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Success Rate ', 'FontSize', 28,'fontname', 'arial');
        title(['Success Rate to First Stim, n=', num2str(size(success_rate_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_success_rate)])
        
        %save figure  
cd(path_output)
filename='Success Rate';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 

%% Failures
failures_Y= -1*(spikes_stat(stim_num).success_rate-1)';
failures_X(1,:)=ones(1,size(failures_Y,2));
failures_X(2,:)=2*ones(1,size(failures_Y,2));
E = std(failures_Y,0,2);
linex=[1;2];
my=max(max(failures_Y))*1.1; 
liney=[my;my];
if spikes_stat.wilcoxon_p_success_rate>0.05 
    asterisk='n.s.';
else if spikes_stat.wilcoxon_p_success_rate<0.05 && spikes_stat.wilcoxon_p_success_rate>0.01
    asterisk='*';
    else if spikes_stat.wilcoxon_p_success_rate<0.01 && spikes_stat.wilcoxon_p_success_rate>0.001
            asterisk='**';
    else if spikes_stat.wilcoxon_p_success_rate<0.001
             asterisk='***';
        end
        end
    end
end
figure
hold on
line(failures_X,failures_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(failures_X(:,1), mean(failures_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
hold off

        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits = [-0.1 1.1];
        y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',{'NB-','NB+'} ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Failures ', 'FontSize', 28,'fontname', 'arial');
        title(['Failures to First Stim, n=', num2str(size(failures_Y,2)), ', p=', num2str(spikes_stat.wilcoxon_p_success_rate)])
        
        %save figure  
cd(path_output) 
filename='Failures';
saveas(gcf,filename,'fig'); 
print(gcf,filename,'-dpng','-r600','-opengl') 