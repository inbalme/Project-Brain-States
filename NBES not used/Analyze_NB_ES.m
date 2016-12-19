%% Electrical Stimulation of Nucleus Basalis
clear all
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
% cd 'D:\Inbal M.Sc\Data PhD\DoubleNB\Extracted Data';
% load 2015-12-31-009_h1-4
% load 2015-12-28-006_h1-6
% load 2015-07-29-002_h1-6
%  load Set2-2015-04-02-001_h3-12
 load Set2-2015-05-06-003_h2-7
%  load Set2-2015-06-01-003_h1-6
%%
 orig_header=7;
 galvano_flag=0;
 header=find(headers==orig_header);
 stim2_X=[];
 sf = 4000;
 dt=1/sf;
 sf_galvano = 800;
 dt_galvano = 1/sf_galvano;
 color_table = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];
 
 mean_trace =  mean(data(header).x_value(1).data1(:,:),2);
 std_trace =  std(data(header).x_value(1).data1(:,:),0,2);
 galvano2_trace = data(header).x_value(1).galvano2(:,2);
 time_axis_galvano =length(galvano2_trace)*dt_galvano;
 galvano_vec = zeros(length(galvano2_trace),1);
 galvano_sub_mean = galvano2_trace - mean(galvano2_trace);
%  galvano_threshold=0.5;
%  galvano_vec(galvano2_trace > galvano_threshold)=1;   %turning the galvano trace into binary
 galvano_threshold = abs(mean(galvano_sub_mean)) + 4.*abs(std(galvano_sub_mean));
 galvano_vec(abs(galvano_sub_mean) > galvano_threshold)=1;   %turning the galvano trace into binary
 
%                  if isempty(find(galvano_vec(:,1)==1, 1))~=1 %condition will be fulfilled if there was galvano activation.
                if galvano_flag==1;
                    galvano_vec_shifted(:,1) = galvano_vec(2:end,1)-galvano_vec(1:end-1, 1);
                    galvano_begin(:,1) = find(galvano_vec_shifted(:,1)==1); %find the locations where galvano pulse starts
                    galvano_begin(:,1) = galvano_begin(:,1)+1; %correction for the shift
                    galvano_end(:,1) = find(galvano_vec_shifted(2:end,1)==-1); %find the locations where galvano pulse ends
                        if length(galvano_begin(:,1))>length(galvano_end(:,1)) %if galvano_begin is larger than galvano_end, take only the points in galvano_begin which has a match in galvano_end
                            galvano_begin_trunc(:,1)=galvano_begin(1:length(galvano_end(:,1)),1);
                            galvano_begin(:,1)=[];
                            galvano_begin=galvano_begin_trunc(:,1);
                        end
                    locations_x_galvano(1,:) = galvano_begin(:,1); %arranging the galvano begin and end locations in one variable for plotting
                    locations_x_galvano(2,:) = galvano_end(:,1);         
                    
                     if isempty(locations_x_galvano)
                            stim2_X = [];
                        else
                            stim2_X = locations_x_galvano(:,:);
                        end
                end
 %%
 figure
 subplot(2,1,1)
    plot((1:size(mean_trace,1)).*dt,mean_trace(:,:), 'LineWidth',0.2,'color', color_table(1,:));
 subplot(2,1,2)
 plot((1:size(std_trace,1)).*dt,std_trace(:,:), 'LineWidth',0.2,'color', color_table(1,:));
 
figure
 subplot(3,1,1)
 hold on
    plot((1:size(mean_trace(1.9*sf:2.2*sf,:),1)).*dt,mean_trace(1.9*sf:2.2*sf,:), 'LineWidth',0.2,'color', color_table(1,:));
    plot((1:size(mean_trace(11.9*sf:12.2*sf,:),1)).*dt,mean_trace(11.9*sf:12.2*sf,:), 'LineWidth',0.2,'color', color_table(2,:));
    hold off
 subplot(3,1,2)
 hold on
 plot((1:size(std_trace(1.9*sf:2.2*sf,:),1)).*dt,std_trace(1.9*sf:2.2*sf,:,:), 'LineWidth',0.2,'color', color_table(1,:)); 
 plot((1:size(std_trace(11.9*sf:12.2*sf,:),1)).*dt,std_trace(11.9*sf:12.2*sf,:), 'LineWidth',0.2,'color', color_table(2,:)); 
 hold off
 subplot(3,1,3)
 hold on
 plot((1:size(galvano2_trace(1.9*sf_galvano:2.2*sf_galvano,:),1)).*dt_galvano,galvano2_trace(1.9*sf_galvano:2.2*sf_galvano,:), 'LineWidth',0.2,'color', color_table(1,:)); 
 plot((1:size(galvano2_trace(11.9*sf_galvano:12.2*sf_galvano,:),1)).*dt_galvano,galvano2_trace(11.9*sf_galvano:12.2*sf_galvano,:), 'LineWidth',0.2,'color', color_table(2,:)); 
 hold off
 
 figure
 subplot(3,1,1)
 hold on
    plot((1:size(mean_trace(1.5*sf:3.5*sf,:),1)).*dt,mean_trace(1.5*sf:3.5*sf,:), 'LineWidth',0.2,'color', color_table(1,:));
    plot((1:size(mean_trace(11.5*sf:13.5*sf,:),1)).*dt,mean_trace(11.5*sf:13.5*sf,:), 'LineWidth',0.2,'color', color_table(2,:));
    hold off
 subplot(3,1,2)
 hold on
 plot((1:size(std_trace(1.5*sf:3.5*sf,:),1)).*dt,std_trace(1.5*sf:3.5*sf,:), 'LineWidth',0.2,'color', color_table(1,:)); 
 plot((1:size(std_trace(11.5*sf:13.5*sf,:),1)).*dt,std_trace(11.5*sf:13.5*sf,:), 'LineWidth',0.2,'color', color_table(2,:)); 
 hold off
  subplot(3,1,3)
 hold on
 plot((1:size(galvano2_trace(1.5*sf_galvano:3.5*sf_galvano,:),1)).*dt_galvano,galvano2_trace(1.5*sf_galvano:3.5*sf_galvano,:,:), 'LineWidth',0.2,'color', color_table(1,:)); 
 plot((1:size(galvano2_trace(11.5*sf_galvano:13.5*sf_galvano,:),1)).*dt_galvano,galvano2_trace(11.5*sf_galvano:13.5*sf_galvano,:,:), 'LineWidth',0.2,'color', color_table(2,:)); 
 hold off

%% Variance plot - one x_value at a time, for flexible data
       x_value = 1;
       orig_header=2;
       header=find(headers==orig_header);
        plot_data = data(header).x_value(1).data1(:,1:15); 
        plot_data_mean = mean(plot_data,2);
        plot_data_std =  std(plot_data,0,2);
        plot_data_CV = (-1).*plot_data_std./plot_data_mean; 

       trace_ind = [1:size(plot_data,2)] ; %1:size(data_BP(:,:,laser_flag==1),2); %trace_ind is the numbers of traces to be plotted
       l=size(plot_data(:,trace_ind, x_value),1);
        l=l/2-1;
        DC=20.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);
        trace_to_plot = plot_data(:,trace_ind, x_value)+DC ;
        
            [Fig1,h1]= fn_Plot_Trace_v2(trace_to_plot, dt, dt_galvano, stim2_X);
            [Fig2,h2]= fn_Plot_Trace_v2(plot_data_std(:,x_value), dt, dt_galvano, stim2_X);
            [Fig3,h3] = fn_Plot_Trace_v2(plot_data_mean(:,x_value), dt, dt_galvano, stim2_X);


        ax1 = get(Fig1, 'children');
        pos1 = [0.08 , 0.7 , 0.8 , 0.2];
        top1 = pos1(1,2)+pos1(1,4);

        ax2 = get(Fig2, 'children');
        pos2 = [0.08 , 0.4 , 0.8 , 0.2];
        top2 = pos2(1,2)+pos2(1,4);

        ax3 = get(Fig3, 'children');
        pos3 = [0.08 , 0.1 , 0.8 , 0.2];
        top3 = pos3(1,2)+pos3(1,4);

        F = figure;
        set(gcf,'color','w');
        set(gcf,'DefaultAxesFontSize',12);
        set(gcf,'DefaultAxesFontName','helvetica');
        set(gcf, 'PaperType', 'A4');
        set(gcf,'PaperUnits','centimeters','PaperPosition',[1.2 1.2 20 27]); %[left, bottom, width, height] 
        set(gcf,'PaperOrientation','portrait');
        set(gcf,'Units','centimeters','Position',get(gcf,'paperPosition')+[0 0 0 0]);

        ax_copy1 = copyobj(ax1,F); % ax1 to new fig
        set(ax_copy1(1),'position',pos1(1,:)) % Set its position  

        ax_copy2 = copyobj(ax2,F); % ax2 to new fig
        set(ax_copy2(1),'position',pos2(1,:)) % Set its position  

        ax_copy3 = copyobj(ax3,F); % ax3 to new fig
        set(ax_copy3(1),'position',pos3(1,:)) % Set its position  
 %% Plot power spectrum for channel data1
orig_header=[]; header=[];
figure
 orig_header=[2,7]; %[2,4:7]
 for count=1:length(orig_header)
 header(count)=find(headers==orig_header(count));
 
start_time=[7,12]; %[sec]
duration = 2; %[sec]
x_value = 1;
Y_abs = []; f = [];

  for t=1:length(start_time);
DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = []; interval = []; 

start_sample = start_time(t).*sf;
if start_time(t)==0
    start_sample = 1;
end
end_sample = start_sample+duration.*sf-1;
interval = start_sample:end_sample;
spec_mat = data(header(count)).x_value(1).data2(interval,:);
DC= mean(spec_mat,1);
l=size(spec_mat,1);
l=l/2-1;
DC=wextend('addrow', 'per', DC, l);
spec_mat_noDC = spec_mat-DC;
spec_mat_mean = mean(spec_mat,2);

[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf);
  end
  
%  % Plot power spectrum with shaded error bar
%  for t=1:length(start_time);
%    figure(gcf)
%      ax(count)=subplot(2,1,count);
%       hold on
%    Fig1= fn_shadedErrorBar(f(:,t),mean(Y_abs(:,:,t),2),std(Y_abs(:,:,t),0,2),{'LineWidth',1,'color', color_table(t,:)},1);
%     xlim([0 20]); ylim([-0.1 0.9])
%     set(gca,'fontsize',14,'fontname','myriad','fontweight','bold')
%  end
 
  % Plot power spectrum without shaded error bar

     figure(gcf)
     ax(count)=subplot(2,1,count);
for t=1:length(start_time);
      hold on
        plot(f(:,t),mean(Y_abs(:,:,t),2),'color', color_table(t,:)) 
      hold off
        xlim([0 20]); ylim([-0.01 1])
%         title('')
        set(gca,'fontsize',20, 'fontname', 'myriad', 'box', 'off','linewidth',2)
        xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
        ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
 end
 
 end
 
%         ax1 = get(gcf, 'children');
%         pos1 = [0.08 , 0.7 , 0.8 , 0.2];
%         F=figure
%         ax_copy1 = copyobj(ax1,F); % ax1 to new fig
%         set(ax_copy1(1),'position',pos1(1,:)) % Set its position  
  %% Plot power spectrum for channel EEG
orig_header=[]; header=[]; sf_EEG=sf_galvano;
figure
 orig_header=[3:12]; %[2,4:7]
 for count=1:length(orig_header)
 header(count)=find(headers==orig_header(count));
 
start_time=[7,12]; %[sec]
duration = 2; %[sec]
x_value = 1;
Y_abs = []; f = [];

  for t=1:length(start_time);
DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = []; interval = []; 

start_sample = start_time(t).*sf_EEG;
if start_time(t)==0
    start_sample = 1;
end
end_sample = start_sample+duration.*sf_EEG-1;
interval = start_sample:end_sample;
spec_mat = data(header(count)).x_value(1).EEG(interval,:); 
DC= mean(spec_mat,1);
l=size(spec_mat,1);
l=l/2-1;
DC=wextend('addrow', 'per', DC, l);
spec_mat_noDC = spec_mat-DC;
spec_mat_mean = mean(spec_mat,2);

[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf);
  end
  
  
 for t=1:length(start_time);
     figure(gcf)
     ax(count)=subplot(5,2,count);
      hold on
    fn_shadedErrorBar(f(:,t),mean(Y_abs(:,:,t),2),std(Y_abs(:,:,t),0,2),{'LineWidth',1,'color', color_table(t,:)},1);
    xlim([0 35]); %ylim([-0.05 0.05])
 end
 end
%% Plot single traces of LFP
x_value=1; gain=20; count=0; figure
for header=[6,8,10];
    count=count+1;
 plot_data = data(header).x_value(1).data1(:,:); 
%  plot_data(9.95*sf:10.55*sf,:) = nan; %plot_data(9.95*sf,:); %putting nan instead of the ES artifact
 
  trace_ind = [1,3,6]; %[1:size(plot_data,2)] ;  %trace_ind is the numbers of traces to be plotted
       l=size(plot_data(:,1),1);
        l=l/2-1;
        DC=20.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l); 
        trace_to_plot = (plot_data(:,trace_ind)+DC)./gain ;
      
        figure(gcf)
        subplot(3,1,(count))
        for i=1:length(trace_ind)         
            hold on    
                plot((1:size(trace_to_plot,1))'.*dt,trace_to_plot(:,i),'k', 'linewidth',1)
            hold off
%             xlabel('Time [sec]' ,'FontSize', 14);
%             ylabel('Vm [mV]', 'FontSize', 14);
            set(gca,'ylim', [-1 3])
        end
end
%% 
filename = '20150402-001_LFP raw data traces_depth_summary_selection_+ES'; 
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\LFP depth';

print (1, '-depsc2', filename)
%%
figure(7)
set(gca,'position',[0.2,0.2,0.5,0.4])
set(gca,'fontsize',14,'fontname','myriad','fontweight','bold')