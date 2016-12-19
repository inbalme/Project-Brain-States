%% Electrical Stimulation of Nucleus Basalis
clear all; close all

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
 load Set2-2015-04-02-001_h3-12
 
%%
 galvano_flag=0;
 stim2_X=[];
 sf = 4000;
 dt=1/sf;
 sf_galvano = 800;
 dt_galvano = 1/sf_galvano;
 color_table = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];
 %% Plot power spectrum for channel data1
orig_header=[]; header=[];
figure
 orig_header=[6,10]; %[2,4:7]
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

[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf,0,0);
  end 
 
  % Plot power spectrum without shaded error bar

     figure(gcf)
     ax(count)=subplot(2,1,count);
for t=1:length(start_time);
      hold on
        plot(f(:,t),mean(Y_abs(:,:,t),2),'color', color_table(t,:)) 
      hold off
        xlim([0 20]); ylim([-0.01 0.5])
%         title('')
        set(gca,'fontsize',20, 'fontname', 'myriad', 'box', 'off','linewidth',2)
        xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
        ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
 end
 
 end

%% Plot single traces of LFP
x_value=1; gain=20; count=0; figure
for header=[6,10];
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
% filename = '20150402-001_LFP raw data traces_depth_summary_selection_+ES'; 
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\LFP depth';
% 
% print (1, '-depsc2', filename)