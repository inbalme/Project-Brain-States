%% Analyze extracellular spikes

clear all
close all
global Exp
cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';

% load Set2-2013-01-09-001  %headers 1; 1 cluster
% load Set2-2013-02-27-002    %headers 3-7, SNR not good enough for clustering...
% load Set2-2013-02-27-003  %headers 1-11; 3 clusters h4,11-2Hz,h8-6Hz,h9-10Hz. choose cluster, light_freq 13
load Set2-2013-03-20-002  %headers 1-12; 2 clusters h2-2Hz, h5-5Hz,h6-6Hz,h10-10Hz,h12-15Hz choose cluster, light_freq 26 
% load Set2-2013-04-03-002  %headers 1-9; 1 clusters h2-2Hz,h3-5Hz, h5-10Hz, light_freq 13

cd 'C:\Users\Inbal\Dropbox\Inbal M.Sc\MATLAB\Project Brain States'
% cd 'C:\Users\inbalme\Dropbox\Inbal M.Sc\MATLAB\Project Brain States'
% cd 'D:\Inbal M.Sc\MATLAB\Project Brain States'

for header = 10;
x_value = 1;
channel = 1; % V1 - 1, Ext1 - 2
clusters = 2; %according to the number of neurons to be sorted

exp_type = Param.header(header).general.exp_type;
exp_type_logical = strcmp('Standby', exp_type);
% if exp_type_logical
%     x_value = 1;
% end
light_freq = -1.*Param.header(header).stim.facade(26); %[Hz] if there was airpuff or galvano then facade(26)
switch channel

    case 1
        raw_data = data.header(header).x_value(x_value).Vm;
        sf = Param.header(1,header).stim.sf;
        

    case 2
        raw_data = -1.*data.header(header).x_value(x_value).Ext1; %for extracellular recording, spikes are upside-down
        sf = Param.header(1,header).stim.sf_Ext1;

        end
        
        dt=1/sf; %[sec]
        time_axis = (1:size(raw_data,1))*dt;
        sf_laser = Param.header(header).stim.sf_laser; %[1/sec]
        dt_laser = 1/sf_laser;
        time_axis_laser = (1:size(data.header(header).x_value(x_value).laser,1))*dt_laser; 
        laser_vec = zeros(length(data.header(header).x_value(x_value).laser(:,1)),1);
        laser_vec(data.header(header).x_value(x_value).laser(:,1) > 2)=1;   %turning the laser trace into binary
        
        laser_vec_shifted = laser_vec(2:end)-laser_vec(1:end-1);
        laser_begin = find(laser_vec_shifted==1); %find the locations where laser pulse starts, ACCORDING TO LASER SF
        laser_begin = laser_begin+1; %correction for the shift
        laser_end = find(laser_vec_shifted==-1); %find the locations where laser pulse ends
        locations_x_laser(1,:) = laser_begin; %arranging the laser begin and end locations in one variable for plotting
        locations_x_laser(2,:) = laser_end;
      

high_pass_freq = 300; %[Hz]
low_pass_freq = 3000; %[Hz]
BaselineTime = [];
robustness = 3;
peaks = [];
interpeak_time = 1; %[ms]

if isempty(clusters)
    clusters = input('enter the number of clusters according to which the spikes will be grouped');
end

color_table = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];

% parameters for preparing vectors at the same length of intervals around
% the peaks:
begin_peak_time = 0.5; %[ms]
end_peak_time = 1; %[ms]
interval_begin = ceil(begin_peak_time./1000.*sf);
interval_end = ceil(end_peak_time./1000.*sf);
interval_around_peak = interval_begin + interval_end +1;
peak_counter = 0;

for i = 1:clusters
    raster(:,:,i) = zeros(size(raw_data));
end

[data_HP] = fn_High_Pass (raw_data, sf, high_pass_freq);
[data_HP] = fn_Low_Pass (data_HP, sf, low_pass_freq);
   
%% Plotting a raw data trace
header_to_plot = header;
trace_to_plot = [10 15];
stim_X = locations_x_laser;

[Fig] = fn_Plot_Trace(data_HP, header_to_plot, trace_to_plot, dt, dt_laser, stim_X);

%% 
for trace = 1:size(data_HP,2)
    data_HP_trace = data_HP(:,trace);
    
    threshold_peak(1,trace) = fn_Threshold_Peak (data_HP_trace, sf, BaselineTime, robustness);
%     threshold_peak(1,trace)= 0.2;
    
%find peaks that are above the threshold and with minimal distance between
%them specified by "threshold_peak" and "interpeak_time":
    [peaks(1,trace).Value,peaks(1,trace).Location] = fn_Detect_Spike_HPF ...
        (data_HP_trace, sf, interpeak_time, threshold_peak(1,trace));
    
    % take an interval around each peak of 0.5ms before and 1ms after the peak,
% for PCA.
    if isempty(peaks(1,trace).Location)
        continue
    else
        peaks(1,trace).Begin = peaks(1,trace).Location - interval_begin;
        peaks(1,trace).End = peaks(1,trace).Location + interval_end;
        
            for peak_ind = 1:length(peaks(1,trace).Value)
                if peaks(1,trace).Begin(peak_ind) <= 0 || peaks(1,trace).End(peak_ind) > length(data_HP_trace)...
                        peaks(1,trace).Peak_Mat(peak_ind,1:interval_around_peak) = NaN;
                else
                    peaks(1,trace).Peak_Mat(peak_ind,:) = data_HP_trace(peaks(1,trace).Begin(peak_ind):peaks(1,trace).End(peak_ind));
                    
                    peak_counter = peak_counter+1;
                    peaks_mat(peak_counter,:) = peaks(1,trace).Peak_Mat(peak_ind,:); %a matrix, each row is a peak
                    peaks_trace(peak_counter, 1) = trace;
                    peaks_number_in_trace(peak_counter, 1) = peak_ind;
                    peaks_loc(peak_counter, 1) = peaks(1,trace).Location(peak_ind);
                    peaks_value(peak_counter, 1) = peaks(1,trace).Value(peak_ind);
                    peaks_begin(peak_counter, 1) = peaks(1,trace).Begin(peak_ind);
                    peaks_end(peak_counter, 1) = peaks(1,trace).End(peak_ind);
%                 
                    
                end
            end
    end
end    
% %% Spike sorting using PCA, fitting gaussians and clustering.
% 
% [U,S,V] = svds(peaks_mat,2);
% 
% obj  = gmdistribution.fit(U,clusters, 'Replicates',500); %fitting the scattered dots into gaussians.
% idx = cluster(obj,U); %returns a vector same length as U, with the group number for each entry according to the clustering

%% PCA
[COEFF,SCORE] = princomp(peaks_mat);
U = SCORE(:,1:clusters);
obj  = gmdistribution.fit(U,clusters, 'Replicates',500); %fitting the scattered dots into gaussians.
idx = cluster(obj,U); %returns a vector same length as U, with the group number for each entry according to the clustering
%% Plotting the PCA and the sorted spikes

if clusters~=1;
    
figure
set(gcf,'color','w')
              
subplot(2,1,1)
     set( gca, 'ticklength', [0.010 0.010],'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('PC1' ,'FontSize', 12);
        ylabel('PC2', 'FontSize', 12);
hold on
for sorted = 1:clusters
    scatter(U(idx==sorted,1),U(idx==sorted,2), 3, color_table(sorted,:),'filled');
end
hold off

subplot(2,1,2)
     set( gca, 'ticklength', [0.010 0.010],'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );       
        xlabel('Time [ms]' ,'FontSize', 12);
        ylabel('V [mV]', 'FontSize', 12);
hold on
for sorted = 1:clusters
    plot((1:size(peaks_mat,2))./sf*1000,peaks_mat(idx==sorted,:),'color', color_table(sorted,:))
end
hold off

else
end
%% making rasters for the sorted spikes
% raster is a 3D matrix where rasters (x and y axes) are stored for
% different sorted neurons (z axis). the matrices are sparse (all 0,
% except for 1 at spike times.

for peak = 1:length(peaks_loc)
    for rast = 1:clusters
        if idx(peak) == rast
            raster(peaks_loc(peak), peaks_trace(peak),rast) = 1;
        end
    end
end

%% Plotting the rasters and PSTHs
bin_time = 10; %[ms]
bin_size = ceil(bin_time./1000.*sf);

for rast_ind = 1:clusters
    A = raster(:,:,rast_ind)'; % convert from sparse to full
    % Plot a line on each spike location 
    [M, N] = size(A);
    [X,Y] = meshgrid(1:N,0:M-1);
    locations_X(1,:) = X(:);
    locations_X(2,:) = X(:);
    locations_Y(1,:) = [Y(:)+1].*A(:);
    locations_Y(2,:) = [Y(:)+1.5].*A(:);
    indxs = find(locations_Y(1,:) ~= 0);
    locations_x = locations_X(:,indxs);
    locations_y = locations_Y(:,indxs);
    locations_y_laser_raster = ones(size(locations_x_laser)).*(max(locations_Y(2,:))+6);   
        clear locations_X locations_Y
        
    PSTH(:,:,rast_ind) = fn_PSTH(bin_time, bin_size, raster(:,:,rast_ind));
    locations_y_laser_PSTH = ones(size(locations_x_laser)).*(max(PSTH(:,:,rast_ind))+4);   
%%
    figure
    set(gcf,'color','w')
    
    position_raster = [0.1 , 0.42 , 0.8 , 0.3];
    position_PSTH = [0.1 0.1 0.8 0.3];
        
        x1limits = [0 size(A,2).*dt];
        x1ticks = [];
        y1limits = [0 size(A,1)+8];
        y1ticks = [];

        axes('position', position_raster);
        set( gca, 'xlim', x1limits, 'ylim', y1limits,'xtick', x1ticks,...
        'ticklength', [0.010 0.010],'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
%         title(['file', num2str(Param.name)] ,'FontSize', 14);
%         xlabel('Time [sec]' ,'FontSize', 12);
        ylabel('Trial no.', 'FontSize', 12);

    hold on
    line(locations_x.*dt,locations_y,'LineWidth',2,'color', color_table(rast_ind,:))    %'Color','k')
    line(locations_x_laser.*dt_laser,locations_y_laser_raster,'LineWidth',6,'Color','c')
    hold off
    
        x2limits = [0 size(A,2).*dt];
        x2ticks = x2limits; %[0:size(A,2).*dt/10:size(A,2).*dt];
        y2limits = [0 max(PSTH(:,:,rast_ind))+6];
        y2ticks = [];
   
   axes('position', position_PSTH);
   set(gca, 'xlim', x2limits, 'ylim', y2limits,'xtick', x2ticks,'XMinorTick','on', ...
        'ticklength', [0.010 0.010],'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );      
% title(['file', num2str(Param.name), ' header', num2str(data.header(header).x_value(x_value).original_header), ' PSTH Spikes ', ...
%           'bin size [ms]=', num2str(bin_time)] ,'FontSize', 14);
xlabel('Time [sec]' ,'FontSize', 12);
ylabel('#Spikes/sec', 'FontSize', 12);
hold on
bar(time_axis(bin_size:bin_size:end),PSTH(:,:,rast_ind), 'k')
% line(locations_x_laser.*dt_laser,locations_y_laser_PSTH,'LineWidth',6,'Color','c')
hold off
end

%% plot zoom-in near pulse
% for pulse_num = 1:length(laser_begin)
pulse_num = 5;
time_of_pulse = laser_begin.*dt_laser.*1000; %[ms]
pre_pulse_time = time_of_pulse-4; %time in ms
post_pulse_time = time_of_pulse+16; %time in ms
zoom_duration = post_pulse_time(pulse_num)-pre_pulse_time(pulse_num);
location_pulse = ceil(time_of_pulse.*sf./1000); %laser locations IN Vm SF!
zoom_interval_begin = ceil(pre_pulse_time.*sf./1000); %[samples]
zoom_interval_end = ceil(post_pulse_time.*sf./1000); %[samples]


for rast_ind = 1:clusters
A_zoom = raster(zoom_interval_begin(pulse_num):zoom_interval_end(pulse_num),:, rast_ind)'; % convert from sparse to full
    % Plot a line on each spike location 
    [M, N] = size(A_zoom);
    [X,Y] = meshgrid(1:N,0:M-1);
    locations_X(1,:) = X(:);
    locations_X(2,:) = X(:);
    locations_Y(1,:) = [Y(:)+1].*A_zoom(:);
    locations_Y(2,:) = [Y(:)+1.5].*A_zoom(:);
    indxs = find(locations_Y(1,:) ~= 0);
    locations_x_zoom = locations_X(:,indxs);
    locations_x_zoom = locations_x_zoom + zoom_interval_begin(pulse_num);
    locations_y_zoom = locations_Y(:,indxs);
    locations_y_laser_raster = ones(size(locations_x_laser(:,pulse_num))).*(max(locations_Y(2,:))+4); 
        clear locations_X locations_Y
figure
    set(gcf,'color','w')
        position_raster = [0.1 , 0.1 , 0.8 , 0.3];
        
        x3limits = [zoom_interval_begin(pulse_num).*dt zoom_interval_end(pulse_num).*dt];
        x3ticks = x3limits; %[0 zoom_duration];
        x1ticklabels = [0 zoom_duration];
        y3limits = [0 size(A_zoom,1)+8];
        y3ticks = [];   
        
    axes('position', position_raster);
        set( gca, 'xlim', x3limits, 'ylim', y3limits,'xtick', x3ticks, 'XTickLabel', x1ticklabels, 'XMinorTick','on',...
        'ticklength', [0.010 0.010],'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
%         title(['file', num2str(Param.name)] ,'FontSize', 14);
        xlabel('Time [ms]' ,'FontSize', 12);
        ylabel('Trial no.', 'FontSize', 12);

    hold on
    line(locations_x_zoom.*dt,locations_y_zoom,'LineWidth',2,'color', color_table(rast_ind,:))    %'Color','k')
    line(locations_x_laser(:,pulse_num).*dt_laser,locations_y_laser_raster,'LineWidth',3,'Color','c')
    hold off

end

%% Calculating latency and variance for each light pulse in a given frequency

for rast_ind = 1:clusters
    for pulse = 1:length(laser_begin);
        range = location_pulse+20.*sf./1000;
        raster_pulse(:,:,rast_ind) = raster(location_pulse(pulse):range(pulse),:, rast_ind);
        for trace = 1:size(data_HP,2)
            if isempty(find(raster_pulse(:,trace, rast_ind)==1))
                 cluster(rast_ind).pulse(1,pulse).latency(trace) = 0;
            else
                cluster(rast_ind).pulse(1,pulse).latency(trace) = find(raster_pulse(:,trace, rast_ind)==1,1,'first'); 
            end
        end
        cluster(rast_ind).pulse(1,pulse).latency = cluster(rast_ind).pulse(1,pulse).latency(cluster(rast_ind).pulse(1,pulse).latency~=0);
        cluster(rast_ind).pulse(1,pulse).latency_mean = mean(cluster(rast_ind).pulse(1,pulse).latency);
        cluster(rast_ind).pulse(1,pulse).latency_std = std(cluster(rast_ind).pulse(1,pulse).latency);
    end
end
 %%
% cell(1).freq(1).pulse = cluster(1).pulse;
% cell(1).file_name = Param.name;
% cell(1).sf = sf;
%     fname = Param.name;
%     name = [fname(1:(strfind(fname,'.')-1)),'h',num2str(header),'freq',num2str(light_freq)];
% cd 'C:\Users\Inbal\Dropbox\Inbal M.Sc\Data PhD\ChAT Extracellular'
%     save( name, 'cell', 'dt')  


end
%%
cd 'C:\Users\Inbal\Dropbox\Inbal M.Sc\Data PhD\ChAT Extracellular'
load Set2-2013-02-27-003h4freq2.mat
cell_lat(1) = cell;
load Set2-2013-02-27-003h7freq5.mat
cell_lat(1).freq(2) = cell.freq;
load Set2-2013-02-27-003h8freq6.mat
cell_lat(1).freq(3) = cell.freq;
load Set2-2013-02-27-003h9freq10.mat
cell_lat(1).freq(4) = cell.freq;
load Set2-2013-03-20-002h2freq2.mat
cell_lat(2) = cell;
load Set2-2013-03-20-002h5freq5.mat
cell_lat(2).freq(2) = cell.freq;
load Set2-2013-03-20-002h6freq6.mat
cell_lat(2).freq(3) = cell.freq;
load Set2-2013-03-20-002h10freq10.mat
cell_lat(2).freq(4) = cell.freq;
load Set2-2013-03-20-002h12freq15.mat
cell_lat(2).freq(5) = cell.freq;
load Set2-2013-04-03-002h2freq2.mat
cell_lat(3) = cell;
load Set2-2013-04-03-002h3freq5.mat
cell_lat(3).freq(2) = cell.freq;
load Set2-2013-04-03-002h5freq10.mat
cell_lat(3).freq(3) = cell.freq;

cell_lat(1).frequencies = [2 5 6 10];
cell_lat(2).frequencies = [2 5 6 10 15];
cell_lat(3).frequencies = [2 5 10];
%% Latencies and STDs
% important: the latencies and variance and std are not in units of time,
% but number of samples! need to change...
for cell_ind = 1:length(cell_lat);
  for freq = [2,5,10,15]
   if isempty(find(cell_lat(cell_ind).frequencies==freq))
       continue
   else
    for pulse = 1:length(cell_lat(cell_ind).freq(find(cell_lat(cell_ind).frequencies==freq)).pulse)
        mean_latency{cell_ind,freq}(pulse) = cell_lat(cell_ind).freq(find(cell_lat(cell_ind).frequencies==freq)).pulse(pulse).latency_mean;
        std_latency{cell_ind,freq}(pulse) = cell_lat(cell_ind).freq(find(cell_lat(cell_ind).frequencies==freq)).pulse(pulse).latency_std;
        var_latency{cell_ind,freq}(pulse) = var(cell_lat(cell_ind).freq(find(cell_lat(cell_ind).frequencies==freq)).pulse(pulse).latency);
    end
   end
  end
end

%% Ploting Latencies:
max_latency(1) = max([max(mean_latency{1,2}) max(mean_latency{2,2}) max(mean_latency{3,2}(1:2))]);
max_latency(1) = ceil(max_latency(1)./10).*10;
min_latency(1) = min([min(mean_latency{1,2}) min(mean_latency{2,2}) min(mean_latency{3,2}(1:2))]);
min_latency(1) = floor(min_latency(1)./10).*10;
x1limits = [1 2];
 x1ticks = [1:2];
 y1limits = [min_latency(1) max_latency(1)];
 y1ticks = y1limits;
 
figure
set( gca, 'xlim', x1limits, 'xtick', x1ticks,'ylim', y1limits, 'ytick', y1ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],...
        'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('Pulse No.' ,'FontSize', 12);
        ylabel('Mean Latency', 'FontSize', 12);
hold on
    plot([1:2],mean_latency{1,2},...
        '-s','markersize',5,'MarkerFaceColor',[1 0 0],'color',[1 0 0],'linewidth',2)
    plot([1:2],mean_latency{2,2},...
        '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
    plot([1:2],mean_latency{3,2}(1:2),'-s','markersize',5,...
        'MarkerFaceColor',[0 0 1],'color',[0 0 1],'linewidth',2)
hold off

max_latency(2) = max([max(mean_latency{1,5}(1:5)) max(mean_latency{2,5}(1:5)) max(mean_latency{3,5}(1:5))]);
max_latency(2) = ceil(max_latency(2)./10).*10;
min_latency(2) = min([min(mean_latency{1,5}(1:5)) min(mean_latency{2,5}(1:5)) min(mean_latency{3,5}(1:5))]);
min_latency(2) = floor(min_latency(2)./10).*10;
 x2limits = [1 5];
 x2ticks = [1:5];
 y2limits = [min_latency(2) max_latency(2)];
 y2ticks = y2limits;
 
figure
set( gca, 'xlim', x2limits, 'xtick', x2ticks,'ylim', y2limits, 'ytick', y2ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],...
        'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('Pulse No.' ,'FontSize', 12);
        ylabel('Mean Latency', 'FontSize', 12);
hold on
plot([1:5],mean_latency{1,5}(1:5),...
    '-s','markersize',5,'MarkerFaceColor',[1 0 0],'color',[1 0 0],'linewidth',2)
plot([1:5],mean_latency{2,5}(1:5),...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
plot([1:5],mean_latency{3,5}(1:5),'-s','markersize',5,...
    'MarkerFaceColor',[0 0 1],'color',[0 0 1],'linewidth',2)
hold off

max_latency(3) = max([max(mean_latency{1,10}(1:9)) max(mean_latency{2,10}(1:9)) max(mean_latency{3,10}(1:9))]);
max_latency(3) = ceil(max_latency(3)./10).*10;
min_latency(3) = min([min(mean_latency{1,10}(1:9)) min(mean_latency{2,10}(1:9)) min(mean_latency{3,10}(1:9))]);
min_latency(3) = floor(min_latency(3)./10).*10;
 x3limits = [1 9];
 x3ticks = [1:9];
 y3limits = [min_latency(3) max_latency(3)];
 y3ticks = y3limits;
 
figure
set( gca, 'xlim', x3limits, 'xtick', x3ticks,'ylim', y3limits, 'ytick', y3ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],...
        'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('Pulse No.' ,'FontSize', 12);
        ylabel('Mean Latency', 'FontSize', 12);
hold on
plot([1:9],mean_latency{1,10}(1:9),...
    '-s','markersize',5,'MarkerFaceColor',[1 0 0],'color',[1 0 0],'linewidth',2)
plot([1:9],mean_latency{2,10}(1:9),...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
plot([1:9],mean_latency{3,10}(1:9),'-s','markersize',5,...
    'MarkerFaceColor',[0 0 1],'color',[0 0 1],'linewidth',2)
hold off

max_latency(4) = max(mean_latency{2,15}(1:13));
max_latency(4) = ceil(max_latency(4)./10).*10;
min_latency(4) = min(mean_latency{2,15}(1:13));
min_latency(4) = floor(min_latency(4)./10).*10;
 x4limits = [1 13];
 x4ticks = [1:13];
 y4limits = [min_latency(4) max_latency(4)];
 y4ticks = y4limits;
 
figure
plot([1:13],mean_latency{2,15}(1:13),...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
set( gca, 'xlim', x4limits, 'xtick', x4ticks,'ylim', y4limits, 'ytick', y4ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],...
        'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
    xlabel('Pulse No.' ,'FontSize', 12);
    ylabel('Mean Latency', 'FontSize', 12);        
%% Plotting Standard Deviations:

x1limits = [1 2];
 x1ticks = [1 2];
 y1limits = [];
 y1ticks = [];
 
figure
set( gca, 'xlim', x1limits, 'xtick', x1ticks,...
        'YMinorTick','on','ticklength', [0.010 0.010],...
        'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('Pulse No.' ,'FontSize', 12);
        ylabel('Variance', 'FontSize', 12);
hold on
plot([1:2],std_latency{1,2},...
    '-s','markersize',5,'MarkerFaceColor',[1 0 0],'color',[1 0 0],'linewidth',2)
plot([1:2],std_latency{2,2},...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
plot([1:2],std_latency{3,2}(1:2),'-s','markersize',5,...
    'MarkerFaceColor',[0 0 1],'color',[0 0 1],'linewidth',2)
hold off

figure
hold on
plot([1:5],std_latency{1,5}(1:5),...
    '-s','markersize',5,'MarkerFaceColor',[1 0 0],'color',[1 0 0],'linewidth',2)
plot([1:5],std_latency{2,5}(1:5),...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
plot([1:5],std_latency{3,5}(1:5),'-s','markersize',5,...
    'MarkerFaceColor',[0 0 1],'color',[0 0 1],'linewidth',2)
hold off

figure
hold on
plot([1:9],std_latency{1,10}(1:9),...
    '-s','markersize',5,'MarkerFaceColor',[1 0 0],'color',[1 0 0],'linewidth',2)
plot([1:9],std_latency{2,10}(1:9),...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
plot([1:9],std_latency{3,10}(1:9),'-s','markersize',5,...
    'MarkerFaceColor',[0 0 1],'color',[0 0 1],'linewidth',2)
hold off

figure
plot([1:13],std_latency{2,15}(1:13),...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)

%% Plotting variance:

max_var(1) = max([max(var_latency{1,2}) max(var_latency{2,2}) max(var_latency{3,2}(1:2))]);
max_var(1) = ceil(max_var(1)./10).*10;
min_var(1) = min([min(var_latency{1,2}) min(var_latency{2,2}) min(var_latency{3,2}(1:2))]);
min_var(1) = floor(min_var(1)./10).*10;
x1limits = [1 2];
 x1ticks = [1:2];
 y1limits = [min_var(1) max_var(1)];
 y1ticks = y1limits;
 
figure
set( gca, 'xlim', x1limits, 'xtick', x1ticks,'ylim', y1limits, 'ytick', y1ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],...
        'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('Pulse No.' ,'FontSize', 12);
        ylabel('Variance', 'FontSize', 12);
        
hold on
plot([1:2],var_latency{1,2},...
    '-s','markersize',5,'MarkerFaceColor',[1 0 0],'color',[1 0 0],'linewidth',2)
plot([1:2],var_latency{2,2},...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
plot([1:2],var_latency{3,2}(1:2),'-s','markersize',5,...
    'MarkerFaceColor',[0 0 1],'color',[0 0 1],'linewidth',2)
hold off

max_var(2) = max([max(var_latency{1,5}(1:5)) max(var_latency{2,5}(1:5)) max(var_latency{3,5}(1:5))]);
max_var(2) = ceil(max_var(2)./10).*10;
min_var(2) = min([min(var_latency{1,5}(1:5)) min(var_latency{2,5}(1:5)) min(var_latency{3,5}(1:5))]);
min_var(2) = floor(min_var(2)./10).*10;
 x2limits = [1 5];
 x2ticks = [1:5];
 y2limits = [min_var(2) max_var(2)];
 y2ticks = y2limits;
 
figure
set( gca, 'xlim', x2limits, 'xtick', x2ticks,'ylim', y2limits, 'ytick', y2ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],...
        'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('Pulse No.' ,'FontSize', 12);
        ylabel('Variance', 'FontSize', 12);
hold on
plot([1:5],var_latency{1,5}(1:5),...
    '-s','markersize',5,'MarkerFaceColor',[1 0 0],'color',[1 0 0],'linewidth',2)
plot([1:5],var_latency{2,5}(1:5),...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
plot([1:5],var_latency{3,5}(1:5),'-s','markersize',5,...
    'MarkerFaceColor',[0 0 1],'color',[0 0 1],'linewidth',2)
hold off

max_var(3) = max([max(var_latency{1,10}(1:9)) max(var_latency{2,10}(1:9)) max(var_latency{3,10}(1:9))]);
max_var(3) = ceil(max_var(3)./10).*10;
min_var(3) = min([min(var_latency{1,10}(1:9)) min(var_latency{2,10}(1:9)) min(var_latency{3,10}(1:9))]);
min_var(3) = floor(min_var(3)./10).*10;
 x3limits = [1 9];
 x3ticks = [1:9];
 y3limits = [min_var(3) max_var(3)];
 y3ticks = y3limits;
 
figure
set( gca, 'xlim', x3limits, 'xtick', x3ticks,'ylim', y3limits, 'ytick', y3ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],...
        'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('Pulse No.' ,'FontSize', 12);
        ylabel('Variance', 'FontSize', 12);
hold on
plot([1:9],var_latency{1,10}(1:9),...
    '-s','markersize',5,'MarkerFaceColor',[1 0 0],'color',[1 0 0],'linewidth',2)
plot([1:9],var_latency{2,10}(1:9),...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
plot([1:9],var_latency{3,10}(1:9),'-s','markersize',5,...
    'MarkerFaceColor',[0 0 1],'color',[0 0 1],'linewidth',2)
hold off

max_var(4) = max(var_latency{2,15}(1:13));
max_var(4) = ceil(max_var(4)./10).*10;
min_var(4) = min(var_latency{2,15}(1:13));
min_var(4) = floor(min_var(4)./10).*10;
 x4limits = [1 13];
 x4ticks = [1:13];
 y4limits = [min_var(4) max_var(4)];
 y4ticks = y4limits;
 
figure 
plot([1:13],var_latency{2,15}(1:13),...
    '-s','markersize',5,'MarkerFaceColor',[0 1 0],'color',[0 1 0],'linewidth',2)
set( gca, 'xlim', x4limits, 'xtick', x4ticks,'ylim', y4limits, 'ytick', y4ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],...
        'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
    xlabel('Pulse No.' ,'FontSize', 12);
    ylabel('Variance', 'FontSize', 12); 
%%
% % plotting the rasters:
% %multiplication of each trace in the number of the trace so that in the
% %plot, the spikes from different taces will be one above the other.
% %this is done by multiplication by a diagonal matrix.
% for rast = 1:clusters
%     D(rast,:) = [1:size(raster(:,:,rast),2)];
%     A(:,:,rast) = diag(D(rast,:));
%     raster_transpose(:,:,rast) = A(:,:,rast)*raster(:,:,rast)';
%     raster_plot(:,:,rast) =  raster_transpose(:,:,rast)';
%         
%     figure
%     set(gcf,'color','w')
%         
%         x1limits = [0 length(time_axis).*dt];
%         x1ticks = [];
% %         y1limits = ;
% %         y1ticks = ;
% 
%         set( gca, 'xlim', x1limits, ...
%         'ticklength', [0.010 0.010],'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
%         title(['file', num2str(Param.name)] ,'FontSize', 14);
%         xlabel('Time [sec]' ,'FontSize', 12);
%         ylabel('Trial no.', 'FontSize', 12);
%             hold on
%                 for i=1:size(raster(:,:,rast),2);
%                     scatter(time_axis(raster(:,i,rast)==1), i.*ones(sum(raster(:,i,rast)),1),...
%                         3, color_table(rast,:),'filled')
%                 end
%                     plot(time_axis_laser(laser_vec~=0), laser_vec(laser_vec~=0), '.c')
%             hold off    
%                
% end

% %% Making PSTHs:
% bin_time = 10; %[ms]
% bin_size = ceil(bin_time./1000.*sf);
% 
% for sorted = 1:clusters
% PSTH(:,:,sorted) = fn_PSTH(bin_time, bin_size, raster(:,:,sorted));
% 
% figure
% title(['file', num2str(Param.name), ' header', num2str(data.header(header).x_value(x_value).original_header), ' PSTH Spikes ', ...
%           'bin size [ms]=', num2str(bin_time)] ,'FontSize', 14);
% xlabel('Time [sec]' ,'FontSize', 12);
% ylabel('#Spikes/sec', 'FontSize', 12);
% hold on
% bar(time_axis(bin_size:bin_size:end),PSTH(:,:,sorted), 'k')
% plot(time_axis_laser(laser_vec~=0), laser_vec(laser_vec~=0), '.c')
% hold off
% 
% end
% 
