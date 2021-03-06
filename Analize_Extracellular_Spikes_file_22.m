%% Analyze extracellular spikes

clear all
close all
global Exp
cd 'D:\Inbal M.Sc\Data\Extracted_Data';
load Set2-2013-03-20-002 
cd 'D:\Inbal M.Sc\MATLAB\Project Brain States'

header = 2;
channel = 1; % V1 - 1, Ext1 - 2

switch channel
    case 1
        raw_data = data.header(header).Vm;
        sf = Param.header(1,header).stim.sf;
    
    case 2
        raw_data = -1.*data.header(header).Ext1; %for extracellular recording, spikes are upside-down
        sf = Param.header(1,header).stim.sf_Ext1;
end

dt=1/sf; %[sec]
time_axis=(1:size(raw_data,1))*dt;

high_pass_freq = 300; %[Hz]
low_pass_freq = 3000; %[Hz]
BaselineTime = [];
robustness = 3;
peaks = [];
interpeak_time = 3; %[ms]
clusters = 1; %according to the number of neurons to be sorted

% parameters for preparing vectors at the same length of intervals around
% the peaks:
begin_peak_time = 0.5; %[ms]
end_peak_time = 1; %[ms]
interval_begin = ceil(begin_peak_time./1000.*sf);
interval_end = ceil(end_peak_time./1000.*sf);
interval_around_peak = interval_begin + interval_end +1;
peak_counter = 0;

raster_1 = zeros(size(raw_data));
raster_2 = zeros(size(raw_data));
raster_3 = zeros(size(raw_data));

[data_HP] = fn_High_Pass (raw_data, sf, high_pass_freq);
[data_HP] = fn_Low_Pass (data_HP, sf, low_pass_freq);
   


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
%% Spike sorting using PCA, fitting gaussians and clustering.
% clusters = 1;

[U,S,V] = svds(peaks_mat,2);

obj  = gmdistribution.fit(U,clusters, 'Replicates',500); %fitting the scattered dots into 2 gaussians.
idx = cluster(obj,U); %returns a vector same length as U, with the group number for each entry according to the clustering

figure
subplot(2,1,1)
hold on
scatter(U(idx==1,1),U(idx==1,2),3,'b','filled');
scatter(U(idx==2,1),U(idx==2,2),3,'r','filled');
% scatter(U(idx==3,1),U(idx==3,2),3,'g','filled');
% scatter(U(idx==4,1),U(idx==4,2),5,'y','filled');
hold off


subplot(2,1,2)
hold on
plot([1:size(peaks_mat,2)],peaks_mat(idx==1,:),'b')
plot([1:size(peaks_mat,2)],peaks_mat(idx==2,:),'r')
% plot([1:size(peaks_mat,2)],peaks_mat(idx==3,:),'g')
% plot([1:size(peaks_mat,2)],peaks_mat(idx==4,:),'y')
hold off

%% making rasters for the sorted spikes

for peak = 1:length(peaks_loc)
    if idx(peak) == 1
        raster_1(peaks_loc(peak), peaks_trace(peak)) = 1;
    else if idx(peak) == 2;
            raster_2(peaks_loc(peak), peaks_trace(peak)) = 1;
         else if idx(peak) == 3;
                raster_3(peaks_loc(peak), peaks_trace(peak)) = 1;
             end
        end
    end
end

% plotting the rasters:
D1 = [1:size(raster_1',1)];
A1 = diag(D1);
raster_1_plot = A1*raster_1';
raster_1_plot = raster_1_plot';

D2 = [1:size(raster_2',1)];
A2 = diag(D2);
raster_2_plot = A2*raster_2';
raster_2_plot = raster_2_plot';

% D3 = [1:size(raster_3',1)];
% A3 = diag(D3);
% raster_3_plot = A3*raster_3';
% raster_3_plot = raster_3_plot';
% 
            figure
            title(['file', num2str(Param.name), ' header', num2str(data.header(header).original_header)] ,'FontSize', 14);
            xlabel('Time [sec]' ,'FontSize', 12);
            ylabel('Trial no.', 'FontSize', 12);
            hold on
                for i=1:size(raster_1,2);
                    plot(time_axis(raster_1(:,i)==1), raster_1_plot(raster_1(:,i)==1,i), 'b.')
                end
            hold off

            figure
            title(['file', num2str(Param.name), ' header', num2str(data.header(header).original_header)] ,'FontSize', 14);
            xlabel('Time [sec]' ,'FontSize', 12);
            ylabel('Trial no.', 'FontSize', 12);
            hold on
                for i=1:size(raster_2,2);
                    plot(time_axis(raster_2(:,i)==1), raster_2_plot(raster_2(:,i)==1,i), 'r.')
                end
            hold off
            
%             figure
%             title(['file', num2str(Param.name), ' header', num2str(data.header(header).original_header)] ,'FontSize', 14);
%             xlabel('Time [sec]' ,'FontSize', 12);
%             ylabel('Trial no.', 'FontSize', 12);
%             hold on
%                 for i=1:size(raster_3,2);
%                     plot(time_axis(raster_3(:,i)==1), raster_3_plot(raster_3(:,i)==1,i), 'g.')
%                 end
%             hold off

%% Making PSTHs:
bin_time = 10; %[ms]
bin_interval = ceil(bin_time./1000.*sf);

PSTH_1 = fn_PSTH(bin_interval, raster_1);
PSTH_2 = fn_PSTH(bin_interval, raster_2);
% PSTH_3 = fn_PSTH(bin_interval, raster_3);

figure
plot(PSTH_1, 'b')

figure
plot(PSTH_2, 'r')

% figure
% plot(PSTH_3, 'g')
