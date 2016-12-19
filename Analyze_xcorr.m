% 
clear all
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load f16_workspace

%% cross-correlation Vm-EEG - spontaneous activity
clear c lags c_mean lags_new start_time duration start_time start_sample end_sample interval interval_mat
x_value=1; coeffs=[];
data_downsamp = data_no_spikes(1:5:end,:,:);
 start_time = [0,4]; %[sec] %[0,5]
 duration = 3; %[sec] 
for t=1:length(start_time);
 start_sample(:,t) = start_time(t).*sf_EEG;
if start_time(t)==0
    start_sample(:,t) = 1;
end
end_sample(:,t) = start_sample(:,t)+duration.*sf_EEG-1;
interval(:,t) = round(start_sample(:,t):end_sample(:,t));
% interval_mat1(:,:,t) = data_downsamp(interval(:,t),:,x_value);
% interval_mat2(:,:,t) = data_EEG(interval(:,t),:,x_value);

     for trace = 1:size(data_no_spikes,2)
         x=data_downsamp(interval(:,t),trace,x_value)-mean(data_downsamp(interval(:,t),trace,x_value));
         y=data_EEG(interval(:,t),trace,x_value)-mean(data_EEG(interval(:,t),trace,x_value));
        [c{t}(:,trace),lags{t}(:,1)] = xcorr(x,y,'coeff') ;
        
        [r,p] = corrcoef(x,y);
        coeffs{t}(1,trace)= r(1,2);
        pvals{t}(1,trace)= p(1,2);
     end      
c_mean(:,t) = mean(c{t},2);
coeffs_mean(1,t)=mean(coeffs{t});
pvals_mean(1,t)=mean(pvals{t});

end

%% plot crosscorrelation
figure
hold on
plot( lags{1,1},c_mean(:,1), 'b-')
plot( lags{1,1},c_mean(:,2), 'r-')
hold off
%%
trace=5;
% plotting one trace of data(downsampled) and EEG against each other -
% before ES
figure
subplot(3,1,1)
hold on
plot(interval(:,1),data_downsamp(interval(:,1),trace,x_value), 'k-')
plot(interval(:,1),data_EEG(interval(:,1),trace,x_value).*10-50,'g-')
hold off

% after ES
subplot(3,1,2)
hold on
plot(interval(:,2),data_downsamp(interval(:,2),trace,x_value), 'k-')
plot(interval(:,2),data_EEG(interval(:,2),trace,x_value).*10-50,'g-')
hold off

%plotting the crosscorrelation for a single trace
subplot(3,1,3)
hold on
plot( lags{1,1},c{1}(:,trace), 'b-')
plot( lags{1,1},c{2}(:,trace), 'r-')
hold off