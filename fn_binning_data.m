%% averages over bins of Vm and compares them to the average baseline
%accepts a matrix in which each column is a trace

function [binned_data_mean binned_data_std] = fn_binning_data(raw_data,dt_data, bin_win);

% bin_win = 0.01; %[sec]
bin_size = bin_win./dt_data; %#data points

 bin_num = size(raw_data,1)./bin_size;
 row_num = bin_size;
 col_num = size(raw_data, 2).*bin_num;
 reshaped_data = reshape(raw_data, row_num, col_num);
 reshaped_data_mean = mean(reshaped_data,1);
 reshaped_data_std = std(reshaped_data,0,1);
 binned_data_mean = reshape(reshaped_data_mean, bin_num, size(raw_data, 2));
 binned_data_std = reshape(reshaped_data_std, bin_num, size(raw_data, 2));

 
end