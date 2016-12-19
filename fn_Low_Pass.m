%% Low-Pass Filtering using butterworth
% fn_low_Pass is low-pass filtering using butterworth filter.
% inputs:
% data is a vector or a matrix in which columns are traces.
% sf is the scanning frequency of data
% low_pass_freq is in Hz.

function [data_HP] = fn_Low_Pass (data, sf, low_pass_freq)

    data_HP = zeros(size(data));

    WnHigh=low_pass_freq/(0.5*sf);
    [filterCoeff1 filterCoeff2]=butter(2,WnHigh,'low');
    
    for trace=1:size(data,2);
          data_HP(:,trace)=filtfilt(filterCoeff1,filterCoeff2,data(:,trace));
    end
    
end