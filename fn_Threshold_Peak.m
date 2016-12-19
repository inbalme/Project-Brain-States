% fn_Threshold_Peak is calculating peak threshold for each trace, based on
% a calculation of mean baseline-std. the mean baseline is taken from the
% beginning of the trace. For this reason, it is only suitable for
% high-pass filtered data. 

% need 4 inputs:
% 1. data is a vector or a matrix where columns are traces
% 2. sf is the scanning frequency
% 3. BaselineTime is the time (in ms) at the beggining of each trace that will
%    be averaged to yield the meanBaseline. default is 10ms.
%4. robustness is the number of std that will be used to calculate the threshold.
%   default is 3 std.

function [threshold_peak] = fn_Threshold_Peak (data, sf, BaselineTime, robustness)

% if nargin == 3
%     robustness = 3;
% else if nargin == 2
%         robustness = 3;
%         BaselineTime = 10; %[ms]
%     end
% end

if isempty(BaselineTime) || BaselineTime==0
    BaselineInterval = length(data);
else if isempty(robustness) || robustness==0
        robustness = 3;
    else
        BaselineTime = BaselineTime/1000; %[sec]
        BaselineInterval = floor(BaselineTime.*sf);
    end
end
%         meanBaseline = mean(data(1:BaselineInterval,:));
        medianBaseline = median(data(1:BaselineInterval,:));
        stdBaseline = std(data(1:BaselineInterval,:));
        threshold_peak(1,:) = medianBaseline+robustness.*stdBaseline;
%         threshold_peak(1,:) = robustness.*stdBaseline;
end