%% Detect Spikes for High-Pass-Filtered data: using the function findpeaks.
% need 4 inputs:
% 1. data is a high-pass-filtered vector (one trace).
% 2. sf is the scanning frequency
% 3. interpeak_time [ms] is the minimum time interval between two adjacent
%    peaks.
% 4. threshold_peak is the threshold for peak detection. can be provided by
%    fn_Threshold_Peak.
% fn_Detect_Spike_HPF returns two column vectors: peaksValue and peakLoc -
% amplitude and location of peaks.
 
function [peaksValue,peaksLoc] = fn_Detect_Spike_HPF(data, sf, interpeak_time, threshold_peak)

   interpeak_distance = ceil(interpeak_time./1000.*sf);
%    data = data.*(-1); % for extracellular recording, spikes are up-side down
   
   [peaksValue, peaksLoc] = findpeaks(data, ...
      'minpeakheight', threshold_peak, 'minpeakdistance',interpeak_distance);

        
end
