%% PSTH
% This function only works if the number of data points in each trace (the
% number of rows) is devided by the binsize without leftovers...

function [PSTH] = fn_PSTH(bin_time, bin_size, raster)

nbin = size(raster, 1)/bin_size;
 nrow = bin_size;
 ncol = size(raster, 2).*nbin;
 rasterReshaped = reshape(raster, nrow, ncol);
 rasterReshapedSum = sum(rasterReshaped);
 binSum = reshape(rasterReshapedSum, nbin, size(raster, 2));
 PSTH = mean(binSum, 2)./bin_time.*1000; %[spikes/sec]
 
end