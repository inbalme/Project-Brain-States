% Removes spikes, and puts spline interpolation instead.
% If unsuccessful (there is not enough data, or there are too many spikes)
% an error (code < 0) will be reported. Otherise code is the number of
% spikes removed.
% spike_positions is a vector with the starting points of spikes.
% Parameters:
%   voltages - A row vector, describing the voltange in the cell over time.
%   dt - the time step (in seconds).
%edited on 20/12/2016 by Inbal Meir. added the variable 'spike_positions'

function [code, result, spike_positions] = fn_RemoveSpikesMO2(voltages, dt)
    % RemoveSpikes2 does one direction, so we use it twice.
    [code, r, spike_positions] = RemoveSpikes2(voltages, dt);
    r = -1 * r;
    [d, r, s] = RemoveSpikes2(r, dt);
    result = -1 * r;
    
%     figure
%     plot(voltages, 'b');
%     hold on;
%     plot(result, 'r');

    
function [code, result, spike_positions] = RemoveSpikes2(voltages, dt,spike_derivative_threshold)
    if nargin<3;
        spike_derivative_threshold = 3e3;
    end
    ffff = spike_derivative_threshold;
    if (length(voltages) < 1e3)
        code = -1;
        result = [];
        spike_positions = -1;
        return;
    end;
    Dvoltages = zeros(1, length(voltages));
    % Dvoltages(2:length(voltages)) = diff(voltages)/dt;
    Dvoltages(2:length(voltages)) = diff(sgolayfilt(voltages, 3, 7))/dt;
    Dvoltages(1) = Dvoltages(2);
  
    %spike_derivative_threshold = 30e3;
   
        
    min_calm_interval = 5e-3/dt; % 5 mseconds
    
    current_pos = 10;
    result = voltages;
    code = 0;
    spike_positions = [];
    
    spikes = find(Dvoltages(current_pos:end) > spike_derivative_threshold);
    if length(spikes) == 0 % no more spikes
        return;
    else
        current_pos = current_pos + min(spikes) - 1;
        code = 1;
    end;
    
    while 1
        % We want to find where the spike / train of spikes ends.
        % For this we need to find an interval of length min_calm_interval with
        % small derivatives.
        step = 5;
        if current_pos+step+50 > length(voltages)
            return;
        end;
        spikes = find(Dvoltages(current_pos+step:current_pos+step+min_calm_interval) > spike_derivative_threshold);
        while length(spikes) > 0 
            step = step + 5;  
            if current_pos+step+min_calm_interval > length(voltages)
                return;
            end;
            spikes = find(Dvoltages(current_pos+step:current_pos+step+min_calm_interval) > spike_derivative_threshold);
        end;
        % our spike / train of spikes is in the 
        % current_pos...current_pos+step
        % to find the precise end-point we look for the time where
        % the voltage goes below the spike starting point.
        
        tmp = find(voltages(current_pos+step:current_pos+step+min_calm_interval) < voltages(max(1,current_pos-1)));
        if length(tmp) > 0
            step = step + min(tmp);
        else
            step = step + 50;
        end;
       
        x = [];
        y = [];
        points_num = 5;
        
        if current_pos + step + points_num - 1 > length(voltages)
            step = step - (current_pos + step + points_num - 1 - length(voltages));
        end;
               
        x(1:points_num) = [-points_num+1:0];
        x(points_num+1:2*points_num) = [1+step:step + points_num];
        y(1:points_num) = voltages(current_pos - points_num : current_pos - 1);
        y(points_num+1:2*points_num) = voltages(current_pos + step : current_pos + step + points_num - 1);
    
        % That's the basic thing. 
        % Now, if it is not a single spike but several, we will find
        % the minimums between spikes and pass the spline through them.
                
        tmp_counter = 1;
        for i = 0 : step-1
            if Dvoltages(current_pos + i) <= 0 & Dvoltages(current_pos + i + 1) > 0
                x(2*points_num+tmp_counter) = i + 1;
                y(2*points_num+tmp_counter) = voltages(current_pos + i);
                tmp_counter = tmp_counter + 1;
                code = code + 1; % counts the number of spikes
                spike_positions=[spike_positions, current_pos];
            end;
        end;
        
        pp = csapi(x,y); %spline fit
        
        for i = 0 : step
            result(i + current_pos) = fnval(pp, i + 1);
        end;
        
        % disp(['--> from ', int2str(current_pos), ' to ', int2str(current_pos+step)]);
                
        current_pos = current_pos + step; % We move on...
        spikes = find(abs(Dvoltages(current_pos:end)) > spike_derivative_threshold);
        if length(spikes) == 0 % no more spikes
            return;
        else
            current_pos = current_pos + min(spikes) - 1;
        end;
                
    end;
    
    