% Returns a list of events detected in the trace (starting point,
% amplitude, amplitude position, half width of each event).
% Inputs: voltages - trace of recorded data, typically voltages (in mV),
% but actually it can be anything (e.g., currents, ECoG, etc.).
%                   dt - length (in seconds)  of single step in voltages
%                   finalAmp_Thres (optional) - minimum required amplitude of an event.
%                   negative values are interpreted in terms of standard  deviation
%                   to be considered an event (in mV.)
%                   doPlot (optional) - if true a plot with detected events
%                   marked is presented.



function [voltages, starting, amplitude, ampPos, halfWidth, halfWidthS, halfWidthE] = EventDetector_v2(voltages, dt, finalAmp_Thres, doPlot,I_temp)
starting = [];
amplitude = [];
ampPos = [];
halfWidth = []; halfWidthS = []; halfWidthE = [];

dt = dt*1e3; % now it's in msec.
% transpose the vector if nessecery
if size(voltages,1) > size(voltages,2)
    voltages=voltages';
end

% By medfilt1 uses a matrix of size length(voltages)*81, which
% can be way more than the available memory.
% The BLKSZ parameter is for 256MB buffer computer.
%voltages = medfilt1(voltages, 81, floor(256e6/(64*81))); % smoothed voltages

smooth_voltages = sgolayfilt(voltages, 3, 81);
Dvoltages = diff(smooth_voltages) / dt; % In mV per msec.
Dvoltages = [Dvoltages Dvoltages(end)]; % so it'll be the same length with voltages
DvoltagesSTD = std(Dvoltages);
D_Thres = DvoltagesSTD;
interval_length = 10; % in msec.


%D_Thres = 1;% mV per msec

%Amp_Thres = 1; %mV
Amp_Thres = std(voltages)/5; % <--- this way it can also work for detection of events in currents, EEG, etc.

if nargin < 3
    finalAmp_Thres = std(voltages); %mV
elseif finalAmp_Thres < 0
    finalAmp_Thres = abs(finalAmp_Thres) * std(voltages); %mV
end;
if nargin < 4
    doPlot = false;
end;

if Amp_Thres ~= 0 || D_Thres ~= 0
    convVector = zeros(1,ceil(Amp_Thres / D_Thres)/dt); 
    convVector(1) = +1;
    convVector(end) = -1; %convolving will make a derivative but with a large dt
    
    broadDvoltages = conv(convVector,smooth_voltages) / ceil(abs(Amp_Thres / D_Thres)); % In mV per msec.
    broadDvoltages = broadDvoltages(1:length(voltages));
    smoothDvoltages = medfilt1(Dvoltages, 9);
    eventsApprox1 = find(broadDvoltages > D_Thres);
    
    
    while ~isempty(eventsApprox1)
        %length(eventsApprox1)
        
        % here we try to find the real onset, by looking on the "real"
        % derivative.
        tmp = find(smoothDvoltages(eventsApprox1(1):-1:1) <= 0);
        %tmp = [];
        
        if isempty(tmp)
            %eventsApprox1(1)
            starting_point  = eventsApprox1(1);
        else
            starting_point = eventsApprox1(1) - tmp(1) + 1;
        end;
  %finding the peak position (the point after the onset where the derivative becomes negative):   
        tmp = find(smoothDvoltages(eventsApprox1(1):end-2) <= 0 & ...
            smoothDvoltages(eventsApprox1(1)+1:end-1) <= 0  & ...
            smoothDvoltages(eventsApprox1(1)+2:end) <= 0);
        if isempty(tmp)
            break;
        end;
        [maxval tmp] = max(voltages(eventsApprox1(1):eventsApprox1(1)+tmp(1)-1));
        
        tmp = tmp + eventsApprox1(1) - 1;
        amplitude_current =  voltages(tmp) - voltages(starting_point);
        amplitude_pos = tmp;
        
        
        
        while true
            %Half width of the events:
            t1 = find(voltages(starting_point:amplitude_pos-1) < voltages(starting_point)+amplitude_current/2 & ...
                voltages(starting_point+1:amplitude_pos) >= voltages(starting_point)+amplitude_current/2);
            t2 = find(voltages(amplitude_pos:end-1) >= voltages(starting_point)+amplitude_current/2 & ...
                voltages(amplitude_pos+1:end) < voltages(starting_point)+amplitude_current/2);
            
            if ~isempty(t1) & ~isempty(t2)
                t1 = t1(1) + starting_point;
                t2 = t2(1) + amplitude_pos-1;
                
                % Update the max. point:
                [maxval tmp] = max(voltages(t1:t2));
                tmp = tmp + t1 - 1;
                amplitude_current =  voltages(tmp) - voltages(starting_point);
                if abs(tmp - amplitude_pos) < 0.5 / dt % nothing significant changed
                    break;
                else
                    amplitude_pos = tmp;
                end;
            else
                break;
            end;
            
        end;
        
        % Here we remove all the event points belonging to the same event (i.e.
        % up to the max pos).
        if eventsApprox1(1) >= amplitude_pos
            eventsApprox1 = eventsApprox1(2:end);
            continue; % This means the whole event is before the first approx. ==> we ignore it.
        elseif ~isempty(ampPos) & starting_point <= ampPos(end) % This means this event begins before the peak of the previous event
            eventsApprox1 = eventsApprox1(2:end);
            continue;
        else
            eventsApprox1 = eventsApprox1(find(eventsApprox1 > amplitude_pos));
        end;
        %     tmp = find(broadDvoltages(eventsApprox1(1):end) <= 0);
        %     if isempty(tmp)
        %         break;
        %     end;
        %     % Else we remove all the other events between the last starting point
        %     % and where the derivative becomes <= 0
        %     tmp = tmp + eventsApprox1(1) - 1;
        %     eventsApprox1 = eventsApprox1(find(eventsApprox1 > tmp(1) ));
        
        
        
        % Some sanity checks
        if amplitude_current > finalAmp_Thres & starting_point < amplitude_pos & ~isempty(t1) & ~isempty(t2)
            starting = [starting starting_point];
            amplitude = [amplitude amplitude_current];
            ampPos = [ampPos amplitude_pos];
            halfWidth = [halfWidth (t2 - t1)*dt*1e-3]; %dt is in msec
            halfWidthS = [halfWidthS t1];
            halfWidthE = [halfWidthE t2];
        end;
    end;
    if doPlot
        if I_temp >0
            
            voltages = -1*voltages;
            I_temp;
        end
        figure(4); clf;
        plot(voltages, 'b');
        hold on;
        plot(starting,voltages(starting), 'r*');
        plot(ampPos,voltages(ampPos), 'k*');
        
        %text(ampPos,voltages(ampPos), num2str(halfWidth'));
        for i = 1:length(halfWidthS)
            plot([halfWidthS(i) halfWidthE(i)], [voltages(halfWidthS(i)), voltages(halfWidthE(i))], 'r');
        end;
    end;
end




end
% while current_pos + interval_length < length(voltages)
%     tmp_min = min(voltages(current_pos:current_pos + interval_length));
%     tmp_max = max(voltages(current_pos:current_pos + interval_length));
%
%     % We check if they are indeed local minima / maxima
%     if tmp_min <= min(voltages(current_pos:current_pos + round(interval_length*0.3))) & ...
%             tmp_min <= min(voltages(current_pos + round(interval_length*0.7):current_pos + interval_length)) & ...
%             tmp_min == min(voltages(current_pos + round(interval_length*0.3):current_pos + round(interval_length*0.7)))
%         local_extrema_x(counter) = current_pos + min(find(voltages(current_pos:current_pos + interval_length) <= tmp_min));
%         local_extrema_y(counter) = voltages(local_extrema_x(counter));
%         counter = counter + 1;
%     end;
%
%     if tmp_max >= max(voltages(current_pos:current_pos + round(interval_length*0.3))) & ...
%             tmp_max >= max(voltages(current_pos + round(interval_length*0.7):current_pos + interval_length)) & ...
%             tmp_max == max(voltages(current_pos + round(interval_length*0.3):current_pos + round(interval_length*0.7)))
%         local_extrema_x(counter) = current_pos + min(find(voltages(current_pos:current_pos + interval_length) >= tmp_max));
%         local_extrema_y(counter) = voltages(local_extrema_x(counter));
%         counter = counter + 1;
%     end;
%
%     % disp(['---']);
%     current_pos = current_pos + round(0.3 * interval_length);
% end;
%
%
% % We order the result into an ascending order of the x values.
% % And remove duplicates
% tmp = [local_extrema_x; local_extrema_y];
% tmp = unique(sortrows(tmp'), 'rows');
% local_extrema_x = tmp(:,1);
% local_extrema_y = tmp(:,2);

