function res = fn_CalculateCapacitResist(I, voltages, frequency, doPlot, v_baseline)

% CalculateCapacitResist - claculate the Capacitance and Resistance of the
% cell and electrode.
% Inputs:
%           1)I - the amount of current injected during the current step (in nA).
%           2)voltages - a trace of the voltages starting from the point the
%                       current step started (in mV).
%           3)frequency - sampling frequency of data (in KHz, e.g., 10,15,20).
%           4)doPlot    - boolean variable that controls if a plot of  the fitting is
%                       presented.
%           5)v_baseline - if provided, the level to which the exponent descends
%
% Output:   result(1) = R1, 
%           result(2) = C1 - resistance and capacitance of the electrode
%           result(3) = R2, 
%           result(4) = C2 - resistance and capacitance of the cell
%           result(5) - the max distance between actual data and the fit (in the
%                       interval which is used for fitting).
%           result(6) - the error in L2 metric (sum of square roots of errors at each
%                       point).

%% var checks
% What follows below assumes that the voltages curve is decreasing, so
% if the case the opposite, we simply reverse it.
if I > 0 & voltages(3) > voltages(1) & voltages(6) > voltages(3)
    voltages = -1*voltages;
    I = -1*I;
    if exist('v_baseline')
        v_baseline = -1*v_baseline;
    end;
end;

DO_PLOTS = false;
if nargin >= 4
    DO_PLOTS_FINAL_ONLY = doPlot;
else
    DO_PLOTS_FINAL_ONLY = false;
end;
frequency = frequency * 1e3;
I = I*1e-9;
voltages = voltages*1e-3; % the voltage decay response at end of injection
if exist('v_baseline')
    v_baseline = v_baseline*1e-3;
end;

if size(voltages,1) ~= 1 % not a vector
    error('wrong voltages format');
elseif length(voltages) < frequency*0.1
    % We don't have 100 ms
    error('the trace is not long enough');
end;

% And we make sure the beginning is at 0.
tmp = voltages(1);
voltages = voltages - tmp;
if exist('v_baseline')
    v_baseline = v_baseline - tmp;
end;

if DO_PLOTS | DO_PLOTS_FINAL_ONLY % plot the voltage decay
    figure;
    hold on;
    plot(voltages(1:(length(voltages/2))))
end;

%%
% We have to determine the length of the fitting interval.
% It cannot be set to some constant value in advance (e.g. 10 msecs)
% because the time constant varies quite a lot
% (QX slows it down by almost an order of magnitude).

% Dvoltages = diff(sgolayfilt(voltages(1:frequency*0.1),3,51));
Dvoltages = medfilt1(diff(medfilt1(voltages(1:frequency*0.1),30)),30);
Dvoltages = Dvoltages(20:length(Dvoltages)); % because of the medfilt1ering the first values are 0s.
if Dvoltages(1) < 0
    FITTING_DURATION = min(find(Dvoltages > 0))/frequency;
else
    FITTING_DURATION = min(find(Dvoltages < 0))/frequency;
end;
if isempty(FITTING_DURATION)
    FITTING_DURATION = min(find(Dvoltages == 0))/frequency;
end;
if exist('v_baseline') & ~isempty(find(voltages < v_baseline))
    FITTING_DURATION = min(find(voltages < v_baseline))/frequency;
end;
if isempty(FITTING_DURATION) | (FITTING_DURATION < 0.15 & ceil(frequency*(FITTING_DURATION+0.01)) > length(voltages))
    error('Problem with finding FITTING_DURATION');
    return;
    % FITTING_DURATION = length(voltages)/frequency - 0.015;
end;
if FITTING_DURATION > 0.2
    FITTING_DURATION = 0.2;
end;
% FITTING_DURATION = max(FITTING_DURATION,0.010); % 10 mseconds

%%

% we will have several fittings from which the best will be chosen.
num_of_results = 0;

% According to the voltage after FITTING_DURATION mseconds
%check that voltage is steady ( back to base line)
if max(voltages(ceil(frequency*FITTING_DURATION):ceil(frequency*(FITTING_DURATION+0.01)))) ...
        - min(voltages(ceil(frequency*FITTING_DURATION):ceil(frequency*(FITTING_DURATION+0.01)))) < 0.002
    final_voltage = mean(voltages(ceil(frequency*FITTING_DURATION):ceil(frequency*(FITTING_DURATION+0.01))));
else
    % the variablity is still too high - the cell does not behave
    % nicely
    final_voltage = min(voltages(ceil(frequency*FITTING_DURATION):ceil(frequency*(FITTING_DURATION+0.01))));
end;

R_total = final_voltage/I;
R_total = abs(R_total)


%% First we do PEELING
time = (0:1:length(voltages)-1);
pealing_trial = 1;

while 1 % Calculate the initial guess values with polyfit
    if pealing_trial > 10 % we won't do more than 10 retries
        error('Problem with initial guess');
        return;
    end;
    % The slower exponent (of the cell).
    p2 = polyfit((ceil(frequency*15/1e4):1:ceil(50*frequency/1e4))/frequency,log(voltages(ceil(frequency*15/1e4):ceil(50*frequency/1e4))-final_voltage),1);
    x2 = p2(1);
    alpha2 = exp(p2(2));

    % The faster exponent (of the electrode).
    %if isempty(find((voltages(1:5)-final_voltage-alpha2*exp(x2*(1:1:5)/frequency)) <= 0))
    p1 = polyfit((1:1:ceil(frequency*5/1e4))/frequency,log(voltages(1:ceil(frequency*5/1e4))-final_voltage-alpha2*exp(x2*(1:1:ceil(frequency*5/1e4))/frequency)),1);
    x1 = p1(1);
    alpha1 = exp(p1(2));
    %end;



    % At the beginning voltage derivative is a good approximation for I/C1,
    % thus:
    v_derivative_1 = (voltages(2)-voltages(1))*frequency;
    v_derivative_2 = (voltages(3)-voltages(2))*frequency;
    v_derivative_3 = (voltages(4)-voltages(3))*frequency;
    C1 = I * 1/((v_derivative_1 + v_derivative_2 + v_derivative_3)/3);

    %The other initial guesses:
    R1 = 1/(-1*x1*C1);
    R2 = R_total - R1;
    C2 = 1/(-1*x2*R2);

    if C1 <= 0 | C2 <= 0 | R1 <= 0 | R2 <= 0
        pealing_trial = pealing_trial + 1;
        final_voltage = final_voltage - 2e-4; % 0.2 mV
    else
        break;
    end;
end; % while
%% initial values

tau1=R1*C1;
tau2=R2*C2;
rho= R2/R1;

num_of_results = num_of_results + 1;
result(num_of_results,1) = R1;
result(num_of_results,2) = C1;
result(num_of_results,3) = R2;
result(num_of_results,4) = C2;
result(num_of_results,5) = checkFitting(R1,C1,R2,C2,false);
%%
% Nested function - used for fitting all parameters
    function diff = requiredVoltage_all(beta)
        R1 = abs(beta(1));
        C1 = abs(beta(2));
        R2 = abs(beta(3));
        C2 = abs(beta(4));

        if (R1 + R2 > R_total)
            R2 = max(1e7,R_total - R1);
        end;


        rho= R2/R1;
        tau1 = C1*R1;
        tau2 = C2*R2;

        x1 = -((1+rho)*tau1+tau2)-sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x1 = x1/(2*tau1*tau2);
        x2 = -((1+rho)*tau1+tau2)+sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x2 = x2/(2*tau1*tau2);
        alpha2 = I*(1+(R1+R2)*C1*x1)/(C1*x2 - C1*x1);
        alpha1 = -I*(R1+R2)-alpha2;

        time = (0:1:ceil(FITTING_DURATION*frequency));
        v = alpha1*exp(x1*time/frequency) + alpha2*exp(x2*time/frequency)+I*(R1+R2);


        diff = norm(v - voltages(1:1+ceil(FITTING_DURATION*frequency)),FITTING_NORM);
    end % function requiredVoltage_all


%%
    % Nested function - used for fitting the taus (the two time
    % constants).
    function diff = requiredVoltage_taus(beta)
        tau1 = abs(beta(1));
        tau2 = abs(beta(2));

        R1 = tau1/C1;
        R2 = abs(R_total - R1); % normally it will be positive anyway
        rho= R2/R1;
        C2 = tau2/R2;

        x1 = -((1+rho)*tau1+tau2)-sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x1 = x1/(2*tau1*tau2);
        x2 = -((1+rho)*tau1+tau2)+sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x2 = x2/(2*tau1*tau2);
        alpha2 = I*(1+(R1+R2)*C1*x1)/(C1*x2 - C1*x1);
        alpha1 = -I*(R1+R2)-alpha2;

        time = (0:1:ceil(FITTING_DURATION*frequency));
        v = alpha1*exp(x1*time/frequency) + alpha2*exp(x2*time/frequency)+I*(R1+R2);


        diff = norm(v - voltages(1:1+ceil(FITTING_DURATION*frequency)),FITTING_NORM);
    end % function requiredVoltage_taus
%%
    % Nested function - used for fitting C1  (electrode capacitance).
    function diff = requiredVoltage_C1(beta)
        C1 = abs(beta);

        R1 = tau1/C1;
        R2 = abs(R_total - R1); % normally it will be positive anyway
        rho= R2/R1;
        C2 = tau2/R2;

        x1 = -((1+rho)*tau1+tau2)-sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x1 = x1/(2*tau1*tau2);
        x2 = -((1+rho)*tau1+tau2)+sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x2 = x2/(2*tau1*tau2);
        alpha2 = I*(1+(R1+R2)*C1*x1)/(C1*x2 - C1*x1);
        alpha1 = -I*(R1+R2)-alpha2;

        time = (0:1:ceil(FITTING_DURATION*frequency));

        v = alpha1*exp(x1*time/frequency) + alpha2*exp(x2*time/frequency)+I*(R1+R2);
        diff = norm(v - voltages(1:1+ceil(FITTING_DURATION*frequency)),FITTING_NORM);
    end % function requiredVoltage_C1
%%
    % Nested function - used for fitting the two resistances.
    function diff = requiredVoltage_R1R2(beta)
        R1 = abs(beta(1));
        R2 = abs(beta(2));

        tau1=R1*C1;
        rho= R2/R1;
        tau2 = R2*C2;

        x1 = -((1+rho)*tau1+tau2)-sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x1 = x1/(2*tau1*tau2);
        x2 = -((1+rho)*tau1+tau2)+sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x2 = x2/(2*tau1*tau2);
        alpha2 = I*(1+(R1+R2)*C1*x1)/(C1*x2 - C1*x1);
        alpha1 = -I*(R1+R2)-alpha2;

        time = (0:1:ceil(FITTING_DURATION*frequency));
        v = alpha1*exp(x1*time/frequency) + alpha2*exp(x2*time/frequency)+I*(R1+R2);

        diff = norm(v - voltages(1:1+ceil(FITTING_DURATION*frequency)),FITTING_NORM);
    end % function requiredVoltage_R1R2

%%
    % dist is  the distance of the fit from the data (in max metric)
    % r2dist is  the distance of the fit from the data (in L2 metric)
    % also plots if needToPlot == true.
    function [dist r2dist] = checkFitting(R1,C1,R2,C2,needToPlot)
        if C1 <= 0 | C2 <= 0 | R1 <= 0 | R2 <= 0
            error('Negative arguments to checkFitting()');
            return;
        end;
        tau1=R1*C1;
        rho= R2/R1;
        tau2 = R2*C2;

        x1 = -((1+rho)*tau1+tau2)-sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x1 = x1/(2*tau1*tau2);
        x2 = -((1+rho)*tau1+tau2)+sqrt(((1+rho)*tau1+tau2)^2 - 4*tau1*tau2);
        x2 = x2/(2*tau1*tau2);
        alpha2 = I*(1+(R1+R2)*C1*x1)/(C1*x2 - C1*x1);
        alpha1 = -I*(R1+R2)-alpha2;

        time = (0:1:ceil(10*FITTING_DURATION*frequency));
        v = alpha1*exp(x1*time/frequency)+ alpha2*exp(x2*time/frequency)+I*(R1+R2);



        dist = max(abs(v(1:ceil(FITTING_DURATION*frequency))-voltages(1:ceil(FITTING_DURATION*frequency))));
        if needToPlot
            index = find(abs(v(1:ceil(FITTING_DURATION*frequency))-voltages(1:ceil(FITTING_DURATION*frequency))) >= dist);
            plot(voltages(1:ceil(FITTING_DURATION*frequency)),'*');
            plot(v(1:floor(2*FITTING_DURATION*frequency)), 'r')
            plot((1:1:floor(2*FITTING_DURATION*frequency)),final_voltage,'k');
            plot(index, v(index), '+r');
        end;
        r2dist = norm(v(1:ceil(FITTING_DURATION*frequency)) -  voltages(1:ceil(FITTING_DURATION*frequency)));
    end;
%%
    % And now... to fitting!

    % The code is somewhat twisted, since the values (namely R1,C1,R2,C2) are changed
    % inside the functions called by fminsearch. Since it seems like fminsearch
    % converges to the best value (and thus the final calls are the best
    % ones) it works...


    % 1st fitting - few parameters each time.

    for FITTING_NORM = 1:2

        fminsearch(@requiredVoltage_taus, [tau1 tau2], optimset('MaxFunEvals',1000,'Display', 'off'));
        fminsearch(@requiredVoltage_C1, [C1], optimset('MaxFunEvals',1000, 'Display', 'off'));

        temp_res = +Inf;
        while temp_res - checkFitting(R1,C1,R2,C2,false) > 1e-13  % otherwise the program can get stuck for a very long time trying to
            % get very small improvements
            temp_res = checkFitting(R1,C1,R2,C2,false);
            fminsearch(@requiredVoltage_R1R2, [R1 R2],optimset('Display', 'off'));
            fminsearch(@requiredVoltage_C1, [C1],optimset('Display', 'off'));
            fminsearch(@requiredVoltage_taus, [tau1 tau2],optimset('Display', 'off'));
        end;

        temp_res = +Inf;
        while temp_res - checkFitting(R1,C1,R2,C2,false) > 1e-13
            temp_res = checkFitting(R1,C1,R2,C2,false);
            fminsearch(@requiredVoltage_R1R2, [R1 R2],optimset('Display', 'off'));
        end;

        num_of_results = num_of_results + 1;
        result(num_of_results,1) = abs(R1);
        result(num_of_results,2) = abs(C1);
        result(num_of_results,3) = abs(R2);
        result(num_of_results,4) = abs(C2);
        [result(num_of_results,5)  result(num_of_results,6)]= checkFitting(R1,C1,R2,C2,DO_PLOTS);

    end; % for FITTING_NORM

    % 2nd fitting - all together (with 2 norms)
    FITTING_NORM = 1;
    temp_res = +Inf;
    while temp_res - checkFitting(R1,C1,R2,C2,false) > 1e-13
        temp_res = checkFitting(R1,C1,R2,C2,false);
        y = fminsearch(@requiredVoltage_all, [R1 C1 R2 C2],optimset('Display', 'off'));
    end;
    num_of_results = num_of_results + 1;
    result(num_of_results,1) = abs(R1);
    result(num_of_results,2) = abs(C1);
    result(num_of_results,3) = abs(R2);
    result(num_of_results,4) = abs(C2);
    [result(num_of_results,5) result(num_of_results,6)] = checkFitting(R1,C1,R2,C2,DO_PLOTS);
    FITTING_NORM = 2;
    temp_res = +Inf;
    while temp_res - checkFitting(R1,C1,R2,C2,false) > 1e-13
        temp_res = checkFitting(R1,C1,R2,C2,false);
        y = fminsearch(@requiredVoltage_all, [R1 C1 R2 C2],optimset('Display', 'off'));
    end;
    num_of_results = num_of_results + 1;
    result(num_of_results,1) = abs(R1);
    result(num_of_results,2) = abs(C1);
    result(num_of_results,3) = abs(R2);
    result(num_of_results,4) = abs(C2);
    [result(num_of_results,5) result(num_of_results,6)] = checkFitting(R1,C1,R2,C2,DO_PLOTS);


    [min_val min_index] = min(result(:,5));
    res = result(min_index,:);

    % To plot the final fit only:
    checkFitting(res(1),res(2),res(3),res(4),DO_PLOTS_FINAL_ONLY);



end



