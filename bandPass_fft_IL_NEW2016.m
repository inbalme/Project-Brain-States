function y = bandPass_fft_IL_NEW2016(x,dt,f1,f2,bandstop,ploton);
%bandPass(x,dt,f1,f2);
% filter 1d array (x) using dt and frequencies 1 and 2; Everything in
% between these two frequencies will be filterout.
% filtering.
%For low pass filtering use f1 negative.
%for high pass filtering use f2 that is negative.

%for f1 <0 low pass
% for f2 <0 high pass


%If f2 is larger than the nyquest frequency the data f1 will used as cuttoff for high passs

%if nothing said data between f1 and f2 will be preserved. the rest will be killed. If f2 is larger than
% than the nyquest freq it will be ignored.
%if bandstop == 1, meaning cleaning in the range of f1 to f2,
%if bandstop == 0; meaning filterour below f1 and above f2

if nargin <5;
    bandstop = 0;
end;



meandata = mean(x);
x = x - meandata;

fmax = 1/dt;
fx = fft(x);
fxo = fx; %for testing
df = 1/(dt*length(x));

f1i = round(f1/df);
f2i = round(f2/df);

l = length(fx);

display(strcat('nyquest    ',num2str(l)))
if f1<0 & f2<0% do not change the data;
    y = x;
end;

if f1<0 & f2 >0 % low pass filtering.
    l = length(fx);
    fx(f2i:ceil(l/2)) = 0;
    fx(ceil(l/2):end-f2i) = 0;
    display(strcat('perfoming:     ','lowpass'))
end;

if f1>0 & f2 <0; % high pass filtering
    l = length(fx);
    fx(1:f1i) = 0;
    fx(end-f1i:end) = 0;
    display(strcat('perfoming:     ','highpass'))
end;

if f1>0 & f2 >0; % band pass filtering
    l = length(fx);
    if f1i <l & f2i <l ;% both below nyquest
        if bandstop;
            l = length(fx);
            fx(f1i:f2i) = 0;
            fx(end-f2i:end-f1i) = 0;
            display(strcat('perfoming:     ','bandstop'))
            
        else
           
            %performing 
            fx(1:f1i) = 0; %clean below f1
            fx(end-f1i:end) = 0;
            
            fx(f2i:ceil(l/2)) = 0; %clean below f1 and above f2;
            fx(ceil(l/2):end-f2i) = 0;
            display(strcat('perfoming:     ','bandpass'))
        end;
    end;
end
    
    
    
    y = real(ifft(fx));
    
    
    %hold;
    if nargin>4;
        
        if ploton
            close all
            figure
            plot([1:1:length(fx)]*df,abs(fxo));
            
            hold on;
            
            plot([1:1:length(fx)]*df,abs(fx),'r');
            legend('original','filtered')
            
            
            figure
            plot(x,'b');
            hold on;
            plot(y,'r')
            
            legend('original','filtered')
            
        end
    end
    
    y = y+meandata;
    
    
    
    
