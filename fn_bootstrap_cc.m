%This function calculates the mean cross-correlation between two signals 
%The input matrices should be after DC-subtraction (subtracting the mean of
%the trace)
% inputs:
%1. data1,data2: Two matrices with each column a trial.
%2. iterations: number of bootstrap iterations
%optional: DC_sub. if DC_sub==0 then the mean was already subtracted from the traces. If DC_sub==1,
%then this function will first subtract the mean from each trace and then calculate the
%cross-correlations. The default is DC_sub==0;

function [mean_cc_shuffled] = fn_bootstrap_cc(data1, data2,iterations,DC_sub)
        trials=size(data1,2);
        if nargin==4;
            if DC_sub==1   
                mean_traces_data1= mean(data1,1);
                l=size(data1,1);
                l=l/2-1;
                DC1=wextend('addrow', 'per', mean_traces_data1, l);
                data1 = data1-DC1;
 %subtract DC from LFP matrix
        mean_traces_data2= mean(data2,1);
        l=size(data2,1);
        l=l/2-1;
        DC2=wextend('addrow', 'per', mean_traces_data2, l);
        data2 = data2-DC2;
            end
        end
      
        for n = 1:iterations
                f = 1;
                while f
                    Indexs1 = randperm(trials);
                    Indexs2 = randperm(trials);
                    f = sum(Indexs1==Indexs2); %the loop will run until two sets of shuffled traces with no overlaps are generated
                end
                for j = 1:length(Indexs1) %running on all traces
                    w1 = Indexs1(j);
                    w2 = Indexs2(j);               
                    cc_shuffled_temp(:,j) = xcorr(data1(:,w1),data2(:,w2),'coeff');                                
                end
                cc_shuffled_temp2(:,n)=mean(cc_shuffled_temp,2); %average across traces
        end
        mean_cc_shuffled=mean(cc_shuffled_temp2,2); %average across bootstrap repetitions