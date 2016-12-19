%This function takes two simultaneously recorded signals and shuffle them using bootstrap
%The input matrices should be after DC-subtraction (subtracting the mean of
%the trace)
% inputs:
%1. data1,data2: Two matrices with each column a trial.
%2. iterations: number of bootstrap iterations
%optional: DC_sub. if DC_sub==0 then the mean was already subtracted from the traces. If DC_sub==1,
%then this function will first subtract the mean from each trace and then return two matrices with corresponding shuffled traces.
%The default is DC_sub==0;

function [data1_shuff, data2_shuff] = fn_bootstrap_shuff_trace_couples(data1, data2,iterations,DC_sub)
        trials=size(data1,2);
        data1_shuff=[]; data2_shuff=[];
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
                    data1_shuff_temp(:,j) = data1(:,w1); 
                    data2_shuff_temp(:,j) = data2(:,w2);                                
                end
               data1_shuff = [data1_shuff, data1_shuff_temp]; 
               data2_shuff = [data2_shuff, data2_shuff_temp]; 
        end