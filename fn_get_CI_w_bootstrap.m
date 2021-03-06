% Computing Confidence Interval with bootstrap method

%dataset is a row vector or a matrix. the number of elements in the dataset
%is the number of rows. 
%k is the number of samples to pick from the dataset *with substitution*.
%if k==0 the function will use the number of columns as the number of
%samples.
%lower_bound and upper_bound are the percentiles for the CI. the default is
%2.5 and 97.5
%if the number of iterations is not provided, the default is 5000
%type can take 'mean' or 'median'. type speifies whether to take the mean or the median of the dataset. the
%default is the mean.

function [prcntile1, prcntile2]=fn_get_CI_w_bootstrap(dataset,k,iterations, lower_bound, upper_bound,type)

if k==0;
    k= size(dataset,2);     
end
if nargin<5
    lower_bound=2.5;
    upper_bound=97.5;
end
if nargin<3
    iterations=5000;
end
if nargin<6
    type='mean';
else if sum(strcmpi({'mean','median'},{type}))==0
   error('wrong string input')     
    end
end
 for bootstrap_ind=1:iterations;
     resampled_cells=datasample([1:k],k); %with replacement is the default
         tmp_set(:,:)= dataset(:,resampled_cells);
         if strcmpi('mean',type)
             bs_set(:,bootstrap_ind)=mean(tmp_set,2);
         else if strcmpi('median',type)
                 bs_set(:,bootstrap_ind)=median(tmp_set,2);
             end
         end
 end
 prcntile1 = prctile(bs_set,lower_bound,2);
 prcntile2 = prctile(bs_set,upper_bound,2);
 
end