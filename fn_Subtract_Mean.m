function [ mat_no_mean ] =  fn_Subtract_Mean( mat, interval )
%This function subtracts the mean of each trace from all the data points in that trace.
%   mat is a 2D matrix with traces as column vectors. for a 3D matrix,
%   somehow the wextend function works only on z<=3, i.e., only 3 x-values
%if only one input argument is provided - the function will calculate the
%mean based on the whole trace and then subtract it. if "interval" input is
%provided - the function will calculate the mean over this interval and
%then subtract it from the whole trace.
if nargin==1
    mat_mean= mean(mat,1);
else if nargin==2
        mat_mean=mean(mat(interval,:),1);
    end
end
l=size(mat,1);
if mod(l,2)==0 % if l is not odd then add another row if you get odd number of rows after the extension.
l=ceil(l/2-1);
mat_mean=wextend('addrow', 'per', mat_mean, l);
else % if l is odd
    l=ceil(l/2-1);
mat_mean=wextend('addrow', 'ppd', mat_mean, l); % do not add another row
end

mat_no_mean =mat-mat_mean;


end
