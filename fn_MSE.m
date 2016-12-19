% This function gets two vectors and computes the mean squared error between the elements.
% If the inputs are matrices, the function treats the columns as vectors
%optional: DC_sub. if DC_sub==0 then the mean was already subtracted from the traces. 
% If DC_sub then the function will first subtract the mean from each trace.
%The default is DC_sub==0;

function mse = fn_MSE(a,b,DC_sub)
assert(isequal(size(a),size(b)), 'inputs should be the same size')
if DC_sub==1;
    a=a-mean(a);
    b=b-mean(b);
end
errors=b-a;
sq_errors=errors.^2;
mse=mean(sq_errors);
end
