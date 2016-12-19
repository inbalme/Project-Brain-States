
% Does Y Granger Cause X?
% [F,c_v] = granger_cause(x,y,alpha,max_lag)
alpha=0.05;
max_lag=1000;
fileind=1;
trace=2;
t=1;
x=sig2{fileind,t}(:,trace);
y=sig1{fileind,t}(:,trace);

[F,c_v] = granger_cause(x,y,alpha,max_lag);