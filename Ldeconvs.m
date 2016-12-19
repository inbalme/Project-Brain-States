% Prediction of a transfer function between signal from x to y
% Done in the frequency domain.
% Inputs: x,y - the two signals (vectors of equal length)
%           L - length of the kernel we are looking for (by default L = numel(x))
%        wlev - "water level" for regularization, by default 5%  
% Output: G - the kernel, such that wconv1(G,x) is y's estimate 

function G = Ldeconvs(x, y, L, wlev,reg_flag)

assert(numel(x) == numel(y), 'wrong inputs')

if nargin < 3
  L = numel(x);
end
if nargin < 4
  wlev = 0.05;
end
if nargin < 5
    reg_flag = 1;
end

% We break the signals into L parts, for each of which we will do the
% computation separately, and then average.
x = x(1:floor(length(x)/L)*L);
y = y(1:floor(length(y)/L)*L);

x = x-mean(x);
y = y - mean(y);

x = reshape(x, L, []);
y = reshape(y, L, []);


f_x = fft(x, [], 1);
f_y = fft(y, [], 1);

% Implementation of the WATER LEVEL REGULARIZATION, e.g. see chapter 8.4 of
% "Parameter Estimation and Inverse Problems" book
wlev = wlev*max(abs(f_x(:)));
f_xw = f_x;
I = abs(f_x) < wlev;
f_xw(I) = wlev*(f_x(I)./abs(f_x(I)));
f_xw(f_x == 0) = wlev;

if reg_flag==1;
    G = ifft(mean(f_y ./ f_xw, 2));
else
    G = ifft(mean(f_y ./ f_x, 2));
end

% G = ifft(f_y ./ f_xw, [], 1);
% if ismatrix(G) % i.e. L is such that we had to break into several pieces) 
%   G = mean(G, 2);
% end

% Finally we "half-rotate" G
G = circshift(torow(G)', numel(G)/2)';
% The reason is the following. The full solution for something convolved with x (of length L) that gives y (of length L) is a kernel
% of length 2*L-1. The fft thing we've done only returns something of length L, therefore the interesting central part starts not in the centre
% but already at the edge (and then loops over) 


return


%% Example of using the function

N = 211;
t=linspace(-5,100,N);
sig=1;
g = (exp(-((t(1:N-1)-8).^2/(sig^2*2)))+0.5*exp(-((t(1:N-1)-20).^2/(sig^2*2)))+0.25*exp(-((t(1:N-1)-30).^2/(sig^2*2))))';
g=g/max(g);

sig=2;
% mtrue = (exp(-((t(1:N-1)-8).^2/(sig^2*2)))+0.5*exp(-((t(1:N-1)-25).^2/(sig^2*2))))';
% mtrue=mtrue/max(mtrue);

mtrue = randn(1e4, 1);
mtrue = removeFrequencies(mtrue, 100, 20, 1e4);

d = conv(g,mtrue);
d = d(length(g)/2:end-length(g)/2);
d = d + 0.05*randn(length(d),1);
G = Ldeconvs([mtrue; zeros(length(d)-length(mtrue), 1)], d, 300, 0.3);
figure, plot(G)
