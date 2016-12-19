% Deconvolution of 2 signals
% Inputs: X,Y - the two signals (row vectors, must be of the same length)
%                   L - the length of the kernel that is sought (1000 by default)
%                   doPLot (optional) - false by default
% Output: G - the optimal linear transfer function from X to Y
%                   cc - the peak of the crosscorrelation betwen  G*X and Y.

function [G cc] = Ldeconv(X,Y, L, doPlot)
if nargin < 4
    doPlot = false;
end;
if nargin < 3
    L = 1e3; % length of the kernel
elseif mod(L,2) == 1
  L = L+1;
end;
intLength = min(1e4, length(X));
G = zeros(floor(length(X)/intLength),L);
for step = 1:floor(length(X)/intLength)
    x = X((step-1)*intLength+1:step*intLength);
    y = Y((step-1)*intLength+1:step*intLength);
    
    x = x - mean(x);
    y = y - mean(y);
    
    A = zeros(L, length(x)-L+1);
    b = zeros(1, length(x)-L+1);
    
    for i = 1:length(x)-L+1
        A(:,i) = x(i:i+L-1)';
        b(i) = y(i+L/2);
    end;
    G(step,:) = b/A;
end;
if step > 1
    G = mean(G);
end;
G = G(end:-1:1);
Ypred = wconv1(G, X);
cc = crosscorr(Ypred(round(L/2):end), Y, round(L/2));
cc = max(cc);
if doPlot
    figure;
    subplot(2,3, [1 4]);
    plot(G);
    title('transfer function');
    subplot(2,3, [2 3]);
    crosscorr(X, Y, round(L/2));
    title('crosscorr(X,Y)');
    subplot(2,3, [5 6]);
    crosscorr(Ypred(round(L/2):end), Y, round(L/2));
    title('crosscorr(GX, Y)');
end;
