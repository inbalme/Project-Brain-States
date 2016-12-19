fftbefore = [];
for kk = 1:20;
    v = before1(:,kk)-mean(before1(:,kk));
    fftbefore(kk,:) = abs(fft(v));
end;

fftduring = [];
for kk = 1:20;
    v = during1(:,kk)-mean(during1(:,kk));
    fftduring(kk,:) = abs(fft(v));
end;

%%

during1 = spec_mat;
size(before1)

ans =

       50000          20

for kk = 1:20;
    v = before1(:,kk)-mean(before1(:,kk));
    fftbefore(kk,:) = fft(v);
end;
fftbefore = [];
for kk = 1:20;
    v = before1(:,kk)-mean(before1(:,kk));
    fftbefore(kk,:) = fft(v);
end;

fftduring = [];
for kk = 1:20;
    v = during1(:,kk)-mean(during1(:,kk));
    fftduring(kk,:) = fft(v);
end;

fftbefore = [];
for kk = 1:20;
    v = before1(:,kk)-mean(before1(:,kk));
    fftbefore(kk,:) = abs(fft(v));
end;

fftduring = [];
for kk = 1:20;
    v = during1(:,kk)-mean(during1(:,kk));
    fftduring(kk,:) = abs(fft(v));
end;
figure
plot(mean
 plot(mean
          |
Error: Expression or statement is incorrect--possibly unbalanced (, {, or [.
 
plot(mean(fftbefore))
set(gca,'xlim',[0 100])
hold on;
plot(mean(fftduring),'r')
dt = 1/10000;
df = 1/(length(v)*dt)

df =

    0.2000

clf
plot(df*[1:1:length(fftduring)],mean(fftduring),'r')
set(gca,'xlim',[0 100])
hold on;
plot(mean(fftbefore))
clf;
plot(df*[1:1:length(fftduring)],mean(fftduring),'r')
hold on;
plot(df*[1:1:length(fftduring)],mean(fftbefore))
set(gca,'xlim',[0 100])
set(gca,'xlim',[0 30])
set(gca,'xlim',[0 50])
