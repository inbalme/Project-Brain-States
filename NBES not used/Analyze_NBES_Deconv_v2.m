% The purpose of this file is to compare the performance of difference
% deconvolution methods for known data and Inbal's data, and to apply the
% method to Inbal's data

cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
%allfiles=[44,46,48,50,52,62,72,75]; these are all the files
allfiles=46; %[44,46,48,50,52,62,72,75];
kernelsize=1002; % 1002 means 50 ms back and forth
intervalsize=20000;
totallength=200000;
type='ev'; % 'ev' for explained variance (R^2), 'mse' for mean-squared error
alltimes=cell(1,totallength/intervalsize-1);

for i=1:totallength/intervalsize-1
    
alltimes{1,i}=[(i)*intervalsize/10000 (i+1)*intervalsize/10000];  

end
% alltimes={[0.5 2.5],[4 6],[6 8],[8 10],[10 12]};

for stim=1

pred=zeros(size(alltimes,2),intervalsize);
preds=zeros(size(alltimes,2),intervalsize);

expvar=zeros(length(allfiles),2,25);
corr=zeros(length(allfiles),2,25);
expvars=zeros(length(allfiles),2,25);
corrs=zeros(length(allfiles),2,25);

allkernels=cell(1,size(alltimes,2));
allkernelss=cell(1,size(alltimes,2));

for file=allfiles
fname = files(file).extracted_name; 
load(fname) 
channel=1;
sf{1} = Param.sf_Vm;
sf{2} = Param.sf_I1;
sf{3} = Param.sf_V2;
sf{4} = Param.sf_I2;
dt=1/sf{channel};
bp_filt=files(file).V2_filter;


trials = size(raw_data{3},2);

for xx=stim
   for trace=1:trials
        jl=raw_data{3}(:,trace,xx);
        jr=bandPass_fft_IL_NEW2016(jl,dt,49,51,1,0);
        Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(jr,dt,bp_filt(1),bp_filt(2),0,0);
        jm=data_no_spikes{channel}(:,trace,xx);
        data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jm,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
        km=data_no_spikes_filt{channel}(:,trace,xx);
        data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,bp_filt(1),bp_filt(2),0,0); %filtering Vm same as LFP
   end
end



% Below is an implementation of cross-validation for our data



for ba=alltimes
ba=ba{1};
for n=1:size(alltimes,2)
if isequal(ba,alltimes{n})
epoch=n;
end
end
kernel=zeros(trials,kernelsize);
kernelc=zeros(trials,kernelsize);


% normal

for xx=stim
        for trace=1:trials
            LFP=Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),trace,xx)/20;
            Vm=data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),trace,xx);
%              LFP=raw_data{3}(ba(1)*(1/dt)+1:(1.2*(1/dt)),trace,xx)/20;
%              Vm=data_no_spikes{channel}(ba(1)*(1/dt)+1:(1.2*(1/dt)),trace,xx);
            kernel(trace,:)=Ldeconvs(LFP,Vm,kernelsize,0.5,1);
%             
%             the problems start here - no semblance of trace although t=1
%             to 2
%             Vpred=conv(kernel(trace,:)-mean(kernel(trace,:)),LFP-mean(LFP))+mean(Vm);
%             figure
%             subplot(3,1,1)
%             hold on;
%             plot(-Vpred(size(kernelc,2)/2:end-size(kernelc,2)/2)','Color','b')
%             plot(-Vm,'Color',[0.6 0.6 0.6])
%             subplot(3,1,2)
%             hold on;
%             plot(LFP','Color','b')
%             xlabel('Time (s)')
%             ylabel('mV') 
%             subplot(3,1,3)
%             plot(kernel(trace,:))

        end
end

for part=1:trials
kernelc=kernel;
kernelc(part,:)=[];
kernelc=mean(kernelc,1);
allkernels{1,epoch}=vertcat(allkernels{1,epoch},kernelc);

predall=conv(kernelc',Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)/20)+mean(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx));
pred(epoch,:)=predall(size(kernelc,2)/2:end-size(kernelc,2)/2);

if part==3
figure
hold on;
if epoch==1
plot(0:0.0001:2-0.0001,pred(epoch,:)','Color','b')
else
plot(0:0.0001:2-0.0001,pred(epoch,:)','Color','r')
end
plot(0:0.0001:2-0.0001,data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx),'Color',[0.6 0.6 0.6])
xlabel('Time (s)')
ylabel('mV')    
end

if isequal(type,'mse')
expvar(file==allfiles,epoch,part)=immse(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx),pred(epoch,:)');

elseif isequal(type,'ev')
expvar(file==allfiles,epoch,part)=corr2(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx),pred(epoch,:)')^2;
end

corr(file==allfiles,epoch,part)=max(xcorr(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)-mean(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)),Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)'./20-mean(Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)'./20),'coeff'));
end


% shuffled

for xx=stim
        alltraces=randperm(trials);
        for trace=1:trials
            Vmtrace=alltraces(trace);
            LFP=Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),trace,xx)/20;
            Vm=data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),Vmtrace,xx);
            kernel(trace,:)=Ldeconvs(LFP,Vm,kernelsize,0.5,1);
            corrs(file==allfiles,epoch,trace)=max(xcorr(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),Vmtrace,xx)-mean(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),Vmtrace,xx)),Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),trace,xx)./20'-mean(Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),trace,xx)'./20),'coeff'));
        end
end

for part=1:trials
kernelc=kernel;
kernelc(part,:)=[];
kernelc=mean(kernelc,1);
allkernelss{1,epoch}=vertcat(allkernelss{1,epoch},kernelc);
predall=conv(kernelc',Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)/20)+mean(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),alltraces(part),xx));
preds(epoch,:)=predall(size(kernelc,2)/2:end-size(kernelc,2)/2);
if isequal(type,'mse')
expvars(file==allfiles,epoch,part)=immse(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),alltraces(part),xx),preds(epoch,:)');
elseif isequal(type,'ev')
expvars(file==allfiles,epoch,part)=corr2(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),alltraces(part),xx),preds(epoch,:)')^2;
end
%corr(1,part)=max(xcorr(data_no_spikes_filt{channel}(ba(1)*(1/dt):(ba(2)*(1/dt)-1),part,xx)-mean(data_no_spikes_filt{channel}(ba(1)*(1/dt):(ba(2)*(1/dt)-1),part,xx)),pred(part,:)'-mean(pred(part,:)'),'coeff'));
end

end



end


% plotting the combined explained variances and xcorrs

sum=cell(1,size(alltimes,2));
for filenum=1:length(allfiles)
    for epoch=1:size(alltimes,2)
sum{1,epoch}=horzcat(sum{1,epoch},mean(expvar(filenum,epoch,expvar(filenum,1,:)~=0),3));
    end
end

figure
% plot(vertcat(sum{1,:}),'LineWidth',1.5)
set(gca,'xtick',1:size(sum,2));
hold on;
xlimit=size(sum,2)+0.2;
xlim([0.8 xlimit])
ylim([0 1])
set(gca,'XTickLabel',1:size(sum,2));
errorbar(mean(vertcat(sum{1,:}),2),std(vertcat(sum{1,:}),0,2)/8,'color',[0 0 0],'LineWidth',3)

sum=cell(1,size(alltimes,2));
for filenum=1:length(allfiles)
    for epoch=1:size(alltimes,2)
sum{1,epoch}=horzcat(sum{1,epoch},mean(expvars(filenum,epoch,expvars(filenum,1,:)~=0),3));
    end
end


set(gca,'xtick',1:size(sum,2));
hold on;
xlimit=size(sum,2)+0.2;
xlim([0.8 xlimit])
set(gca,'XTickLabel',1:size(sum,2));
errorbar(mean(vertcat(sum{1,:}),2),std(vertcat(sum{1,:}),0,2)/sqrt(64),'color',[0.5 0.5 0.5],'LineWidth',3)
ylabel('Explained variance')
xlabel('Epoch')
set(gca,'fontsize',14)




% sum=cell(1,2);
% for file=1:length(allfiles)
%     for epoch=1:size(alltimes,2)
% sum{1,epoch}=horzcat(sum{1,epoch},mean(corr(file,epoch,corr(file,1,:)~=0),3));
%     end
% end
% 
% figure
% plot([sum{1,1}; sum{1,2}],'color',[0 0 0.4],'LineWidth',1)
% set(gca,'xtick',[1 2]);
% hold on;
% xlim([0.8 2.2])
% set(gca,'XTickLabel',{'Before', 'After'});
% errorbar([mean(sum{1,1}); mean(sum{1,2})],[std(sum{1,1})/8; std(sum{1,2})/8],'color',[0 0 0],'LineWidth',3)
% 
% 
% ranksum(sum{1,1},sum{1,2})
% 
% sum=cell(1,2);
% for file=1:length(allfiles)
%     for epoch=1:2
% sum{1,epoch}=horzcat(sum{1,epoch},mean(corrs(file,epoch,corrs(file,1,:)~=0),3));
%     end
% end
% 
% set(gca,'xtick',[1 2]);
% hold on;
% xlim([0.8 2.2])
% set(gca,'XTickLabel',{'Before', 'After'});
% errorbar([mean(sum{1,1}); mean(sum{1,2})],[std(sum{1,1})/8; std(sum{1,2})/8],'color',[0.5 0.5 0.5],'LineWidth',3)
% ylim([0.2 0.6])
% ylabel('Cross correlation peak')
% 
% ranksum(sum{1,1},sum{1,2})






% figure
% hold on;
% all=allkernels{1,1};
% s1=shadedErrorBar(-50:0.1:50,mean(all(:,length(all)/2-500:length(all)/2+500),1),std(all(:,length(all)/2-500:length(all)/2+500))/8,'-b',1);
% plot(-50:0.1:50,mean(all(:,length(all)/2-500:length(all)/2+500),1),'color',[0 0 1 0.6],'LineWidth',3)
% all=allkernels{1,2};
% s2=shadedErrorBar(-50:0.1:50,mean(all(:,length(all)/2-500:length(all)/2+500),1),std(all(:,length(all)/2-500:length(all)/2+500))/8,'-r',1);
% plot(-50:0.1:50,mean(all(:,length(all)/2-500:length(all)/2+500),1),'color',[1 0 0 0.6],'LineWidth',3)
% set(gca,'Ydir','reverse')
% ylim([-50*10^(-3) 20*10^(-3)])
% alpha(s1.patch,0.5)
% alpha(s2.patch,0.5)
% xlabel('Lag (ms)')
% ylabel('Kernel (mV per mV)')
% 
% figure
% hold on;
% all=allkernelss{1,2};
% s2=shadedErrorBar(-50:0.1:50,mean(all(:,length(all)/2-500:length(all)/2+500),1),std(all(:,length(all)/2-500:length(all)/2+500))/8,'-r',1);
% plot(-50:0.1:50,mean(all(:,length(all)/2-500:length(all)/2+500),1),'color',[1 0 0 0.3],'LineWidth',3)
% all=allkernelss{1,1};
% s1=shadedErrorBar(-50:0.1:50,mean(all(:,length(all)/2-500:length(all)/2+500),1),std(all(:,length(all)/2-500:length(all)/2+500))/8,'-b',1);
% plot(-50:0.1:50,mean(all(:,length(all)/2-500:length(all)/2+500),1),'color',[0 0 1 0.3],'LineWidth',3)
% alpha(s1.patch,0.05)
% alpha(s2.patch,0.05)
% set(gca,'Ydir','reverse')
% ylim([-50*10^(-3) 20*10^(-3)])
% xlabel('Lag (ms)')
% ylabel('Kernel (mV per mV)')
% 


% figure
% plot(0:20,5)
% hold on;
% set(gca,'xtick',0:20);
% set(gca,'fontsize',14)
% xlabel('Time (s)')
% xlabel('Epoch')
% box('Off')



% % Below is a comparison of regularized and unregularized filters
% 
% Func=Funcmf(-0.05:0.0001:(0.05-0.0001),[0.01 0]);
% Vmclean=zeros(Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),trace,xx),8,3);
% 
% for xx=1:3
% for trace=1:trials
% A=conv(Func,Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),trace,xx));
% Vmclean(1:length(Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),trace,xx)),trace,xx)=A(length(Func)/2:end-length(Func)/2);
% end
% end
% 
% 
% decon=Ldeconvs([Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),1,1); zeros(length(Ch2_data_filt(0.5*(1/dt):(2.5*(1/dt)-1),1,1))-length(Vmclean(:,1,1)), 1)], Vmclean(:,1,1), 1000, 0.05,1);
% nodecon=Ldeconvs([Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),1,1); zeros(length(Ch2_data_filt(0.5*(1/dt):(2.5*(1/dt)-1),1,1))-length(Vmclean(:,1,1)), 1)], Vmclean(:,1,1), 1000, 0.05,0);
% 
% figure
% hold all;
% plot(-50:0.1:(50-0.1),decon(1:1000),'LineWidth',1.5)
% plot(-0.05:0.0001:(0.05-0.0001),nodecon(1:1000),'LineWidth',1.5)
% plot(-50:0.1:(50-0.1),Func,'LineWidth',1.5)
% box('Off')

% debugging

% 
% 
% ba(1)=6;
% ba(2)=8;
% epoch=1;
% 
% for xx=2
%    for trace=1:trials
%         jl=raw_data{3}(:,trace,xx);
%         jr=bandPass_fft_IL_NEW2016(jl,dt,49,51,1,0);
%         Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(jr,dt,bp_filt(1),bp_filt(2),0,0);
%         jm=data_no_spikes{channel}(:,trace,xx);     
%         data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jm,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
%         km=data_no_spikes_filt{channel}(:,trace,xx);
%         data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,bp_filt(1),bp_filt(2),0,0); %filtering Vm same as LFP
%    end
% end
% 
% 
% 
%         for trace=1:trials
%             LFP=Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),trace,xx)/20;
%             Vm=data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),trace,xx);
%             kernel(trace,:)=Ldeconvs(LFP,Vm,kernelsize,0.5,1);
%         end
% 
% for part=3
% kernelc=kernel;
% kernelc(part,:)=[];
% kernelc=mean(kernelc,1);
% allkernels{1,epoch}=vertcat(allkernels{1,epoch},kernelc);
% predall=conv(kernelc',Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)/20)+mean(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx));
% pred(epoch,:)=predall(size(kernelc,2)/2:end-size(kernelc,2)/2);
% if isequal(type,'mse')
% expvar(file==allfiles,epoch,part)=immse(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx),pred(epoch,:)');
% elseif isequal(type,'ev')
% expvar(file==allfiles,epoch,part)=corr2(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx),pred(epoch,:)')^2;
% end
% corr(file==allfiles,epoch,part)=max(xcorr(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)-mean(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)),Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)'./20-mean(Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)'./20),'coeff'));
% end
% 
% for part=1:trials
% kernelc=kernel;
% kernelc(part,:)=[];
% kernelc=mean(kernelc,1);
% allkernels{1,epoch}=vertcat(allkernels{1,epoch},kernelc);
% 
% predall=conv(kernelc',Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)/20)+mean(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx));
% pred(epoch,:)=predall(size(kernelc,2)/2:end-size(kernelc,2)/2);
% 
% if isequal(type,'mse')
% expvar(file==allfiles,epoch,part)=immse(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx),pred(epoch,:)');
% elseif isequal(type,'ev')
% expvar(file==allfiles,epoch,part)=corr2(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx),pred(epoch,:)')^2;
% end
% corr(file==allfiles,epoch,part)=max(xcorr(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)-mean(data_no_spikes_filt{channel}(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)),Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)'./20-mean(Ch2_data_filt(ba(1)*(1/dt)+1:(ba(2)*(1/dt)),part,xx)'./20),'coeff'));
% end
% 
% 
% 
% figure
% subplot(2,1,1)
% hold all
% plot(Ch2_data_filt(60000:80000,3,2))
% plot([30000 30000],[-10000 10000],'LineStyle','--','color','y')
% plot([35000 35000],[-10000 10000],'LineStyle','--','color','y')
% plot([00000 00000],[-10000 10000],'LineStyle','--','color','r')
% plot([15000 15000],[-10000 10000],'LineStyle','--','color','r')
% ylim([-15 15])
% xlim([00000 20000])
% subplot(2,1,2)
% hold all
% plot(-data_no_spikes_filt{1}(60000:80000,3,2))
% plot([30000 30000],[-10000 10000],'LineStyle','--','color','y')
% plot([35000 35000],[-10000 10000],'LineStyle','--','color','y')
% plot([00000 00000],[-10000 10000],'LineStyle','--','color','r')
% plot([15000 15000],[-10000 10000],'LineStyle','--','color','r')
% ylim([20 100])
% xlim([00000 20000])
% 
% 
% plot(-pred(1,:)','Color','r')
% xlabel('Time (s)')
% ylabel('mV')


end
