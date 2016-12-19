% The purpose of this file is to compare the performance of difference
% deconvolution methods for known data and Inbal's data, and to apply the
% method to Inbal's data
clear all
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
load NBES_Files_v2
files_to_analyze =46;
fileind=1;
fname = files(files_to_analyze).extracted_name; 
load(fname) 
channel=1;
sf{1} = Param.sf_Vm;
sf{2} = Param.sf_I1;
sf{3} = Param.sf_V2;
sf{4} = Param.sf_I2;
dt=1/sf{channel};
start=0.5;
stop=2.5;
bp_filt=files(files_to_analyze(fileind)).V2_filter;

trials = size(raw_data{3},2);
for xx=1:3
   for trace=1:trials
        jl=raw_data{3}(:,trace,xx);
        Ch2_data_filt(:,trace,xx)=bandPass_fft_IL_NEW2016(jl,dt,49,51,1,0);
        jm=data_no_spikes{channel}(:,trace,xx);
        data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(jm,dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm
        km=data_no_spikes_filt{channel}(:,trace,xx);
        data_no_spikes_filt{channel}(:,trace,xx)=bandPass_fft_IL_NEW2016(km,dt,bp_filt(1),bp_filt(2),0,0); %filtering Vm same as LFP
   end
end

% Below is a comparison of regularized and unregularized filters

Func=gaussmf(-0.05:0.0001:(0.05-0.0001),[0.01 0]);
Vmclean=zeros(Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),trace,xx),8,3);

for xx=1:3
for trace=1:trials
A=conv(Func,Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),trace,xx));
Vmclean(1:length(Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),trace,xx)),trace,xx)=A(length(Func)/2:end-length(Func)/2);
end
end


decon=Ldeconvs([Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),1,1); zeros(length(Ch2_data_filt(0.5*(1/dt):(2.5*(1/dt)-1),1,1))-length(Vmclean(:,1,1)), 1)], Vmclean(:,1,1), 1000, 0.05,1);
nodecon=Ldeconvs([Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),1,1); zeros(length(Ch2_data_filt(0.5*(1/dt):(2.5*(1/dt)-1),1,1))-length(Vmclean(:,1,1)), 1)], Vmclean(:,1,1), 1000, 0.05,0);

figure
hold all;
plot(-50:0.1:(50-0.1),decon(1:1000),'LineWidth',1.5)
plot(-0.05:0.0001:(0.05-0.0001),nodecon(1:1000),'LineWidth',1.5)
plot(-50:0.1:(50-0.1),Func,'LineWidth',1.5)
box('Off')



% Below is an implementation of cross-validation for our data

kernel=zeros(8,length(Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),trace,xx)));

for xx=1
        for trace=1:8
    kernel(trace,:)=Ldeconvs(Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),trace,xx),data_no_spikes_filt{channel}(start*(1/dt):(stop*(1/dt)-1),trace,xx),length(Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),trace,xx)),0.05,1);
        end
end

for part=1:8
kernelc=kernel;
kernelc(part,:)=[];
kernelc=mean(kernelc,1);
pred=conv(kernelc',Ch2_data_filt(start*(1/dt):(stop*(1/dt)-1),part,1));
pred=pred(size(kernelc,2)/2:end-size(kernelc,2)/2);
corr2(data_no_spikes_filt{channel}(start*(1/dt):(stop*(1/dt)-1),part,1),pred)^2;
end



% Below is an implementation of the method in Haider et al. 

% sigma=1;
% lambda=0.05;
% 
% for part=1
% startnew=
% stopnew=
%     
% LFPlarge(1,part)=Ch2_data_filt(0.5*(1/dt):(2.5*(1/dt)-1),1,1);
% 
% end

