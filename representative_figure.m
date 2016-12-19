clear all;clc;
cd('D:\Keti\Data\OptoInhibition\cortical gain\Exp files\Database_cortical gain');
load '2013-11-13-001.mat';
trace1 = sgolayfilt(a.avetrail2(4,2000:14000),1,11).*20;
trace2 = sgolayfilt(a.avetrail2(5,2000:14000)-a.avetrail2(2,2000:14000),1,11).*20;
trace1 = trace1 - trace1(3000);
trace2 = trace2 - trace2(3000);
x = (0:1:length(trace1)-1)./10000;
f=figure(1)
clf
set(f,'units','normalized','position',[0.2 0.3 0.4 0.4],'color',[ 1 1 1]);
subplot(2,1,1);
plot(x,trace1,'k','linewidth',1.5);
hold on;
plot(x,trace2,'r','linewidth',1.5);
box off;
xlabel('Time (s)')
ylabel('Current (pA)')
l1 = legend('Total excitation','Thalamic excitation');
set(l1,'EdgeColor',[1 1 1]);
plot(x(1000:11000-1), ones(1,10000).*25,'linewidth',3);
 ylim([-70 30]);
hold on;
for i = 1:10
    line([0.3+(i-1)*0.05,0.3+(i-1)*0.05],[15 19] ,'color',[0 0 0],'linewidth',1.5);
    hold on;
end
subplot(2,1,2)
g= a.thalamic_input;
 bar(g(1,:))
 xlabel('Stimulus number');
 ylabel('Thalamic contribution');
 box off
 
 
 %%
clear all;clc;
path = 'D:\Keti\Data\OptoInhibition\cortical gain';
global Exp;
Exp = [];
fname = '2014-04-01-001.exp2';
expreader(path,fname);
cellData        = {};
Galvo_PW_Data   = [];
Galvo_AW_Data   = [];
PW              = 15;
AW              = 11; 
l               = [32 31 22 9 14];
H               = [3 4 5 6 7];
[s,no_H]        = size(H);
for ii=1:no_H
        last= l(ii);
        for i= 1:last
            blockData(i,:)  = readOneBlockdata(H(ii),i);
            channel1(i,:)   = blockData(i).channel(1).chdata; % Vm values
            %channel10(i,:)  = blockData(i).channel(10).chdata;% Laser
            channelPW(i,:)  = blockData(i).channel(PW).chdata;% Galvo 0 train
            channelAW(i,:)  = blockData(i).channel(AW).chdata;% Galvo 1 train
            cellData{1,ii}(i,:)  = channel1(i,:); %contains cells. in each cell 1 matrix containing x repetion of the same trail defined by the x value
            Galvo_PW_Data(i,:)    = channelPW(i,:);
            Galvo_AW_Data(i,:)    = channelAW(i,:);
            %laserData(Xval,:)        = channel10(i,:);
            spikeData{1,ii}(i,:) = zeros(1,length(channel1(i,:)));
            temp=[];
            temp = SpikeCountKeti(channel1(i,:));
            if ~isempty(temp) && (length(temp)>1 || temp ~= 0)
                spikeData{1,ii}(i,temp)=1;
            end
        end
end


for i = 1:5
    aa=1;
    for j=1:10:length(cellData{1,i}(1,:))
        sumSpike(i,aa) = sum(sum(spikeData{1,i}(:,j:j+9)));
        aa=aa+1;
    end
    sumSpike(i,:)= sumSpike(i,:)./l(i);
end

f=figure(1);
clf
set(f,'units','normalized','position',[0.2 0.3 0.4 0.4],'color',[ 1 1 1]);


for i = 1:2
    subplot(2,1,i);
    bar(sumSpike(i,:),'k');
    if i == 1
        title('Control','fontsize',15);
    else
        title('LED','fontsize',15);
    end
    hold on;
    v=0;
    for ii=1:10
        v=v+50;
        line([v v],[-0.02 0],'color',[1 0 0],'linewidth',3);
    end
     ylim([-0.02 0.2]);
    xlim([0 1000]);
    ylabel('Spikes per stimulus');
    xlabel('Time (ms)');
    box off;
end
line([40 800],[0.2 0.2],'color',[0 0 1],'linewidth',3);

    