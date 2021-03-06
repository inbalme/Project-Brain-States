%% Electrical Stimulation of Nucleus Basalis
clear all; close all
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
% cd 'D:\Inbal M.Sc\Data PhD\DoubleNB\Extracted Data';

%different ES depth. same LFP location in barrel cortex:
% load Set2-2015-03-23-002_h1-7 %data1, different depths
%  load Set2-2015-04-02-001_h3-12 %different depths
 load 2016-10-20-001_h2-34 %different depths data1; [27,28,29,30,31,33,34], gain 10, DC_int=0.4;depth: [5200,5000,4700,4000, 3700,4700,3000, 2300]
%Atropine+Mecamylamine:
% load 2015-12-28-006_h1-6 %data2
% load 2015-07-29-002_h1-6 %data2 ESdel7sec, trace 15sec, cut 1500samples from both sides of the trace
%  load Set2-2015-05-06-003_h2-7 %data1 ESdel10sec, trace 20sec
%  load Set2-2015-06-01-003_h1-6 %data1 ESdel10sec
%GABA blockers:
% load 2015-12-31-009_h1-4 %picrotoxin data2
% load 2016-05-23-001_h3-9 %ch1-V1, ch2-S1
% load 2016-05-23-002_h1-4 %ch1-V1,ch2-S1+GABAb inhibitor
%regular traces:
% load 2015-12-31-002_h1-2 %data2
% load 2015-07-29-001_h2 %data2 ESdel7sec, trace 15sec, cut 1500samples from both sides of the trace
%  load Set2-2015-05-06-001_h2 %data1 ESdel10sec, trace 20sec, multiply DC by 2
% load 2015-12-17-001_h1-2 %data2
% load 2015-12-21-001_h2-3 %data2
% load 2015-12-22-001_h1-3 %data2
% load Set2-2015-05-18-001_h2 %data1
% load Set2-2015-05-24-001_h2
% load Set2-2015-05-25-001_h2
% load Set2-2015-05-27-001_h2
% load Set2-2015-05-29-002_h3-4
% load Set2-2015-06-01-001_h1
% load Set2-2015-06-01-003_h1-6
% load Set2-2015-06-03-001_h2+4
% load Set2-2015-06-05-001_h7+9
% load Set2-2015-06-10-001_h2-3
% load Set2-2015-06-15-001_h5+7
% load Set2-2015-06-18-001_h2-3
% load Set2-2015-06-25-001_h3-4
% load Set2-2015-06-29-001_h3-4
% load Set2-2015-07-07-001_h2+3
% load Set2-2015-07-09-001_h1+2
%  load 2016-10-24-001_h1-8 


%% 2015-05-06-001_LFP_traces_h2_t2  3  4  5  7; 
 orig_header=33; %[27,28,29,30,31,33,34] [1,2,3,4,5,6,7]; 2016-10-20-001_LFP_PSD_1-100Hz_h28_v3
 save_flag=0;
 galvano_flag=0;
 ES_del=10; %10
 dataCh=1;
 gain=10; %20;
 DC_int=0.4;
 depth=3000; %[2500,3000,3500,4000,4300,4600,4100]; 
 trace_ind = [2,3,4,5,6]; %1:size(plot_data,2);  %trace_ind is the numbers of traces to be plotted [2,3,4,5,6];
 header=find(headers==orig_header);
 stim2_X=[];
 sf = 4000;
 dt=1/sf;
 sf_galvano = 800;
 dt_galvano = 1/sf_galvano;
 
stim1_X=[ES_del*sf, (ES_del+0.5)*sf];           

 start_time=[ES_del-3, ES_del+2];
duration = 2.5; %[sec]
x_value = 1;
 for t=1:length(start_time);
start_sample(:,t) = start_time(t).*sf;
if start_time(t)==0
    start_sample(:,t) = 1;
end
end_sample(:,t) = start_sample(:,t)+duration.*sf-1;
interval(:,t) = start_sample(:,t):end_sample(:,t);
 end
%choosing data channel
 if dataCh==1;
    data_mat=data(header).x_value(1).data1;
 else
     data_mat=data(header).x_value(1).data2;
 end
 %dividing by the recording gain:
 data_mat=data_mat./gain;
%bandpass filtering to remove 50Hznoise from LFP channel, and then low pass 200 Hz
%filtering before and after ES separately because the high frequencies insert artifacts when filtered.
f1=1; %[-1];
    for trace= 1:size(data_mat,2)    
            jl=data_mat(:,trace);
            data_mat_filt1(1:stim1_X(1)-1,trace)=bandPass_fft_IL_NEW2016(jl(1:stim1_X(1)-1),dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm           
            data_mat_filt1(stim1_X(1):stim1_X(2),trace)=data_mat(stim1_X(1):stim1_X(2),trace);
            data_mat_filt1(stim1_X(2)+1:length(jl),trace)=bandPass_fft_IL_NEW2016(jl(stim1_X(2)+1:end),dt,49,51,1,0); %filtering out 50Hz noise from LFP and Vm            
            km=data_mat_filt1(:,trace);
             data_mat_filt2(1:stim1_X(1)-1,trace)=bandPass_fft_IL_NEW2016(km(1:stim1_X(1)-1),dt,f1,200,0,0); %low-pass 300Hz
            data_mat_filt2(stim1_X(1):stim1_X(2),trace)=data_mat(stim1_X(1):stim1_X(2),trace);
            data_mat_filt2(stim1_X(2)+1:length(km),trace)=bandPass_fft_IL_NEW2016(km(stim1_X(2)+1:end),dt,f1,200,0,0); %low-pass 300Hz
    end
    
 if galvano_flag==1;
 galvano2_trace = data(header).x_value(1).galvano2(:,2);
 time_axis_galvano =length(galvano2_trace)*dt_galvano;
 galvano_vec = zeros(length(galvano2_trace),1);
 galvano_sub_mean = galvano2_trace - mean(galvano2_trace);
%  galvano_threshold=0.5;
%  galvano_vec(galvano2_trace > galvano_threshold)=1;   %turning the galvano trace into binary
 galvano_threshold = abs(mean(galvano_sub_mean)) + 4.*abs(std(galvano_sub_mean));
 galvano_vec(abs(galvano_sub_mean) > galvano_threshold)=1;   %turning the galvano trace into binary
 end 
%                  if isempty(find(galvano_vec(:,1)==1, 1))~=1 %condition will be fulfilled if there was galvano activation.
                if galvano_flag==1;
                    galvano_vec_shifted(:,1) = galvano_vec(2:end,1)-galvano_vec(1:end-1, 1);
                    galvano_begin(:,1) = find(galvano_vec_shifted(:,1)==1); %find the locations where galvano pulse starts
                    galvano_begin(:,1) = galvano_begin(:,1)+1; %correction for the shift
                    galvano_end(:,1) = find(galvano_vec_shifted(2:end,1)==-1); %find the locations where galvano pulse ends
                        if length(galvano_begin(:,1))>length(galvano_end(:,1)) %if galvano_begin is larger than galvano_end, take only the points in galvano_begin which has a match in galvano_end
                            galvano_begin_trunc(:,1)=galvano_begin(1:length(galvano_end(:,1)),1);
                            galvano_begin(:,1)=[];
                            galvano_begin=galvano_begin_trunc(:,1);
                        end
                    locations_x_galvano(1,:) = galvano_begin(:,1); %arranging the galvano begin and end locations in one variable for plotting
                    locations_x_galvano(2,:) = galvano_end(:,1);         
                    
                     if isempty(locations_x_galvano)
                            stim2_X = [];
                        else
                            stim2_X = locations_x_galvano(:,:);
                        end
                end
                %% Plotting traces
trace_fontsize=12;
scalebar_fontsize=12;
plot_stim_1=[1,0,1];
plot_stim_2=[0,1,1];
x_axis=[];
x_axis_long=[1:size(data_mat_filt2,1)];
x_axis_short(1,:)=[start_sample(1,1):end_sample(1,1)];
x_axis_short(2,:)=[start_sample(1,2):end_sample(1,2)];

clear color_table
    color_table=[0 0 0; [216 22 22]/256; [136 137 138]/256; [255 153 153]/256];
    rectangle_color=[239 239 239]/256;
    
plot_data=data_mat_filt2; 
plot_data(stim1_X(1):stim1_X(2),:)=nan;
plot_data_mean = nanmean(plot_data,2);  
% plot_data_mean(stim1_X(1):stim1_X(2),:)=nan;
plot_data_std =  nanstd(plot_data,0,2); 
% plot_data_std(stim1_X(1):stim1_X(2),:)=nan;
% plot_data_var=var(plot_data,0,2); 
% plot_data_ff = (-1).*plot_data_var./plot_data_mean; 
       l=size(plot_data(:,trace_ind),1);
        l=l/2-1;
        DC=DC_int.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l);        
              
%         x2lim=[stim2_X{x_value(1)}(1,1).*dt-0.5,stim2_X{2}(1,1).*dt+2];
 
%%
    %plot spont traces long:
    trace_to_plot = plot_data(:,trace_ind,1)+DC ; %DC is added just to make space between the traces 
        f2=figure;
        hold on
        rec1=rectangle('Position',[x_axis_long(1,1)*dt,trace_to_plot(x_axis_long(1,1),1),0.1,0.1]);
        rec2=rectangle('Position',[x_axis_long(1,1)*dt,trace_to_plot(x_axis_long(1,1),1),0.1,0.1]);
htrace1=plot(x_axis_long(1500:ES_del*sf)*dt,trace_to_plot(x_axis_long(1500:ES_del*sf),:,:), 'LineWidth',1,'color', color_table(1,:));
htrace2=plot(x_axis_long((ES_del+0.5)*sf:(end-1500))*dt,trace_to_plot(x_axis_long((ES_del+0.5)*sf:(end-1500)),:,:), 'LineWidth',1,'color', color_table(2,:));

% for i=1:length(trace_ind);
%     text(x_axis_long(1)*dt,trace_to_plot(x_axis_long(1),i),[num2str(floor(plot_data(x_axis_long(1),i,1))), ' mV '],'HorizontalAlignment','right','fontsize',trace_fontsize,'fontname','arial')
% end
axis tight
ylim_data=[get(gca,'ylim')]';
xlim_data=get(gca,'xlim');
set(rec1,'Position',[x_axis_short(1,1)*dt,ylim_data(1),(x_axis_short(1,end)-x_axis_short(1,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
set(rec2,'Position',[x_axis_short(2,1)*dt,ylim_data(1),(x_axis_short(2,end)-x_axis_short(2,1))*dt,ylim_data(2)-ylim_data(1)],'FaceColor',rectangle_color, 'edgecolor','none');
rec3=rectangle('Position',[stim1_X(1)*dt-0.05,ylim_data(1),(stim1_X(2)-stim1_X(1))*dt+0.05,ylim_data(2)-ylim_data(1)],'faceColor',color_table(3,:),'edgecolor','none');

%plotting scale bar
% yline_start=ylim_data(1)-0.5; yline_end=yline_start+0.5; %10;
% % yline_start=yline_start./gain; yline_end=yline_end./gain;
% xline_start=x_axis_long(1,1)*dt-1; xline_end=xline_start+1;
% xline=[xline_start xline_start;xline_end xline_start]; yline=[yline_start yline_start; yline_start yline_end];
% stringh=[num2str(xline(2,1)-xline(1,1)), ' s'];
% stringv=[num2str(yline(2,2)-yline(2,1)), ' mV'];
% hline_h(1)=line(xline(:,1),yline(1,:),'color',[0 0 0],'linewidth',2);
% hline_v(2)=line(xline(:,2),yline(:,2),'color',[0 0 0],'linewidth',2);
% htext_h=text(xline(1,1)+(xline(2,1)-xline(1,1))/2,yline(1,1),stringh,'HorizontalAlignment', 'center','VerticalAlignment', 'top','fontsize',scalebar_fontsize); %'color',[0 0 0]
% % htext_v=text(xline(1,1),yline(2,1)+(yline(2,2)-yline(2,1))/2,stringv,'HorizontalAlignment', 'center','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',10);
% htext_v=text(xline(1,1),yline(2,1),stringv,'HorizontalAlignment', 'left','VerticalAlignment', 'bottom', 'rotation',90,'fontsize',scalebar_fontsize);

%plotting scale bar
horiz_vert=1;        lengthh=1;     textit=[num2str(lengthh), ' s'];  c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
 horiz_vert=0;        lengthh=0.5;     textit=[num2str(lengthh),' mV'];  c=[0,0,0];  fonsizes=scalebar_fontsize; perc1=[]; perc2=[];
        [p1,p2] = fn_makeCalibBar2(horiz_vert,lengthh,textit,c,fonsizes,perc1,perc2);
        
hold off
% set(gca,'ylim',[yline_start-0.2, ylim_data(2)],'xlim',[xline_start-0.2, xlim_data(2)]);
% ylabel('Vm [mV]', 'FontSize', 16);  xlabel('Time [sec]' ,'FontSize', 16);
set(gca,'color',[1 1 1],'xticklabel',[],'yticklabel',[],'xtick',[], 'ytick',[])
set(gca, 'visible', 'off') ;
 %% Plot power spectrum for channel data1
% orig_header=[]; header=[];
f1=figure;
%  orig_header=[1]; %[2,4:7]
 for count=1:length(orig_header)
 header(count)=find(headers==orig_header(count));
 
% start_time=[4,9]; %[7,12]; %[sec]

% Y_abs = []; f = [];

  for t=1:length(start_time);
spec_mat_noDC=[]; spec_mat = []; 

 spec_mat = data_mat_filt2(interval(:,t),:);
spec_mat_noDC=fn_Subtract_Mean(spec_mat);

[Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf,0,0);
 Y_abs_sem(:,t)=std(Y_abs(:,:,t),0,2)./sqrt(size(Y_abs(:,:,t),2));
 %calculating the critical t.
alpha=0.05;
df=size(Y_abs(:,:,t),2)-1;
alphaup = 1-alpha/2;
alphalow = alpha/2;
tupp = tinv(alphaup,df);
tlow = tinv(alphalow,df); 
 Y_abs_ci(:,t) = tupp.*Y_abs_sem(:,t);

  
Y_1to10(:,:,t)= Y_abs(f(:,t)<10 & f(:,t)>1,:,t);
Y_1to10_m(:,t)=mean(sum(Y_1to10(:,:,t),1));
Y_1to10_std(:,t)=std(sum(Y_1to10(:,:,t),1));
Y_10to20(:,:,t)= Y_abs(f(:,t)<20 & f(:,t)>10,:,t);
Y_10to20_m(:,t)=mean(sum(Y_10to20(:,:,t),1));
Y_10to20_std(:,t)=std(sum(Y_10to20(:,:,t),1));
Y_30to50(:,:,t)= Y_abs(f(:,t)<50 & f(:,t)>30,:,t);
Y_30to50_m(:,t)=mean(sum(Y_30to50(:,:,t),1));
Y_30to50_std(:,t)=std(sum(Y_30to50(:,:,t),1));

Y_1to10_norm(:,t)=sum(Y_1to10(:,:,t),1)./sum(Y_abs(:,:,t),1);
Y_1to10_norm_m(:,t)=mean(Y_1to10_norm(:,t));
Y_1to10_norm_std(:,t)=std(Y_1to10_norm(:,t));
Y_30to50_norm(:,t)=sum(Y_30to50(:,:,t),1)./sum(Y_abs(:,:,t),1);
Y_30to50_norm_m(:,t)=mean(Y_30to50_norm(:,t));
Y_30to50_norm_std(:,t)=std(Y_30to50_norm(:,t));
  end

  [h1,p1]=ttest(sum(Y_1to10(:,:,1),1),sum(Y_1to10(:,:,2),1));
  [h2,p2]=ttest(sum(Y_30to50(:,:,1),1),sum(Y_30to50(:,:,2),1));
  [h1n,p1n]=ttest(Y_1to10_norm(:,1),Y_1to10_norm(:,2));
  [h2n,p2n]=ttest(Y_30to50_norm(:,1),Y_30to50_norm(:,2));

% %% Plot power spectrum with shaded error bar
 for t=1:length(start_time);
   figure(f1)
%      ax(count)=subplot(2,1,count);
      hold on
   d1(t)= fn_shadedErrorBar(f(f(:,t)<49,t),mean(Y_abs(f(:,t)<49,:,t),2),Y_abs_sem(f(:,t)<49,t),{'LineWidth',1,'color', color_table(t,:)});
%    hpatch = findobj(gcf, 'type', 'patch'); set(hpatch(1),'facecolor',[0 0 0],'edgecolor','none','faceAlpha', 0.3)
%     xlim([0 20]); ylim([-0.1 0.9])    
        hold off
 end
 xlim([1 100]); ylim([0 0.1])
        set(gca,'xscale','log');
        set(gca,'yscale','log');
        set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',1, 'xtick',[1, 10,100],'ytick', [10e-7,10e-6, 10e-5, 10e-4, 10e-3],'yminortick','on') 
         xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
        ylabel('PSD [mV^2/Hz]','fontsize',20, 'fontname', 'arial'); %Power spectral density
l=legend([d1(1).mainLine d1(2).mainLine ],{'NB-','NB+'},'Location','northeast', 'box', 'off');
l.LineWidth=1.5;
l.FontSize=10;
 %% Plot power spectrum without shaded error bar

%      figure(f2)
% %      ax(count)=subplot(2,1,count);
% for t=1:length(start_time);
%       hold on
%         plot(f(:,t),mean(Y_abs(:,:,t),2),'color', color_table(t,:),'LineWidth',2) 
% %         plot(f(:,t),mean(Y_abs(:,:,t),2)./sum(mean(Y_abs(:,:,t),2)),'color', color_table(t,:)) %normalized by the total power
%       hold off
% %         ylim([-0.01 1]);
% %         xlim([0 20]); 
% %         title('')
% xlim([0 1000]); ylim([-0.01 10])
%         set(gca,'xscale','log');
%         set(gca,'yscale','log');
%         set(gca,'fontsize',20, 'fontname', 'arial', 'box', 'off','linewidth',2)
%         xlabel('Frequency [Hz]','fontsize',20, 'fontname', 'arial')
%         ylabel('Power spectral density [mV^2/Hz]','fontsize',20, 'fontname', 'arial')
%  end
 
 end
  
  %% Plot power spectrum for channel EEG
% orig_header=[]; header=[]; sf_EEG=sf_galvano;
% figure
%  orig_header=[3:12]; %[2,4:7]
%  for count=1:length(orig_header)
%  header(count)=find(headers==orig_header(count));
%  
% start_time=[7,12]; %[sec]
% duration = 2; %[sec]
% x_value = 1;
% Y_abs = []; f = [];
% 
%   for t=1:length(start_time);
% DC=[]; spec_mat_noDC=[]; spec_mat = []; spec_mat_mean = []; start_sample = []; end_sample = []; interval = []; 
% 
% start_sample = start_time(t).*sf_EEG;
% if start_time(t)==0
%     start_sample = 1;
% end
% end_sample = start_sample+duration.*sf_EEG-1;
% interval = start_sample:end_sample;
% spec_mat = data(header(count)).x_value(1).EEG(interval,:); 
% DC= mean(spec_mat,1);
% l=size(spec_mat,1);
% l=l/2-1;
% DC=wextend('addrow', 'per', DC, l);
% spec_mat_noDC = spec_mat-DC;
% spec_mat_mean = mean(spec_mat,2);
% 
% [Y_abs(:,:,t),f(:,t)] = fn_Amp_Spectrum(spec_mat_noDC,sf);
%   end
%   
%   
%  for t=1:length(start_time);
%      figure(gcf)
%      ax(count)=subplot(5,2,count);
%       hold on
%     fn_shadedErrorBar(f(:,t),mean(Y_abs(:,:,t),2),std(Y_abs(:,:,t),0,2),{'LineWidth',1,'color', color_table(t,:)},1);
%     xlim([0 35]); %ylim([-0.05 0.05])
%  end
%  end

%% 
if save_flag==1; 
    %for different depths:
filename_PSD = ['2016-10-20-001_LFP_PSD_1-100Hz_h',num2str(orig_header),'_',num2str(depth),'um'];
filename_traces = ['2016-10-20-001_LFP_traces_h',num2str(orig_header),'_t',num2str(trace_ind),'_',num2str(depth),'um'];
% filename_PSD = ['2016-10-20-001_LFP_PSD_1-100Hz_h',num2str(orig_header),'_',num2str(depth(header)),'um'];
% filename_traces = ['20150323-002_LFP_traces_h',num2str(orig_header),'_t',num2str(trace_ind),'_',num2str(depth(header)),'um'];
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\LFP depth';
%for Blockers:
filename_PSD = ['2015-31-12-002_LFP_PSD_1-100Hz_h',num2str(orig_header)];
filename_traces = ['2015-31-12-002_LFP_traces_h',num2str(orig_header),'_t',num2str(trace_ind)];
% cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\LFP+Blockers';
cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\LFP+Blockers\GABA Blockers';


savefig(f1,filename_PSD)
print(f1,filename_PSD,'-dpng','-r600','-opengl')  
savefig(f2,filename_traces)
print(f2,filename_traces,'-dpng','-r600','-opengl')  
% print (1, '-depsc2', filename)
end