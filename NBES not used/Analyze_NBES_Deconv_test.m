%% Deconvolution of the real data          
%deconvolution of sig_1 and sig_2, convolving the kernel with sig_1 to create sig_2_predict.  
% clear all
clearvars -except Ch2_data_filt data_no_spikes_filt interval
fileind=1; t=1; trials=1; kernelsize=1002; w=0.05; interval(:,1)=[4000:28999]; interval(:,2)=[54000:78999]; 
data_type=6;
randcell = randsample(1:10,trials,true);
randLFP = randsample(1:10,trials,true);
randcelltrace=randsample(1:5,trials,true);
randLFPtrace=randsample(1:5,trials,true);
for trace=1:trials;
switch data_type
    case 1 %rand data
        data_LFP{t}=rand(25000,1);
        data_Vm_filt{t}=rand(25000,1);

    case 2 %wgn data
         data_LFP{t}= wgn(25000,1,1); %rand(25000,1);
        data_Vm_filt{t}=wgn(25000,1,1); %rand(25000,1);


    case 3 %cell data
        data_LFP{t}=Ch2_data_filt{randLFP(trace)}(interval(:,t),randLFPtrace(trace),1);
        data_Vm_filt{t}=data_no_spikes_filt{randcell(trace)}(interval(:,t),randcelltrace(trace),1);  
        
    case 4       
     data_LFP{t}= wgn(25000,1,1); %rand(25000,1);   
     data_LFP{t}=sgolayfilt(data_LFP{t}, 1, 55);
     data_Vm_filt{t}=data_no_spikes_filt{randcell(trace)}(interval(:,t),randcelltrace(trace),1);  
     
     case 5
      data_Vm_filt{t}=wgn(25000,1,1); %rand(25000,1); 
      data_LFP{t}=Ch2_data_filt{3}(interval(:,t),1,1);
       case 6 %cell data
        data_LFP{t}=data_no_spikes_filt{randLFP(trace)}(interval(:,t),randLFPtrace(trace),1);
        data_Vm_filt{t}=data_no_spikes_filt{randcell(trace)}(interval(:,t),randcelltrace(trace),1);  
end
%  data_LFP{t}=sgolayfilt(data_LFP{t}, 1, 55);
%   data_Vm_filt{t}=sgolayfilt(data_Vm_filt{t}, 1, 55);
  
 xcorr_L=size(data_LFP{t},1)-1; %0.5*sf{1}; %lag of 0.5sec
   CI_std=2;
  
%   % normalizing:
%                     sig1{fileind,t}=(data_LFP{t}-mean(data_LFP{t}))./std(data_LFP{t});
%                     sig2{fileind,t}=(data_Vm_filt{t}-mean(data_Vm_filt{t}))./std(data_Vm_filt{t});

sig1{fileind,t}(:,trace)=data_LFP{t};
sig2{fileind,t}(:,trace)=data_Vm_filt{t};
%                 for trace = 1:trials;
                    G{fileind,t}(:,trace) = Ldeconvs(sig1{fileind,t}(:,trace),sig2{fileind,t}(:,trace),kernelsize,w,1);
%                 end
                      
                        kernelc{fileind,t}(:,trace) = G{fileind,t}(:,trace);
                    
                                sig2_predict_tmp(:,trace)=conv(kernelc{fileind,t}(:,trace),sig1{fileind,t}(:,trace))+mean(sig2{fileind,t}(:,trace)); %+mean(sig2{fileind,t}(:,trace))
                                sig2_predict{fileind,t}(:,trace)=sig2_predict_tmp(kernelsize/2:end-kernelsize/2,trace);
                                [cc_sig2_predict{fileind,t}(:,trace), cc_sig2_predict_lag{fileind,t}(:,trace), cc_sig2_predict_ci{fileind,t}(:,trace)] =...
                                    crosscorr(sig2_predict{fileind,t}(:,trace), sig2{fileind,t}(:,trace),xcorr_L);      %,xcorr_L,CI_std
                                [r, p] = corr(sig2_predict{fileind,t}(:,trace),sig2{fileind,t}(:,trace));
                                pearson_r_temp(trace,1)=r;
                                pearson_p_temp(trace,1)=p;
                                R_temp(trace,1) =  pearson_r_temp(trace,1)^2;
                                mse(trace,1) = fn_MSE(sig2_predict{fileind,t}(:,trace), sig2{fileind,t}(:,trace),1);
                                 [cc_sig2{fileind,t}(:,trace), cc_sig2_lag{fileind,t}(:,trace), cc_sig2_ci{fileind,t}(:,trace)] =...
                                     crosscorr(sig2{fileind,t}(:,trace), sig2{fileind,t}(:,trace),xcorr_L); %,xcorr_L,CI_std
                                 [cc_sig1sig2{fileind,t}(:,trace), cc_sig1sig2_lag{fileind,t}(:,trace), cc_sig1sig2_ci{fileind,t}(:,trace)] =...
                                     crosscorr(sig1{fileind,t}(:,trace), sig2{fileind,t}(:,trace),xcorr_L); %,xcorr_L,CI_std
                                 [cc_sig1sig2_predict{fileind,t}(:,trace), cc_sig1sig2_predict_lag{fileind,t}(:,trace), cc_sig1sig2_predict_ci{fileind,t}(:,trace)] =...
                                     crosscorr(sig1{fileind,t}(:,trace), sig2_predict{fileind,t}(:,trace),xcorr_L); %,xcorr_L,CI_std
end
                 
              %     Vm-Vm predicted - actual data
            color_table=[0 0 0; [216 22 22]/256; [136 137 138]/256; [255 153 153]/256]; %NB+ in red     
            figure
                        hold on     
                            plot(sig2{fileind,t}(:,:),'LineWidth',1.5,'color', color_table(t+2,:)); 
%                             plot(x_axis,sig1{fileind,t}(:,trace)-60,'LineWidth',2,'color',[0 0 1]); 
                            plot(sig2_predict{fileind,t}(:,:),'LineWidth',1.5,'color', color_table(t,:));   
%                              title(['f', num2str(files_to_analyze(fileind)),'Vm,Vm predicted - actual data']); 
%                             text(sig2_predict{fileind,t}(1,trace),[num2str(floor(sig2_predict{fileind,t}(1,trace))), ' mV '],'HorizontalAlignment','right','fontsize',13,'fontname','arial')
                             axis tight
                              ylim_data=[get(gca,'ylim')]';
                              xlim_data=[get(gca,'xlim')]';

figure
kernelc_m=mean(kernelc{fileind,t}(:,:),2);
hold on
plot(kernelc{fileind,t}(:,:),'LineWidth',1,'color', color_table(t+2,:))
plot(kernelc_m,'LineWidth',2)
hold off

figure
hold on
 plot(mean(cc_sig2{fileind,t}(:,:),2),'LineWidth',2,'color', color_table(t+2,:));
 plot(mean(cc_sig2_predict{fileind,t}(:,:),2),'LineWidth',2,'color', color_table(t,:));
 hold off