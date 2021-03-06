%%
close all
clear all
plot_amp_histogram_flag=1;
save_flag=0;
if plot_amp_histogram_flag==1;
    binsize=0.5; %[mV];
    exp_type=1;
    switch exp_type
        case 1
            fileind=12;
            legend_string={'NB-','NB+'};
            cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\single trial analysis\Spontaneous'
            load('Spontaneous activity.mat')
            color_table=[0 0 0; [216 22 22]/256;  [136 137 138]/256; [255 153 153]/256; [30,75,14]/256; [112,172,90]/256];  
        case 2
            fileind=4;
            legend_string={'Light Off','Light On'};
            cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\single trial analysis\Spontaneous'
            load('Spontaneous activity.mat')
            color_table=[0 0 0; [0 0 204]/256;  [136 137 138]/256; [102, 172,255]./256; [30,75,14]/256; [112,172,90]/256];  
    end

        for t=1:2;
            
              data_vec{t}=event_ongoing(fileind).amplitude{t}(:);  
              nbin(t)=ceil(range(data_vec{t}))./binsize; 
% estimating median ci:
%the standard error is given by: 1.2533*(sigma/sqrt(N))
alpha=0.05;
df=size(data_vec{t}(:,:),1)-1;
alphaup = 1-alpha/2;
alphalow = alpha/2;
tupp = tinv(alphaup,df);
tlow = tinv(alphalow,df); 
ci_median{t}=1.255*tupp.*std(data_vec{t}(:,:),0,1)./sqrt(size(data_vec{t}(:,:),1));
              %bootstrap CI for the median:  
                    dataset=[];
                    dataset=data_vec{t}(:)';
                    k= size(dataset,2);     
                    lower_bound=2.5;
                    upper_bound=97.5;
                    iterations=5000;
                     [prcntile1(t,:), prcntile2(t,:)]=fn_get_CI_w_bootstrap(dataset,0,5000,2.5,97.5,'median');
                     data_vec_median(1,t) = median(data_vec{t}(:)); 
        end
        %for population histogram: take even number of samples from each condition
%             samples=min(length(event_ongoing(fileind).amplitude{1}(:)),length(event_ongoing(fileind).amplitude{2}(:)));
%            data_vec_pop{t}=event_ongoing(fileind).amplitude{t}(1:samples)
%             data_vec_norm{t}=data_vec_pop{t}./max(max(data_vec_pop{1})); %for a population histogram
%              nbin = 30;
 
    Fig1=figure;
            hold on
            [ncounts{1}(1,:),nbins{1}(1,:)]=hist(data_vec{1}(:),nbin(1));
            ncounts{1}(1,:)=ncounts{1}(1,:)./length(data_vec{1}(:));
            [ncounts{2}(1,:),nbins{2}(1,:)]=hist(data_vec{2}(:),nbins{1}(1,:));
            ncounts{2}(1,:)=ncounts{2}(1,:)./length(data_vec{2}(:));
            h1=bar(nbins{1}(1,:),ncounts{1}(1,:));
            h2=bar(nbins{2}(1,:),ncounts{2}(1,:));
            
%             obj = gmdistribution.fit(data_vec(:,1),2);
            set(h2,'FaceColor',color_table(2,:),'EdgeColor','w','faceAlpha', 0.3)
            set(h1,'FaceColor',color_table(1,:),'EdgeColor','w','faceAlpha', 0.7)
            ylim1=get(gca,'ylim');       
            hold off
%             uistack(h1,'top'); uistack(h2,'top'); 
            axis tight
            xlim1=get(gca,'xlim');
            xrim=diff(xlim).*0.1;
            xlim1(1)=xlim1(1)-xrim;
            xlim1(2)=xlim1(2)+xrim;
            set(gca,'xlim',xlim1,'FontSize', 20,'fontname', 'arial')
            xlabel('PSP Amplitude [mV]', 'FontSize', 20,'fontname', 'arial');
            ylabel('Probability', 'FontSize', 20,'fontname', 'arial')

 [l,OBJH,OUTH,OUTM] = legend([h1 h2],legend_string, 'position',[0.8,0.9,0.2,0.1],'fontsize',10, 'box', 'off'); % returns a handle LEGH to the legend axes; a vector OBJH containing handles for the text, lines, and patches in the legend; a vector OUTH of handles to thelines and patches in the plot; and a cell array OUTM containingthe text in the legend.
lpatch=findobj(OBJH,'type','patch');
lpatch(1).Vertices(2:3,2)=lpatch(1).Vertices(1,2)+0.15;
lpatch(2).Vertices(2:3,2)=lpatch(2).Vertices(1,2)+0.15;
lpatch(1).Vertices(:,2)=lpatch(1).Vertices(:,2)+0.06;
lpatch(2).Vertices(:,2)=lpatch(2).Vertices(:,2)+0.08;
lpatch(1).FaceColor=color_table(1,:);
lpatch(1).FaceAlpha=0.7;
lpatch(2).FaceColor=color_table(2,:);
lpatch(2).FaceAlpha=0.3;


    Fig2=figure;
        hold on
            e1=errorbar(1,data_vec_median(1,1),data_vec_median(1,1)-prcntile1(1,1),prcntile2(1,1)-data_vec_median(1,1),'color',color_table(3,:)); %
            e2=errorbar(2,data_vec_median(1,2),data_vec_median(1,2)-prcntile1(2,1),prcntile2(2,1)-data_vec_median(1,2),'color',color_table(4,:)); % 
            set(e1,'Marker','o','MarkerSize',10,'MarkerFaceColor',color_table(1,:),'MarkerEdgeColor',color_table(1,:),'linewidth',2)
            set(e2,'Marker','o','MarkerSize',10,'MarkerFaceColor',color_table(2,:),'MarkerEdgeColor',color_table(2,:),'linewidth',2)
            set(gca,'xtick',[1,2],'xticklabel',legend_string);
            ylabel('Amp. [mV]')

        hold off
        
        Fig3=figure;
        hold on
            e1=errorbar(1,data_vec_median(1,1),ci_median{1},ci_median{1},'color',color_table(3,:)); %
            e2=errorbar(2,data_vec_median(1,2),ci_median{2},ci_median{2},'color',color_table(4,:)); % 
            set(e1,'Marker','o','MarkerSize',10,'MarkerFaceColor',color_table(1,:),'MarkerEdgeColor',color_table(1,:),'linewidth',2)
            set(e2,'Marker','o','MarkerSize',10,'MarkerFaceColor',color_table(2,:),'MarkerEdgeColor',color_table(2,:),'linewidth',2)
            set(gca,'xtick',[1,2],'xticklabel',legend_string);
            ylabel('Amp. (mV)')

        hold off
%             pause
       if save_flag==1;     
            saveas(Fig1,'Event_Amp_Hist','fig') 
            print(Fig1,'Event_Amp_Hist','-dpng','-r600','-opengl')
            saveas(Fig2,'Event_Amp_Median_ci','fig') 
            print(Fig2,'Event_Amp_Median_ci','-dpng','-r600','-opengl')
       end   
end