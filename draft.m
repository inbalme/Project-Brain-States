%% Plot single traces of LFP
x_value=1; gain=20; figure
for header=1:10;
 plot_data = data(header).x_value(1).data1(:,:); 
 plot_data(9.95*sf:10.55*sf,:) = nan; %plot_data(9.95*sf,:); %putting a straight line instead of the ES artifact
 
  trace_ind = [1,3,6]; %[1:size(plot_data,2)] ;  %trace_ind is the numbers of traces to be plotted
       l=size(plot_data(:,1),1);
        l=l/2-1;
        DC=10.*(0:length(trace_ind)-1);
        DC=wextend('addrow', 'per', DC, l); 
        trace_to_plot = (plot_data(:,trace_ind)+DC)./gain ;
      
        figure(gcf)
        subplot(5,2,header)
        for i=1:length(trace_ind)         
            hold on    
                plot((1:size(trace_to_plot,1))'.*dt,trace_to_plot(:,i),'k', 'linewidth',1)
            hold off
%             xlabel('Time [sec]' ,'FontSize', 14);
%             ylabel('Vm [mV]', 'FontSize', 14);
            set(gca,'ylim', [-1 2])
            title('Depth','FontSize', 14)
        end
end