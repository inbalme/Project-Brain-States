function []=fn_plot_sensory_stim(dt, stim2_X,whisker_stim_color)
            figure(gcf)
            patch_xdata=[stim2_X{2}; flipud(stim2_X{2})];
            xlim_data=get(gca,'xlim');
            ylim_data=get(gca,'ylim');
            yex=wextend('1D','sym',ylim_data,1)';
            l=(size(patch_xdata,2)-1)/2;
            patch_ydata=wextend('ac','sym',yex,l);
            patch_cdata=ones(size(patch_xdata));
            patch(patch_xdata.*dt*1000,patch_ydata,patch_cdata,'faceColor',whisker_stim_color(1,:),'edgecolor','none'); %'faceAlpha', 0.3
            set(gca,'children',flipud(get(gca,'children'))) %changes the order of objects plotted so that the bars will appear under the traces.
end