function []=fn_plot_sensory_stim(dt, stim2,whisker_stim_color)
            figure(gcf)
            patch_xdata=[stim2; flipud(stim2)];
            xlim_data=get(gca,'xlim');
            ylim_data=get(gca,'ylim');
            a=ones(size(patch_xdata));
            yex=wextend('1D','sym',ylim_data,1)';
            patch_ydata=bsxfun(@times,a,yex);
            patch_cdata=ones(size(patch_xdata));
            patch(patch_xdata.*dt*1000,patch_ydata,patch_cdata,'faceColor',whisker_stim_color(1,:),'edgecolor','none'); %'faceAlpha', 0.3
            set(gca,'children',flipud(get(gca,'children'))) %changes the order of objects plotted so that the bars will appear under the traces.
end