cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Vm-LFP correlations\LFP_50Hz+BP1-200 Vm_50Hz';
cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Vm-LFP correlations\LFP_50Hz+BP1-200 Vm_50Hz';
figure
hold on
scatter(cc_stat.spont.lag0(:,1),cc_stat.evoked.lag0(:,1),40,'k','fill')
scatter(cc_stat.spont.lag0(:,2),cc_stat.evoked.lag0(:,2),40,'r','fill')
line([0,0.5],[0,0.5])
hold off
[r1,m1,b1] = regression(cc_stat.spont.lag0(:,1)',cc_stat.evoked.lag0(:,1)');
[r2,m2,b2] = regression(cc_stat.spont.lag0(:,2)',cc_stat.evoked.lag0(:,2)');