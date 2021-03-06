%% NBES Rin

clear all
close all

global dt sf dt_galvano sf_galvano data data_no_spikes files Param raw_data current_data Ch2_data stim2_X stim1_X 
 global exp_type
 channel = 1;    % V1 - 1, I1 - 2, V2 - 3, I2 - 4
exp_type=3; %1-NBES, 2-ChAT
% trace_type_input=[3,2]; %1:3
save_flag= 0;
print_flag=0;
norm_flag=0;
no_DC_flag=0;
BP50HzLFP_flag=0; %removing 50Hz noise from LFP signal
BP50HzVm_flag=0; %removing 50Hz noise from Vm signal
BPLFP_flag=0; %filtering LFP. the default filter is the one used to filter LFP in the multiclamp
bp_manual_LFP=[0.1,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)
BPVm_flag=0; %filtering LFP and Vm same as LFP was filtered in the multiclamp
bp_manual_Vm=[0,300]; %if bp_manual=[] the default would be to take bp_filt from Param (the filter used for LFP in the multiclamp)

switch exp_type
    case 1 
        path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Rin';   
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        cd(path_output)       
    case 2
        path_output='D:\Inbal M.Sc\Data PhD\ChAT Data\Figures\Rin';        
        a = exist(path_output,'dir'); %a will be 7 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        cd(path_output)
    case 3 
        path_output='D:\Inbal M.Sc\Data PhD\NB-ES Data\Figures\Rin_VC';   
        a = exist(path_output,'dir'); %a will be 1 if a folder "name" exists and 0 otherwise
        if a==0;
            mkdir(path_output);
        end
        cd(path_output)      
end

switch exp_type
    case 1
        files_to_analyze =[62,72,75,82,84]; %[62,72,75,82,84]; 
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB-', 'NB+'};       

    case 2
        files_to_analyze = [74,76,80,82,84,87,91,93,113,116,117]; %[78,81,83,79];
        cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
        load ChAT_Files_v3
        legend_string={'Light Off', 'Light On'};
        
    case 3
        files_to_analyze =[31,38,42,51,69,71,74];
        cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
        load NBES_Files_v2
        legend_string={'NB+', 'NB-'};    
end
        
% for fileind=1;
    for fileind=1:length(files_to_analyze) ;
%     close all
    fname = files(files_to_analyze(fileind)).extracted_name;        
    path = files(files_to_analyze(fileind)).extracted_path;
    cd(path)
    load(fname) 
    %%
        Ch2_data= raw_data{3}./20; %dividing by the LFP gain           
        if isempty(data_no_spikes)
            current_data=raw_data{channel};
            data_used='raw_data';
        else
        current_data=data_no_spikes{channel}; %raw_data{channel}; 
        data_used='data_no_spikes';
        end
        galvano_nstim = Param.facade(6);
        galvano_freq = Param.facade(7);
        
   data_preprocessing
   
    if ~isempty(current_data_filt)
     current_data=current_data_filt;
    end

%% 
if exp_type~=3
x_value=2:3;
        curr_inj_1=Param.facade(20); %in sec
        curr_inj_dur=Param.facade(22); %in sec
        curr_inj_2=Param.facade(20)+curr_inj_dur+Param.facade(23); %in sec        
        curr_inj_1_end=curr_inj_1+curr_inj_dur;
        curr_inj_2_end=curr_inj_2+curr_inj_dur;
        Iinj=Param.facade(21)/2.5; %[nA]
        
curr_inj_interval_before=[round((curr_inj_1+0.025)*sf{1}):round((curr_inj_1_end-0.025)*sf{1})];
curr_inj_interval_after=[round((curr_inj_2+0.025)*sf{1}):round((curr_inj_2_end-0.025)*sf{1})];
baseline_before=round((curr_inj_1-5/1000)*sf{1}:curr_inj_1*sf{1});
baseline_after=round((curr_inj_2-5/1000)*sf{1}:curr_inj_2*sf{1});
Vm_x2_curr_inj(1,1)=mean(mean(current_data(curr_inj_interval_before,:,2)));
Vm_x2_curr_inj(1,2)=mean(mean(current_data(curr_inj_interval_after,:,2)));
Vm_x2_curr_inj_amp(1,1)=Vm_x2_curr_inj(1,1)-mean(mean(current_data(baseline_before,:,2)));
Vm_x2_curr_inj_amp(1,2)=Vm_x2_curr_inj(1,2)-mean(mean(current_data(baseline_after,:,2)));
Rin_x2_curr_inj(1,1)=Vm_x2_curr_inj_amp(1,1)./Iinj;
Rin_x2_curr_inj(1,2)=Vm_x2_curr_inj_amp(1,2)./Iinj;
Vm_x3_curr_inj(1,1)=mean(mean(current_data(curr_inj_interval_before,:,3)));
Vm_x3_curr_inj(1,2)=mean(mean(current_data(curr_inj_interval_after,:,3)));
Vm_x3_curr_inj_amp(1,1)=Vm_x3_curr_inj(1,1)-mean(mean(current_data(baseline_before,:,3)));
Vm_x3_curr_inj_amp(1,2)=Vm_x3_curr_inj(1,2)-mean(mean(current_data(baseline_after,:,3)));
Rin_x3_curr_inj(1,1)=Vm_x3_curr_inj_amp(1,1)./Iinj;
Rin_x3_curr_inj(1,2)=Vm_x3_curr_inj_amp(1,2)./Iinj;
Rin(fileind,:)=[Rin_x2_curr_inj(1,1),Rin_x2_curr_inj(1,2),Rin_x3_curr_inj(1,1),Rin_x3_curr_inj(1,2)];
%here I didn't subtract the value of Vm at beginning of current inj. from
%the value it reached during the curr. inj.
%Need also to get the Rin from each trace and test for each cell if it is
%different with ttest.

voltages1=mean(raw_data{1}(curr_inj_1*sf{1}:curr_inj_1_end*sf{1},:,3),2)';
voltages2=mean(raw_data{1}(curr_inj_2*sf{1}:curr_inj_2_end*sf{1},:,3),2)';
before_ES  = fn_CalculateCapacitResist(Iinj, voltages1, sf{1}/1000, 1,Vm_x3_curr_inj(1,1));
after_ES  = fn_CalculateCapacitResist(Iinj, voltages2, sf{1}/1000, 1,Vm_x3_curr_inj(1,2));
% for trace=1:7; %this is a mistake. need to put the mean trace in the function
% curr_inj_del=[Param.facade(20), Param.facade(20)+Param.facade(22)+Param.facade(23)];
% curr_inj_duration=Param.facade(22);
% voltages1=raw_data{1}(curr_inj_del(1)*sf{1}:(curr_inj_del(1)+curr_inj_duration)*sf{1},trace,3)';
% voltages2=raw_data{1}(curr_inj_del(2)*sf{1}:(curr_inj_del(2)+curr_inj_duration)*sf{1},trace,3)';
% before_ES(trace,:)  = fn_CalculateCapacitResist(Iinj, voltages1, sf{1}/1000, 1,Vm_x3_curr_inj(1,1));
% after_ES(trace,:)  = fn_CalculateCapacitResist(Iinj, voltages2, sf{1}/1000, 1,Vm_x3_curr_inj(1,2));
%     end
R_electrode(fileind,:)=[before_ES(:,1), after_ES(:,1)];
C_electrode(fileind,:)=[before_ES(:,2), after_ES(:,2)];
R_cell(fileind,:)=[before_ES(:,3), after_ES(:,3)];
C_cell(fileind,:)=[before_ES(:,4), after_ES(:,4)];
Error1(fileind,:)=[before_ES(:,5), after_ES(:,5)];
Error2(fileind,:)=[before_ES(:,6), after_ES(:,6)];
cell_props.fileind(fileind).fname=files_to_analyze(fileind);
% cell_props(fileind).trace=trace;
cell_props.fileind(fileind).R_electrode=[before_ES(:,1), after_ES(:,1)];
cell_props.fileind(fileind).C_electrode=[before_ES(:,2), after_ES(:,2)];
cell_props.fileind(fileind).R_cell=[before_ES(:,3), after_ES(:,3)];
% cell_props(fileind).R_cell_m=nanmean(cell_props(fileind).R_cell);
% cell_props(fileind).R_cell_std=nanstd(cell_props(fileind).R_cell);
cell_props.fileind(fileind).C_cell=[before_ES(:,4), after_ES(:,4)];
cell_props.fileind(fileind).Error1=[before_ES(:,5), after_ES(:,5)];
cell_props.fileind(fileind).Error2=[before_ES(:,6), after_ES(:,6)];
cell_props.fileind(fileind).R_total=[before_ES(:,1)+before_ES(:,3), after_ES(:,1)+after_ES(:,3)];
cell_props.fileind(fileind).R_total_NoExpFit=Rin(fileind,:);
end

if exp_type==3
    x_value=2;
        curr_inj_1=1.2; %in sec
        curr_inj_dur=0.5; %in sec
        curr_inj_2=6.2; %in sec        
        curr_inj_1_end=curr_inj_1+curr_inj_dur;
        curr_inj_2_end=curr_inj_2+curr_inj_dur;

curr_inj_interval_before=[round((curr_inj_1)*sf{1}):round((curr_inj_1_end)*sf{1})];
curr_inj_interval_after=[round((curr_inj_2)*sf{1}):round((curr_inj_2_end)*sf{1})];
baseline_before=round(1:0.1*sf{1});
baseline_after=round(4.6*sf{1}:4.7*sf{1});
Im_x2_curr_inj(1,1)=mean(mean(current_data(curr_inj_interval_before,:,2)));
Im_x2_curr_inj(1,2)=mean(mean(current_data(curr_inj_interval_after,:,2)));
Im_x2_curr_inj_amp(1,1)=Im_x2_curr_inj(1,1)-mean(mean(current_data(baseline_before,:,2)));
if files_to_analyze(fileind)==74;
    Im_x2_curr_inj_amp(1,1)=Im_x2_curr_inj(1,1)-0; %the cell had a lot of spontaneous activity but in fact he baseline is 0
end
Im_x2_curr_inj_amp(1,2)=Im_x2_curr_inj(1,2)-mean(mean(current_data(baseline_after,:,2)));
Im_x2_curr_inj_amp=Im_x2_curr_inj_amp./100; %turn from pA to nA
% Vc=10.*mean(mean(raw_data{2}(round(curr_inj_interval_before./5),:,2)));
Vc=100.*Param.facade(2,1);
Rin_x2_curr_inj(1,1)=(-1).*Vc./Im_x2_curr_inj_amp(1,1);
Rin_x2_curr_inj(1,2)=(-1).*Vc./Im_x2_curr_inj_amp(1,2);
% Vm_x3_curr_inj(1,1)=mean(mean(current_data(curr_inj_interval_before,:,3)));
% Vm_x3_curr_inj(1,2)=mean(mean(current_data(curr_inj_interval_after,:,3)));
% Vm_x3_curr_inj_amp(1,1)=Vm_x3_curr_inj(1,1)-mean(mean(current_data(baseline_before,:,3)));
% Vm_x3_curr_inj_amp(1,2)=Vm_x3_curr_inj(1,2)-mean(mean(current_data(baseline_after,:,3)));
% Rin_x3_curr_inj(1,1)=Vm_x3_curr_inj_amp(1,1)./Iinj;
% Rin_x3_curr_inj(1,2)=Vm_x3_curr_inj_amp(1,2)./Iinj;
Rin(fileind,:)=[Rin_x2_curr_inj(1,1),Rin_x2_curr_inj(1,2)];
end
    end
    if exp_type~=3
        cell_props.R_electrode=R_electrode;
        cell_props.C_electrode=C_electrode;
        cell_props.R_cell=R_cell;
        cell_props.C_cell=C_cell;
        cell_props.Error1=Error1;
        cell_props.Error2=Error2;
        cell_props.R_total=R_electrode+R_cell;
        cell_props.R_total_NoExpFit=Rin;
        
        [cell_props.C_cell_lillietest_h,cell_props.C_cell_lillietest_p]=lillietest(cell_props.C_cell(:,2)-cell_props.C_cell(:,1));
        [cell_props.C_cell_ttest_h,cell_props.C_cell_ttest_p]=ttest(cell_props.C_cell(:,1),cell_props.C_cell(:,2));
        [cell_props.C_cell_wilcoxon_p,cell_props.C_cell_wilcoxon_h]=signrank(cell_props.C_cell(:,1),cell_props.C_cell(:,2));
    end

    if exp_type==3
        cell_props.R_cell=Rin;
    end
    
[cell_props.R_cell_lillietest_h,cell_props.R_cell_lillietest_p]=lillietest(cell_props.R_cell(:,2)-cell_props.R_cell(:,1));
[cell_props.R_cell_ttest_h,cell_props.R_cell_ttest_p]=ttest(cell_props.R_cell(:,1),cell_props.R_cell(:,2));
[cell_props.R_cell_wilcoxon_p,cell_props.R_cell_wilcoxon_h]=signrank(cell_props.R_cell(:,1),cell_props.R_cell(:,2));


%%
if any(cell_props.R_cell(:,1)>1000000)
    tmp_Y= [cell_props.R_cell./1000000]';
else
    tmp_Y= cell_props.R_cell';
end
tmp_X(1,:)=ones(1,size(tmp_Y,2));
tmp_X(2,:)=2*ones(1,size(tmp_Y,2));
E = std(tmp_Y');
if cell_props.R_cell_wilcoxon_p >0.05 
    asterisk_sp='n.s.'; a1_fontsize=13;
else if cell_props.R_cell_wilcoxon_p<0.05 && cell_props.R_cell_wilcoxon_p>0.01
    asterisk_sp='*'; a1_fontsize=17;
    else if cell_props.R_cell_wilcoxon_p<0.01 && cell_props.R_cell_wilcoxon_p>0.001
            asterisk_sp='**'; a1_fontsize=17;
    else if cell_props.R_cell_wilcoxon_p<0.001
             asterisk_sp='***'; a1_fontsize=17;
        end
        end
    end
end
linex=[1;2];
my=max(max(tmp_Y)).*1.05; 
liney=[my ;my ];

g1=figure;
hold on
line(tmp_X,tmp_Y,'color',[0.7 0.7 0.7],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
errorbar(tmp_X(:,1), nanmean(tmp_Y,2),E,'k','linewidth',2.5,'markersize',10,'markerfacecolor','k')
% line(linex,liney,'color',[0 0 0],'linewidth',1.5,'markersize',10,'markerfacecolor','k')
text(1.5,my,asterisk_sp,'HorizontalAlignment', 'center','verticalAlignment','bottom','fontsize',a1_fontsize)
hold off
        x1limits = [0.75 2.25];
        x1ticks = [1,2];
        y1limits=[0 150]; %get(gca,'ylim');
%         y1limits = [0 1.1];
%         y1ticks = [0,0.5,1];
        set( gca, 'xlim', x1limits,'ylim',[0 150], 'xtick', x1ticks,'fontsize',28,'linewidth',1,...
        'ticklength', [0.010 0.010],'fontname', 'arial','xticklabel',legend_string ,'box', 'off'); %'fontweight', 'bold', 
        ylabel('Rin (M\Omega)', 'FontSize', 28,'fontname', 'arial');
        title(['Rin,  p=' num2str(cell_props.R_cell_wilcoxon_p)] ,'FontSize', 20,'fontname', 'arial');   
        %%
    if save_flag==1
        cd(path_output)
        print(g1,'Rin','-dpng','-r600','-opengl')
        saveas(g1,'Rin','fig') 
        save('Rin', 'files_to_analyze', 'cell_props')
    end
