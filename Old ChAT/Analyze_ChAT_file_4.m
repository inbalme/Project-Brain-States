%% Analyze_ChAT
%calls the functions: 

clear all
close all

cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
load ChAT_Files

for fileind = 4;
    
    fname = files(fileind).extracted_name;        
    path = files(fileind).extracted_path;
    cd(path)
    load(fname) 
    
    cd  'D:\Inbal M.Sc\MATLAB\Project Brain States';
    
        for header = 1:length(files(1, fileind).headers);
            stim_type = Param.header(header).stim.type;
            
            %%  Airpuff/Galvano protocol stimulus:
                sf_Vm = Param.header(header).stim.sf; %[1/sec]
                dt_Vm=1/sf_Vm; %[sec]
                time_axis_Vm=[dt_Vm:dt_Vm:size(data.header(header).x_value(1).Vm,1)*dt_Vm];
                sf_airpuff = Param.header(header).stim.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                time_axis_airpuff = [dt_airpuff:dt_airpuff:size(data.header(header).x_value(1).airpuff,1)*dt_airpuff];
                sf_galvano = Param.header(header).stim.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;
                time_axis_galvano = [dt_galvano:dt_galvano:size(data.header(header).x_value(1).galvano,1)*dt_galvano];

                for x_value = 1:size(data.header(header).x_value,2)
                    
                      data.header(header).x_value(x_value).Vm_mean = mean(data.header(header).x_value(x_value).Vm,2);
                      data.header(header).x_value(x_value).Vm_std = std(data.header(header).x_value(x_value).Vm');
                      data.header(header).x_value(x_value).Ext1_mean = mean(data.header(header).x_value(x_value).Ext1,2);
                      data.header(header).x_value(x_value).Ext1_std = std(data.header(header).x_value(x_value).Ext1');
                end
        end
end
    
%% file 4 

% Depth 784: header 2 - intracellular, protocol airpuff (puff at 1.3sec)
% traces 1-29: 1-20: puff only, 21-29: puff+light
% No response to airpuff

header = 2;

    time_axis_airpuff = [dt_airpuff:dt_airpuff:size(data.header(header).x_value(1).airpuff,1)*dt_airpuff];    
    airpuff_vec = data.header(header).x_value(1).airpuff(:,1)./max(data.header(header).x_value(1).airpuff(:,1));
    sf_laser = Param.header(header).stim.sf_laser; %[1/sec]
    dt_laser = 1/sf_laser;
    time_axis_laser = [dt_laser:dt_laser:size(data.header(header).x_value(2).laser,1)*dt_laser];
    laser_vec = zeros(length(data.header(header).x_value(2).laser(:,1)),1);
    laser_vec(data.header(header).x_value(2).laser(:,1) > 2)=-0.5;
    
    for trace = 1:5
        Puff_alone_trace(:,trace) = data.header(header).x_value(1).Vm(:,trace);
        Puff_light_trace(:,trace) = data.header(header).x_value(2).Vm(:,trace);
    
    
        figure (1)
        set(gcf,'color','w')
        
        x1limits = [0 60000];
        x1ticks = [];
        y1limits = [-60 40];
        y1ticks = [-60:20:40];

        subplot(6,2,2.*trace-1)
        plot(Puff_alone_trace(:,trace))
        set( gca, 'xlim', x1limits, 'xtick', x1ticks, 'ylim', y1limits, 'ytick', y1ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        
        
        subplot(6,2,2.*trace)
        plot(Puff_light_trace(:,trace))
        set( gca, 'xlim', x1limits, 'xtick', x1ticks, 'ylim', y1limits, 'ytick', y1ticks,'yminortick', 'on',...
        'ticklength', [0.010 0.010],'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        
        
        subplot (6,2,11)
        plot(time_axis_airpuff,airpuff_vec, '-k','linewidth', 2)
        set( gca, 'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('Time [sec]' ,'FontSize', 12);
        
        subplot (6,2,12)
        hold on
        plot(time_axis_airpuff,airpuff_vec, '-k','linewidth', 2)
        plot(time_axis_laser(laser_vec~=0), laser_vec(laser_vec~=0), '.c')
        hold off
        set( gca, 'fontname', 'helvetica', 'fontweight', 'bold', 'box', 'off' );
        xlabel('Time [sec]' ,'FontSize', 12);
    end    
 
    mean_puff_only = mean(Puff_alone_trace,2);
    mean_puff_light = mean(Puff_light_trace,2);
    
    figure (2)
    set(gcf,'color','w')
        subplot(2,1,1)
            hold on 
            plot(mean_puff_only, '-b', 'linewidth', 2)
%             plot(LFP_mean_puff_only - std(LFP_mean_puff_only),'k');
            plot(mean_puff_light, '-r', 'linewidth', 1)
            hold off
        
        subplot(2,1,2)
            plot(time_axis_airpuff, airpuff_vec, '-k','linewidth', 2);