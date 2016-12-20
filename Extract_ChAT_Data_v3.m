function Extract_ChAT_Data_v3
% function Extract_ChAT_Data_v3 forms the database 'Extracted Data' used by
% Analyze_ChAT. This is the version that I am using for files in which
% protocols are saved such that every header is one repetition over all the
% traces types (x-values).
%created based on NBES analysis on October 2016

clear all

global Exp

files = Get_ChAT_Files_v3();

for fileind =80:87; 
    clearvars -except  files fileind Exp
no_spikes_flag=1;
    fname = files(fileind).name;
    path = files(fileind).path;
    Exp = [];
    data = [];
    Param = [];
    code=[]; data_no_spikes=[]; ch_gain=1;
    x_value_ind= zeros(20,1);
    readExpMainHeader(fname,path);
    %             file info
    Param.name                        = files(fileind).name;
    Param.path                          = files(fileind).path;
    Param.extracted_name        = files(fileind).extracted_name;
    Param.extracted_path          = files(fileind).extracted_path;
    Param.animal                       = files(fileind).animal;
    Param.anesthesia_initial        = files(fileind).anesthesia_initial;
    Param.anesthesia_maintenance  = files(fileind).anesthesia_maintenance;
    Param.pipette_solution         = files(fileind).pipette_solution;
    Param.recording_type         = files(fileind).recording_type;
    Param.recording_location        = files(fileind).recording_location;
    Param.protocol                    = files(fileind).protocol;
    Param.electrode_resist         = files(fileind).electrode_resist;    
    Param.orig_headers             = files(fileind).headers;
    Param.depth                        = files(fileind).depth;
    Param.channel_Vm = files(fileind).channel_Vm;
    Param.channel_I1 = files(fileind).channel_I1;
    Param.channel_galvano = files(fileind).channel_galvano;
    Param.channel_airpuff = files(fileind).channel_airpuff;
    Param.channel_EKG = files(fileind).channel_EKG;
    Param.channel_EEG = files(fileind).channel_EEG;
    Param.channel_laser = files(fileind).channel_laser;
    Param.EEG_flag = files(fileind).EEG_flag;
    Param.dual_rec_flag = files(fileind).dual_rec_flag;
    Param.channel_V2 = files(fileind).channel_V2;
    Param.channel_I2 = files(fileind).channel_I2;
    if isfield(files,'channel_ES')
        Param.channel_ES = files(fileind).channel_ES;
    end
    if isfield(files,'mode_name')
        Param.mode_name = files(fileind).mode_name;
    end
     if isfield(files,'mode_val')
        Param.mode_val = files(fileind).mode_val;      
        POG_loc=strcmp('PrimaryOutputGain',Param.mode_name);
          if isnan(Param.mode_val{1,POG_loc==1});
            else
                ch_gain=10./Param.mode_val{1,POG_loc==1}; %Ilan's data aquisition system automatically divides by 10;
          end
    end
    Param.pr_whisker                = files(fileind).pr_whisker;
    Param.adj_whisker                = files(fileind).adj_whisker;
%     Param.time                           = block.timesincejan1904;
    Param.galvano                        = files(fileind).galvano;
    Param.offset                        = files(fileind).offset;
    Param.comments                        = files(fileind).comments;
    
    Param.stim_type                             = files(fileind).stim_type;
    Param.stim_type_name                   = {'spont' , 'light alone', 'puff alone', 'puff+light', 'galvano alone', 'galvano+light', 'protocol light only', 'protocol galvano'};  
    
     stim_type = files(fileind).stim_type;
      
    for header = 1:length(files(fileind).headers)       %header is the counter of the headers
        
        run_header = files(fileind).headers(header);  %run_header is the absolut number of the header in the file
        Exp.Header(run_header).headerInfo = readExpOneHeader(run_header);
        %Exp.Header(run_header).headerInfo.subexp_type
        for run_trace = 1:length(Exp.Header(run_header).headerInfo.blockLocations)
            block = readOneBlockdata(run_header,run_trace);
         if isfield(block, 'channel') ==0
             break
         end
         Param.sf_Vm                                = block.scanrate; % the sampling frequency
    Param.dt_Vm                                = 1/(Param.sf_Vm);
    Param.facade                         = block.IncedValues;
    Param.sf_I1                      = Param.sf_Vm/(length(block.channel(1,1).chdata)/length(block.channel(1,Param.channel_I1).chdata));
%    if Param.dual_rec_flag ==1
       Param.sf_V2 = Param.sf_Vm/(length(block.channel(1,1).chdata)/length(block.channel(1,Param.channel_V2).chdata));
        Param.sf_I2 = Param.sf_Vm/(length(block.channel(1,1).chdata)/length(block.channel(1,Param.channel_I2).chdata));
%    end
    Param.sf_galvano                      = Param.sf_Vm/(length(block.channel(1,1).chdata)/length(block.channel(1,Param.channel_galvano).chdata));
    Param.sf_airpuff                     = Param.sf_Vm/(length(block.channel(1,1).chdata)/length(block.channel(1,Param.channel_airpuff).chdata)); 
    Param.sf_laser                     = Param.sf_Vm/(length(block.channel(1,1).chdata)/length(block.channel(1,Param.channel_laser).chdata));      
    Param.sf_EEG                     = Param.sf_Vm/(length(block.channel(1,1).chdata)/length(block.channel(1,Param.channel_EEG).chdata));    
    Param.sf_EKG                     = Param.sf_Vm/(length(block.channel(1,1).chdata)/length(block.channel(1,Param.channel_EKG).chdata));   
    
            if stim_type== 1 || stim_type== 2 || stim_type== 3 || stim_type==4 || stim_type==5 ||stim_type==6
                x_value = 1;
                data.x_value(x_value).Vm(:,run_trace)= block.channel(Param.channel_Vm).chdata;
                data.x_value(x_value).airpuff(:,run_trace) = block.channel(Param.channel_airpuff).chdata;
                data.x_value(x_value).galvano(:,run_trace) = block.channel(Param.channel_galvano).chdata;
               
                 if isempty(block.channel(Param.channel_EKG).chdata)
                 else
                    data.x_value(x_value).EKG(:,run_trace) = block.channel(Param.channel_EKG).chdata;
                 end
                
                if isempty(block.channel(Param.channel_laser).chdata)
                 else
                    data.x_value(x_value).laser(:,run_trace)= block.channel(Param.channel_laser).chdata;
                end
                 
                if isempty(block.channel(Param.channel_EEG).chdata)
                 else
                    data.x_value(x_value).EEG(:,run_trace)= block.channel(Param.channel_EEG).chdata;
                 end
                 
                 if isempty(block.channel(Param.channel_V2).chdata)
                 else
                    data.x_value(x_value).V2(:,run_trace)= block.channel(Param.channel_V2).chdata;
                 end
                   if isempty(block.channel(Param.channel_I1).chdata)
                 else
                    data.x_value(x_value).I1(:,run_trace)= block.channel(Param.channel_I1).chdata;
                 end
                   if isempty(block.channel(Param.channel_I2).chdata)
                 else
                    data.x_value(x_value).I2(:,run_trace)= block.channel(Param.channel_I2).chdata;
                 end
                 
                data.x_value(x_value).original_header(1,run_header) = run_header;
                data.x_value(x_value).original_trace(1,run_trace) = run_trace;
                   
                else if  stim_type== 7 || stim_type== 8;
                
                x_value_ind(block.newXvalue,1) = x_value_ind(block.newXvalue,1)+1; %similar to run_trace index but contains 2 indices for the 2 X values
                        
                data.x_value(block.newXvalue).Vm(:,x_value_ind(block.newXvalue))=  block.channel(Param.channel_Vm).chdata;
                
                if isempty(block.channel(Param.channel_galvano).chdata)
                else
                    data.x_value(block.newXvalue).galvano(:,x_value_ind(block.newXvalue)) = block.channel(Param.channel_galvano).chdata;
                 end 
                if isempty(block.channel(Param.channel_airpuff).chdata)
                else
                    data.x_value(block.newXvalue).airpuff(:,x_value_ind(block.newXvalue)) = block.channel(Param.channel_airpuff).chdata;
                 end 
                
                if isempty(block.channel(Param.channel_laser).chdata)
                else
                    data.x_value(block.newXvalue).laser(:,x_value_ind(block.newXvalue))=  block.channel(Param.channel_laser).chdata;
                 end
                
                 if isempty(block.channel(Param.channel_EKG).chdata)
                else
                    data.x_value(block.newXvalue).EKG(:,x_value_ind(block.newXvalue)) = block.channel(Param.channel_EKG).chdata;
                end
                
                if isempty(block.channel(Param.channel_EEG).chdata)
                else
                    data.x_value(block.newXvalue).EEG(:,x_value_ind(block.newXvalue))=  block.channel(Param.channel_EEG).chdata;
                end
                
                 if isempty(block.channel(Param.channel_V2).chdata)
                 else
                    data.x_value(block.newXvalue).V2(:,x_value_ind(block.newXvalue))=  block.channel(Param.channel_V2).chdata;
                 end
                 
                 if isempty(block.channel(Param.channel_I1).chdata)
                 else
                    data.x_value(block.newXvalue).I1(:,x_value_ind(block.newXvalue))=  block.channel(Param.channel_I1).chdata;
                 end
                 
                    if isempty(block.channel(Param.channel_I2).chdata)
                 else
                    data.x_value(block.newXvalue).I2(:,x_value_ind(block.newXvalue))=  block.channel(Param.channel_I2).chdata;
                 end
                 
                data.x_value(block.newXvalue).original_header(1,run_header) = run_header;
                data.x_value(block.newXvalue).original_trace(1,run_trace) = run_trace;
                      end
            end
         end
%%             block info

    Param.exp_type                  = Exp.Header(run_header).headerInfo.exp_type;
    Param.sub_type                  = Exp.Header(run_header).headerInfo.subexp_type;
   
        if strcmp('Standby', Param.exp_type);
      Param.facade_name             = {''; 'board'; 'channel'; 'waveform'; 'delay'; ''; 'amp'; 'on_duration'; 'cycles'; ''; ''; ''; ''; ''; 'board'; 'channel'; 'waveform'; 'delay'; ''; 'amp'; 'on_duration'; 'cycles'; ''; ''; ''; ''; ''};
            end
    end     
%%    
%   cd 'D:\Dropbox\Inbal M.Sc\MATLAB\Project Brain States'
    
    color_table = [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];
    
laser_flag = zeros(1,size(data.x_value,2)); %for each x-value, laser_flag is1 if there was laser activation
galvano_flag = zeros(1,size(data.x_value,2)); %for each x-value, galvano_flag is1 if there was galvano activation

 locations_x_galvano = [];
 locations_x_ES = [];
    
protocol = Param.protocol;

if isempty(protocol)
    
else
    switch protocol
            
        case 5
             laser_flag(1,:) = [1 0 1 0 1];
             galvano_flag(1,:) = [0 1 1 1 1];
             galvano_del = Param.facade(4);
             galvano_ondur = Param.facade(5);
             galvano_nstim = Param.facade(6);
             galvano_freq = Param.facade(7);
             laser_del = Param.facade(11);
             laser_ondur = Param.facade(12);
             laser_nstim = Param.facade(13);
             laser_freq = Param.facade(14);
        case 6
            if length(laser_flag(1,:))==9
                laser_flag(1,:) = [1 1 1 0 0 0 1 1 1];
                galvano_flag(1,:) = [0 0 0 1 1 1 1 1 1];
            else
             laser_flag(1,:) = [1 1 1 0 0 0 1 1 1 0 0 0 1 1 1];
             galvano_flag(1,:) = [0 0 0 1 1 1 1 1 1 1 1 1 1 1 1];   
            end
        case 8
             laser_flag(1,:) = [1 1 1];
             galvano_flag(1,:) = [0 0 0];  
        case 9
%              laser_flag(1,:) = [0 0 0 0 0 0 0];
%              galvano_flag(1,:) = [0 0 0 0 0 0 0];  
        case 10
            laser_flag(1,:) = [1 0 1];
            galvano_flag(1,:) = [0 1 1 ];  
            galvano_del = Param.facade(16);
            galvano_test_del = Param.facade(17);
            galvano_ondur = Param.facade(5);
            galvano_nstim = Param.facade(6);
            galvano_freq = Param.facade(7);
            laser_del = Param.facade(11);
            laser_ondur = Param.facade(12);
            laser_nstim = Param.facade(13);
            laser_freq = Param.facade(14);
        case 11
            laser_flag(1,:) = [1 1 1 1 1 1];
            galvano_flag(1,:) = [0 0 0 1 1 1 ];  
            galvano_del = [Param.facade(28) Param.facade(31)];
            galvano_test_del = [Param.facade(29) Param.facade(32)];
            galvano_ondur = Param.facade(14);
            galvano_nstim = Param.facade(15);
            galvano_freq = Param.facade(16);
            laser_del = Param.facade(22);
            laser_ondur = Param.facade(23);
            laser_nstim = Param.facade(24);
            laser_freq = Param.facade(25);
            curr_inj_del = [Param.facade(27) Param.facade(30)] ;
            if length(Param.facade)>32
                curr_inj_del = [Param.facade(27) Param.facade(30) Param.facade(33)];
            end
            curr_inj_depo = Param.facade(1);
            curr_inj_hyper = Param.facade(2);
            curr_inj_dur = Param.facade(7);
            
    end
end
                 
    for x_value = 1:size(data.x_value,2)           
                   
        trace = 1:size(data.x_value(x_value).Vm,2);
        
  if protocol==9 && x_value==1
              continue
        else                
 raw_data{1}(:,1:length(trace),x_value) = data.x_value(x_value).Vm(:,trace).*ch_gain; 
 raw_data{2}(:,1:length(trace),x_value) = data.x_value(x_value).I1(:,trace); 
 raw_data{3}(:,1:length(trace),x_value) = data.x_value(x_value).V2(:,trace);
 raw_data{4}(:,1:length(trace),x_value) = data.x_value(x_value).I2(:,trace); 
  end    
          if Param.mode_val{strcmp('VC',Param.mode_name)}==1;
               raw_data{1}(:,1:length(trace),x_value)=-1.* raw_data{1}(:,1:length(trace),x_value);
          end
              
                sf{1} = Param.sf_Vm;
                sf{2} = Param.sf_I1;
                sf{3} = Param.sf_V2;
                sf{4} = Param.sf_I2;
                             
                sf_airpuff = Param.sf_airpuff; %[1/sec]
                dt_airpuff = 1/sf_airpuff;
                sf_galvano = Param.sf_galvano; %[1/sec]
                dt_galvano = 1/sf_galvano;
  
 for ii=1:4
               dt{ii}=1/sf{ii}; %[sec]
                time_axis{ii}(:,x_value) = (1:size(raw_data{ii}(:,:,x_value),1))*dt{ii};
              
                if (strcmp(Param.recording_type,'cell attached') | strcmp(Param.recording_type,'LFP') | strcmp(Param.recording_type,'cell attached-LFP')) %skip this part for extracellular recordings
                    code=[]; data_no_spikes=[];
                    else
                    %Removing spikes
                    if no_spikes_flag==1;
                        for traces=1:size(raw_data{1}(:,:,x_value),2)
                            [code{1}(:,traces,x_value), data_no_spikes{1}(:,traces,x_value)] = fn_RemoveSpikesMO2(raw_data{1}(:,traces,x_value)', dt{1});
                        end
                    end
                end   
 
               
                
%              ES_begin_time = Param. ES_param_val(2);
%              ES_duration_time = Param. ES_param_val(5)*(1/Param. ES_param_val(4)); %number of pulses times 1/f
    if sum(laser_flag)==0
        else
             laser_begin_time = laser_del;
             laser_duration_time = laser_nstim*(1/laser_freq); %number of pulses times 1/f
             if laser_nstim==1;
                 laser_duration_time=laser_ondur;
             end
             laser_end_time = laser_begin_time+laser_duration_time;
             laser_vec{ii}(:,1) = zeros(length(time_axis{ii}),1);
             laser_vec{ii}(laser_begin_time*sf{ii}:laser_end_time*sf{ii},1)=1;
             
             locations_x_laser{ii}(1,:) = laser_begin_time*sf{ii}; %arranging the galvano begin and end locations in one variable for plotting
             locations_x_laser{ii}(2,:) = laser_end_time*sf{ii}; 
    end
                  
 end
    

                    if Param.EEG_flag==1;
                    data_EEG(:,1:length(trace),x_value) = data.x_value(x_value).EEG(:,trace);
                         sf_EEG = Param.sf_EEG;
                         dt_EEG = 1/sf_EEG;
                    end

                if Param.sf_airpuff==0
                else
                    time_axis_airpuff(:,x_value) = [1:size(data.x_value(x_value).airpuff,1)]*dt_airpuff;    
                    airpuff_vec(:,x_value) = data.x_value(x_value).airpuff(:,1)./...
                    max(data.x_value(x_value).airpuff(:,1));
                end

                 if Param.sf_galvano==0
                else
                    time_axis_galvano = [1:size(data.x_value(x_value).galvano,1)]*dt_galvano;
%                     galvano_vec(:,x_value) = zeros(length(data.x_value(x_value).galvano(:,1)),1);
%                     galvano_sub_mean = data.x_value(x_value).galvano(:,1) - mean(data.x_value(x_value).galvano(:,1));
%                     galvano_threshold = abs(mean(galvano_sub_mean)) + 3.*abs(std(galvano_sub_mean));
%                     galvano_vec(abs(galvano_sub_mean) > galvano_threshold, x_value)=1;   %turning the galvano trace into binary
                 end
             
                    if galvano_flag(1,x_value) ==1; %condition will be fulfilled if there was galvano activation.
                        galvano_ISI = 1/galvano_freq; %[sec]
                        for gg=1:length(galvano_del)
                            gd=galvano_del(gg); galvano_begin_time_temp=[];                           
                            galvano_begin_time_temp = gd(ones(1,galvano_nstim))+([0 1:galvano_nstim-1])*galvano_ISI;
                            galvano_end_time_temp =galvano_begin_time_temp+galvano_ondur(ones(1,galvano_nstim));
                            gt_flag = exist('galvano_test_del','var'); %returns 1 if name exists as a variable
                                if gt_flag==1
                                    galvano_begin_time_temp(end+1) = galvano_test_del(gg);
                                    galvano_end_time_temp(end+1) = galvano_test_del(gg)+galvano_ondur;
                                end  
%                           galvano_begin_time{x_value}(gg,:) = galvano_begin_time_temp;
%                           galvano_end_time{x_value}(gg,:) =galvano_end_time_temp;
                         locations_x_galvano{x_value}(2*gg-1,:) = galvano_begin_time_temp*sf{1}; %arranging the galvano begin and end locations in one variable for plotting
                         locations_x_galvano{x_value}(2*gg,:) = galvano_end_time_temp*sf{1};
                         time_x_galvano{x_value}(2*gg-1,:) = galvano_begin_time_temp; %arranging the galvano begin and end times in one variable for plotting
                         time_x_galvano{x_value}(2*gg,:) = galvano_end_time_temp;            
                        end
                        
                    end      
                
    end
    
     %  Stim vectors
 for x_value = 1:size(data.x_value,2) 
         if isempty(locations_x_laser)
             stim1_X = [];
         else
            stim1_X = locations_x_laser;
         end

        if isempty(locations_x_galvano)
            stim2_X{x_value} = [];
            stim2_time{x_value} = [];
        else
            stim2_X{x_value} = locations_x_galvano{x_value}(:,:);
            stim2_time{x_value} = time_x_galvano{x_value}(:,:);
        end
 end
%%
 name = files(fileind).extracted_name;
    cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data 2016';
    save( name, 'Param','laser_flag','color_table', 'data','raw_data','code', 'data_no_spikes','galvano_flag','stim1_X','stim2_X', 'stim2_time','time_axis')  % for specific parameter from Exp look at GetPatameters

 
end
 