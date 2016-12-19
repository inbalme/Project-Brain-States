function Extract_ChAT_Data
% function Extract_ChAT_Data forms the database 'Extracted Data' used by
% Analyze_ChAT. This is the version that I am using!
clear all

global Exp

files = Get_ChAT_Files();

for fileind = 2:length(files) % analyzing files that were not yet saved.

    fname = files(fileind).name;
    path = files(fileind).path;
    Exp = [];
    data = [];
    Param = [];
    readExpMainHeader(fname,path);
    
    channelVm = files(fileind).channel_Vm;
    channelLFP = files(fileind).channel_LFP;
    channelGalvano = files(fileind).channel_galvano;
    channelAirpuff = files(fileind).channel_airpuff;
    channelEKG = files(fileind).channel_EKG;
    channelEEG = files(fileind).channel_EEG;
    channelLaser = files(fileind).channel_laser;
    
  
    for header = 1:length(files(fileind).headers)       %header is the counter of the headers
        x_value_ind= zeros(20,1);
        stim_type = files(fileind).header_stim(header).type;
        run_header = files(fileind).headers(header);  %run_header is the absolut number of the header in the file
        Exp.Header(run_header).headerInfo = readExpOneHeader(run_header);
        %Exp.Header(run_header).headerInfo.subexp_type
        for run_trace = 1:length(Exp.Header(run_header).headerInfo.blockLocations)
            block = readOneBlockdata(run_header,run_trace);
            
            if stim_type== 1 || stim_type== 2 || stim_type== 3 || stim_type==4 || stim_type==5 ||stim_type==6
                x_value = 1;
                data.header(header).x_value(x_value).Vm(:,run_trace)= block.channel(channelVm).chdata;
                data.header(header).x_value(x_value).LFP(:,run_trace)= block.channel(channelLFP).chdata;
                data.header(header).x_value(x_value).airpuff(:,run_trace) = block.channel(channelAirpuff).chdata;
                data.header(header).x_value(x_value).galvano(:,run_trace) = block.channel(channelGalvano).chdata;
                data.header(header).x_value(x_value).laser(:,run_trace)= block.channel(channelLaser).chdata;
                data.header(header).x_value(x_value).EKG(:,run_trace) = block.channel(channelEKG).chdata;
                data.header(header).x_value(x_value).EEG(:,run_trace)= block.channel(channelEEG).chdata;
                data.header(header).x_value(x_value).original_header = run_header;
                   
                else if  stim_type== 7 || stim_type== 8;
                
                x_value_ind(block.newXvalue,1) = x_value_ind(block.newXvalue,1)+1; %similar to run_trace index but contains 2 indices for the 2 X values
                        
                data.header(header).x_value(block.newXvalue).Vm(:,x_value_ind(block.newXvalue))=  block.channel(channelVm).chdata;
                data.header(header).x_value(block.newXvalue).LFP(:,x_value_ind(block.newXvalue))=  block.channel(channelLFP).chdata;
                data.header(header).x_value(block.newXvalue).airpuff(:,x_value_ind(block.newXvalue)) = block.channel(channelAirpuff).chdata;
                data.header(header).x_value(block.newXvalue).galvano(:,x_value_ind(block.newXvalue)) = block.channel(channelGalvano).chdata;
                data.header(header).x_value(block.newXvalue).laser(:,x_value_ind(block.newXvalue))=  block.channel(channelLaser).chdata;
                data.header(header).x_value(block.newXvalue).EKG(:,x_value_ind(block.newXvalue)) = block.channel(channelEKG).chdata;
                data.header(header).x_value(block.newXvalue).EEG(:,x_value_ind(block.newXvalue))=  block.channel(channelEEG).chdata;
                data.header(header).original_header = run_header;
%                 data.header(header).stim_direction = files(fileind).header_stim(header).stim_direction;           %1=rost-caud, 2=dors-vent, 3= caud-rost, 4=vent-dors
                      end
            end
%             block info
        end
%             header info 
    Param.header(header).general.exp_type                  = Exp.Header(run_header).headerInfo.exp_type;
    Param.header(header).general.sub_type                  = Exp.Header(run_header).headerInfo.subexp_type;
    Param.header(header).stim.type                             = files(fileind).header_stim(header).type;
    Param.header(header).stim.type_name                   = {'spont alone', 'spont+laser', 'puff alone', 'puff+laser', 'galvano alone', 'galvano+laser', 'protocol airpuff', 'protocol galvano'};
    Param.header(header).stim.sf                                 = block.scanrate; % the sampling frequency
    Param.header(header).stim.dt                                = 1/(Param.header(header).stim.sf);
    Param.header(header).stim.facade                         = block.IncedValues;
   
        if strcmp('Standby', Param.header(header).general.exp_type);
      Param.header(header).stim.facade_name             = {''; 'board'; 'channel'; 'waveform'; 'delay'; ''; 'amp'; 'on_duration'; 'cycles'; ''; ''; ''; ''; ''; 'board'; 'channel'; 'waveform'; 'delay'; ''; 'amp'; 'on_duration'; 'cycles'; ''; ''; ''; ''; ''};
        else if strcmp('Whisker_Stim', Param.header(header).general.exp_type); 
           Param.header(header).stim.facade_name        = {'low_amp'; 'mid_amp'; 'high_amp'; 'channel'; 'board'; 'waveform'; 'delay'; 'on_duration'; 'off_duration'; 'cycles'};
            end
       end
    Param.header(header).stim.sf_galvano                      = Param.header(header).stim.sf/(length(block.channel(1,1).chdata)/length(block.channel(1,channelGalvano).chdata));
    Param.header(header).stim.sf_airpuff                     = Param.header(header).stim.sf/(length(block.channel(1,1).chdata)/length(block.channel(1,channelAirpuff).chdata)); 
    Param.header(header).stim.sf_laser                     = Param.header(header).stim.sf/(length(block.channel(1,1).chdata)/length(block.channel(1,channelLaser).chdata));      
    Param.header(header).stim.sf_LFP                     = Param.header(header).stim.sf/(length(block.channel(1,1).chdata)/length(block.channel(1,channelLFP).chdata));           
    Param.header(header).stim.sf_EEG                     = Param.header(header).stim.sf/(length(block.channel(1,1).chdata)/length(block.channel(1,channelEEG).chdata));    
    Param.header(header).stim.sf_EKG                     = Param.header(header).stim.sf/(length(block.channel(1,1).chdata)/length(block.channel(1,channelEKG).chdata));       

    
    if isfield(files(fileind).header_stim(header), 'stim_direction')
        Param.header(header).stim.direction                      = files(fileind).header_stim(header).stim_direction;
        Param.header(header).stim.direction_name            = {'rost-caud'; 'dors-vent'; 'caud-rost'; 'vent-dors'};
    else
        Param.header(header).stim.direction                      = [];
    end
    
    end
%             file info
    

    Param.name                        = files(fileind).name;
    Param.path                          = files(fileind).path;
    Param.exp                           = files(fileind).exp;
    Param.orig_headers                     = files(fileind).headers;
    Param.depth                        = files(fileind).depth;
    Param.pr_whisker                = files(fileind).pr_whisker;
    Param.adj_whisker                = files(fileind).adj_whisker;
%     Param.pr_direction               = files(fileind).pr_direction;
    Param.time                           = block.timesincejan1904;
    Param.electrode_resist         = files(fileind).electrode_resist;
    
%    
%     
%     if isfield(Exp.Header(1).headerInfo,'positionsof4motors')
%         Param.stimulation.locationIn  = Exp.Header(run_head).headerInfo.positionsof4motors;
%         Param.stimulation.locationEx  = Exp.Header(run_head).headerInfo.positionsof4motorsB;
%     else
%         Param.stimulation.locationIn    = [ 0 ; 0 ; 0 ;0];
%         Param.stimulation.locationEx    = [ 0 ; 0 ; 0 ;0];
%     end

    name = fname(1:(strfind(fname,'.')-1));
    cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
    save( name, 'Exp', 'data', 'Param')  % for specific parameter from Exp look at GetPatameters
    files(fileind).extracted_name = name;
    files(fileind).extracted_path = 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
    
end
    name = 'ChAT_Files';
    cd 'D:\Inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
    save( name,'files')