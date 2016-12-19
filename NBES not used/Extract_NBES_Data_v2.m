function Extract_NBES_Data_v2
% function Extract_NBES_Data_v2 forms the database 'Extracted Data' used by
% Analyze_NBES. This is the version that I am using for files in which
% protocols are saved such that every header is one repetition over all the
% traces types (x-values)

clear all

global Exp

files = Get_NBES_Files_v2();

for fileind =5; 

    fname = files(fileind).name;
    path = files(fileind).path;
    Exp = [];
    data = [];
    Param = [];
    x_value_ind= zeros(20,1);
    readExpMainHeader(fname,path);
    
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
%             block info
        end

    Param.exp_type                  = Exp.Header(run_header).headerInfo.exp_type;
    Param.sub_type                  = Exp.Header(run_header).headerInfo.subexp_type;
   
        if strcmp('Standby', Param.exp_type);
      Param.facade_name             = {''; 'board'; 'channel'; 'waveform'; 'delay'; ''; 'amp'; 'on_duration'; 'cycles'; ''; ''; ''; ''; ''; 'board'; 'channel'; 'waveform'; 'delay'; ''; 'amp'; 'on_duration'; 'cycles'; ''; ''; ''; ''; ''};
            end
        end
    
     
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
    Param. ES_param_name      = files(fileind).ES_param_name;
    Param. ES_param_val      = files(fileind).ES_param_val;
    Param.electrode_resist         = files(fileind).electrode_resist;    
    Param.orig_headers             = files(fileind).headers;
    Param.depth                        = files(fileind).depth;
    if isfield(files,'mode_name')
        Param.mode_name = files(fileind).mode_name;
    end
     if isfield(files,'mode_val')
        Param.mode_val = files(fileind).mode_val;
    end
    Param.pr_whisker                = files(fileind).pr_whisker;
    Param.adj_whisker                = files(fileind).adj_whisker;
%     Param.time                           = block.timesincejan1904;
    Param.galvano                        = files(fileind).galvano;
    Param.offset                        = files(fileind).offset;
    Param.comments                        = files(fileind).comments;
    
    Param.stim_type                             = files(fileind).stim_type;
    Param.stim_type_name                   = {'spont' , 'light alone', 'puff alone', 'puff+light', 'galvano alone', 'galvano+light', 'protocol light only', 'protocol galvano'};
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
    
    
%    
%     
%     if isfield(Exp.Header(1).headerInfo,'positionsof4motors')
%         Param.stimulation.locationIn  = Exp.Header(run_head).headerInfo.positionsof4motors;
%         Param.stimulation.locationEx  = Exp.Header(run_head).headerInfo.positionsof4motorsB;
%     else
%         Param.stimulation.locationIn    = [ 0 ; 0 ; 0 ;0];
%         Param.stimulation.locationEx    = [ 0 ; 0 ; 0 ;0];
%     end

    name = files(fileind).extracted_name;
    cd(files(fileind).extracted_path);
%     cd 'D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
    save( name, 'Exp', 'data', 'Param')  % for specific parameter from Exp look at GetPatameters
 
end
 