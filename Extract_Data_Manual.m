%function [Exp TotalHeader] =Extract_Data ();
%this function is used when extracting data manually, file by file
% close all; 
clear all; clc;
set=1; %1 for set 1, 2 for set 2
traces_to_remove = []; %1:12; %1:length(Exp.Header(run_header).headerInfo.blockLocations);
ex
disp('after opening file of interest, press any key to continue');
pause

EXPGLOBALS;

global Exp;

protocol_type = 2; %protocol_type is 1 if each header is a whole protocol and 2 if each header is a single repeat of the protocol or standby
headers=[1:2]; %headers to extract
switch set
    case 1
        ch_Data1=1;
        ch_Data2=3; 
        ch_EEG=13;
        ch_EKG=5;
        ch_Puff=6; 
        ch_Galvo1=7;
        ch_Galvo2=8;
        ch_Laser=10;
        
    case 2
        ch_Data1=1;
        ch_Data2=3; 
        ch_EEG=5;
        ch_EKG=6;
        ch_Puff=13; 
        ch_Galvo1=14;
        ch_Galvo2=15;
        ch_Laser=10;
end

switch protocol_type
    case 1
 for header = 1:length(headers)       %header is the counter of the headers
   x_value_ind= zeros(20,1);
   data = [];
   
    run_header = headers(header);  %run_header is the absolut number of the header in the file
        Exp.Header(run_header).headerInfo = readExpOneHeader(run_header);

        for run_trace=1:length(Exp.Header(run_header).headerInfo.blockLocations);
            block = readOneBlockdata(run_header,run_trace);
         if isfield(block, 'channel') ==0
             break
         end
                 x_value_ind(block.newXvalue,1) = x_value_ind(block.newXvalue,1)+1; %similar to run_trace index but contains several indices for all the X values
                        
                data(header).x_value(block.newXvalue).data1(:,x_value_ind(block.newXvalue))=  block.channel(ch_Data1).chdata;
                if isempty(block.channel(ch_Data2).chdata)
                else
                    data(header).x_value(block.newXvalue).data2(:,x_value_ind(block.newXvalue))=  block.channel(ch_Data2).chdata;
                 end
                
                if isempty(block.channel(ch_Galvo1).chdata)
                else
                    data(header).x_value(block.newXvalue).galvano1(:,x_value_ind(block.newXvalue)) = block.channel(ch_Galvo1).chdata;
                end 
                  if isempty(block.channel(ch_Galvo2).chdata)
                else
                    data(header).x_value(block.newXvalue).galvano2(:,x_value_ind(block.newXvalue)) = block.channel(ch_Galvo2).chdata;
                 end 
                if isempty(block.channel(ch_Puff).chdata)
                else
                    data(header).x_value(block.newXvalue).airpuff(:,x_value_ind(block.newXvalue)) = block.channel(ch_Puff).chdata;
                 end 
                
                if isempty(block.channel(ch_Laser).chdata)
                else
                    data(header).x_value(block.newXvalue).laser(:,x_value_ind(block.newXvalue))=  block.channel(ch_Laser).chdata;
                 end
                
                 if isempty(block.channel(ch_EKG).chdata)
                else
                    data(header).x_value(block.newXvalue).EKG(:,x_value_ind(block.newXvalue)) = block.channel(ch_EKG).chdata;
                end
                
                if isempty(block.channel(ch_EEG).chdata)
                else
                    data(header).x_value(block.newXvalue).EEG(:,x_value_ind(block.newXvalue))=  block.channel(ch_EEG).chdata;
                end
                             
                data(header).x_value(block.newXvalue).original_header(1,run_header) = run_header;
                data(header).x_value(block.newXvalue).original_trace(1,run_trace) = run_trace;
        end
                      
     path='D:\Inbal M.Sc\Data PhD\SSA Data\Extracted Data';
     fname = Exp.namepath.name;
     name = fname(1:(findstr(fname,'.')-1));
     name = [name , '_header', num2str(headers(header))];
     cd(path);
     save( name, 'data');
     cd(Exp.namepath.path);
end
    case 2
 x_value_ind= zeros(20,1);
%%
  for header = 1:length(headers)       %header is the counter of the headers
        
        run_header = headers(header);  %run_header is the absolut number of the header in the file
        Exp.Header(run_header).headerInfo = readExpOneHeader(run_header);
        
        for run_trace=1:length(Exp.Header(run_header).headerInfo.blockLocations);
            block = readOneBlockdata(run_header,run_trace);
         if isfield(block, 'channel') ==0
             break
         end
            if strcmp('Standby', Exp.Header(run_header).headerInfo.exp_type)  
                x_value = 1;
                data(header).x_value(x_value).data1(:,run_trace)= block.channel(ch_Data1).chdata;
                 if isempty(block.channel(ch_Data2).chdata)
                 else
                    data(header).x_value(x_value).data2(:,run_trace)= block.channel(ch_Data2).chdata;
                 end
                if isempty(block.channel(ch_Puff).chdata)
                else
                data(header).x_value(x_value).airpuff(:,run_trace) = block.channel(ch_Puff).chdata;
                end
                if isempty(block.channel(ch_Galvo1).chdata)
                else
                data(header).x_value(x_value).galvano1(:,run_trace) = block.channel(ch_Galvo1).chdata;
                end
                if isempty(block.channel(ch_Galvo2).chdata)
                else
                data(header).x_value(x_value).galvano2(:,run_trace) = block.channel(ch_Galvo2).chdata;
                end
                 if isempty(block.channel(ch_EKG).chdata)
                 else
                    data(header).x_value(x_value).EKG(:,run_trace) = block.channel(ch_EKG).chdata;
                 end
                
                if isempty(block.channel(ch_Laser).chdata)
                 else
                    data(header).x_value(x_value).laser(:,run_trace)= block.channel(ch_Laser).chdata;
                end
                 
                if isempty(block.channel(ch_EEG).chdata)
                 else
                    data(header).x_value(x_value).EEG(:,run_trace)= block.channel(ch_EEG).chdata;
                end    
                 
                data(header).x_value(x_value).original_header(1,run_header) = run_header;
                data(header).x_value(x_value).original_trace(1,run_trace) = run_trace;
                
            else %if it is from a protocol
                
                x_value_ind(block.newXvalue,1) = x_value_ind(block.newXvalue,1)+1; %similar to run_trace index but contains several indices for all the X values
                        
                data(1).x_value(block.newXvalue).data1(:,x_value_ind(block.newXvalue))=  block.channel(ch_Data1).chdata;
                if isempty(block.channel(ch_Data2).chdata)
                else
                    data(1).x_value(block.newXvalue).data2(:,x_value_ind(block.newXvalue))=  block.channel(ch_Data2).chdata;
                 end
                
                if isempty(block.channel(ch_Galvo1).chdata)
                else
                    data(1).x_value(block.newXvalue).galvano1(:,x_value_ind(block.newXvalue)) = block.channel(ch_Galvo1).chdata;
                end 
                  if isempty(block.channel(ch_Galvo2).chdata)
                else
                    data(1).x_value(block.newXvalue).galvano2(:,x_value_ind(block.newXvalue)) = block.channel(ch_Galvo2).chdata;
                 end 
                if isempty(block.channel(ch_Puff).chdata)
                else
                    data(1).x_value(block.newXvalue).airpuff(:,x_value_ind(block.newXvalue)) = block.channel(ch_Puff).chdata;
                 end 
                
                if isempty(block.channel(ch_Laser).chdata)
                else
                    data(1).x_value(block.newXvalue).laser(:,x_value_ind(block.newXvalue))=  block.channel(ch_Laser).chdata;
                 end
                
                 if isempty(block.channel(ch_EKG).chdata)
                else
                    data(1).x_value(block.newXvalue).EKG(:,x_value_ind(block.newXvalue)) = block.channel(ch_EKG).chdata;
                end
                
                if isempty(block.channel(ch_EEG).chdata)
                else
                    data(1).x_value(block.newXvalue).EEG(:,x_value_ind(block.newXvalue))=  block.channel(ch_EEG).chdata;
                end
                             
                data(1).x_value(block.newXvalue).original_header(1,run_header) = run_header;
                data(1).x_value(block.newXvalue).original_trace(1,run_trace) = run_trace;
                      end
        end
        data(1).x_value.data1(:,traces_to_remove)=[];
        data(1).x_value.data2(:,traces_to_remove)=[];
        data(1).x_value.airpuff(:,traces_to_remove)=[];
        data(1).x_value.galvano1(:,traces_to_remove)=[];
        data(1).x_value.galvano2(:,traces_to_remove)=[];
        data(1).x_value.EKG(:,traces_to_remove)=[];
        data(1).x_value.laser(:,traces_to_remove)=[];
        data(1).x_value.original_trace(:,traces_to_remove)=[];
  end   
  
    
                

path='D:\Inbal M.Sc\Data PhD\NB-ES Data\Extracted Data';
fname = Exp.namepath.name;
 name = [fname(1:(findstr(fname,'.')-1)) '_h1-2'];
 cd(path);
 save( name, 'data', 'headers');
cd(Exp.namepath.path);

end
