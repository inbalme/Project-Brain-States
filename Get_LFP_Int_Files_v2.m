function files = Get_LFP_Int_Files_v2()
% list of file names Headers and directories 
%cd 'D:\Inbal M.Sc\Data PhD\Awake LFP Data\Extracted Data'
% called by Extract_Awake_Data_v2
% bulid to deal with files in which every header is one repetition over the
% protocol traces, such that there will be an equal number of repetitions
% from each trace type.
% protocols: 1- light only 1 pulse 20ms, 2- light only 1 pulse 1sec, 3- light only 20 pulses 20 Hz 20ms.   
%                 4 - light+galvano (3 x-values): light 1 pulse 1 sec, galvo 10 Hz 10 pulses
%                 5 - light+galvano(5 x-values): light  1 pulse, galvo 1 pulse, light train, galvo train.
%                 6 - light+galvo+current injections (9 x-values)
%                 7 - ChAT: Stim+light (yonatan's protocol). 4 x-values:
%                 light train, single galvo, light train+galvo+depo, single
%                 galvo+depo
%                 8 - ChAT intracellular: light only (3 x-values): 1) train 20 Hz, 10ms ondur, 50 pulses. 2) single pulse 10ms ondur. 3) single pulse 5 sec
%                 9 - IC: I-V curve
%


files(1).name = 'Set2-2014-07-14-001.exp2';
files(1).path = 'D:\Inbal M.Sc\Data PhD\Awake LFP Data\Raw Data';
files(1).extracted_name = 'File_1';
files(1).extracted_path = 'D:\Inbal M.Sc\Data PhD\Awake LFP Data\Extracted Data';
files(1).animal = 'ChAT-ChR2(+/+) BD18.5.14 not responding to light';   
files(1).anesthesia_initial = 'ket+very low xyl';
files(1).anesthesia_maintenance = 'halothane';
files(1).pipette_solution = 'patch+bio+QX';
files(1).recording_type = 'whole cell';        %LFP, Cell attached, whole cell, extracellular
files(1).recording_location = 'Barrel cortex';
files(1).electrode_resist = [];
files(1).depth = 650;
files(1).depth_2 = 350;
files(1).pr_whisker = [];
files(1).adj_whisker = [];                                            
files(1).headers = 2;
files(1).protocol = [];
files(1).stim_type = 1;                        % 1=spont , 2=light alone, 3=puff alone, 4=puff+light, 5=galvano alone, 6=galvano+light, 7=protocol light only, 8=protocol galvano
files(1).EEG_flag = 0;
files(1).dual_rec_flag = 1;
files(1).channel_Vm = 1; 
files(1).channel_V2 = 3; 
files(1).channel_I1 = 2; 
files(1).channel_galvano = 15;
files(1).channel_airpuff = 13;
files(1).channel_EKG = 6;
files(1).channel_EEG = 5;
files(1).channel_laser = 10;                %optogenetic stimulation
files(1).galvano = 2;
files(1).comments = [];
files(1).offset = 150;                                                   %need to subtract the offset

files(2).name = 'Set2-2014-07-14-001.exp2';
files(2).path = 'D:\Inbal M.Sc\Data PhD\Awake LFP Data\Raw Data';
files(2).extracted_name = 'File_2';
files(2).extracted_path = 'D:\Inbal M.Sc\Data PhD\Awake LFP Data\Extracted Data';
files(2).animal = 'ChAT-ChR2(+/+) BD18.5.14 not responding to light';   
files(2).anesthesia_initial = 'ket+very low xyl';
files(2).anesthesia_maintenance = 'halothane';
files(2).pipette_solution = 'patch+bio+QX';
files(2).recording_type = 'whole cell';        %LFP, Cell attached, whole cell, extracellular
files(2).recording_location = 'Barrel cortex';
files(2).electrode_resist = [];
files(2).depth = 650;
files(2).depth_2 = 350;
files(2).pr_whisker = [];
files(2).adj_whisker = [];                                            
files(2).headers = 3:7;
files(2).protocol = 9;
files(2).stim_type = 7;                        % 1=spont , 2=light alone, 3=puff alone, 4=puff+light, 5=galvano alone, 6=galvano+light, 7=protocol light only, 8=protocol galvano
files(2).EEG_flag = 0;
files(2).dual_rec_flag = 1;
files(2).channel_Vm = 1; 
files(2).channel_V2 = 3; 
files(2).channel_I1 = 2; 
files(2).channel_galvano = 15;
files(2).channel_airpuff = 13;
files(2).channel_EKG = 6;
files(2).channel_EEG = 5;
files(2).channel_laser = 10;                %optogenetic stimulation
files(2).galvano = 2;
files(2).comments = [];
files(2).offset = 150;                                                   %need to subtract the offset

files(3).name = 'Set2-2014-07-14-001.exp2';
files(3).path = 'D:\Inbal M.Sc\Data PhD\Awake LFP Data\Raw Data';
files(3).extracted_name = 'File_2';
files(3).extracted_path = 'D:\Inbal M.Sc\Data PhD\Awake LFP Data\Extracted Data';
files(3).animal = 'ChAT-ChR2(+/+) BD18.5.14 not responding to light';   
files(3).anesthesia_initial = 'ket+very low xyl';
files(3).anesthesia_maintenance = 'halothane';
files(3).pipette_solution = 'patch+bio+QX';
files(3).recording_type = 'whole cell';        %LFP, Cell attached, whole cell, extracellular
files(3).recording_location = 'Barrel cortex';
files(3).electrode_resist = [];
files(3).depth = 650;
files(3).depth_2 = 350;
files(3).pr_whisker = [];
files(3).adj_whisker = [];                                            
files(3).headers = 3:7;
files(3).protocol = 9;
files(3).stim_type = 7;                        % 1=spont , 2=light alone, 3=puff alone, 4=puff+light, 5=galvano alone, 6=galvano+light, 7=protocol light only, 8=protocol galvano
files(3).EEG_flag = 0;
files(3).dual_rec_flag = 1;
files(3).channel_Vm = 1; 
files(3).channel_V2 = 3; 
files(3).channel_I1 = 2; 
files(3).channel_galvano = 15;
files(3).channel_airpuff = 13;
files(3).channel_EKG = 6;
files(3).channel_EEG = 5;
files(3).channel_laser = 10;                %optogenetic stimulation
files(3).galvano = 2;
files(3).comments = [];
files(3).offset = 150;                                                   %need to subtract the offset


cd 'D:\inbal M.Sc\Data PhD\ChAT Data\Extracted Data';
 
save('Awake LFP_Files_v2')