%% repeated measures ANOVA for a response parameter (e.g. amplitude) over the train stim.
% rmANOVA with matlab function
function[stat]=fn_sensory_train_response_2_way_rmANOVA(data_factorA1, data_factorA2)
%data_factorA1 is a matrix where columns are the repeated measures and rows
%are the value for each subject (values along the row are within-subject),
%where factor A =1. data_factorA2 is the same but factor A=2. for examle,
%before and after electrical stimulation.
d1=num2cell(data_factorA1);
d2=num2cell(data_factorA2);
 S_1(:,1)=1:size(data_factorA1,1); %S_1 is the numbers of subjects
 s=num2cell(S_1);
 
ta_vector_names={'data_factorA1','data_factorA2'};
varnames{1}='cells';
for i=2:2*size(data_factorA1,2)+1;
    varnames{i}=['measure',num2str(i-1)];
end

d=[s,d1,d2];
ta=cell2table(d);
ta.Properties.VariableNames=varnames;
factorNames = {'NB','Time'};
for j=1:size(data_factorA1,2);
    factorA1{j,1}='N';
    factorA2{j,1}='Y';
    factorRM{j,1}=num2str(j);
end
factors=[[factorA1; factorA2] [factorRM;factorRM]];
within=table([factorA1; factorA2], [factorRM;factorRM],'VariableNames',factorNames);
% fit the repeated measures model
rm_model_string=[varnames{2},'-',varnames{end},'~1'];
rm = fitrm(ta,rm_model_string,'WithinDesign',within);
% rm = fitrm(ta,'Y1-Y22~1','WithinDesign',within);
%test for sphericity
mauchly_sphericity=mauchly(rm);
% run my repeated measures anova here
[ranovatbl] = ranova(rm,'withinmodel','NB*Time');
%multiple comparisons:
NBbyTime = multcompare(rm,'NB','By','Time');
TimebyNB = multcompare(rm,'Time','By','NB');
eps_amp = epsilon(rm);
eps=table2cell(eps_amp(1,2));
if eps{1}<0.75
    rmANOVA_NB_effect_p= table2cell(ranovatbl(3,6));
    rmANOVA_Time_effect_p= table2cell(ranovatbl(5,6));
    rmANOVA_interaction_effect_p= table2cell(ranovatbl(7,6));
else
     rmANOVA_NB_effect_p= table2cell(ranovatbl(3,7));
     rmANOVA_Time_effect_p= table2cell(ranovatbl(5,7));
     rmANOVA_interaction_effect_p= table2cell(ranovatbl(7,7));
end
table_column1=table2cell(NBbyTime(:,1));

for stim=1:size(data_factorA1,2) 
    row= find(strcmp(num2str(stim),table_column1),1);
    p_stim{1}(stim,:)=table2cell(NBbyTime(row,6));
end
    
stat.table=ta;
stat.table_data_vecs=ta_vector_names;
stat.within_design=within;
stat.rm=rm;
stat.mauchly=mauchly_sphericity;
stat.eps_all=eps_amp;
stat.eps=eps;
stat.ANOVA=ranovatbl;
stat.multcomp_NBbyTime=NBbyTime;
stat.multcomp_TimebyNB=TimebyNB;
stat.rmANOVA_NB_effect_p= rmANOVA_NB_effect_p;
stat.rmANOVA_Time_effect_p= rmANOVA_Time_effect_p;
stat.rmANOVA_interaction_effect_p= rmANOVA_interaction_effect_p;
stat.p_stim=p_stim{1};
