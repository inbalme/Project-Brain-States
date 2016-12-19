%% SSA 

%for short ISI

 load Set2-2015-01-28-001_header7

 dt=1/10000;
 color_table = [0 0 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 0 0 0];
 
 mean_trace_standard_1 =  mean(data.x_value(1).data1(16001:19500,:),2);
 mean_trace_deviant_1 =  mean(data.x_value(2).data1(20001:23500,:),2);
 mean_trace = [mean_trace_standard_1 mean_trace_deviant_1];
%  [Fig,h] = fn_Plot_Trace_v2(mean_trace, 1/10000);
% title(gca, ['header ',num2str(data.x_value(1).original_header(end)),' whisker 1']);
figure;
subplot (2,2,1)
hold on
for i = 1:size(mean_trace,2);
   plot((1:size(mean_trace,1)).*dt,mean_trace(:,i), 'LineWidth',0.2,'color', color_table(i,:));
end
hold off

 mean_trace_standard_2 =  mean(data.x_value(2).data1(16001:19500,:),2);
 mean_trace_deviant_2 =  mean(data.x_value(1).data1(20001:23500,:),2);
 mean_trace = [mean_trace_standard_2 mean_trace_deviant_2];
%  [Fig,h] = fn_Plot_Trace_v2(mean_trace, 1/10000);
% title(gca, ['header ',num2str(data.x_value(1).original_header(end)),' whisker 2']);

gcf
subplot (2,2,2)
hold on
for i = 1:size(mean_trace,2);
   plot((1:size(mean_trace,1)).*dt,mean_trace(:,i), 'LineWidth',0.2,'color', color_table(i,:));
end
hold off

 
%% for long ISI

load Set2-2015-01-28-001_header9

 mean_trace_standard_1 =  mean(data.x_value(1).data1(60001:64500,:),2);
 mean_trace_deviant_1 =  mean(data.x_value(2).data1(75001:79500,:),2);
 mean_trace = [mean_trace_standard_1 mean_trace_deviant_1];
%  [Fig,h] = fn_Plot_Trace_v2(mean_trace, 1/10000);
%  title(gca, ['header ',num2str(data.x_value(1).original_header(end)),' whisker 1']);

gcf
subplot (2,2,3)
hold on
for i = 1:size(mean_trace,2);
   plot((1:size(mean_trace,1)).*dt,mean_trace(:,i), 'LineWidth',0.2,'color', color_table(i,:));
end
hold off

 mean_trace_standard_2 =  mean(data.x_value(2).data1(60001:64500,:),2);
 mean_trace_deviant_2 =  mean(data.x_value(1).data1(75001:79500,:),2);
 mean_trace = [mean_trace_standard_2 mean_trace_deviant_2];
%  [Fig,h] = fn_Plot_Trace_v2(mean_trace, 1/10000);
%  title(gca, ['header ',num2str(data.x_value(1).original_header(end)),' whisker 2']);

gcf
subplot (2,2,4)
hold on
for i = 1:size(mean_trace,2);
   plot((1:size(mean_trace,1)).*dt,mean_trace(:,i), 'LineWidth',0.2,'color', color_table(i,:));
end
hold off

 