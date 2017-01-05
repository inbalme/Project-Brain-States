%intervals_to_analyze

 switch trace_type
        case 1
            x_value=[1,1];
        case 2
            x_value=[2:3]; %2:3; %[1,1];
        case 3
            x_value=[2,1]; %for apont. activity takes the "before" from x-value 2 and the "after" from x-value 1. enables taking longer interval
end
clear  duration start_time start_sample end_sample interval interval_mat interval_plot x y patch_xdata patch_ydata yex ylim_data sem_xdata sem_ydata sem_cdata
coeffs=[]; 
switch exp_type
    case 1
        switch trace_type
            case 1
                 start_time = [0.4,5.6]; %[sec] %[0,5]
                 duration = 2.5; %[sec] 
                    for t=1:length(start_time);
                     start_sample(:,t) = ceil(start_time(t).*sf{1});
                        if start_time(t)==0
                            start_sample(:,t) = 1;
                        end
                      end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
                      interval(:,t) = start_sample(:,t):end_sample(:,t);
                    end 
            case 2
                start_sample =stim2_X{x_value(1)}(1,1)-analyze_time_before_train.*sf{1};  %start 100ms before sensory stim
                galvano_nstim = Param.facade(6);
                galvano_freq = Param.facade(7);
                duration = galvano_nstim./galvano_freq+0.1;
                if analyze_train_only_flag==0
                    end_sample = start_sample+duration.*sf{1}-1;
                else
                   end_sample = stim2_X{x_value(1)}(1,end)+0.1.*sf{1}-1;  
                end
                interval(:,1) = round(start_sample:end_sample);
                interval(:,2) = interval(:,1);
            case 3
                 start_time = [0.4,5.6]; %[sec] %[0,5]
                 duration = 3; %[sec] 
                    for t=1:length(start_time);
                     start_sample(:,t) = ceil(start_time(t).*sf{1});
                        if start_time(t)==0
                            start_sample(:,t) = 1;
                        end
                        end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
                        interval(:,t) = start_sample(:,t):end_sample(:,t);
                    end
            end
    case 2
        switch trace_type
            case 1
                 start_time=[0.4, stim1_X{x_value(1)}(1,1).*dt+0.4];
                 duration =(stim1_X{x_value(1)}(2,1)-stim1_X{x_value(1)}(1,1)).*dt;
                    for t=1:length(start_time);
                     start_sample(:,t) = ceil(start_time(t).*sf{1});
                        if start_time(t)==0
                            start_sample(:,t) = 1;
                        end
                      end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
                      interval(:,t) = start_sample(:,t):end_sample(:,t);
                    end 
            case 2
                start_sample =stim2_X{x_value(1)}(1,1)-analyze_time_before_train.*sf{1};  %start 100ms before sensory stim
                galvano_nstim = Param.facade(6);
                galvano_freq = Param.facade(7);
                duration = galvano_nstim./galvano_freq+0.1;
                if analyze_train_only_flag==0
                    end_sample = start_sample+duration.*sf{1}-1;
                else
                   end_sample = stim2_X{x_value(1)}(1,end)+0.1.*sf{1}-1;  
                end
                interval(:,1) = round(start_sample:end_sample);
                interval(:,2) = interval(:,1);
                
           case 3
                 start_time = [0.4,stim1_X{x_value(1)}(1,1).*dt+0.4]; %[sec] %[0,5]
                 duration = 3; %[sec] 
                    for t=1:length(start_time);
                     start_sample(:,t) = ceil(start_time(t).*sf{1});
                        if start_time(t)==0
                            start_sample(:,t) = 1;
                        end
                        end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1})-1;
                        interval(:,t) = start_sample(:,t):end_sample(:,t);
                    end
        end
end