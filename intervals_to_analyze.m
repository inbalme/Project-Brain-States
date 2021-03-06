%intervals_to_analyze
                 switch trace_type
        case 1
            x_value=[1,1];
            if exp_type==3
                x_value=clamp_flag.*[1,1];
            end
        case 2
            x_value=[2:3]; %2:3; %[1,1];
            if exp_type==3
                x_value=(clamp_flag+3).*[1,1];
            end
        case 3
            x_value=[2,1]; %for spont. activity takes the "before" from x-value 2 and the "after" from x-value 1. enables taking longer interval           
                 end

clear  duration start_time start_sample end_sample interval interval_mat interval_plot x y patch_xdata patch_ydata yex ylim_data sem_xdata sem_ydata sem_cdata
coeffs=[]; 
switch exp_type
    case 1
         stim2{1}=stim2_X{2};
         stim2{2}=stim2_X{2};
        switch trace_type
            case 1
                 start_time = [0.4,5.6]; %[sec] %[0,5]
                 duration = 2.5; %[sec] %2.5
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
%                 duration = galvano_nstim./galvano_freq+0.1;
                duration = galvano_nstim./galvano_freq;
                if analyze_train_only_flag==1
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
           stim2{1}=stim2_X{2};
           stim2{2}=stim2_X{2};
        switch trace_type
            case 1
                 start_time=[0.4, stim1_X{x_value(1)}(1,1).*dt+0.4];
                 duration =(stim1_X{x_value(1)}(2,1)-stim1_X{x_value(1)}(1,1)).*dt-0.4;
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
%                 duration = galvano_nstim./galvano_freq+0.1;
                duration = galvano_nstim./galvano_freq;
                if analyze_train_only_flag==1
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
    case 3
           stim2{1}=stim2_X{4}(1:2,:);
           stim2{2}=stim2_X{4}(3:4,:);
           
                curr_inj_del = [Param.facade(27) Param.facade(30)] ;
            if length(Param.facade)>32
                curr_inj_del = [Param.facade(27) Param.facade(30) Param.facade(33)];
            end
            curr_inj_depo = Param.facade(1);
            curr_inj_hyper = Param.facade(2);
            curr_inj_dur = Param.facade(7);
          switch trace_type
            case 1   
                 start_time = curr_inj_del+0.1; %[sec] %[0,5]
        %          start_time(2) = stim2_X{x_value}(1,end).*dt+0.1; %start 0.5 sec after the last sensory stim
                 duration = 1.2; %curr_inj_dur-0.03; %[sec] 
                    for t=1:length(start_time);
                     start_sample(:,t) = ceil(start_time(t).*sf{1});
                      start_sample(mod(start_sample(:,t),2)==1,t)=start_sample(mod(start_sample(:,t),2)==1,t)-1;                 
                        if start_time(t)==0
                            start_sample(:,t) = 1;
                        end
                      end_sample(:,t) = ceil(start_sample(:,t)+duration.*sf{1});
                     interval_temp(:,t)= start_sample(:,t):end_sample(:,t);
                      if mod(size(interval_temp(:,t),1),2)
                           interval(:,t) =interval_temp(1:end-1,t);
                      else
                           interval(:,t) =interval_temp(:,t);
                      end
                    end 
              case 2
                   for t=1:length(x_value);
                    start_sample(:,t) =stim2{t}(1,1)-analyze_time_before_train.*sf{1};  %start 100ms before sensory stim
                    galvano_nstim = Param.facade(15);
                    galvano_freq = Param.facade(16);
%                     duration = galvano_nstim./galvano_freq+0.1;
                    duration = galvano_nstim./galvano_freq;
                        if analyze_train_only_flag==1
                            end_sample(:,t) = start_sample(:,t)+duration.*sf{1}-1;
                        else
                           end_sample(:,t) = stim2{t}(1,end)+0.1.*sf{1}-1;  
                        end
                    interval(:,t) = round(start_sample(:,t):end_sample(:,t));
                   end
          end      
      case 4
           stim2{1}=stim1_X{1}(1:2,:);
           stim2{2}=stim1_X{1}(1:2,:);
        switch trace_type
            case 1
                 start_time=[0.4,4.1]; %[sec]
                 duration = 2.5; %[sec] %2.5
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
%                 duration = galvano_nstim./galvano_freq+0.1;
                duration = galvano_nstim./galvano_freq;
                if analyze_train_only_flag==1
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