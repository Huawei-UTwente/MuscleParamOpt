function fig = computationalTimePlot(homePath,subj, trialNamesOpt, ...
                                     trialType, Phase_list, W2)
    % computationalTimePlot: function to plot computational time of the reflex
    % identification trials. 
    % by: Huawei Wang
    % date: May 16th 2022

    T = length(trialNamesOpt);
    W3 = 0;
    
    Time_list = zeros(100, length(Phase_list));
        
    for  p = 1:length(Phase_list)
        
        tPhase = Phase_list(p);
        
        % optimization saving path
        save_path = sprintf('%s/oriOpt/Subj%02d/WalkComb%d/%dPhases/trial_%s', ...
                            homePath, subj, T, tPhase, trialType);
                        
        folder = sprintf('%s/W2_%d_W3_%d', save_path, W2, W3);
        
        for opt = 1:100
            saving_names = sprintf('%s/optimization_res%02d.mat', folder, opt);
            res = load(saving_names);
            Time_list(opt, p) = res.time;
        end
        
    end
    
    fig = boxplot(Time_list./60);
    ylabel('Minutes')
    xticks(1:length(Phase_list));
    xticklabels(Phase_list);
    xlabel("Phase division")
    title('Computation Time')

end

