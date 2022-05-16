% plot the experimental data
function [h1, h2] = experimental_data_plot(moment, mus_act, lmt, ma)
    
    global T J L M N coordinateNames muscleNames;
    
    % plot joint torques
    h1 = figure();        
    subj_h1 = ceil(sqrt(J));        
    for j = 1:J
        subplot(subj_h1, subj_h1, j)
        for l = 1:L
            for t = 1:T
                plot((l-1)*N + 1:(l-1)*N + N, moment(:, L*J*(t-1) + (l-1)*J + j),...
                    '-', 'linewidth', 2)
                hold on
                if l == 1
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    'k-', 'linewidth', 1.5)
                elseif l == 2
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    '-', 'linewidth', 1.5, 'color', [1, 1, 1]*0.75)
                end

                title(coordinateNames(j))
                ylabel('joint torque, Nm/kg')
            end
        end
    end
    
    % plot muscle activations
    h2 = figure();        
    subj_h2 = ceil(sqrt(M));        
    for m = 1:M
        subplot(subj_h2, subj_h2, m)
        for l = 1:L
            for t = 1:T
                plot((l-1)*N + 1:(l-1)*N + N, mus_act(:, L*M*(t-1) + (l-1)*M + m),...
                    '-', 'linewidth', 2)
                hold on
                if l == 1
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    'k-', 'linewidth', 1.5)
                elseif l == 2
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    '-', 'linewidth', 1.5, 'color', [1, 1, 1]*0.75)
                end
            end
        end
        title(muscleNames(m))
        ylabel('Muscle activation')
    end
    
    % plot muscle length
    h3 = figure();        
    subj_h3 = ceil(sqrt(M));        
    for m = 1:M
        subplot(subj_h3, subj_h3, m)
        for l = 1:L
            for t = 1:T
                plot((l-1)*N + 1:(l-1)*N + N, lmt(:, L*M*(t-1) + (l-1)*M + m),...
                    '-', 'linewidth', 2)
                hold on
                if l == 1
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    'k-', 'linewidth', 1.5)
                elseif l == 2
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    '-', 'linewidth', 1.5, 'color', [1, 1, 1]*0.75)
                end
            end
        end
        title(muscleNames(m))
        ylabel('Muscle length, m')
    end

    % plot muscle moment arms
    for j = 1:J
    figure();        
    subj_h = ceil(sqrt(M));        
    for m = 1:M
        subplot(subj_h, subj_h, m)
        for l = 1:L
            for t = 1:T
                plot((l-1)*N + 1:(l-1)*N + N, ma(:, L*M*(t-1) + (l-1)*M + m),...
                    '-', 'linewidth', 2)
                hold on
                if l == 1
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    'k-', 'linewidth', 1.5)
                elseif l == 2
                    plot([(l-1)*N + 1, (l-1)*N + N], [0, 0],...
                    '-', 'linewidth', 1.5, 'color', [1, 1, 1]*0.75)
                end
            end
        end
        title(muscleNames(m))
        ylabel('Muscle length, m')
    end
    end
end