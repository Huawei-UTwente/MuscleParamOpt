function fig = reflexParPlot(gains_fse, gains_lce, res, tPhase, M, hostMus, refMus)
%UNTITLED Plot the reflex control parameters of one muscle

    fig = figure();
    lw = 1;
    pw = 12;
    lcolor = 'k';
    bcolor = 'k';
    
    fig.Position = [10 10 pw*200 M*200];

    res1_ind = [];
    res2_ind = [];
    
    for m = 1:M

        fse_st = (m - 1)*pw + 1;
        lce_st = (M + m)*pw + 1;
        
        res1_ind = [res1_ind, fse_st + pw - 2];
        res2_ind = [res2_ind, lce_st + pw - 2];
        
        subplot(2*M+1, pw, [fse_st:fse_st+(pw-3)])
        boxplot(gains_fse(:, :, m))
%         hold on
%         plot(1:tPhase, mean(gains_fse(:,:, m), 1), lcolor, 'linewidth', lw);
        ylabel(refMus(m))
        % ylim([-1, 1])
        if m == 1
            title('Force reflex')
        end
        set(gca, 'XTickLabel', [])

        subplot(2*M+1, pw, lce_st:lce_st+(pw-3))
        boxplot(gains_lce(:, :, m))
%         hold on
%         plot(1:tPhase, mean(gains_lce(:, :, m), 1), lcolor, 'linewidth', lw);
        ylabel(refMus(m))
        % ylim([-1, 1])
        if m == M
            xlabel('Gait phase')
        else
            set(gca, 'XTickLabel', [])
        end
        if m == 1
            title('Length reflex')
        end
        
    end
    
    subplot(2*M+1, pw, res1_ind)
    boxplot(res(:, 1))
    % bar(res(1), bcolor);
    ylim([0, 1.5])
    xlabel('fse offset')
    set(gca, 'XTickLabel', [])

    subplot(2*M+1, pw, res2_ind)
    boxplot(res(:, 2))
    % bar(res(2), bcolor);
    ylim([0.5, 1.5])
    xlabel('lce offset')
    set(gca, 'XTickLabel', [])

    subplot(2*M+1, pw, (1:2*M+1)*pw)
    boxplot(res(:, 3))
    % bar(res(3), bcolor)
    ylim([0, 1])
    xlabel('activation offset')
    set(gca, 'XTickLabel', [])
    
    sgtitle(hostMus)
    
end

