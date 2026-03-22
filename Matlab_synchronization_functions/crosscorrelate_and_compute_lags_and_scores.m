function lag = crosscorrelate_and_compute_lags_and_scores(lfp_sync_part, eeg_sync_part, sr_lfp, sr_eeg, plotting)
% CROSSCORRELATE_AND_COMPUTE_LAGS_AND_SCORES
% Computes lag between LFP (Percept) and EEG signals by cross-correlating each channel pair.
% Stores correlations in a cell array to handle different lengths.

    % Step 0: Resample LFP if sampling rates differ
    if sr_lfp ~= sr_eeg
        disp('The sampling rates between EEG and LFP do not match. Resampling LFP...');
        lfp_sync_part = rowwise_resample(lfp_sync_part, sr_lfp, sr_eeg);
    end

    eeg_sync_part1 = eeg_sync_part; % copy for alignment
    n_lfp = size(lfp_sync_part, 1);
    n_eeg = size(eeg_sync_part, 1);

    score = zeros(n_lfp, n_eeg);           % max correlation per channel pair
    scores = cell(n_lfp, n_eeg);           % store full correlation vector per pair
    lag_vals = zeros(n_lfp, n_eeg);        % lag in samples

    % Step 1: Compute normalized cross-correlations
    for it_lfp = 1:n_lfp
        a = lfp_sync_part(it_lfp, :);
        a = a / norm(a);

        for it_eeg = 1:n_eeg
            b = eeg_sync_part1(it_eeg, :);
            b = b / norm(b);

            correlation = xcorr(a, b);                     % full cross-correlation
            lags = -(length(b)-1):(length(a)-1);          % corresponding lags

            % Save results
            [score_max, idx_max] = max(correlation);
            score(it_lfp, it_eeg) = score_max;
            scores{it_lfp, it_eeg} = correlation;        % use cell array
            lag_vals(it_lfp, it_eeg) = lags(idx_max);
        end
    end

    % Step 2: Find best channel pair
    [~, ind_linear] = max(score(:));
    [best_lfp_chan, best_eeg_chan] = ind2sub(size(score), ind_linear);
    fprintf('The best LFP and EEG channels are %d and %d\n', best_lfp_chan, best_eeg_chan);

    lag = -1 * lag_vals(best_lfp_chan, best_eeg_chan); % lag definition
    fprintf('Lag: %d samples, maximum correlation: %.4f\n', lag, score(best_lfp_chan, best_eeg_chan));

    % Step 3: Align signals based on lag safely
    if lag > 0
        eeg_sync_part_new = eeg_sync_part1(:, lag+1:end);
        min_len = min(size(lfp_sync_part,2), size(eeg_sync_part_new,2));
        lfp_sync_part_new = lfp_sync_part(:, 1:min_len);
        eeg_sync_part_new = eeg_sync_part_new(:, 1:min_len);
    elseif lag < 0
        lfp_sync_part_new = lfp_sync_part(:, -lag+1:end);
        min_len = min(size(lfp_sync_part_new,2), size(eeg_sync_part1,2));
        eeg_sync_part_new = eeg_sync_part1(:, 1:min_len);
        lfp_sync_part_new = lfp_sync_part_new(:, 1:min_len);
    else
        min_len = min(size(lfp_sync_part,2), size(eeg_sync_part1,2));
        lfp_sync_part_new = lfp_sync_part(:, 1:min_len);
        eeg_sync_part_new = eeg_sync_part1(:, 1:min_len);
    end
    % Step 4: Optional plotting
    if plotting
        
        figure;
        plot(zscore(eeg_sync_part1(best_eeg_chan, :)), 'b', 'LineWidth', 1.2); hold on;
        plot(zscore(lfp_sync_part(best_lfp_chan, :)), 'r', 'LineWidth', 1.2);
        title('EEG and LFP before aligning');
        legend('EEG', 'LFP');

        figure;
        plot(zscore(eeg_sync_part_new(best_eeg_chan, :)), 'b', 'LineWidth', 1.2); hold on;
        plot(zscore(lfp_sync_part_new(best_lfp_chan, :)), 'r', 'LineWidth', 1.2);
        title('EEG and LFP after aligning');
        legend('EEG', 'LFP');

        % Plot cross-correlation scores of the best channel pair
        best_corr = scores{best_lfp_chan, best_eeg_chan};
        len_a = length(lfp_sync_part(best_lfp_chan,:));
        len_b = length(eeg_sync_part1(best_eeg_chan,:));
        
        % Compute lag vector
        
        lags= -(len_b-1):(len_a-1);
        plot_lags=-lags
        figure;
        plot(plot_lags, best_corr, 'k', 'LineWidth', 1.2);
        xlabel('Lag (samples)');
        ylabel('Correlation');
        title('Cross-correlation scores of best LFP and EEG channel pair');
        grid on;

% Highlight 0 lag
hold on;
xline(0,'r--','LineWidth',1.5);
        figure;
        histogram(scores{best_lfp_chan, best_eeg_chan}, 30, 'FaceColor', [0.53 0.81 0.98]);
        hold on;
        xline(mean(scores{best_lfp_chan, best_eeg_chan}), 'r', 'LineWidth', 2);
        xline(mean(scores{best_lfp_chan, best_eeg_chan}) + 3*std(scores{best_lfp_chan, best_eeg_chan}), 'k--', 'LineWidth', 1.5);
        title('Distribution of cross-correlation scores of best channel pair');
    end
end