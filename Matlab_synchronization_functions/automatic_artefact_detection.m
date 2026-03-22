function [filtered_data, arr] = automatic_artefact_detection(data, sr, start_freq, end_freq, sig_val, rel_chan, len_art)
% AUTOMATIC_ARTEFACT_DETECTION Automatically detects artifacts in LFP/EEG data
%
% INPUTS:
%   data      : channels x samples
%   sr        : sampling rate
%   start_freq: lower bound of artifact frequency
%   end_freq  : upper bound of artifact frequency
%   sig_val   : threshold in std deviations (default = 1)
%   rel_chan  : channel index for detection (default = 1)
%   len_art   : sliding window length in seconds (default = 0.1)
%
% OUTPUTS:
%   filtered_data : filtered data in the specified frequency band
%   arr           : binary vector indicating artifact points (1 = artifact)

    % Step 1: Filter data in the specified frequency range
    filtered_data = filter_frequency(data, sr, start_freq, end_freq);

    % Step 2: Compute global mean and std on the relevant channel
    channel_data = filtered_data(rel_chan, :);
    ges_mean = mean(channel_data);
    ges_std  = std(channel_data);

    % Step 3: Sliding window detection
    N_win = round(len_art * sr);  % window length in samples
    arr = zeros(1, length(channel_data));

    for i = 1:length(channel_data) - N_win
        current_part = channel_data(i:i+N_win-1);
        current_mean = mean(current_part);
        if current_mean > (ges_mean + sig_val * ges_std)
            arr(i) = 1;
        end
    end

    % Step 4: Plot z-scored filtered data and artifacts
    figure;
    plot(zscore(channel_data), 'b'); hold on;
    plot(arr * max(zscore(channel_data)), 'r'); % scale arr for visibility
    xlabel('Samples');
    ylabel('Amplitude (z-scored)');
    title(['Automatic Artifact Detection - Channel ', num2str(rel_chan)]);
    legend('Filtered Data (z-score)', 'Detected Artifacts');

end