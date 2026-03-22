function upsample_data = rowwise_resample(data, sr_old, sr_new)
% ROWWISE_RESAMPLE Resamples input data row-wise to a given sampling rate
%
%   data: 1D vector or 2D matrix (rows = channels, columns = samples)
%   sr_old: original sampling rate
%   sr_new: new sampling rate
%
%   Returns:
%       upsample_data: resampled data (same shape as input but new sample length)

    if isvector(data)
        % 1D data
        N_old = length(data);
        N_new = round(N_old * sr_new / sr_old);
        upsample_data = resample(data, N_new, N_old);

    elseif ismatrix(data)
        % 2D data (rows = channels)
        n_channels = size(data,1);
        N_old = size(data,2);
        N_new = round(N_old * sr_new / sr_old);
        upsample_data = zeros(n_channels, N_new);

        for ch = 1:n_channels
            upsample_data(ch,:) = resample(data(ch,:), N_new, N_old);
        end

    else
        error('rowwise_resample: input must be 1D or 2D numeric array');
    end
end