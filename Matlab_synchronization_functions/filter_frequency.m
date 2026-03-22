
function data_filtered=filter_frequency(data, sr, low_freq,high_freq)
    [b,a] = butter(3, [low_freq high_freq] / (sr/2), 'bandpass');
    filtered = filtfilt(b,a,data')';   % filter along time
    data_filtered = abs(hilbert(filtered')');
end