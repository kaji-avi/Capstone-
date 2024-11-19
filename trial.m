
seq_length = 64;         % Length of ZC sequence (prime number)
root = 1;               % Root of ZC sequence
SNR_dB = 10;             % Signal-to-noise ratio (in dB)
cp_length = 16;          % Cyclic prefix length
CFO = 20;              % Carrier Frequency Offset (in Hz)
fs = 10000;              % Sampling frequency (in Hz)

tx_signal = generate_signal(seq_length, root, cp_length);  % Generate transmitted signal (ZC sequence with CP and data)
rx_signal = AWGN_with_CFO(tx_signal, SNR_dB, CFO, fs);   % Generate recieved signal with CFO, delay, and noise

% Call synchronization method to get the time and frequency estimate along
% with the correlation output
[time_estimate, corr_output, cfo_estimate] = synchronization_with_CFO(rx_signal, tx_signal, cp_length, seq_length);
cfo_estimate_hz = cfo_estimate * fs;

% Plot Transmitted, Received Signal, and Correlation
my_plot(tx_signal, rx_signal, corr_output);
% Print start position of detected signal
fprintf('Start of received signal detected at sample index: %d\n', time_estimate)

fprintf('Estimated CFO in Hz: %.2f Hz\n', cfo_estimate_hz);

