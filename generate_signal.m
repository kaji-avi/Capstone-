function tx_signal_plain = generate_signal(seq_length, root, cp_length) % Function to generate ZC sequence with CP and data
    % Generate Zadoff-Chu (ZC) sequence
    n = 0:seq_length-1;
    zc_sequence = exp(-1i*pi*root*n.*n/seq_length);  % ZC sequence formula
    % Add Cyclic Prefix (CP) to ZC sequence
    zc_with_cp = [zc_sequence(end-cp_length+1:end), zc_sequence];
    data_val = randi([1, 100]);
    tx_signal_plain = [zc_with_cp, zeros([1, data_val])]; % ZC with cp and random length zero data
end


