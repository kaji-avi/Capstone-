function [c, lags] = corr(rx_signal, tx_signal)
    % Get lengths of signals
    len_rx = length(rx_signal);
    len_tx = length(tx_signal);
    
    % Determine range of lags
    maxLag = len_rx + len_tx - 1;
    lags = -(len_tx-1):(len_rx-1);
    
    % Pre-allocate correlation vector
    c = zeros(1, maxLag);
    
    % Compute correlation for each lag
    for i = 1:maxLag
        lag = lags(i);
        if lag < 0
            c(i) = sum(conj(rx_signal(1:len_tx+lag)).*tx_signal(-lag+1:end));
        else
            c(i) = sum(conj(rx_signal(lag+1:len_tx)).*tx_signal(1:len_tx-lag));
        end
    end
    
    % Normalize correlation
    c = c / sqrt(sum(abs(rx_signal).^2) * sum(abs(tx_signal).^2));
end
