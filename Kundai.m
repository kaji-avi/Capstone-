clc
clear
close all
%% Transmitter

% Set a random seed for reproducibility
rng(42); 

% Parameters
N = 64;  % Number of OFDM subcarriers
M = 16;  % 16-QAM
numBits = N * log2(M);  % Total number of bits
bitsPerSymbol = log2(M);  % Bits per symbol (4 bits for 16-QAM)
cyclicPrefixLength = N * 0.25;


% % Roberto' Parameters
SNR = 40; 
to = 0;
h = 1;


dataBits = randi([0, 1], numBits, 1);  % Random bit stream

% QAM Mapping (Gray Mapping for 16-QAM)
numSymbols = numBits / bitsPerSymbol;  % Number of QAM symbols
qamSymbols = zeros(numSymbols, 1);  % Initialize QAM symbols
for i = 1:numSymbols
    bits = dataBits((i-1)*bitsPerSymbol + (1:bitsPerSymbol))';  % Extract 4 bits
    grayIndex = bi2de(bits, 'left-msb');  
    qamSymbols(i) = gray_to_qam(grayIndex, M); 
end

% OFDM Modulation (TX): Serial to Parallel, IFFT
parallelData = reshape(qamSymbols, N, []);  % Reshape to N subcarriers x symbols
timeDomainSymbols = ifft_function(parallelData, N);
timeDomainSymbols = timeDomainSymbols* sqrt(N);

% Add Cyclic Prefix (CP)
timeDomainSymbolsWithCP = [timeDomainSymbols(end-cyclicPrefixLength+1:end, :); timeDomainSymbols];
x = timeDomainSymbolsWithCP(:);  % Convert to serial for transmission


% Pass through the channel (provided function)
c = channelEmulation(x, SNR, to, h);

%% Receiver

% Remove CP, Serial to Parallel, FF
receivedParallelData = reshape(c, N + cyclicPrefixLength, []);  % Serial to parallel conversion
receivedSymbolsWithoutCP = receivedParallelData(cyclicPrefixLength+1:end, :);  % Remove CP
receivedFrequencyDomainSymbols = fft_function(receivedSymbolsWithoutCP, N);  % FFT
receivedFrequencyDomainSymbols= receivedFrequencyDomainSymbols/ sqrt(N); % normalization

% Equalization Using Known Channel (h_hat)
h_hat = fft(h, N) / sqrt(N);  % Take FFT with normalization
signalPower = mean(abs(timeDomainSymbols(:)).^2);  % Signal power
noiseVariance = signalPower / (10^(SNR / 10));  % Noise variance
equalizedSymbols = mmse_equalization(receivedFrequencyDomainSymbols, h_hat, noiseVariance, signalPower);

% QAM Demapping: Convert Equalized Symbols Back to Bits
receivedBits = zeros(numBits, 1);
for i = 1:numSymbols
    receivedBits((i-1)*bitsPerSymbol + (1:bitsPerSymbol)) = qam_to_gray(equalizedSymbols(i), M);
end

% Calculate Bit Error Rate (BER)
numErrors = sum(dataBits ~= receivedBits);  % Count bit errors
ber = numErrors / numBits;  % Compute BER
fprintf('Bit Error Rate (BER): %f\n', ber);

% Plot transmitted and received constellations
scatterplot(qamSymbols);  % Plot original transmitted symbols
title('Transmitted QAM Constellation');

scatterplot(equalizedSymbols(:));  % Plot received symbols after equalization
title('Received Constellation after Equalization');

% Debugging
% After OFDM Modulation
txPowerBeforeCP = mean(abs(timeDomainSymbols(:)).^2);
fprintf('Power before CP: %f\n', txPowerBeforeCP);

% After Adding CP
txPowerAfterCP = mean(abs(timeDomainSymbolsWithCP(:)).^2);
fprintf('Power after CP: %f\n', txPowerAfterCP);

% Received Signal Power (Before Equalization)
rxPowerBeforeEQ = mean(abs(receivedFrequencyDomainSymbols(:)).^2);
fprintf('Received Power before EQ: %f\n', rxPowerBeforeEQ);

% Received Signal Power (After Equalization)
rxPowerAfterEQ = mean(abs(equalizedSymbols(:)).^2);
fprintf('Received Power after EQ: %f\n', rxPowerAfterEQ);




function c = channelEmulation(x,SNR,to,h)

% c = zeros(1,to); % create some offset
% c = [c conv(h,x)];
% 
% % scale with SNR and multiply CFO
% c = sqrt(snr).*c;
% 
% % add awgn
% c = c+sqrt(1/2).*(randn(size(c))+1i.*randn(size(c)));

noiseVariance = 1 / (10^(SNR / 10));
c = conv(h, x) + sqrt(noiseVariance / 2) * (randn(size(x)) + 1i * randn(size(x)));

end

% Function for MMSE equalization
function equalizedSymbols = mmse_equalization(receivedSymbols, H, noiseVariance, signalPower)
    mmseFilter = conj(H) ./ (abs(H).^2 + noiseVariance / signalPower);
    equalizedSymbols = receivedSymbols .* mmseFilter;
    equalizedSymbols = equalizedSymbols / sqrt(mean(abs(mmseFilter).^2));
end


% Function for IFFT calculation
function timeDomainSymbols = ifft_function(parallelData, Nsubcarriers)
    [Nrows, Ncols] = size(parallelData);
    if Nrows ~= Nsubcarriers
        error('Number of rows in parallelData must match Nsubcarriers');
    end
    timeDomainSymbols = zeros(Nsubcarriers, Ncols);
    for col = 1:Ncols
        X = parallelData(:, col);
        N = length(X);
        x = zeros(N, 1);
        for n = 0:N-1
            for k = 0:N-1
                x(n+1) = x(n+1) + X(k+1) * exp(1i * 2 * pi * k * n / N);
            end
        end
        timeDomainSymbols(:, col) = x / N;
    end
    % Normalization below ensures that the power of the signal remains 
    % consistent during the transformation from the frequency domain to the time domain.
    timeDomainSymbols = timeDomainSymbols * sqrt(Nsubcarriers);
end

% Function for FFT calculation
function receivedFrequencyDomainSymbols = fft_function(receivedSymbolsWithoutCP, Nsubcarriers)
    [Nrows, Ncols] = size(receivedSymbolsWithoutCP);
    if Nrows ~= Nsubcarriers
        error('Number of rows in receivedSymbolsWithoutCP must match Nsubcarriers');
    end
    receivedFrequencyDomainSymbols = zeros(Nsubcarriers, Ncols);
    for col = 1:Ncols
        x = receivedSymbolsWithoutCP(:, col);
        N = length(x);
        X = zeros(N, 1);
        for k = 0:N-1
            for n = 0:N-1
                X(k+1) = X(k+1) + x(n+1) * exp(-1i * 2 * pi * k * n / N);
            end
        end
        receivedFrequencyDomainSymbols(:, col) = X;
    end
end

% Function for Gray to QAM Mapping (16-QAM)
function qamSymbol = gray_to_qam(grayIndex, M)
    if M == 16
        % Predefined QAM constellation points for 16-QAM
        constellation = [-3 -1 3 1] / sqrt(10);
        realPart = constellation(mod(grayIndex, 4) + 1);
        imagPart = constellation(floor(grayIndex / 4) + 1);
        qamSymbol = realPart + 1i * imagPart;
    else
        error('Mapping not implemented for this value of M');
    end
end

% Function for QAM to Gray Demapping (16-QAM)
function bits = qam_to_gray(qamSymbol, M)
    if M == 16
        qamSymbol = qamSymbol * sqrt(10);  % Rescale for comparison
        constellation = [-3 -1 1 3];  % Original constellation points
        [~, realIdx] = min(abs(real(qamSymbol) - constellation));
        [~, imagIdx] = min(abs(imag(qamSymbol) - constellation));
        grayIndex = (imagIdx - 1) * 4 + (realIdx - 1);
        bits = de2bi(grayIndex, log2(M), 'left-msb')';
    else
        error('Demapping not implemented for this value of M');
    end
end
