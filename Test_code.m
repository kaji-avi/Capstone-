<<<<<<< HEAD
clc
clear
close all
=======
clc;
clear;
>>>>>>> e2b27fdbd4fa405e18faecadf10dfa9b41394c11
%% Transmitter

% Set a random seed for reproducibility
rng(42); 

% Parameters
N = 64;  % Number of OFDM subcarriers
M = 16;  % 16-QAM
numBits = N * log2(M);  % Total number of bits
bitsPerSymbol = log2(M);  % Bits per symbol (4 bits for 16-QAM)
<<<<<<< HEAD
cyclicPrefixLength = 48;
=======
cyclicPrefixLength = 8;
>>>>>>> e2b27fdbd4fa405e18faecadf10dfa9b41394c11

% Roberto's Parameters
SNR_dB = 20;  % SNR in dB
to = 0;  % No time offset
h = 1;  % Channel impulse response

% Convert SNR from dB to linear scale
SNR = 10^(SNR_dB / 10);  

<<<<<<< HEAD

=======
% Generate random data bits
>>>>>>> e2b27fdbd4fa405e18faecadf10dfa9b41394c11
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

% Add Cyclic Prefix (CP)
timeDomainSymbolsWithCP = [timeDomainSymbols(end-cyclicPrefixLength+1:end, :); timeDomainSymbols];
x = timeDomainSymbolsWithCP(:);  % Convert to serial for transmission

% Pass through the channel (provided function)
y = channelEmulation(x, SNR_dB, to, h);

%% Receiver

% Remove CP, Serial to Parallel, FFT
receivedParallelData = reshape(y(1:length(x)), N + cyclicPrefixLength, []);  % Serial to parallel conversion
receivedSymbolsWithoutCP = receivedParallelData(cyclicPrefixLength+1:end, :);  % Remove CP
receivedFrequencyDomainSymbols = fft_function(receivedSymbolsWithoutCP, N);  % FFT

% Equalization Using Known Channel (h_hat)
h_hat = fft(h, N);  % Assuming h is given in time domain, take its FFT

% Calculate the actual noise variance
signalPower = mean(abs(timeDomainSymbols(:)).^2);  % Signal power
noiseVariance = signalPower / SNR;  % Noise variance

% MMSE Equalization
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

%% Function for MMSE equalization
function equalizedSymbols = mmse_equalization(receivedSymbols, H, noiseVariance, signalPower)
    mmseFilter = conj(H) ./ (abs(H).^2 + noiseVariance / signalPower);
  
    equalizedSymbols = mmseFilter .* receivedSymbols;
end

%% Function for IFFT calculation
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
    timeDomainSymbols = timeDomainSymbols * sqrt(Nsubcarriers);
end

%% Function for FFT calculation
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

%% Function for Gray to QAM Mapping
function qamSymbol = gray_to_qam(grayIndex, M)
    if M == 16
<<<<<<< HEAD
        realMap = [-3 -1 1 3];
        imagMap = [-3 -1 1 3];
=======
        % Define the real and imaginary maps for 16-QAM
        realMap = [-3 -1 1 3];
        imagMap = [-3 -1 1 3];

        % Ensure Gray index is within bounds
        assert(all(grayIndex >= 0 & grayIndex < M), 'Gray index out of range');

        % Calculate real and imaginary indices (1-based)
>>>>>>> e2b27fdbd4fa405e18faecadf10dfa9b41394c11
        realIdx = mod(grayIndex, 4) + 1;
        imagIdx = floor(grayIndex / 4) + 1;

        % Extract corresponding real and imaginary parts
        realPart = realMap(realIdx);
        imagPart = imagMap(imagIdx);

        % Normalize power
        qamSymbol = (realPart + 1i * imagPart) / sqrt(10);
    else
        error('Mapping not implemented for this value of M');
    end
end

%% Function for QAM to Gray Demapping
function bits = qam_to_gray(qamSymbol, M)
    if M == 16
        % Define the real and imaginary maps for 16-QAM
        realMap = [-3 -1 1 3];
        imagMap = [-3 -1 1 3];

        % Scale back to original values
        qamSymbol = qamSymbol * sqrt(10);

        % Quantize to nearest constellation points
        realPart = quantize_to_nearest(real(qamSymbol), realMap);
        imagPart = quantize_to_nearest(imag(qamSymbol), imagMap);

        % Find corresponding indices
        realIdx = find(realMap == realPart);
        imagIdx = find(imagMap == imagPart);

        % Calculate Gray-coded index
        grayIndex = (imagIdx - 1) * 4 + (realIdx - 1);

        % Convert index to binary bits
        bits = de2bi(grayIndex, log2(M), 'left-msb')';
    else
        error('Demapping not implemented for this value of M');
    end
end

<<<<<<< HEAD
=======
% Helper Function to Quantize to Nearest Value
function nearestVal = quantize_to_nearest(value, map)
    [~, idx] = min(abs(value - map));
    nearestVal = map(idx);
end
>>>>>>> e2b27fdbd4fa405e18faecadf10dfa9b41394c11
