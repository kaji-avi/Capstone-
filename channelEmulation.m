function y = channelEmulation(x,snr,to,h) 

y = zeros(1,to); % create some offset
y = [y conv(h,x)];

% scale with SNR and multiply CFO
y = sqrt(snr).*y;

% add awgn
y = y+sqrt(1/2).*(randn(size(y))+1i.*randn(size(y)));

end