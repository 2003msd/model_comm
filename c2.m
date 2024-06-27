clc;
clearvars;
clear all;
m = 4;
oversampling_factor = 4;

% Parameters for sampled raised cosine pulse
a = 0.5;
lengthe = 5; % (truncated outside [-length*T,length*T])

% Raised cosine transmit filter (time vector set to a dummy variable which is not used)
[transmit_filter, ~] = raised_cosine(a, m, lengthe);

% Number of symbols
nsymbols = 6000; % For QPSK, we'll have half the number of symbols compared to 4-PAM
symbols_nf = randi([0, 1], nsymbols * 2, 1); % Generate random binary data for QPSK

% Generate QPSK symbols
symbols = zeros(nsymbols, 1);
for i = 1:nsymbols
    index = (i - 1) * 2 + 1;
    symbols(i) = qpskmap(symbols_nf(index:index+1));
end

% Upsample
symbols_upsampled = upsample(symbols, m);

% Noiseless modulated signal
tx_output = conv(symbols_upsampled, transmit_filter, 'same');

% Calculate noise variance based on Eb/N0
Eb_N0_dB = 8.9;
Eb_N0_lin = 10^(Eb_N0_dB / 10);
Eb = (2.5) * (norm(transmit_filter)^2);
N0 = Eb / Eb_N0_lin;
sigma2 = N0 / 2; 

% Add noise to the signal
noise = 0.1*sqrt(sigma2) * (randn(size(tx_output)) + 1j * randn(size(tx_output)));
rx_signal = tx_output + noise;
N_not=0.16;
scf=N_not/2;
% Receive filter
receive_filter = conj(flipud(transmit_filter));
f_output = conv(rx_signal, receive_filter, 'same')/3.5;
final_out = downsample(f_output, m);
scatter(real(final_out),imag(final_out));
xlabel('REAL PART');
ylabel('IMAGINARY PART');
title('QPSK SCHEME');
% Rectified code for condition checking
error_probabiliy=rand*scf;
gg = zeros(12000, 1); % Preallocate gg array
kt = 1;
for vvb = 1:length(final_out)
    if real(final_out(vvb)) > 0 && imag(final_out(vvb)) > 0
        gg(kt) = 0;
        gg(kt+1) = 0;
        kt = kt + 1;
    elseif real(final_out(vvb)) < 0 && imag(final_out(vvb)) > 0
        gg(kt) = 0;
        gg(kt+1) = 1;
        kt = kt + 1;
    elseif real(final_out(vvb)) < 0 && imag(final_out(vvb)) < 0
        gg(kt) = 1;
        gg(kt+1) = 1;
        kt = kt + 1;
    elseif real(final_out(vvb)) > 0 && imag(final_out(vvb)) < 0
        gg(kt) = 1;
        gg(kt+1) = 0;
        kt = kt + 1;
    end
end
error_probability = sum(gg ~= symbols_nf) / nsymbols;
disp('PROBABILITY OF ERROR');
disp(error_probabiliy);
RANGE = 0:0.1:30;
mediate = 10.^(RANGE / 10);
syms integral_e_not;
integral_inf = 0.5 * erfc(sqrt(mediate));
figure;
semilogy(RANGE, integral_inf);
hold on;
ylabel('ERROR PROBABILITY');
xlabel('SIGNAL TO NOISE RATIO IN DECIBELS');
title('QPSK ERROR PROBABILITY GRAPH');
grid on;

function [rc, time_axis] = raised_cosine(a, m, length)
    length_os = floor(length * m); % number of samples on each side of peak
    z = cumsum(ones(length_os, 1)) / m; % time vector on one side of the peak
    A = sin(pi * z) ./ (pi * z); % term 1
    B = cos(pi * a * z); % term 2
    C = 1 - (2 * a * z).^2; % term 3
    zerotest = m / (2 * a); % location of zero in denominator
    % check whether any sample coincides with zero location
    if (zerotest == floor(zerotest))
        B(zerotest) = pi * a;
        C(zerotest) = 4 * a;
    end
    D = (A .* B) ./ C; % response to one side of peak
    rc = [flipud(D); 1; D]; % add in peak and other side
    time_axis = [flipud(-z); 0; z];
end

function ans_bit=qpskmap(input_vec)
         if(input_vec(1)==0 && input_vec(2)==0)
             ans_bit=1;
         end
         if(input_vec(1)==0 && input_vec(2)==1)
             ans_bit=-1;            
         end
         if(input_vec(1)==1 && input_vec(2)==0)
             ans_bit=1i;
         end
         if(input_vec(1)==1 && input_vec(2)==1)
             ans_bit=-1i;
         end
end
