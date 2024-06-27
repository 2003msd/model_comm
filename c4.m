clc;
clearvars;
close all;
m = 4;
oversampling_factor = 1000;
a = 0.5;
lengthe = 5;
[transmit_filter, ~] = raised_cosine(a, m, lengthe);
nsymbols = 16; 
symbols_nf = randi([0, 1], nsymbols * 2, 1);
symbols = zeros(nsymbols, 1);
for i = 1:nsymbols
    index = (i - 1) * 2 + 1;
    symbols(i) = qpskmap(symbols_nf(index:index+1));
end
symbols_upsampled = upsample(symbols, 1000);
tx_output = conv(symbols_upsampled, transmit_filter, 'same');
receive_filter = conj(flipud(transmit_filter));
Ts = 0.001;
symtime = 1;
fs = 1 / Ts; % Symbol rate (sampling frequency)
t = linspace(0, (nsymbols) * symtime, length(tx_output)); % Time axis from 0 to n*Ts
figure;
subplot(2,1,1);
plot(t, real(tx_output));
title('real');
xlabel('Time');
ylabel('Amplitude');
subplot(2,1,2);
plot(t, imag(tx_output));
title('imag');
xlabel('Time');
ylabel('Amplitude');
%modulation block
mul1 = cos(2*pi*100*t);
p1 = real(tx_output).*mul1';
mul2 = sin(2*pi*100*t);
p2 = imag(tx_output).*mul2';
%noise adder block
p3=p1+p2;
noise_inside_channel=0*randn(size(p3));
p3_with=p3+noise_inside_channel;
figure;
subplot(2,1,1);
plot(t,p3);
subplot(2,1,2);
plot(t,p3_with);
p3_with1=p3_with.*mul1';
p3_with2=p3_with.*mul2';
cutoff_frequency = 100;
tx_output_real=lowpass(p3_with1,50,fs,ImpulseResponse="iir",Steepness=0.95);
tx_output_imag=lowpass(p3_with2,50,fs,ImpulseResponse="iir",Steepness=0.95);
figure;
subplot(2,1,1);
plot(t, tx_output_real);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Real Part of tx\_output');
grid on;
subplot(2,1,2);
plot(t, tx_output_imag);
xlabel('Time (s)');
ylabel('Amplitude');
title('Recovered Imaginary Part of tx\_output');
grid on;
tx_output_reconstructed = tx_output_real + 1i * tx_output_imag;
tx_output_filtered = conv(tx_output_reconstructed, receive_filter, 'same');
figure;
subplot(2,1,1);
plot(t, real(tx_output_filtered));
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed and Filtered tx\_output');
grid on;
subplot(2,1,2);
plot(t, imag(tx_output_filtered));
xlabel('Time (s)');
ylabel('Amplitude');
title('Reconstructed and Filtered tx\_output');
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
