clc;
clearvars;
close all;
[wavdata, Fs] = audioread('project.wav');
% quantise the interval -1 to 1 in 128 levels and take closest quantised
% level for wavdata(1) to wavdata(410624) and scale by 128 and store it in
% bits_vector as a 8 bit number bits_vector=zeros(410624*8,1)
 wavnew=zeros(410624,1);
 wavedata_cpy=wavdata;
 for kk=1:410624
     wavnew(kk)=wavdata(kk);
end
temp=wavnew+1;
symbols_tx=round(temp*64);
snew=dec2bin(symbols_tx) - '0';
tri=snew';
bits_fuf=zeros(410624*8,1);
for gv=1:length(tri)
    bits_fuf(gv)=tri(gv);
end
%%
sd = 6;
m = 4;
oversampling_factor = 17;
a = 0.5;
lengthe = 5;
symbols = zeros(length(bits_fuf)/2, 1);
[transmit_filter, ~] = raised_cosine(a, oversampling_factor, lengthe);
for i = 1:length(symbols)
    index = (i - 1) * 2 + 1;
    symbols(i) = qpskmap(bits_fuf(index:index+1));
end
symbols_upsampled = upsample(symbols, oversampling_factor);
tx_output = conv(symbols_upsampled, transmit_filter, 'same');
receive_filter = conj(flipud(transmit_filter));
fc=1000000;
Ts =1/(17*fc);
fs = 1 / Ts; % Symbol rate (sampling frequency)
symtime =0.000001;
t = linspace(0, (length(symbols)) * symtime, length(tx_output)); % Time axis from 0 to n*Ts

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
 sd = sd/10;
 mul1 = cos(2*pi*fc*t);
 p1 = real(tx_output).*mul1';
 mul2 = sin(2*pi*fc*t);
 p2 = imag(tx_output).*mul2';
 %noise adder block
 p3=p1+p2;
 noise_inside_channel=sd*randn(size(p3));
 p3_with=p3+noise_inside_channel;
 figure;
 subplot(2,1,1);
 plot(t,p3);
 subplot(2,1,2);
 plot(t,p3_with);
 p3_with1=p3_with.*mul1';
 p3_with2=p3_with.*mul2';
 tx_output_real=lowpass(p3_with1,1000000,fs,ImpulseResponse="iir",Steepness=0.95);
 tx_output_imag=lowpass(p3_with2,1000000,fs,ImpulseResponse="iir",Steepness=0.95);
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
 dec_sym=zeros(length(symbols),1);
 rb=1;
 tx_output_filtered2=downsample(tx_output_filtered,oversampling_factor);
 for yy=1:length(tx_output_filtered2)
 if abs(real(tx_output_filtered2(yy)))>abs(imag(tx_output_filtered2(yy))) & real(tx_output_filtered2(yy))>0
     dec_sym(rb)=1;
     rb=rb+1;
 end
 if abs(real(tx_output_filtered2(yy)))>abs(imag(tx_output_filtered2(yy))) & real(tx_output_filtered2(yy))<0
     dec_sym(rb)=-1;
     rb=rb+1;
 end
 if abs(real(tx_output_filtered2(yy)))<abs(imag(tx_output_filtered2(yy))) & imag(tx_output_filtered2(yy))>0
     dec_sym(rb)=1j;
     rb=rb+1;
 end
 if abs(real(tx_output_filtered2(yy)))<abs(imag(tx_output_filtered2(yy))) & imag(tx_output_filtered2(yy))<0
 dec_sym(rb)=-1j;
     rb=rb+1;
 end
 end
 ber_rate=sum(symbols~=dec_sym)/length(symbols);
 disp(ber_rate);
 rvr=zeros(2*length(dec_sym),1);
 rf=1;
 for h=1:length(dec_sym)
     if dec_sym(h)==1
         rvr(rf)=0;
         rvr(rf+1)=0;
         rf=rf+2;
     end
     if dec_sym(h)==-1
         rvr(rf)=0;
         rvr(rf+1)=1;
         rf=rf+2;
     end
     if dec_sym(h)==1j
         rvr(rf)=1;
         rvr(rf+1)=0;
         rf=rf+2;
     end
     if dec_sym(h)==-1j
         rvr(rf)=1;
         rvr(rf+1)=1;
         rf=rf+2;
     end
 end
 ber_fin=sum(rvr~=bits_fuf)/length(bits_fuf);
 disp(ber_fin);
 %keep on taking 8 bits from rvr convert them to decimal and store in an
 %array
 % Calculate the number of decimal values to extract
 num_values = length(rvr) / 8;
 ratf = zeros(num_values, 1);
 dj=1;
 for j = 1:num_values
       curr_value=128*rvr(dj)+64*rvr(dj+1)+32*rvr(dj+2)+16*rvr(dj+3);
       curr_value=curr_value+8*rvr(dj+4)+4*rvr(dj+5)+2*rvr(dj+6)+rvr(dj+7);
     ratf(j) =curr_value;
     dj=dj+8;
 end
 ratf=ratf/64;
 ratf=ratf-1;
 for kk2=1:410624
      wavedata_cpy(kk2)=ratf(kk2);
 end
 sound(wavedata_cpy,Fs); 


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
