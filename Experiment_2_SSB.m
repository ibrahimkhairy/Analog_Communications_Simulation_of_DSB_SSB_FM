clc
clear
close all 

[y,~] = audioread('test.mp3');
%sound(y,fs);
y=resample(y,5000,441);
fs=500000;
t=linspace(0,10,10*fs);
message = 0.5*(y(:,1) + y(:,2));  %getting channels average
message = message(1:10*fs);         %10 seconds are enough
plot(linspace(0,10,length(message)),message);
title('Original message')
xlabel('Time (Sec)')
%-------------------
%(1)
message_length = length (message);
original_spectrum= fft(message,message_length); %getting fourrier transform 
f_spacing=fs/(message_length-1);                  %frequency spacing  
f_scale=1/f_spacing;                                %frequency scale , to be able to visualize it    
f= (-fs/2:f_spacing:fs/2);%dividing into samples which has spacing equal to 1/(message length in seconds) to have spectrum in hertz
original_spectrum_centered=fftshift(original_spectrum);%centering spectrum
figure;
% plot (f/f_scale,original_spectrum_centered);
plot (f,abs(original_spectrum_centered));
title('Original message spectrum')
xlabel('Frequency (Hz)')
%-------------------
%(2)
num_of_seconds=message_length/fs;
filtering_frequency=4000;
M2 = ones(round(2*filtering_frequency*f_scale),1);
M1= zeros(round((message_length-2*filtering_frequency*f_scale)/2),1);
M = [M1 ;M2; zeros(message_length-length(M1)-length(M2),1)];%rectangular signal
filtered_spectrum = original_spectrum_centered .* M;
figure; plot (f,abs(filtered_spectrum));
title('Filtered message spectrum')
xlabel('Frequency (Hz)')
%-------------------
%(3)
 filtered_message = ifft(ifftshift(filtered_spectrum),message_length,'symmetric');%used ifftshift because ifft have to be applied on the noncentered spectrum
 figure;  plot (linspace(0,10,length(message)),filtered_message);
 title('Filtered message')
xlabel('Time (Sec)')
%  sound (resample(filtered_message,441,5000),fs*441/5000); %to be able to hear it i resambled into 44100

 %-------------------
%(4)
[psnr,mse,maxerr,l] = measerr(message,filtered_message);
 Error_in_filtered_baseband=mse
 %-------------------
%(5)
fc=100000*f_scale;
t_mod=linspace(1,length(filtered_message)/fs,length(filtered_message)) ; 
 SC_modulated = (filtered_message).*cos(2*pi*fc*t_mod)';%SC_MOD

SC_modulated_spectrum=fftshift(fft(SC_modulated,message_length));

figure
plot (f,abs(SC_modulated_spectrum));
title('SC modulated spectrum')
xlabel('Frequency (Hz)')
 %-------------------
%(6)
M=[zeros(1500000-1,1); ones(2000002,1) ;zeros(1500000-1,1)];
LSB_SC_spectrum= SC_modulated_spectrum.*M ;
figure
plot (f,abs(LSB_SC_spectrum))
title('LSB-SC spectrum using IBPF')
xlabel('Frequency (Hz)')
 %-------------------
%(7)
LSB_SC=ifft(ifftshift(LSB_SC_spectrum),message_length,'symmetric');
unfiltered_LSB_SC= LSB_SC.*cos(2*pi*fc*t_mod)'; %The coherent detector.
ideal_filtered_LSB_SC=fftshift(fft(unfiltered_LSB_SC)).*[zeros(2460000,1);ones(80000,1);zeros(2460000,1)];%low pass filtering
figure
demod_LSB_SC_message_coherent_no_noise= ifft(ifftshift(ideal_filtered_LSB_SC),message_length,'symmetric');
plot (linspace(0,10,length(demod_LSB_SC_message_coherent_no_noise)),demod_LSB_SC_message_coherent_no_noise)
title('demod LSB-TC message using IBPF with no noise')
xlabel('Time (S)')

 [~,mse,~,~] = measerr(message,demod_LSB_SC_message_coherent_no_noise);
 Error_in_filtered_coherent_LSB_SC_no_noise=mse %printing ERROR 
 
%sound (resample(demod_LSBmessage_coherent_no_noise,48,500),48000); %playing the audio 
%-------------------
%(8)
% 8.a trimming upper side band
[b,a] = butter(3,[0.3 0.7],'bandpass');  %getting filter coefficiants for the baseband filtering
y=abs(freqz(b,a));%getting filter spectrum
bandpass_butterworth_filter=resample(y,5000000,length(y));%upsampling filter spectrum

LSB_SC_spectrum= SC_modulated_spectrum.*bandpass_butterworth_filter ;
figure

plot (f,abs(LSB_SC_spectrum))
title('LSB-TC spectrum using practical BPF')
xlabel('Frequency (Hz)')

figure
plot (f,abs(bandpass_butterworth_filter));
title('Filter spectrum')
xlabel('Frequency (Hz)')

% 8.B demodulating ,trimming upper components
LSB_SC=ifft(ifftshift(LSB_SC_spectrum),message_length,'symmetric');
unfiltered_practical_LSB_SC= LSB_SC.*cos(2*pi*fc*t_mod)'; %The coherent detector.

[b,a] = butter(3,[0.492 0.508],'bandpass');  %getting filter coefficiants for the baseband filtering
y=abs(freqz(b,a));%getting filter spectrum
butterworth_filter=resample(y,5000000,length(y));%upsampling filter spectrum

practical_filtered_LSB_SC=fftshift(fft(unfiltered_practical_LSB_SC)).*butterworth_filter;%practical low pass filtering
practical_filtered_LSB_SC(2500000+1)=0;%DC block
practical_demod_LSBmessage_coherent_no_noise= ifft(ifftshift(practical_filtered_LSB_SC),message_length,'symmetric');
figure
plot (linspace(0,10,length(practical_demod_LSBmessage_coherent_no_noise)),practical_demod_LSBmessage_coherent_no_noise)
title('demod LSB-SC message using practical BPF with no noise')
xlabel('Time (S)')


  [~,mse,~,~] = measerr(message,practical_demod_LSBmessage_coherent_no_noise);
 Error_in_filtered_coherent_LSB_SC_no_noise_practical=mse 
%sound (resample(practical_demod_LSBmessage_coherent_no_noise,48,500),48000); %playing the audio 
%-------------------
%(9)
SNR = [0 10 20];

   noisy_signal(:,1) = awgn(LSB_SC,SNR(1));
unfiltered_noisy_LSB_SC(:,1)=  noisy_signal(:,1).*cos(2*pi*fc*t_mod)'; %The coherent detector.
practical_noisy_filtered_LSB_SC(:,1)=fftshift(fft(unfiltered_noisy_LSB_SC(:,1))).*butterworth_filter;%practical low pass filtering
practical_noisy_filtered_LSB_SC(2500000+1,1)=0;%DC block
practical_demod_LSBmessage_coherent_noisy(:,1)= ifft(ifftshift(practical_noisy_filtered_LSB_SC(:,1)),message_length,'symmetric');

   noisy_signal(:,2) = awgn(LSB_SC,SNR(2));
unfiltered_noisy_LSB_SC(:,2)=  noisy_signal(:,2).*cos(2*pi*fc*t_mod)'; %The coherent detector.
practical_noisy_filtered_LSB_SC(:,2)=fftshift(fft(unfiltered_noisy_LSB_SC(:,2))).*butterworth_filter;%practical low pass filtering
practical_noisy_filtered_LSB_SC(2500000+1,2)=0;%DC block
practical_demod_LSBmessage_coherent_noisy(:,2)= ifft(ifftshift(practical_noisy_filtered_LSB_SC(:,2)),message_length,'symmetric');

   noisy_signal(:,3) = awgn(LSB_SC,SNR(3));
unfiltered_noisy_LSB_SC(:,3)=  noisy_signal(:,3).*cos(2*pi*fc*t_mod)'; %The coherent detector.
practical_noisy_filtered_LSB_SC(:,3)=fftshift(fft(unfiltered_noisy_LSB_SC(:,3))).*butterworth_filter;%practical low pass filtering
practical_noisy_filtered_LSB_SC(2500000+1,3)=0;%DC block
practical_demod_LSBmessage_coherent_noisy(:,3)= ifft(ifftshift(practical_noisy_filtered_LSB_SC(:,3)),message_length,'symmetric');

 
  
  [~,mse,~,~] = measerr(message,practical_demod_LSBmessage_coherent_noisy(:,1));
 Error_in_0_SNR=mse %printing ERROR at 0 SNR 
 
   [~,mse,~,~] = measerr(message,practical_demod_LSBmessage_coherent_noisy(:,2));
 Error_in_10_SNR=mse %printing ERROR at 10 SNR 
 
   [~,mse,~,~] = measerr(message,practical_demod_LSBmessage_coherent_noisy(:,3));
 Error_in_20_SNR=mse %printing ERROR at 20 SNR 
 
figure
plot (linspace(0,10,length(practical_demod_LSBmessage_coherent_noisy(:,1))),practical_demod_LSBmessage_coherent_noisy(:,1))
title('demod LSB-SC message using practical BPF with SNR = 0')
xlabel('Time (S)')

figure
plot (linspace(0,10,length(practical_demod_LSBmessage_coherent_noisy(:,2))),practical_demod_LSBmessage_coherent_noisy(:,2))
title('demod LSB-SC message using practical BPF with SNR = 10')
xlabel('Time (S)')

figure
plot (linspace(0,10,length(practical_demod_LSBmessage_coherent_noisy(:,3))),practical_demod_LSBmessage_coherent_noisy(:,3))
title('demod LSB-SC message using practical BPF with SNR = 20')
xlabel('Time (S)')

%sound (resample(practical_demod_LSBmessage_coherent_noisy(:,1),48,500),48000); %playing the audio 
%sound (resample(practical_demod_LSBmessage_coherent_noisy(:,2),48,500),48000); %playing the audio 
%sound (resample(practical_demod_LSBmessage_coherent_noisy(:,3),48,500),48000); %playing the audio 
%-------------------
%(10)
% A)
unfiltered_noisy_LSB_SC(:,2)=  noisy_signal(:,2).*cos(2*pi*fc*(1.0001/1)*t_mod)'; %The coherent detector.
practical_noisy_filtered_LSB_SC_freq_error(:,2)=fftshift(fft(unfiltered_noisy_LSB_SC(:,2))).*butterworth_filter;%practical low pass filtering
practical_demod_LSBmessage_coherent_noisy_freq_error(:,2)= ifft(ifftshift(practical_noisy_filtered_LSB_SC_freq_error(:,2)),message_length,'symmetric');
figure
plot (f,abs(practical_noisy_filtered_LSB_SC_freq_error(:,2)))
title('LSB-TC spectrum with frequency mismatch')
xlabel('Frequency (Hz)')
%sound (resample(practical_demod_LSBmessage_coherent_noisy_freq_error(:,2),48,500),48000); %playing the audio 

% B)
unfiltered_noisy_LSB_SC(:,2)=  noisy_signal(:,2).*cos(2*pi*fc*t_mod+1*180/pi)'; %The coherent detector.
practical_noisy_filtered_LSB_SC_phase_error(:,2)=fftshift(fft(unfiltered_noisy_LSB_SC(:,2))).*butterworth_filter;%practical low pass filtering
practical_demod_LSBmessage_coherent_noisy_phase_error(:,2)= ifft(ifftshift(practical_noisy_filtered_LSB_SC_phase_error(:,2)),message_length,'symmetric');
figure
plot (f,abs(practical_noisy_filtered_LSB_SC_phase_error(:,2)))
title('LSB-TC spectrum with phase mismatch')
xlabel('Frequency (Hz)')
%sound (resample(practical_demod_LSBmessage_coherent_noisy_phase_error(:,2),48,500),48000); %playing the audio 

%-------------------
%(11)


TC_modulated = (filtered_message-2*min(min(filtered_message))).*cos(2*pi*fc*t_mod)';%TC_MOD
TC_modulated_spectrum=fftshift(fft(TC_modulated,message_length));
LSB_TC_modulated_spectrum= TC_modulated_spectrum.*M ;
figure
plot (f,abs(LSB_TC_modulated_spectrum))
title('LSB-TC spectrum using IBPF')
xlabel('Frequency (Hz)')

mod_LSB_SC_time_domain = ifft(ifftshift(LSB_SC_spectrum),message_length,'symmetric');
Envelope_Detected_LSB_TC = abs(hilbert(mod_LSB_SC_time_domain));%The envelope detector.


figure
plot (linspace(0,10,length(Envelope_Detected_LSB_TC)),Envelope_Detected_LSB_TC)
title('Envelope detected LSB-TC signal')
xlabel('Time (S)')

   [~,mse,~,~] = measerr(message,Envelope_Detected_LSB_TC);
Error_in_Envelope_Detected_LSB_TC=mse %printing ERROR

%sound (resample(Envelope_Detected_LSB_TC,48,500),48000); %playing the audio 

%-------------------
%(12)
%Disadvantages of SSB:
% Complex circuitry is required for high order filter
% oscillator frequencies involved are very critical and must be stabilized otherwise distortion will occur
% bit expansive in comparison to DSB
%-------------------
%(13)
LSB_using_helbert= filtered_message.*cos(2*pi*fc*t_mod)'-imag(hilbert(filtered_message)).*sin(2*pi*fc*t_mod)';

LSB_using_helbert_spectrum=fftshift(fft(LSB_using_helbert,message_length));
figure
plot (f,abs(LSB_using_helbert_spectrum));
title('LSB spectrum -using helbert')
xlabel('Frequency (Hz)')