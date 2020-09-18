clc
clear all
close all

[y,fs] = audioread('eric.wav');
m = y(1:5*fs,1);                                 % message equals first 5 seconds of the audio file
ts = 1/fs;                                       % sampling time
tf = length(m)/fs;                               % time limit of the message
t = linspace(0,tf,length(m));                    % time vector of the sampled message

L = length(m);                                   
N = 2^nextpow2(L);
M = fft(m,N);                                    % fourier transform of the message
f = (-N/2:N/2-1)/(N*ts);                         % frequency vector of the sampled message                        
MS = fftshift(M);                                % message centered around 0HZ
BW = fs/2;                                       % band width of the message 

%-------------------------------------------------------------------------%
%              REMOVING ALL FREQUENCIES GREATER THAN 4KHZ                 %
%-------------------------------------------------------------------------%
wn = 4000/BW;                                    % normalized cut off frequency
a = fir1(300,wn);                                % numerator coefficients of filter transfer function with order=300 & normalized cut off frequency=4KHZ/BW
b = 1;
m_fil = filter(a,b,m);                           % filter message at 4KHZ
M_fil = fftshift(fft(m_fil,N));                  % filtered message at frequency domain

[psnr,mse,maxerr,l] = measerr(m,m_fil);          % calculate Mean-Square-Error between the filtered message & the original one

                    %###### DOUBLE SIDEBAND MODULATION TC ########%
%-------------------------------------------------------------------------%
% MODULATION AT THE DESIRED SAMPLING FREQUENCY & DESIRED CARRIER FREQUENCY%
%-------------------------------------------------------------------------%
[p,q] = rat(500000/fs);
m_up = resample(m_fil,p,q);                      % resample message to the desired sampling frequency (500KHZ)
v = var(m_up);                                   % variance-total power- of the message
ts_up = 1/500000;                                % sampling time of the resampled message
dc = 2*abs(min(m_up));                           % DC offset to ensure envelope detector can be used
m_plus_dc = m_up + dc;                           
fc = 100000;                                     % carrier frequency (100KHZ)
t_up = linspace(0,tf,length(m_up));              % time vector of the resampled message
c = cos(2*pi*fc*t_up);                           % modulation carrier
s = 2*m_plus_dc'.*c;                             % modulated message in time domain
LN = length(m_plus_dc);                          
NN = 2^nextpow2(LN);
S = fftshift(fft(s,NN));                         % modulated message in frequency domain
ff = (-NN/2:NN/2-1)/(NN*ts_up);

%-------------------------------------------------------------------------%
%        DETECTION USING ENVELOPE DETECTOR WITHOUT NOISE CONTRIBUTION     %
%-------------------------------------------------------------------------%
envel = abs(hilbert(s)) - dc;                    % detected message using envelope detector demodulation method
envel = envel';
envel_f = fftshift(fft(envel,NN));               % detected message in frequency domain

D = abs(envel - m_up).^2;
error = sum(D(:))/numel(envel);                  % Mean-Square-Error between the detected message and the transmitted one

% sound(m,fs);
% [w,r] = rat(fs/500000);
% envel_down = resample(envel,w,r);
% sound(envel_down,fs);

%-------------------------------------------------------------------------%
%                    DEMODULATION IN PRESENCE OF NOISE                    %
%-------------------------------------------------------------------------%
SNR = [0 10 20];
noisy_signal(:,1) = awgn(s,SNR(:,1));                      % modulated message plus noise at SNR=0dB
envel_noisy(1,:) = abs(hilbert(noisy_signal(:,1))) - dc;   % detcted message in presence of noise

noisy_signal(:,2) = awgn(s,SNR(:,2));                      % modulated message plus noise at SNR=10dB
envel_noisy(2,:) = abs(hilbert(noisy_signal(:,2))) - dc;   % detcted message in presence of noise

noisy_signal(:,3) = awgn(s,SNR(:,3));                      % modulated message plus noise at SNR=20dB
envel_noisy(3,:) = abs(hilbert(noisy_signal(:,3))) - dc;   % detcted message in presence of noise
envel_noisy = envel_noisy';

% sound(m,fs);
% [w,r] = rat(fs/500000);
% envel_down_noisy = resample(envel_noisy,w,r);
% sound(envel_down_noisy,fs);

efficiency = 100 * (v / (v+dc^2));               % efficiency = m^2 / m^2 + A^2   where A is the dc offset

%-------------------------------------------------------------------------%
%                      PLOTING SIGNALS IN TIME DOMAIN                     %
%-------------------------------------------------------------------------%
figure(1)
subplot(2,2,1)
plot(m)
title('Original Message In Time Domain')
grid on

subplot(2,2,2)
plot(t,m_fil)
title('Filtered Message In Time Domain')
grid on

subplot(2,2,3)
plot(t_up,s)
title('Modulated Message in Time Domain')
grid on

subplot(2,2,4)
plot(t_up,envel)
title('Detected Message In Time Domain')
grid on

clear y M m_plus_dc D

%-------------------------------------------------------------------------%
%                      PLOTING SIGNALS IN FREQUENCY DOMAIN                %
%-------------------------------------------------------------------------%
figure(2)
subplot(2,2,1)
plot(f,real(MS))
title('Original Message In Frequency Domain')
grid on

subplot(2,2,2)
plot(f,M_fil)
title('Filtered Message In Frequency Domain')
grid on

subplot(2,2,3)
plot(ff,S)
title('Modulated Message in Frequency Domain')
axis([-1.5e5 1.5e5 -6000 6000])
grid on

subplot(2,2,4)
plot(ff,envel_f)
title('Detected Message In Frequency Domain')
axis([-4000 4000 -6000 6000])
grid on
                 %###### DOUBLE SIDEBAND MODULATION SC ########%
%-------------------------------------------------------------------------%
% MODULATION AT THE DESIRED SAMPLING FREQUENCY & DESIRED CARRIER FREQUENCY%
%-------------------------------------------------------------------------%
[y_sc,fs_sc] = audioread('eric.wav'); 
m_sc = y_sc(1:5*fs_sc,1);
ts_sc = 1/fs_sc;
tf_sc = length(m_sc)/fs_sc;
plot(m_sc);
title(' Original Message ')
L_sc = length(m_sc);
N_sc = 2^nextpow2(L_sc);
M_sc = fft(m_sc,N_sc);
f_sc = (-N_sc/2:N_sc/2-1)/(N_sc*ts_sc);
MS_sc = fftshift(M_sc);
figure
plot(f_sc,real(MS_sc))
title('Original Message In Frequency Domain')
grid on

a_sc = fir1(300,0.1814);

m_fil_sc = filter(a_sc,b,m_sc); % filtering the massage 
M_fil_sc = fftshift(fft(m_fil_sc,N_sc));
figure
plot(f,real(M_fil_sc))
title('Original Message After Filtiring')
grid on

m_up_sc = resample(m_fil_sc,p,q);                                 
t_sc = linspace(0,tf_sc,length(m_up_sc));                % time vector of the resampled message
c_sc = cos(2*pi*fc*t_sc);                           % modulation carrier
s_sc = 2*m_up_sc'.*c_sc;                             % modulated message in time domain
LN_sc = length(m_up_sc);                          
NN_sc = 2^nextpow2(LN_sc);
S_sc = fftshift(fft(s_sc,NN_sc));                         % modulated message in frequency domain
ff_sc = (-NN_sc/2:NN_sc/2-1)/(NN_sc*ts_up);
figure
plot(ff_sc,S_sc)
xlabel('Spectrum of Message Signal After cohernt Modualtion ');
grid on
%-------------------------------------------------------------------------%
%        DETECTION USING ENVELOPE DETECTOR WITHOUT NOISE CONTRIBUTION     %
%-------------------------------------------------------------------------%
envelope=abs(hilbert(s_sc)); % detected message using envelope detector demodulation method
figure();
plot(envelope)
xlabel('The Message Signal After ENVELOPE DETECTOR WITHOUT NOISE CONTRIBUTION ');
grid on                   

%D_sc = abs(envelope - m_up).^2;
%error_SC = sum(D_sc(:))/numel(envelope);    % Mean-Square-Error between the detected message and the transmitted one
%-------------------------------------------------------------------------%
%        DETECTION USING COHERANT DETECTOR WITHOUT NOISE CONTRIBUTION     %
%-------------------------------------------------------------------------%
r=s_sc.*c_sc;
m_fil_r = filter(a_sc,b,r);
M_fil_r = fftshift(fft(m_fil_r,N_sc));
figure();
plot(m_fil_r)
xlabel('The Message Signal After COHERANT DETECTOR WITHOUT NOISE CONTRIBUTION After Filtiring in Time Domain ');
grid on
figure
plot(f,real(M_fil_r))
xlabel('The Message Signal After COHERANT DETECTOR WITHOUT NOISE CONTRIBUTION After Filtiring in Frequancy Domain ');
grid on
%-------------------------------------------------------------------------%
%                    DEMODULATION IN PRESENCE OF NOISE                    %
%-------------------------------------------------------------------------%
noisy_signal_sc=awgn(s_sc,10);
r_n=noisy_signal_sc.*c_sc;
m_fil_r_n = filter(a,b,r_n);
M_fil_r_n = fftshift(fft(m_fil_r_n,N));
figure();
plot(m_fil_r_n)
xlabel('The Message Signal After COHERANT DETECTOR WITH NOISE CONTRIBUTION After Filtiring in Time Domain ');
grid on
[psnr,mse,maxerr,l_n] = measerr(s_sc,m_fil_r_n);%\
%-------------------------------------------------------------------------%
%   DETECTION USING COHERANT DETECTOR WITH Frequnacy Error                %
%-------------------------------------------------------------------------%
c_f= cos(2*pi*(fs_sc+0.001)*t_sc);
r_n1=noisy_signal_sc.*c_f;
m_fil_r_n1 = filter(a,b,r_n1);
M_fil_r_n1 = fftshift(fft(m_fil_r_n1,N));
figure();
plot(m_fil_r_n1)
xlabel('The Message Signal After COHERANT DETECTOR WITH Frequnacy Error  After Filtiring in Time Domain ');
grid on
[psnr,mse,maxerr,l_n1] = measerr(s,m_fil_r_n1);
%-------------------------------------------------------------------------%
%   DETECTION USING COHERANT DETECTOR WITH Phase Error                    %
%-------------------------------------------------------------------------%

c_p= cos((2*pi*fs_sc*t_sc)+10);
r_n2=noisy_signal_sc.*c_p;
m_fil_r_n2 = filter(a,b,r_n2);
M_fil_r_n2 = fftshift(fft(m_fil_r_n2,N));
figure();
plot(m_fil_r_n2)
xlabel('The Message Signal After COHERANT DETECTOR WITH Phase Error  After Filtiring in Time Domain ');

grid on
[psnr,mse,maxerr,l_n2] = measerr(s,m_fil_r_n2);
%%%%%%%%%%%%%%%%% power efficiency%%%%%%%%%%%%
v_sc = var(m_up);
power_efficiency = 100 * (v_sc / v_sc);

