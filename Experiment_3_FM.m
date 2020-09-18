close all; % to close all windows.
clear all; %to clear all values.
%% 
fs = 48000; %sampling frequency.
y = audiorecorder(fs, 8, 1); % recording an audio with the given sampling frequency fs, 8 bits sample size and 1 channel.
recordblocking(y, 10); % recording 10 seconds.
A = getaudiodata(y); % saving the recorded file data in a matrix form.
time = 0:1/fs:(10-(1/fs));
figure; plot(time,A); title('The Audio Message In Time Domain');% plotting the audio in time domain.
s = size (A); % size of audio :  [480000  1].
f= -fs./2 : fs./(s-1): fs./2; % establishing the frequency axis.
z= fftshift(fft(A)); % getting the spectrum of the audio msg.
figure; plot (f,abs(z)); title('Spectrum Of The Message');% plotting the spectrum of the audio msg.
%-----------------
%% 
% Low Pass Filter
%-----------------
M2 = ones(8000*10,1); %rect of ones = ( BW = 8KH * time = 10s )
M1= zeros((48000*10-8000*10)/2,1); % zeros = (All the samples - ones )/2;
 M = [M1 ;M2; M1]; %concatinating zeros with ones to form the rect.
 H = z .* M; %Getting the msg after the filter.
 figure; plot (f,abs(H));title('The Spectrum Of The Message After LPF'); %plotting the msg after the filter.
 t = ifft(H); % getting the fourier inverse to get the msg in time domain.
 H1 = real (t); % getting the absolute of the msg in time domain to be plotted.
 figure;  plot (time,H1); title('The Time Domain Of The Message After LPF');% plotting the msg in time domain after LPF.
%sound ( abs(t),48000,8); % here you can hear the msg after LPF and you'll find it the same with low noise.
%--------------------
%% 
% Calculating the SME
%--------------------
for i=1:480000;
     SME1 = 0; %initializing the SME1.
     N1 = ( real(H(i)) - A(i) )^2; 
     SME1 = SME1 + N1; 
 end
 SME1 = SME1/480000; %getting the error. 
 %% Generating the NBFM
 fc = 100000; % carrier frequency.
 fs_NBFM = 500000; %sampling frequency of NBFM. 
 wc = 2*pi*fc; 
 t1= 0 : 1/fs_NBFM :(10- (1/fs_NBFM)); %time sampling.
 H_re = resample(H1,500,48); % resampling the signal after LPF.
 NBFM =  cos (wc*t1) -  (cumsum(H_re).'.* sin(wc*t1)); %NBFM equation
 
 s1 = size (NBFM); % size of audio :  [480000  1].
fn= -fs_NBFM./2 : fs_NBFM./(s1-1): fs_NBFM./2; % establishing the frequency axis.
O = fftshift(fft(NBFM)); % getting the spectrum of the audio msg.
%figure; plot (fn,abs(NBFM)); title('Spectrum Of The NBFM');% plotting the spectrum of the audio msg.
 
 %figure; plot (fn,abs(NBFM)); title('Spectrum Of NBFM'); 
ED = abs(hilbert(NBFM)); %The envelope detector.
D =diff(ED); %the diff of ED.
msg = [D 0]; % concatinating the diff with zero to maintain same matrix size.
msg_re = resample(msg,48,500); %resampling to get clear audio.
sound (10*abs(msg_re),48000,8); %playing the audio with a gain 10 ..
%% Getting SME of NBFM
for i= 1:480000;
    SME2 = 0;
    N2 = ( real(msg_re(i)) - A(i))^2;
    SME2 = SME2 + N2;
end
SME2 = SME2/480000;
%% Adding The Noise.
NBFM1 = awgn (NBFM,0); %adding noise of 0 dB.
ED1 = abs(hilbert(NBFM1));
D1 =diff(ED1);
msg1 = [D1 0];
msg_re1 = resample(msg1,48,500);
% getting SME for msg with noise of 0 dB.
for i= 1:480000;
    SME3 = 0;
    N3 = ( real(msg_re1(i)) - A(i))^2;
    SME3 = SME3 + N3;
end
SME3 = SME3/480000;
% ----
%% Adding noise of 10 dB.
NBFM2 = awgn (NBFM,10); 
ED2 = abs(hilbert(NBFM2));
D2 =diff(ED2);
msg2 = [D2 0];
msg_re2 = resample(msg2,48,500);
% getting SME for msg with noise of 0 dB.
for i= 1:480000;
    SME4 = 0;
    N4 = ( real(msg_re2(i)) - A(i))^2;
    SME4 = SME4 + N4;
end
SME4 = SME4/480000;
%------
%% Adding noise of 20 dB.
 NBFM3 = awgn (NBFM,20);
ED3 = abs(hilbert(NBFM3));
D3 =diff(ED3);
msg3 = [D3 0];
msg_re3 = resample(msg3,48,500);
 % getting SME for msg with noise of 0 dB.
for i= 1:480000;
    SME5 = 0;
    N5 = ( real(msg_re3(i)) - A(i))^2;
    SME5 = SME5 + N5;
end
SME5 = SME5/480000;
