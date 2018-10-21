%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description   Exercice 3: transform coding
%
% Author        A. Llave - CentraleSupélec
%
% Date          19 October 2018
%
% Version       1.0
%
%
%   History
% Date      Author  Version     Content
% 18/10/19   ALl     1.0        Initialization
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% PATH
addpath('../resources')

%% Parameters

% WOLA
win_size_n  = 255; % window size (sample)
hop         = floor(win_size_n/2); % overlap (sample)
nb_FFT_n    = win_size_n; % FFT sample number
nb_freq_n   = floor(nb_FFT_n/2)+1; % frequency number
win_analysis    = sqrt(hann(win_size_n))'; % Analysis window
win_synthesis   = sqrt(hann(win_size_n))'; % Synthesis window

if mod(nb_FFT_n, 2)
    idx_end_n = nb_freq_n;
else
    idx_end_n = nb_freq_n-1;
end

% Quantizer
R0_n    = 4;

% PLOT
plot_frame_b    = true;

%% LOAD AUDIO FILE
[voiceOrig_v, Fs]   = audioread('parole_16k.wav');
voiceOrig_v         = voiceOrig_v(1:Fs*5, 1);

numSample_n         = numel( voiceOrig_v );
R_init_n            = 16;

freq_v  = linspace(0,Fs/2,nb_freq_n);

%% WOLA
numFrame_n      = floor(numSample_n/(win_size_n-hop)) - 1;

axis_4plot_v   = zeros(1,numFrame_n);
FlatCoef_v     = zeros(1,numFrame_n);
sigQFMean_v    = zeros(1,numel(voiceOrig_v));
sigQFopt_v     = zeros(1,numel(voiceOrig_v));
sigQFhsm_v     = zeros(1,numel(voiceOrig_v));
sigQFgrd_v     = zeros(1,numel(voiceOrig_v));

exe_time_nve_v = zeros(1,numFrame_n);
exe_time_opt_v = zeros(1,numFrame_n);
exe_time_grd_v = zeros(1,numFrame_n);
exe_time_hsm_v = zeros(1,numFrame_n);

for frm_id = 1:numFrame_n
    if mod(frm_id, 50) == 0
        disp('--------------------------')
        fprintf('Frame %.0f/%.0f\n', frm_id, numFrame_n);
    end
    % DEFINE START AND END SAMPLES OF THE FRAME
    in_n    = (frm_id-1)*hop+1;
    out_n   = (frm_id-1)*hop+win_size_n;
    axis_4plot_v(frm_id) = in_n;
    
    % SELECT FRAME FROM ORIGINAL SIGNAL AND WINDOW IT
    sig_v   = voiceOrig_v(in_n:out_n)' .* win_analysis;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Compute PSD, Flatness coeff
    sigTF_v             = fft( sig_v, nb_FFT_n ) ./ sqrt(win_size_n/2);
    sigTF_v             = sigTF_v(1:nb_freq_n);
    FlatCoef_v(frm_id)  = geomean(abs(sigTF_v).^2) / mean(abs(sigTF_v).^2);
    
    coeff_v             = [real(sigTF_v) imag(sigTF_v)];
    PowerCoeff_v        = coeff_v.^2;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- FREQUENCY DOMAIN QUANTIZING - NAIVE ALLOCATION
    tic % Start execution time counting
    coeffQFMean_v   = myQuantize2( coeff_v, R0_n );
    exe_time_nve_v(frm_id) = toc*1000; % store execution time
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % ---- FREQUENCY DOMAIN QUANTIZING - OPTIMAL ALLOCATION FLOORED
    tic % Start execution time counting
    % -------------------------------
    % ---- INSERT YOUR CODE HERE ----
    % -------------------------------
    exe_time_opt_v(frm_id) = toc*1000; % store execution time

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % ---- FREQUENCY DOMAIN QUANTIZING - modified Huang-Schulteiss algorithm
    tic % Start execution time counting
    % -------------------------------
    % ---- INSERT YOUR CODE HERE ----
    % -------------------------------
    exe_time_hsm_v(frm_id) = toc*1000; % store execution time
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % ---- FREQUENCY DOMAIN QUANTIZING - greedy algorithm
    tic % Start execution time counting
    % -------------------------------
    % ---- INSERT YOUR CODE HERE ----
    % -------------------------------
    exe_time_grd_v(frm_id) = toc*1000; % store execution time
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    % ---- Coeff to spectrum
    spcQFMean_v = coeffQFMean_v(1:nb_freq_n) + 1i * coeffQFMean_v(nb_freq_n+1:end);
    % -- Insert your code here (optimal floored) --
    % -- Insert your code here (Huang-Schulteiss modified) --
    % -- Insert your code here (Greedy algo) --
    % ---- Symmetrize spectrum
    fftQFMean_v = [ spcQFMean_v fliplr(conj(spcQFMean_v(2:idx_end_n))) ];
    % -- Insert your code here (optimal floored) --
    % -- Insert your code here (Huang-Schulteiss modified) --
    % -- Insert your code here (Greedy algo) --
    % ---- Freq to Time domain
    frmQFMean_v = ifft( fftQFMean_v .* sqrt(win_size_n/2), nb_FFT_n );
    % -- Insert your code here (optimal floored) --
    % -- Insert your code here (Huang-Schulteiss modified) --
    % -- Insert your code here (Greedy algo) --
    % ---- Overlap-add signals
    sigQFMean_v(in_n:out_n)   = sigQFMean_v(in_n:out_n) + win_synthesis .* frmQFMean_v(1:win_size_n);
    % -- Insert your code here (optimal floored) --
    % -- Insert your code here (Huang-Schulteiss modified) --
    % -- Insert your code here (Greedy algo) --
    
    % Plot
    if plot_frame_b && mod(frm_id,10) == 0
        numRow_n = 1;
        numCol_n = 2;
        
        figure(1)
        hold off
        suptitle(sprintf('Wiener entropy: %0.2f', FlatCoef_v(frm_id)))
        subplot(numRow_n,numCol_n,1) %%%%%%%%% 1
        plot( sig_v, 'displayname', 'Original signal' )
        hold on
        plot( frmQFMean_v, 'displayname', 'naive freq Q' )
        legend show
        xlim([1 win_size_n])
        ylim(0.5*[-1 1])
        xlabel('Sample')
        ylabel('Magnitude')
        title('Signal')
        subplot(numRow_n,numCol_n,2) %%%%%%%%% 2
        plot( freq_v, abs(sigTF_v).^2, 'DisplayName', 'Inst. power spectrum' )
        grid on
        xlim([freq_v(1) freq_v(end)])
        ylim(0.5*[0 1])
        xlabel('Frequency')
        ylabel('Magnitude')
        title('Power Spectrum Density')

        pause(0.01)
    end % end
end % ll

%% TIME DOMAIN QUANTIZING
sigQT_v = myQuantize2( voiceOrig_v, R0_n );

%% SNR
noiseT_v        = voiceOrig_v' - sigQT_v;
noiseFMean_v    = voiceOrig_v' - sigQFMean_v;
% -- Insert your code here (optimal floored) --
% -- Insert your code here (Huang-Schulteiss modified) --
% -- Insert your code here (Greedy algo) --

SNR_T_f     = mag2db( rms( voiceOrig_v ) / rms( noiseT_v ) );
SNR_FMean_f = mag2db( rms( voiceOrig_v ) / rms( noiseFMean_v ) );
% -- Insert your code here (optimal floored) --
% -- Insert your code here (Huang-Schulteiss modified) --
% -- Insert your code here (Greedy algo) --

disp('------------------------ SNR ----------------------------')
fprintf('SNR (time domain quantizing): %.2f dB\n', SNR_T_f);
fprintf('SNR (naive freq domain quantizing): %.2f dB\n', SNR_FMean_f);
% -- Insert your code here (optimal floored) --
% -- Insert your code here (Huang-Schulteiss modified) --
% -- Insert your code here (Greedy algo) --

%% EXECUTION TIME
disp('------------------------ Exe. time ----------------------------')
fprintf('Exe time (naive freq domain quantizing):\n');
fprintf('   mean: %.2f ms\n', mean( exe_time_nve_v ));
fprintf('   standard deviation: %.2f ms\n', std( exe_time_nve_v ));
% -- Insert your code here (optimal floored) --
% -- Insert your code here (Huang-Schulteiss modified) --
% -- Insert your code here (Greedy algo) --

%% PLOT
figure(2)
hold off
plot( noiseT_v, 'DisplayName', 'Time domain Q noise' )
hold on
plot( noiseFMean_v, 'DisplayName', 'Naive Freq domain Q noise')
legend show
xlabel('Time (sample)')
ylabel('Magnitude')
title('Quantizing noise')


figure(3)
hold off
plot( voiceOrig_v, 'DisplayName', 'Original' )
hold on
plot( sigQT_v, 'DisplayName', 'Time domain quantizing' )
plot( sigQFMean_v, 'DisplayName', 'Naive freq domain quantizing' )
plot(axis_4plot_v, FlatCoef_v, 'displayname', 'Flatness coefficient' )
legend show
title('Signals')

figure(4)
hold off
plot( axis_4plot_v, exe_time_nve_v, 'displayname', 'Naive algo exe. time' )
hold on
plot( axis_4plot_v, exe_time_opt_v, 'displayname', 'Optimal floor algo exe. time' )
legend show
xlabel('Time (sample)')
ylabel('Execution time (ms)')
title('Execution time for various algo')

%% Sound
sound(voiceOrig_v,Fs)
pause(numel(voiceOrig_v)/Fs+0.5)
sound(sigQT_v,Fs)
pause(numel(voiceOrig_v)/Fs+0.5)
sound(sigQFMean_v,Fs)
% -- Insert your code here (optimal floored) --
% -- Insert your code here (Huang-Schulteiss modified) --
% -- Insert your code here (Greedy algo) --
