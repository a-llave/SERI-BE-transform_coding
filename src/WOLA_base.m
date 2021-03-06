clear all
close all
clc

%%
addpath('../resources')

%% Param

% WOLA
win_size_n  = 255; % window size (sample)
hop         = floor(win_size_n/2); % overlap (sample)
nb_FFT_n    = win_size_n; % FFT sample number
nb_freq_n   = nb_FFT_n/2+1; % frequency number
win_analysis    = sqrt(hann(win_size_n))'; % Analysis window
win_synthesis   = sqrt(hann(win_size_n))'; % Synthesis window

% Quantizer
R_targ_n    = 4;

% PLOT
plot_frame_b = true;

%% LOAD AUDIO FILE
% [voiceOrig_v, Fs]   = audioread('sweep_16k.wav');
[voiceOrig_v, Fs]   = audioread('parole_16k.wav');
voiceOrig_v         = voiceOrig_v(1:Fs*5, 1);

numSample_n         = numel( voiceOrig_v );
R_init_n            = 16;

freq_v  = linspace(0,Fs/2,nb_freq_n);

%% WOLA
numFrame_n      = floor(numSample_n/(win_size_n-hop))-1;

axis_4plot_v    = zeros(1,numFrame_n);
FlatCoef_v      = zeros(1,numFrame_n);
voice_v         = zeros(1,numel(voiceOrig_v));

for frm_id = 1:numFrame_n
    % DEFINE START AND END SAMPLES OF THE FRAME
    in_n    = (frm_id-1)*hop+1;
    out_n   = (frm_id-1)*hop+win_size_n;
    axis_4plot_v(frm_id) = in_n;
    
    % SELECT FRAME FROM ORIGINAL SIGNAL
    sig_v   = voiceOrig_v(in_n:out_n)' .* win_analysis;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- Compute PSD, Flatness coeff and power of PSD
    sigTF_v             = fft( sig_v, nb_FFT_n ) ./ sqrt(win_size_n/2);
    sigTF_v             = sigTF_v(1:nb_freq_n);
    FlatCoef_v(frm_id)  = 0;
    
    % Reconstruct total signals quant
    voice_v(in_n:out_n)   = voice_v(in_n:out_n) + win_synthesis .* sig_v;
    
    % Plot
    if plot_frame_b && mod(frm_id,10) == 0
        numRow_n = 1;
        numCol_n = 2;
        
        figure(1)
        suptitle(sprintf('Wiener entropy: %0.2f', FlatCoef_v(frm_id)))
        subplot(numRow_n,numCol_n,1) %%%%%%%%% 1
        plot( sig_v, 'displayname', 'Original signal' )
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

%% PLOT
figure(2)
hold off
plot( voiceOrig_v )
xlabel('Time (sample)')
title('Original signal')

%% Sound
sound(voiceOrig_v,Fs)









