clear all
close all
clc

%%
addpath('../resources')

%% Param

% WOLA
win_size_n  = 256; % window size (sample)
hop         = floor(win_size_n/2); % overlap (sample)
nb_FFT_n    = 256; % FFT sample number
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
plot_frame_b    = false;
plot_greedy_b   = false;

%% LOAD AUDIO FILE
[voiceOrig_v, Fs]   = audioread('artaud_16k.wav');
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

exe_time_opt_v = zeros(1,numFrame_n);
exe_time_grd_v = zeros(1,numFrame_n);
exe_time_hsm_v = zeros(1,numFrame_n);

iterLimit_n = 100;

for frm_id = 1:numFrame_n
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
    coeffQFMean_v   = myQuantize2( coeff_v, R0_n );
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % ---- FREQUENCY DOMAIN QUANTIZING - OPTIMAL ALLOCATION
    tic % Start execution time counting
    nbBit_v         = zeros(1, nb_freq_n*2);
    maskZero_v      = true(1, nb_freq_n*2);
    maskZero_v(PowerCoeff_v == 0) = false;
    
    for ii = 1:iterLimit_n
        % COMPUTE OPTIMAL BITS ALLOCATION
        nbBit_v(maskZero_v)     = R0_n + 0.5 * log2( PowerCoeff_v(maskZero_v) / geomean(PowerCoeff_v(maskZero_v)) ); % Optimal bit allocation
          
        % Check if there is negative value
        if sum(nbBit_v < 0) ~= 0
            maskZero_v( nbBit_v < 0 ) = false;
            nbBit_v( nbBit_v < 0 )  = 0; % check non-negative value
        else
            break;
        end
        
%         figure(100)
%         hold off
%         plot( nbBit_v )
%         hold on
%         plot( PowerCoeff_v ./ max(PowerCoeff_v) * R_targ_n )
%         plot( maskZero_v * R_targ_n )
%         pause
        
        if ii == iterLimit_n
            warning(sprintf('The condition was not achieved under %i iterations',iterLimit_n));
        end
    end % iterLimit_n
    
    nbBit_v         = floor( nbBit_v );
    coeffQFopt_v    = zeros(1, numel(coeff_v));
    for jj = 1:nb_freq_n
        coeffQFopt_v(jj)             = myQuantize2( real(sigTF_v(jj)), nbBit_v(jj) );
        coeffQFopt_v(jj+nb_freq_n)   = myQuantize2( imag(sigTF_v(jj)), nbBit_v(jj+nb_freq_n) );
    end
    exe_time_opt_v(frm_id) = toc*1000; % store execution time

    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % ---- FREQUENCY DOMAIN QUANTIZING - modified Huang-Schulteiss algorithm
    tic % Start execution time counting
    nbBit_v         = zeros(1, nb_freq_n*2);
    maskZero_v      = true(1, nb_freq_n*2);
    maskZero_v(PowerCoeff_v == 0) = false;
    
    for ii = 1:iterLimit_n
        % COMPUTE OPTIMAL BITS ALLOCATION
        nbBit_v(maskZero_v)     = R0_n + 0.5 * log2( PowerCoeff_v(maskZero_v) / geomean(PowerCoeff_v(maskZero_v)) ); % Optimal bit allocation
          
        % Check if there is negative value
        if sum(nbBit_v < 0) ~= 0
            maskZero_v( nbBit_v < 0 ) = false;
            nbBit_v( nbBit_v < 0 )  = 0; % check non-negative value
        else
            break;
        end
        
%         figure(100)
%         hold off
%         plot( nbBit_v )
%         hold on
%         plot( PowerCoeff_v ./ max(PowerCoeff_v) * R_targ_n )
%         plot( maskZero_v * R_targ_n )
%         pause
        
        if ii == iterLimit_n
            warning(sprintf('The condition was not achieved under %i iterations',iterLimit_n));
        end
    end % iterLimit_n
    
    nbBit_old_v = nbBit_v;
    nbBit_v     = floor( nbBit_v );
    
    coeffQFhsm_v    = zeros(1, numel(coeff_v));
    for jj = 1:nb_freq_n
        coeffQFhsm_v(jj)             = myQuantize2( real(sigTF_v(jj)), nbBit_v(jj), 'with', [-1 1] );
        coeffQFhsm_v(jj+nb_freq_n)   = myQuantize2( imag(sigTF_v(jj)), nbBit_v(jj+nb_freq_n), 'with', [-1 1] );
    end
    
%     nb_bit_2alloc_n = R_targ_n * nb_freq_n * 2 - sum( nbBit_v );
%     coeffQFhsm_v    = zeros(1, numel(coeff_v));
%     for ii = 1:nb_bit_2alloc_n
%         for jj = 1:nb_freq_n
%             if nbBit_v(jj) == 0
%                 coeffQFhsm_v(jj)             = 0;
%             else
%                 coeffQFhsm_v(jj)             = myQuantize2( real(sigTF_v(jj)), nbBit_v(jj), 'with', [-1 1] );
%             end
%             if nbBit_v(jj+nb_freq_n) == 0
%                 coeffQFhsm_v(jj+nb_freq_n)   = 0;
%             else
%                 coeffQFhsm_v(jj+nb_freq_n)   = myQuantize2( imag(sigTF_v(jj)), nbBit_v(jj+nb_freq_n), 'with', [-1 1] );
%             end
%         end
%         if max( nbBit_v ) >= R_init_n
%             break;
%         end
%         mask_coeff_v        = coeff_v ~= 0;
%         err_bitalloc_v      = zeros(1, numel( coeff_v ));
%         err_bitalloc_v(mask_coeff_v)      = abs(coeff_v(mask_coeff_v) - coeffQFhsm_v(mask_coeff_v));
%         [~, idx_n]          = max( err_bitalloc_v );
%         nbBit_v( idx_n )    = nbBit_v( idx_n ) + 1;
%         
% %         if mod(frm_id,1) == 0
% %             figure(140)
% %             subplot(1,3,1)
% %             hold off
% %             plot( err_bitalloc_v )
% %             hold on
% %             plot( [idx_n idx_n], [0 max(err_bitalloc_v)] )
% %             title( num2str( R_targ_n * nb_freq_n * 2 - sum( nbBit_v ) ) )
% %             subplot(1,3,2)
% %             hold off
% %             plot( nbBit_v )
% %             subplot(1,3,3)
% %             hold off
% %             plot( [real(sigTF_v) imag(sigTF_v)] )
% %             hold on
% %             plot( spcQFMean_v )
% %             pause(0.1)
% %         end
% 
%     end % ii
    exe_time_hsm_v(frm_id) = toc*1000; % store execution time

%     coeffQFhsm_v = coeffQFopt_v;
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- FREQUENCY DOMAIN QUANTIZING - greedy algorithm
    tic % Start execution time counting
    % ---- Init
    nbBit_v = zeros(1, numel( PowerCoeff_v ));
    D_v     = PowerCoeff_v;
    % ---- Process
    for ii = 1:iterLimit_n
        nbBitReste_n    = nb_freq_n * 2 * R0_n - sum( nbBit_v );
        if nbBitReste_n <= 0
            break;
        end
        [ ~, idx_max_n ]    = max( D_v );
        nbBit_v(idx_max_n)  = nbBit_v(idx_max_n) + 1;
        D_v( idx_max_n )    = D_v(idx_max_n) / 4;
        if ii == iterLimit_n
            warning(sprintf('The condition was not achieved under %i iterations',iterLimit_n));
            nbBit_v( nbBit_v < 0 )      = 0;
        end
        if plot_greedy_b && sum( abs(sig_v > 0.1) > 0.1 ) ~= 0
            figure(1)
            subplot(2,1,1)
            hold off
            xlim([1 numel(nbBit_v)])
            yyaxis left
            plot( D_v )
            ylim([0 0.2])
            yyaxis right
            plot( nbBit_v )
            hold on
            plot( PowerCoeff_v, 'g' )
            ylim([0 1])
            subplot(2,1,2)
            hold off
            plot( sig_v )
            ylim([-0.6 0.6])
            xlim([1 win_size_n])
            pause(0.001)
        end
        
    end % ii
    
    coeffQFgrd_v = zeros(1,nb_freq_n*2);
    for kk = 1:nb_freq_n*2
        coeffQFgrd_v(kk)   = myQuantize2( coeff_v(kk), nbBit_v(kk), 'with', [-1 1] );
    end % kk
    exe_time_grd_v(frm_id) = toc*1000; % store execution time
    
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    % >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    % ---- Coeff to spectrum
    spcQFMean_v = coeffQFMean_v(1:nb_freq_n) + 1i * coeffQFMean_v(nb_freq_n+1:end);
    spcQFopt_v  = coeffQFopt_v(1:nb_freq_n) + 1i * coeffQFopt_v(nb_freq_n+1:end);
    spcQFhsm_v  = coeffQFhsm_v(1:nb_freq_n) + 1i * coeffQFhsm_v(nb_freq_n+1:end);
    spcQFgrd_v  = coeffQFgrd_v(1:nb_freq_n) + 1i * coeffQFgrd_v(nb_freq_n+1:end);
    % ---- Symmetrize spectrum
    fftQFMean_v = [ spcQFMean_v fliplr(conj(spcQFMean_v(2:idx_end_n))) ];
    fftQFopt_v  = [ spcQFopt_v fliplr(conj(spcQFopt_v(2:idx_end_n))) ];
    fftQFhsm_v  = [ spcQFhsm_v fliplr(conj(spcQFhsm_v(2:idx_end_n))) ];
    fftQFgrd_v  = [ spcQFgrd_v fliplr(conj(spcQFgrd_v(2:idx_end_n))) ];
    % ---- Freq to Time domain
    frmQFMean_v = ifft( fftQFMean_v .* sqrt(win_size_n/2), nb_FFT_n );
    frmQFopt_v  = ifft( fftQFopt_v .* sqrt(win_size_n/2), nb_FFT_n );
    frmQFhsm_v  = ifft( fftQFhsm_v .* sqrt(win_size_n/2), nb_FFT_n );
    frmQFgrd_v  = ifft( fftQFgrd_v .* sqrt(win_size_n/2), nb_FFT_n );
    % ---- Overlap-add signals
    sigQFMean_v(in_n:out_n)   = sigQFMean_v(in_n:out_n) + win_synthesis .* frmQFMean_v(1:win_size_n);
    sigQFopt_v(in_n:out_n)    = sigQFopt_v(in_n:out_n)  + win_synthesis .* frmQFopt_v(1:win_size_n);
    sigQFhsm_v(in_n:out_n)    = sigQFhsm_v(in_n:out_n)  + win_synthesis .* frmQFhsm_v(1:win_size_n);
    sigQFgrd_v(in_n:out_n)    = sigQFgrd_v(in_n:out_n)  + win_synthesis .* frmQFgrd_v(1:win_size_n);
    
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
        plot( frmQFopt_v, 'displayname', 'Opt floor freq Q' )
        plot( frmQFgrd_v, 'displayname', 'Greedy freq Q' )
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
noiseFopt_v     = voiceOrig_v' - sigQFopt_v;
noiseFhsm_v     = voiceOrig_v' - sigQFhsm_v;
noiseFgrd_v     = voiceOrig_v' - sigQFgrd_v;

SNR_T_f     = mag2db( rms( voiceOrig_v ) / rms( noiseT_v ) );
SNR_FMean_f = mag2db( rms( voiceOrig_v ) / rms( noiseFMean_v ) );
SNR_Fopt_f  = mag2db( rms( voiceOrig_v ) / rms( noiseFopt_v ) );
SNR_Fhsm_f  = mag2db( rms( voiceOrig_v ) / rms( noiseFhsm_v ) );
SNR_Fgrd_f  = mag2db( rms( voiceOrig_v ) / rms( noiseFgrd_v ) );

disp('------------------------ SNR ----------------------------')
fprintf('SNR (time domain quantizing): %.2f dB\n', SNR_T_f);
fprintf('SNR (naive freq domain quantizing): %.2f dB\n', SNR_FMean_f);
fprintf('SNR (optimal freq domain quantizing): %.2f dB\n', SNR_Fopt_f);
fprintf('SNR (HSM freq domain quantizing): %.2f dB\n', SNR_Fopt_f);
fprintf('SNR (greedy freq domain quantizing): %.2f dB\n', SNR_Fgrd_f);

%% EXECUTION TIME
disp('------------------------ Exe. time ----------------------------')
fprintf('Exe time (optimal freq domain quantizing):\n');
fprintf('   mean: %.2f ms\n', mean( exe_time_opt_v ));
fprintf('   standard deviation: %.2f ms\n', std( exe_time_opt_v ));
fprintf('Exe time (greedy freq domain quantizing):\n');
fprintf('   mean: %.2f ms\n', mean( exe_time_grd_v ));
fprintf('   standard deviation: %.2f ms\n', std( exe_time_grd_v ));

%% PLOT
figure(2)
hold off
plot( noiseT_v, 'DisplayName', 'Time domain Q noise' )
hold on
plot( noiseFMean_v, 'DisplayName', 'Naive Freq domain Q noise')
plot( noiseFopt_v, 'DisplayName', 'Optimal Freq domain Q noise')
plot( noiseFgrd_v, 'DisplayName', 'Greedy Freq domain Q noise')
legend show
xlabel('Time (sample)')
ylabel('Magnitude')
title('Quantizing noise')


figure(3)
subplot(2,1,1)
hold off
plot( voiceOrig_v, 'DisplayName', 'Original' )
hold on
plot( sigQT_v, 'DisplayName', 'Time domain quantizing' )
plot( sigQFMean_v, 'DisplayName', 'Naive freq domain quantizing' )
plot( sigQFopt_v, 'DisplayName', 'Optimal floor freq domain quantizing' )
plot( sigQFgrd_v, 'DisplayName', 'Greedy freq domain quantizing' )
plot(axis_4plot_v, FlatCoef_v, 'displayname', 'Flatness coefficient' )
legend show
title('Signals')
subplot(2,1,2)
hold off
plot( axis_4plot_v, exe_time_opt_v, 'displayname', 'Optimal floor algo exe. time' )
hold on
plot( axis_4plot_v, exe_time_grd_v, 'displayname', 'Greedy algo exe. time' )
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
pause(numel(voiceOrig_v)/Fs+0.5)
sound(sigQFopt_v,Fs)
pause(numel(voiceOrig_v)/Fs+0.5)
sound(sigQFhsm_v,Fs)
pause(numel(voiceOrig_v)/Fs+0.5)
sound(sigQFgrd_v,Fs)


