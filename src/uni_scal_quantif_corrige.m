% close all
clear all

%% PARAMETERS
R_init_v = 16;
R_targ_n = 4;

rounding_method_s = 'round';

%% SIGNAL TEST GENERATION
rampe_v     = linspace(-1, 1, 2^(R_init_v));

value_nz_v  = linspace(-1, 1, 2^R_targ_n);
value_z_v   = linspace(-1, 1-2/2^R_targ_n, 2^R_targ_n);

value_qz_v = myQuantize2( rampe_v, R_targ_n, 'with', [-1 1], rounding_method_s );
value_qnz_v = myQuantize2( rampe_v, R_targ_n, 'without', [-1 1], rounding_method_s );

figure,
subplot(1,2,1)
plot( rampe_v, rampe_v, 'black')
hold on
plot( rampe_v, value_qz_v, 'b')
plot( rampe_v, value_qnz_v, 'r' )
plot( value_z_v, value_z_v, 'bo' )
plot( value_nz_v, value_nz_v, 'ro' )
grid on
xlabel('Input')
ylabel('Output')
title('Quantizer response')
subplot(1,2,2)
hold on
plot( rampe_v, abs(value_qz_v-rampe_v), 'b' )
plot( rampe_v, abs(value_qnz_v-rampe_v), 'r' )
grid on
xlabel('Input')
ylabel('Error')
title('Quantification error')
suptitle([ 'rounding method: ' rounding_method_s])

%% LOAD SIGNAL
duration = 5;
[voiceOrig_v, Fs]   = audioread('artaud_16k.wav');
numSample_n         = Fs * duration;
voiceOrig_v         = voiceOrig_v(1:numSample_n,1);

%% QUANTIZE
sig_qz_v = myQuantize2( voiceOrig_v, R_targ_n, 'with', [-1 1], 'round' );
sig_qnz_v = myQuantize2( voiceOrig_v, R_targ_n, 'without', [-1 1], 'round' );

%% PLOT

figure,
plot( voiceOrig_v, 'displayname', 'Original' )
hold on
plot( sig_qz_v, 'displayname', 'Quantized with zero' )
plot( sig_qnz_v, 'displayname', 'Quantized without zero' )
legend show

%% PLAY
sound( voiceOrig_v, Fs)
pause(duration+0.5)
sound( sig_qz_v, Fs)
pause(duration+0.5)
sound( sig_qnz_v, Fs)
pause(duration+0.5)
