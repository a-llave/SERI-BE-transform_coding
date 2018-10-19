% close all
clear all

%% PARAMETERS
R_init_v = 16;
R_targ_n = 3;

rounding_method_s = 'round';

%% SIGNAL TEST GENERATION
rampe_v     = linspace(-1, 1, 2^(R_init_v));

value_nz_v  = linspace(-1, 1, 2^R_targ_n);
value_z_v   = linspace(-1, 1-2/2^R_targ_n, 2^R_targ_n);

value_qz_v = myQuantize2( rampe_v, R_targ_n, 'with', [-1 1], rounding_method_s );
value_qnz_v = myQuantize2( rampe_v, R_targ_n, 'without', [-1 1], rounding_method_s );

figure(1)
hold off
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

%% LOAD SIGNAL

%% QUANTIZE

%% PLOT

%% PLAY
% sound( sig_original_v, Fs)
% pause(duration_f+0.5)

