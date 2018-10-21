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

value_qz_v  = myQuantize2( rampe_v, R_targ_n, 'with', [-1 1], rounding_method_s );
value_qnz_v = myQuantize2( rampe_v, R_targ_n, 'without', [-1 1], rounding_method_s );

figure(1)
hold off
plot( rampe_v, rampe_v, 'black', 'DisplayName', 'Original value')
hold on
plot( rampe_v, value_qz_v, 'b', 'DisplayName', 'Quantizer output with zero')
plot( rampe_v, value_qnz_v, 'r', 'DisplayName', 'Quantizer output without zero' )
plot( value_z_v, value_z_v, 'bo' )
plot( value_nz_v, value_nz_v, 'ro' )
grid on
legend show
xlabel('Input')
ylabel('Output')
title('Quantizer response')

%% LOAD SIGNAL
% audioread

%% QUANTIZE

%% PLOT

%% PLAY
% sound( sig_v, Fs)
% pause(numel(sig_v)/Fs)

