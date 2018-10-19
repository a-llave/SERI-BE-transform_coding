%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Description   Easy uniform scalar quantizer
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
    %----------------------------------------------------------------------
    %INPUTS
    % sigIn_m:              <mat> Input signal
    % R_targ_n:             <int> Number of bits per sample target
    % zero_s:               <string> (optional: default 'with') 'with' or 'without' zero state quantizer
    % range_v:              <vector 1x2> (optional: default [-1 1]) possible value range
    % rounding_method_s:    <string> (optional: default 'round') 'round', 'floor' and 'ceil'
    %
    %OUTPUTS
    % sigOut_m:             <mat> quantized signal
    %
    %EXAMPLE
    % [ sigOut_m ] = myQuantize2( sigIn_m, 8 )
    % [ sigOut_m ] = myQuantize2( sigIn_m, 8, 'with', [-1 1], 'round' )
    %
    %----------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ sigOut_m ] = myQuantize2( sigIn_m, R_targ_n, zero_s, range_v, rounding_method_s )
%% Input check
if nargin < 5
    rounding_method_s = 'round';
end
if nargin < 4
    range_v = [-1 1];
end
if nargin < 3
    zero_s = 'with';
end
if nargin < 2
    error('Need at least 2 input arguments');
end

assert( R_targ_n == int8(R_targ_n), 'ERROR: Second input argument have to be an integer' )

if R_targ_n <= 0
    sigOut_m = zeros( size(sigIn_m) );
    return;
end

%% Prepare
switch zero_s
    case 'with'
        value_v   = linspace(range_v(1), range_v(2)-2/2^R_targ_n, 2^R_targ_n);
    case 'without'
        value_v  = linspace(range_v(1), range_v(2), 2^R_targ_n);
    otherwise
        value_v   = linspace(range_v(1), range_v(2)-2/2^R_targ_n, 2^R_targ_n);
end

%% Process
switch rounding_method_s
    case 'round'
        TMP         = bsxfun(@(x,y) abs(x-y), sigIn_m(:), reshape(value_v,1,[]));
    case 'floor'
        TMP         = bsxfun(@(x,y) abs(x-y), sigIn_m(:)-1/2^R_targ_n, reshape(value_v,1,[]));
    case 'ceil'
        TMP         = bsxfun(@(x,y) abs(x-y), sigIn_m(:)+1/2^R_targ_n, reshape(value_v,1,[]));
    otherwise
        TMP         = bsxfun(@(x,y) abs(x-y), sigIn_m(:)+1/2^R_targ_n, reshape(value_v,1,[]));
end
[~, idxB]   = min(TMP,[],2);
sigOut_m    = value_v(idxB);

end % function



