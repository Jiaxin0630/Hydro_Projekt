function W = gaussian(L,cap)

% GAUSSIAN delivers the spherical harmonic coefficients of a gaussian
% smoothing filter. The coefficients are calculates according to Wahr et.
% al. (1998) equation (34) and Swenson and Wahr equation (34)
%
% How:      Wl = gaussian(L,cap)
%
% Input:    L    [1 x 1]    maximum degree
%           cap  [1 x 1]    half width of Gaussian smoothing function [km]
%
% Output:   Wl   [L+1 x 1]  smoothing coefficients
%
%--------------------------------------------------------------------------
% Weigelt, GI Stuttgart                                            15.02.07
%--------------------------------------------------------------------------

% Input checking
% -------------------------
error(nargchk(2,2,nargin));
if length(L)   ~= 1, error('Degree must be scalar.'); end
if L < 2, error('Maximum degree must be higher than 2.'); end
if length(cap) ~= 1, error('Cap size must be scalar.'); end

% Calculation
% --------------------------
W = zeros(L+1,1);
b = log(2)/(1-cos(cap/6371));

% recursive calculation of the weighting coefficients
W(1) = 1;
W(2) = 10^(log10( (1+exp(-2*b))/(1-exp(-2*b)) - 1/b ));
for i = 2:L
    W(i+1) = W(i-1) - (2*(i-1)+1)/b * W(i);
    if W(i+1) > W(i) || W(i+1) < 0, W(i+1) = 0; end
end

    

