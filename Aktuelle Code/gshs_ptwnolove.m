function f = gshs_ptw(field, lamRAD, phiRAD, r, a_E, varargin)
% GSHS_PTW calculates a pointwise global spherical harmonic synthesis. 
%
% f = gshs_ptw(field, lamRAD, phiRAD, r, a_E);
%
% IN:
%    field ....... gravity field in |c\s| or /s|c\ format
%    lamRAD ...... longitude [rad]                                    [n, 1]
%    phiRAD ...... latitude  [rad]                                    [n, 1]
%    r ........... radius    [m]                                      [1, 1]
%    a_E ......... semi major axis of the Earth rotational ellipsoid  [1, 1]
% OPTIONAL:
%    'max_lm' .... maximum degree/order (default: determined from field)
%    'quant' ..... field quantity; possible values:
%                  - 'potential' ... (default), potential [m^2/s^2], needs 'GM'
%                  - 'tr' .......... gravity disturbance [mGal], needs 'GM'
%                  - 'trr' ......... 2nd rad. derivative [E], needs 'GM'
%                  - 'none' ........ coefficients define the output
%                  - 'geoid' ....... geoid height [m]
%                  - 'dg', 'gravity' gravity anomaly [mGal], needs 'GM'
%                  - 'slope' ....... slope [arcsec]
%                  - 'water' ....... equivalent water thickness [m]
%                  - 'smd' ......... surface mass density [kg/m^2]
%    'sub_WGS84' . if set, subtracts WGS84 reference ellipsoid (default: true)
%    'GM' ........ geocentric gravitational constant GM
%    'legendre' .. (default: 'plm') Legendre functions algorithm:
%                  - 'plm' ... unstable (for d/o > ~1800) PLM
%                  - 'mex' ... X-number stabilized LEGENDRE_MEX
%    'waitbar' ... if set, shows a waitbar (default: false)
%    'chunklength' sets the length of phiRAD chuncks which are fed to plm
%                  (default: tries to calculate an optimal value)
% OUT:
%    f ......... field quantity                                      [n, 1]
%
% EXAMPLE: see SHBUNDLE/example/example_gshsptw.m
%
% USES:
%    normalklm, cs2sc, plm, Legendre_mex, lovenr, 
%    UBERALL/checkcoor, UBERALL/twaitbar
%
% SEE ALSO:
%    GSHS_, GSHS_GRID, REGIONALSHS, GSHA


% -------------------------------------------------------------------------
% project: SHBundle 
% -------------------------------------------------------------------------
% authors:
%    Matthias WEIGELT (MW), DoGE, UofC
%    Markus ANTONI (MA), GI, Uni Stuttgart 
%    Matthias ROTH (MR), GI, Uni Stuttgart
% -------------------------------------------------------------------------
% revision history:
%    2018-11-27: MA, extra warning if a reference field is subtracted
%                    (request of external users)
%    2015-09-04: MR, remove the autorotate code in plm, revert ghsh_ptw
%    2015-09-02: MR, remove bug of plm.m: for length(theRAD) == 1, result is 
%                    the transposed of what we need
%    2015-05-21: MR, brush up comments and help text
%    2015-05-20: MR, change waitbar to ASCII version (switched off by default)
%    2015-02-04: MR, add waitbar (can be switched off)
%    2014-11-05: MR, reprogram parameter interface (hence, rename function), 
%                    revise help text
%    2014-01-15: MR, revise help text and code
%    2013-02-13: MR, change function names, brush up comments
%    2013-01-30: MA, comments/removing of smoothing option
%    2013-01-23: MA, input in radian
%    2009-03-04: MW, change input co-latitude -> latitude
%    2005-10-25: MW, initial version
% -------------------------------------------------------------------------
% license:
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the  Free  Software  Foundation; either version 3 of the License, or
%    (at your option) any later version.
%  
%    This  program is distributed in the hope that it will be useful, but 
%    WITHOUT   ANY   WARRANTY;  without  even  the  implied  warranty  of 
%    MERCHANTABILITY  or  FITNESS  FOR  A  PARTICULAR  PURPOSE.  See  the
%    GNU General Public License for more details.
%  
%    You  should  have  received a copy of the GNU General Public License
%    along with Octave; see the file COPYING.  
%    If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------

%% define defaults and check optional parameters
defaultparams = {'max_lm', inf;
                 'quant', 'potential';
                 'sub_wgs84', true;
                 'gm', nan;
                 'legendre', 'plm';
                 'waitbar', false;
                 'chunklength', -1}; 
params = getopt(defaultparams, false, varargin);  

% check Legendre function algorithm
plm_func = false;
switch params.legendre
    case 'mex'
        if exist('Legendre_mex', 'file') == 3; % check if compiled Legendre_mex exists
            plm_func = @Legendre_mex;
        else
            warning('Legendre_mex is not compiled! Falling back to plm...');
            plm_func = @plm;
        end
    case 'plm'
        plm_func = @plm;
    otherwise
        error('Unknown Legendre function.');
end

% field size determination, rearrange field and subtract reference field
[row, col] = size(field);

if col == row              % field in |C\S|-format 
    field = cs2sc(field);  % convert to /S|C\-format
elseif col ~= 2 * row - 1
   error('Input "field" not in cs or sc format');
% else: field is already in /S|C\-format
end

if params.max_lm > (row - 1)     % desired max_lm is bigger than what the field provides
    params.max_lm = row - 1;     % adjust max_lm 
elseif params.max_lm < (row - 1) % if max_lm is smaller than what the field provides
    field = field(1:(params.max_lm + 1), (row - params.max_lm):(row + params.max_lm)); % adjust field
    row = size(field, 1); % update row
% else: everything is ok
end

if params.sub_wgs84
    field = field - cs2sc(full(normalklm(params.max_lm, 'wgs84')));
    if params.display_warning == true
        warning('A reference field (WGS84) is removed from your coefficients')
    end
end

% prepare the coordinates
[lamRAD, phiRAD, r] = checkcoor(lamRAD, phiRAD, r, a_E, 'pointwise');
theRAD              = (pi/2 - phiRAD);
num_coord           = numel(lamRAD);
didx                = 1:num_coord;

% prepare l and m
m       = (0:params.max_lm);
l       = m';
l_trans = l';

% try to calculate optimal chunklength (if not set)
if params.chunklength == -1
    params.chunklength = fix(120000 / params.max_lm);
end
if params.chunklength < 1 % catch a miscalculation (or user stupidity)
    params.chunklength = 1;
end

% apply transfer function
switch lower(params.quant)
    case 'none'
        tk = ones(size(l));                  % []
        tf = ones(size(l));         
        tl = l_trans;
    case 'potential' % default
        if isnan(params.gm), error('GM not defined!'); end;
        tk = params.gm ./ a_E;	             % [m^2/s^2]
        tf = ones(size(l));	
        tl = l_trans + 1;
    case 'geoid'
        tk = a_E;                            % [m]
        tf = ones(size(l));
        tl = l_trans - 1;
    case {'gravity', 'dg'}
        if isnan(params.gm), error('GM not defined!'); end;
        tk = params.gm / a_E / a_E * 1e5;    % [mGal]
        tf = l - 1;
        tl = l_trans + 2;
    case 'tr'
        if isnan(params.gm), error('GM not defined!'); end;
        tk = -params.gm / a_E / a_E * 1e5;   % [mGal]    
%         tf = l + 1;	
%         tl = l_trans + 2;
        tl = l_trans;
        tf=l+1;
    case 'trr'
        if isnan(params.gm), error('GM not defined!'); end;
        tk = params.gm / a_E / a_E / a_E * 1e9; % [E]
        tf = (l + 1) .* (l + 2);	
        tl = l_trans + 3;
    case 'slope'
        tk = ones(size(l));                  % [rad]
        tf = sqrt(l .* (l + 1));				
        tl = l_trans;
    case 'water'
        tk = a_E * 5.517 / 3;                % [m]
        tf = (2 .* l + 1) ./ (1 + lovenr(l));
        tl = l_trans- 1;
    case 'smd'
        tk = a_E * 5517 / 3;                 % [kg/m^2]
        tf = (2 .* l + 1) ./ (1 + lovenr(l));
        tl = l_trans - 1;
    otherwise
        error('Requested functional QUANT not available.')
end

field = field .* ((tk .* tf(:)) * ones(1, 2 * params.max_lm + 1));

%% calculation
f  = NaN(size(lamRAD));

% prepare index
idx    = all(~isnan([lamRAD theRAD r]), 2);
didx   = didx(idx);
lamRAD = lamRAD(idx);
theRAD = theRAD(idx);
r      = r(idx);

% calculation loop
parts = ceil(num_coord / params.chunklength);
if params.waitbar
    WAIT = twaitbar('init', [], 'gshs_ptw: calculating');
end;

for ridx = 1:parts
    idx     = (ridx - 1) * params.chunklength + 1:min(ridx * params.chunklength, num_coord);
    idx_len = length(idx);
     
    % prepare ratio Earth radius over radius
    TF = bsxfun(@power, a_E ./ r(idx), tl); % before tl'
    TA = []; TA(idx_len, row) = 0; % faster than TA = zeros(length(idx), row);
    TB = []; TB(idx_len, row) = 0;
    
    % unwrapping to avoid if ... else for order m = 0
    Cnm = field(:, row);                     % get Cnm coefficients for order 0
    
    TFP = TF .* plm_func(l, 0, theRAD(idx)); % multiply TF with fully normalized Legendre Polynoms
    TA(:, 1) = TFP * Cnm;                    % for m = 0: Snm = 0, hence TB also = 0
    for m = 1:(row - 1) % rest of order m
        Cnm = field(:, row + m);             % get Cnm coefficients for order m
        Snm = field(:, row - m);             % get Snm coefficients for order m
        TFP = TF .* plm_func(l, m, theRAD(idx)); % multiply TF with fully normalized Legendre Polynoms
        TA(:, m + 1) = TFP * Cnm;
        TB(:, m + 1) = TFP * Snm;
    end
    
    % now do the final summation
    mlam = lamRAD(idx) * l_trans;         % matrix of size(length(lam), length(m))
    f(didx(idx)) = sum(TA .* cos(mlam) + TB .* sin(mlam), 2);
   
    if params.waitbar; WAIT = twaitbar(ridx / parts, WAIT); end;
end

if params.waitbar; twaitbar('close', WAIT); end;




