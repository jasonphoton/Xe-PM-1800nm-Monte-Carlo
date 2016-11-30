% function that returns an enevelope and a related time grid 
% ... E A I are also given out
% n - pulse duration in number of cycles gauss fwhm, sin2 something

function [Env E I A ALPHA BETA v r t_grid E_fh] = fct_get_EnvEIAAlphBetvr_sin2_nN(dt_au,wvlm,IWcm,CEP,n,N)

% derives from inputs laser
c0    = 299792458./(2.1876912633e6);                        % [a.u.]
wvl   = wvlm/(5.2917720859e-11);                            % [au]
omega = (2*pi*c0)/wvl;                                      % [au], angular frequency
period= 2*pi/omega;                                         % [atomic units], one period
dt    = dt_au;
E0    = sqrt(IWcm/3.509e16);

%envelope
%N               = 1; % range is not needed
t_grid          = -(n)*period:dt:(N+n)*period;
t_grid_interp   = t_grid + dt/2;
    
E0_fh = @(t) (sin(omega/(4*n)*t-pi/2).^2.*and((t >= -n*period),(t < 0))) ...
        +(1.*and((t >= 0),(t < (N)*period))) ...
        +(sin(omega/(4*n)*((N+n)*period-t)).^2.*and((t >= (N)*period),(t < (N+n)*period)));

% function handles
E_fh    = @(t)   E0_fh(t).*E0.*cos(omega.*t + CEP); 
Eabs_fh = @(t)   sqrt(E_fh(t).^2);
I_fh    = @(t)   ((E0_fh(t).*E0).^2);

% discrete fields
Env  = E0_fh(t_grid).*E0;
E    = E_fh(t_grid);
I    = I_fh(t_grid);
Eabs = Eabs_fh(t_grid);

% vector potential
A = -cumsum(E_fh(t_grid_interp))*dt;
A_interp = [0 A];
 for i = 1:length(A)
      A_interp(i) = (A_interp(i)+A_interp(i+1))/2;
 end
A_interp = A_interp(1:1:end-1);
A        = [0 A(1:1:end-1)];

% integration vector potential from -infinity to t
ALPHA = cumsum(A_interp(:))*dt;
ALPHA = [0 ; ALPHA(1:1:end-1)];

%integration (vector potential)² from -infinity to t
BETA = cumsum(A_interp.^2)*dt;
BETA = [0 BETA(1:1:end-1)];

%trajectorie of electron released at time tr
v = @(tr_index) A - A(tr_index);
r = @(tr_index) ALPHA - ALPHA(tr_index)-A(tr_index).*(t_grid'-t_grid(tr_index));
