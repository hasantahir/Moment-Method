function [rcs_te, rcs_tm,  j_te, j_tm] = CylinderRCSByMode(radius, freq, nphi, nang, nmode)

% Compute the circumferential induced electric surface current density and
% the radar cross section of a perfectly conducting, infinitely long
% cylinder. Follows the treatment in Chapter 4 of 
%
% Ruck, et. al. "Radar Cross Section Handbook", Plenum Press, 1970.
%  
% The incident electric fields are TMz and TEz polarized. The incident wave
% is assumed to be at an angle of phi = 0. 
% 
% Inputs:
%   radius: Radius of the cylinder (meters)
%   frequency: Operating frequency (Hz)
%   nphi: Number of scattering angles (between 0 and 2*pi).
%   nang: Number of circumferential angles (between 0 and 2*pi) to compute
%       the surface currents.
%   nmode: Maximum mode number for computing Bessel functions
%
% Outputs:
%   rcs_te: Forward-scatter TEz-polarized RCS between 0 and 2*pi
%   rcs_tm: Forward-scatter TMz-polarized RCS between 0 and 2*pi
%   j_te: Azimuthal (TEz-excited) surface currents between 0 and 2*pi
%   j_tm: Axial (TMz-excited) surface currents between 0 and 2*pi
% 
%   Author: Walton C. Gibson, Tripoint Industries, Inc.
%   email: kalla@tripoint.org

c = 299792458;
mu = 4.0e-7*pi;
epsilon = 1.0/c/mu/c;
w = 2.0*pi*freq;
k = w*sqrt(mu*epsilon);
eta = sqrt(mu/epsilon);

ang = linspace(0.0, 2.0*pi, nang);
 
% Compute the circumferential surface currents Kz and Kphi via Ruck, et. al.
% (4.1-15, TMz polarization) and (4.1-6, TEz polarization)

j_te = zeros(1,nang);
j_tm = zeros(1,nang);

for n = 0:nmode
    
    an = cos(n*ang)/besselh(n,1,k*radius);

    if n == 0
        val = ((-i)^n)*cos(n*ang);
    else
        val = 2*((-i)^n)*cos(n*ang);
    end

        % Hankel Function
    h = besselh(n,1,k*radius);
        % derivative of Hankel Function
    hd = n*besselh(n,1,k*radius)/(k*radius) - besselh(n+1,1,k*radius);

    % Ruck, et. al. (4.1-16, TEz polarization)
    j_te = j_te + val/hd;
    % Ruck, et. al. (4.1-15, TMz polarization)
    j_tm = j_tm + val/h;
    
end

j_te = -i*j_te*2/(pi*k*radius*eta);
j_tm = j_tm*2/(pi*k*radius*eta);

% Compute the scattered far electric field via Ruck, et. al. via
% (4.1-5, TMz polarization) and (4.1-6, TEz polarization)

phi = linspace(0.0, 2.0*pi, nphi);
es_te = zeros(1, nphi);
es_tm = zeros(1, nphi);

for n = 0:nmode
    
    if n == 0
        val = ((-1)^n)*cos(n*phi);
    else
        val = 2*((-1)^n)*cos(n*phi);
    end
    
        % regular Bessel and Hankel functions (Ruck, et. al. 4.1-13)
    an = -besselj(n,k*radius)/besselh(n,1,k*radius);
       
        % derivatives of regular Bessel and Hankel functions
    jd = n*besselj(n,k*radius)/(k*radius) - besselj(n+1,k*radius);
    hd = n*besselh(n,1,k*radius)/(k*radius) - besselh(n+1,1,k*radius);
        % (Ruck, et. al. 4.1-14)
    bn = -jd/hd;
    
        % (Ruck, et. al. 4.1-5)
    es_tm = es_tm + val*an;
    
        % (Ruck, et. al. 4.1-6)
    es_te = es_te + val*bn;
     
end


es_tm = es_tm*sqrt(2/pi)/sqrt(k);
es_te = es_te*sqrt(2/pi)/sqrt(k);

% we do not include the term exp(i*(k*r - pi/4))/sqrt(r) as it does not
% factor into the RCS
rcs_tm = 20.0*log10(sqrt(2.0*pi)*abs(es_tm));
rcs_te = 20.0*log10(sqrt(2.0*pi)*abs(es_te));