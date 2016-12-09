% MoM code for thin wires, using collocation and
% piecewise sinusoisdal basis functions, to analyze a six-element Yagi antenna
% Developed from the original by David B Davidson in Aug 2011. This 
% application was suggested by Matthys M Botha.
% Geometry is hard-coded for the Yagi in [1,Section 5.3]
% Note: magnetic frill not currently implemented.

% Reference
% [1] D B Davidson, "Computational Electromagnetics for RF and Microwave
% Engineering", 2nd Edition, Cambridge Univ Press

clear;
clf;
c0= 2.997925e8; % m/s
N_seg=input('Enter number of segments per element (must be even): ');
if isempty(N_seg)
    N_seg = 20 % sensible default
end
N_freq = 51;
f_start = 275e6; %Start frequency in Hz
f_end   = 325e6; %Stop frequency in Hz
if N_freq > 1 % Avoid div by 0 error
    del_f   = (f_end - f_start)/(N_freq-1);
else
    del_f = 0 % Arbitrary value
end
L(1) = 0.482; % Reflector, all dimensions in wavelengths
L(2) = 0.475; % Driven element
L(3) = 0.428; % Director 1
L(4) = 0.420; % Director 2
L(5) = 0.420; % Director 3
L(6) = 0.428; % Director 4
x(1) = -0.2;  % Reflector position
x(2) = 0;     % Driven element assumed at x=0
x(3)  = 1*0.25;%
x(4)  = 2*0.25;%
x(5)  = 3*0.25;%
x(6)  = 4*0.25;%
% All wire radii the same
a     = 0.00425;
lambda_c = 1;
L = L*lambda_c; % Scale by wavelength
a = a*lambda_c;
b = 3.5*a % for a 75 ohm system.

dz = L/N_seg;

for l_freq=1:N_freq
    freq(l_freq) = f_start+del_f*(l_freq-1);
    lambda = c0/freq(l_freq);
    B = 2*pi/lambda; % This is beta.
    for m_elem=1:6 % loop over number of Yagi elements (tests)
        for n_elem=1:6 % loop over number of Yagi elements (sources)
            Bdz=B*dz(n_elem);
            N=N_seg-1; % Number of degrees of freedom on each wire is one less than segments.
            for m=1:N
                z_m = -L(m_elem)/2+m*dz(m_elem);
                x_m = x(m_elem);
                for n=1:N
                    z_np1 = -L(n_elem)/2+(n+1)*dz(n_elem);
                    z_n   = -L(n_elem)/2+n*dz(n_elem);
                    z_nm1 = -L(n_elem)/2+(n-1)*dz(n_elem);
                    x_n = x(n_elem);
                    
                    Rnp1 = sqrt((x_m-x_n+a)^2+(z_np1-z_m)^2); % Rn+1. Uses [1,Eq (4.30)] for z_n and z_m, with extensions for elements
                    Rn   = sqrt((x_m-x_n+a)^2+(z_n-z_m)^2);
                    Rnm1 = sqrt((x_m-x_n+a)^2+(z_nm1-z_m)^2);
                    A = exp(-1i*B*Rnp1)/Rnp1;
                    BB = exp (-1i*B*Rn)*sin(2*Bdz)/(Rn*sin(Bdz));
                    C = exp(-1i*B*Rnm1)/Rnm1;
                    Z((m_elem-1)*N+m,(n_elem-1)*N+n) = -1i*30*(A-BB+C)/sin(Bdz);
                end
            end
        end
    end
        
    % Delta gap feed model
    CentreSeg = ceil(N/2)+N; % located in middle of 2nd wire
    deltaZ=dz(2);
    for n = (1:6*N);
        V(n) = 0;
    end;
    V(CentreSeg) = -1/deltaZ;
    I_delta_gap = Z \ V';
    Zin_delta_gap(l_freq)=1/I_delta_gap(CentreSeg); %  this is the dipole impedance
end

% compute reflection coefficients in 50 Ohm system (magnitude, dB)
gamma_delta = 20*log10(abs((Zin_delta_gap-50)./(Zin_delta_gap+50)));

% yagi_data % load data computed using NEC2 and FEKO, as in [1,Fig 5.7]
plot(freq/1e6,gamma_delta,'--')
grid on
axis([275 325 -20 0])
xlabel('Freq [MHz]')
ylabel('S_{11} [dB]')
title('Comparison of MoM simulations of 6 element Yagi, \Delta \approx \lambda_0/40 in all cases')
%legend('FEKO','NEC2',['This code ',num2str(N_seg),' segments per wire'],0)
% legend('FEKO','NEC2','This code',0)
% print -deps yagi_cmp_s11
% print -dpng yagi_cmp_s11