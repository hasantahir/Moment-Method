clear; close all
lambda = 1;
load em_constants.mat
mu0 = mu_0;
ep0 = epsilon_0;
ep1 = 1;
ep2 = (4);
eta0 = sqrt(mu0/ep0);
omega = 2*pi*c/lambda;
K0 = 2*pi/lambda;
K1 = 2*pi/lambda*sqrt(ep1);
w = lambda/2;
L = 6*w;
len = 4/2;
M = 400; % Sections on the strip
deltax = L / M;
gamma = exp(vpa(eulergamma));
e = exp(1);
euler = vpa(eulergamma);
num_theta = 361;
theta = linspace(0,pi,num_theta);
% x = linspace(0 + deltax , L - deltax, M);
x = linspace(0  , L , M);
% x = linspace(-0.5*L , 0.5*L , M);
Z = zeros(1,M);
Y = zeros(1,M);
%%



for j = 1 : 1 % Only one iteration is needed as impedance matrix can be constructed through its Toeplitz property
    for k = 1 : M
        
        fun1 = @(xp) besselh(0,2,K0*abs(x(j) - xp)) + ...
            besselh(0,2,K1*abs(x(j) - xp));% Make a symbolic Function of x_prime
        fun2 = @(xp) (ep1*besselh(0,2,K0*abs(x(j) - xp)) + ...
            ep2*besselh(0,2,K1*abs(x(j) - xp)) + ...
            ep1*besselh(2,2,K0*abs(x(j) - xp)) + ...
            ep2*besselh(2,2,K1*abs(x(j) - xp)));
        % This calculation is based on the reference cited in the code
        % introduction.
        % Mathematical Modeling of Eq. 4.63
        xp_upper = x(k) + deltax; % Upper limit of integration
        xp_lower = x(k); % Lower limit of integration
        
        %% Numerical Integration Using Gauss_Quadratures
        % Self Terms
        if j == k
            Z(k) = deltax*(1 - 1i*2/pi*log(gamma*K0*deltax/(4*e))) + ...
                deltax*(1 - 1i*2/pi*log(gamma*K1*deltax/(4*e)));
            eZ = ep1*deltax*(1 - 1i*2/pi*log(gamma*K0*deltax/(4*e))) + ...
                ep2*deltax*(1 - 1i*2/pi*log(gamma*K1*deltax/(4*e)));
            Y1 = deltax + K0^2*deltax^3/96 - 1i*deltax/pi*...
                (2*log(gamma*K0*deltax/(4*e)) + 16/(K0*deltax)^2);
            Y2 = deltax +  K1^2*deltax^3/96 - 1i*deltax/pi*...
                (2*log(gamma*K1*deltax/(4*e)) + 16/(K1*deltax)^2);
            Y(k) = (eZ + ep1*Y1 + ep2*Y2);
            
        else
            
            int_part1 =  integral(fun1,xp_lower,xp_upper,'RelTol',1e-12,'AbsTol',1e-12);
            int_part2 =  integral(fun2,xp_lower,xp_upper,'RelTol',1e-12,'AbsTol',1e-12);
            %
            Z(j,k) =  int_part1; % Pocklington Integral for j not equal to k
            Y(j,k) =  int_part2; % Pocklington Integral for j not equal to k
            
        end
    end
end
Z = (toeplitz(real(Z(:))) + 1i*toeplitz(imag(Z(:)))); % Make a Toeplitz matrix out of a row vector
Y = .5*(toeplitz(real(Y(:))) + 1i*toeplitz(imag(Y(:)))); % Make a Toeplitz matrix out of a row vector
% V = exp(1i*K*x'*cos_th);

% Angle of Incidence
phi = 1*pi/2;

% V = 4\(omega*mu0)*exp(1i*K*x'*cos(phi));
% V = zeros(M,1); % Initialize Source (RHS)
% V(floor((M-2)/2)+1) = -1i*4*pi*omega*eps0*(1.0/deltax);   % Delta Source

% Define Plane wave
% V = 4\(omega*mu0)*exp(1i*K*x'*cos(phi));
Ve = exp(1i*K0*x'*cos(phi));
Vm = -sin(phi)*exp(1i*K0*x'*cos(phi))/eta0;

% Find current for Pulse Basis Functions and Delta Point-Matching
Ie = Z\Ve;
Im = Y\Vm;


% Find the Scattered Field
for i = 1 : num_theta
    cosTheta = cos(theta(i));
    sinTheta = sin(theta(i));
    E_rad(i) = 0;
%     H_rad(i) = 0;
    % numerical integration
    for m = 1:M
        x_m = x(m) + deltax;
        E_rad(i) = E_rad(i) + deltax*(-Ie(m)*eta0 - Im(m)*sinTheta)*exp(1i*K0*x_m*cosTheta); % Summation representation of the integral
%         H_rad(i) = H_rad(i) + deltax*Im(m)*exp(1i*K0*x_m*cosTheta); % Summation representation of the integral
    end
%     E_rad(i) = mu0/(4.0*pi)*E_rad(i);
%     H_rad(i) = ep0/(4.0*pi)*H_rad(i);
%     E_rad(i) = -sinTheta*E_rad(i) - eta0*H_rad(i);
end

% Normalize the field
E_db = (abs(E_rad)./max(abs(E_rad)));
% H_db = mag2db(abs(H_rad)./max(abs(H_rad)));
E_db = mag2db(abs(E_rad));
% First Polar Plot - E-field 
figure(1)
h3 = plot(theta, E_db);
ax = gca;
h3.Color = 'black';
h3.LineWidth = 1.4;
title(['$\theta$ Pattern of a Dielctric Plate of length \t',int2str(len), '$\lambda$'],'Interpreter','latex')
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
% xlim([ pi/2 pi]);
ax.XTick = [0 pi/4 pi/2 3*pi/4 pi];
ax.XTickLabel = { '0', '\pi/4','\pi/2','3\pi/4','\pi'};
xlabel('$\phi \mathrm{(rad)}$','interpreter','latex')
ylabel('$\sigma_{\theta} \mathrm{(dB)}$','interpreter','latex')
grid on
cleanfigure();
% matlab-2tikz('filename',sprintf('dielectric_e_pat.tex'));

%
figure(2) % Polar Plot
%
h1 = polar(theta,abs(E_rad)./max(abs(E_rad)));
h1.Color = 'black';
h1.LineWidth = 1.4;
% title(['Radiation Pattern of a PEC Box at f =  ',int2str(finc/1e6), ' MHz'],'Interpreter','latex')
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
% cleanfigure();
% matlab2tikz('filename',sprintf('ECEN637_HW9_Polar_plot.tex'));
%
%
% % First Polar Plot - E-field 
% figure(2)
% h3 = plot(theta, H_db);
% ax = gca;
% h3.Color = 'black';
% h3.LineWidth = 1.4;
% title(['$\phi$ Pattern of a Dielctric Plate of length \t',int2str(len), '$\lambda$'],'Interpreter','latex')
% set(gcf,'Color','white'); % Set background color to white
% set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
% xlim([ 0 pi]);
% ax.XTick = [0 pi/4 pi/2 3*pi/4 pi];
% ax.XTickLabel = { '0', '\pi/4','\pi/2','3\pi/4','\pi'};
% xlabel('$\phi \mathrm{(rad)}$','interpreter','latex')
% ylabel('$\sigma_{\phi} \mathrm{(dB)}$','interpreter','latex')
% grid on
% cleanfigure();
% matlab2tikz('filename',sprintf('dielectric_h_pat.tex'));

