clear; close all
lambda = 1;
load em_constants.mat
mu0 = mu_0;
ep0 = epsilon_0;
omega = 2*pi*c/lambda;
K = 2*pi/lambda;
w = lambda/2;
L = 2*w;
len = 3;
M = 300; % Sections on the strip
deltax = L / M;
gamma = exp(vpa(eulergamma));
e = exp(1);
euler = vpa(eulergamma);
num_theta = 360*2;
theta = linspace(0,2*pi,num_theta);
% x = linspace(0 + deltax , L - deltax, M);
x = linspace(0  , L , M);
% x = linspace(-0.5*L , 0.5*L , M);
Z = zeros(1,M);
%%



for j = 1 : 1 % Only one iteration is needed as impedance matrix can be constructed through its Toeplitz property
    for k = 1 : M
        
        fun = @(xp) besselh(0,2,K*abs(x(j) - xp));% Make a symbolic Function of x_prime
        % This calculation is based on the reference cited in the code
        % introduction.
        % Mathematical Modeling of Eq. 4.63
        xp_upper = x(k) + deltax; % Upper limit of integration
        xp_lower = x(k); % Lower limit of integration
        
        %% Numerical Integration Using Gauss_Quadratures
        if abs( j - k) >= 1
            int_part =  integral(fun,xp_lower,xp_upper,'RelTol',1e-12,'AbsTol',1e-12);
%             
            Z(j,k) =  int_part; % Pocklington Integral for j not equal to k
            
        elseif j == k
            Z(j,k) = deltax*(1 - 1i*2/pi*log(gamma*K*deltax/(4*e)));
        end
    end
end

Z = toeplitz(real(Z(1,:))) + 1i*toeplitz(imag(Z(1,:))); % Make a Toeplitz matrix out of a row vector
% Z = toeplitz(Z);
%
% cos_th = cos(theta);
% cos_th = pi/2;
% sin_th = sin(theta);
% % V = 4\(omega*mu0)*exp(1i*K*x'*cos(phi));
% V = exp(1i*K*x'*cos_th);
phi = 1*pi/2;
% V = 4\(omega*mu0)*exp(1i*K*x'*cos(phi));
V = exp(-1i*K*x'*cos(phi));
% V = zeros(M,1); % Initialize Source (RHS)
% V(floor((M-2)/2)+1) = -1i*4*pi*omega*eps0*(1.0/deltax);   % Delta Source

I = Z\V;

for i = 1 : num_theta
    cosTheta = cos(theta(i));
    sinTheta = sin(theta(i));
    E_rad(i) = 0;
    % numerical integration
    for m = 1:M
        x_m = x(m) + deltax;
        E_rad(i) = E_rad(i) + deltax*I(m)*exp(1i*K*x_m*cosTheta); % Summation representation of the integral 
    end
    E_rad(i) = mu0/(4.0*pi)*E_rad(i);
end


E_db = mag2db(abs(E_rad)./max(abs(E_rad)));
% E_db = mag2db(abs(E_rad));
% First Polar Plot
figure(1)
h3 = plot(theta, E_db);
ax = gca;
h3.Color = 'black';
h3.LineWidth = 1.4;
title(['RCS of a PEC Plate of length \t',int2str(len), '$\lambda$'],'Interpreter','latex')
set(gcf,'Color','white'); % Set background color to white
set (gca,'FontName','times new roman') % Set axes fonts to Times New Roman
% xlim([ 0 pi]);
% ax.XTick = [0 pi/4 pi/2 3*pi/4 pi];
xlabel('$\phi \mathrm{[rad]}$','interpreter','latex')
ylabel('$\sigma_{2D} \mathrm{[dB]}$','interpreter','latex')
% ax.XTickLabel = { '0', '\pi/4','\pi/2','3\pi/4','\pi'};
xlim([ 0 pi]);
ax.XTick = [0 pi/4 pi/2 3*pi/4 pi];
ax.XTickLabel = { '0', '\pi/4','\pi/2','3\pi/4','\pi'};
xlabel('$\phi \mathrm{(rad)}$','interpreter','latex')
ylabel('$\sigma_{\theta} \mathrm{(dB)}$','interpreter','latex')
grid on
cleanfigure();
% matlab2tikz('filename',sprintf('ECEN637_HW10_Polar_plot_h_lambda.tex'));

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
