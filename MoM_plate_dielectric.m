clear; close all
lambda = 1;
load em_constants.mat
mu0 = mu_0;
ep0 = epsilon_0;
ep1 = 1;
ep2 = (4);
eta0 = sqrt(mu0/ep0);
omega = 2*pi*c/lambda;
K0 = omega*sqrt(mu0*ep0*ep1);
K1 = omega*sqrt(mu0*ep0*ep2);
w = lambda/2;
L = 7.52*w;
len = 1/2;
M = 400; % Sections on the strip
deltax = L / M;
gamma = exp(vpa(eulergamma));
e = exp(1);
euler = vpa(eulergamma);
% x = linspace(0 + deltax , L - deltax, M);
x = linspace(-0.5*L , 0.5*L , M);
Z = zeros(1,M);
Y = zeros(1,M);
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
        if abs( j - k) >= 1
            int_part1 =  integral(fun1,xp_lower,xp_upper,'RelTol',1e-12,'AbsTol',1e-12);
            int_part2 =  integral(fun2,xp_lower,xp_upper,'RelTol',1e-12,'AbsTol',1e-12);
            %
            Z(k) =  int_part1; % Pocklington Integral for j not equal to k
            Y(k) =  int_part2; % Pocklington Integral for j not equal to k
            
            
            
        elseif j == k
            Z1 = deltax*(1 - 1i*2/pi*log(gamma*K0*deltax/(4*e)));
            Z2 = deltax*(1 - 1i*2/pi*log(gamma*K1*deltax/(4*e)));
            Z(k) = Z1 + Z2;
            Y1 = deltax + K0^2*deltax^3/96 - 1i*deltax/pi*...
                (2*log(gamma*K0*deltax/(4*e)) + 16/(K0*deltax)^2);
            Y2 = deltax +  K1^2*deltax^3/96 - 1i*deltax/pi*...
                (2*log(gamma*K1*deltax/(4*e)) + 16/(K1*deltax)^2);
            Y(k) = ep1*Z1 + ep2*Z2 + ep1*Y1 + ep2*Y2;
        end
    end
end

Z = (toeplitz(real(Z(:))) + 1i*toeplitz(imag(Z(:)))); % Make a Toeplitz matrix out of a row vector
% Z = toeplitz(real(Z(2,M-1)),imag(Z(2,M-1))); % Make a Toeplitz matrix out of a row vector
% Z = padarray(Z,[1 1]);
Y = (toeplitz(real(Y(:))) + 1i*toeplitz(imag(Y(:)))); % Make a Toeplitz matrix out of a row vector
% Y = padarray(Y,[1 1], 1e8);
% Z = toeplitz(Z);
%
phi = 1*pi/2;
% V = 4\(omega*mu0)*exp(1i*K*x'*cos(phi));
Ve = exp(1i*K0*x'*cos(phi));
Vm = -sin(phi)*exp(1i*K0*x'*cos(phi))/eta0;

Ie = Z\Ve;
% I1(2:M-1) =Z(2:M-1,2:M-1)\V1;
Im = Y\Vm;
% I2(2:M-1) =Y(2:M-1,2:M-1)\V2;

%% PLOTTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First plot the displacement

% Electric Current Plot
figure(1)
% x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
H = plot(x, real(Ie),x, imag(Ie));
ax = gca;
H(1).Color = 'black';
H(1).LineWidth = 1.4;
H(2).Color = 'black';
H(2).LineWidth = 1.4;
H(2).LineStyle = '--';
% title(['Current on the wire of half-length $ h = .35\lambda$ at f =  ',int2str(f/1e6), ' MHz'],'Interpreter','latex')
set(gcf,'Color','white'); % Set background color to white
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% ax.XTick = [-.3498 -0.2625  -0.1750 -0.0875 0 0.0875 0.1750 0.2625 0.3498];
% ax.XTickLabel = { '-h','-.75h','-.5h','-.25h' , '0' ,'.25h', '5h', '.75h', 'h'};
% ax.YTick = [-2e-3   -1e-3 0 1e-3 2e-3];
% ax.YTickLabel = { '-.002','-.001','0' , '.001', '.002'};
% axis([ -.3498 .3498 -2.5e-3 2.5e-3]);
hold on
title(['Current Distribution on a Dielectric plate of length ',int2str(len), '$\lambda$ at $\phi_i = \pi/2$'],'Interpreter','latex')
xlabel('$\frac{x}{\lambda}$','interpreter','latex')
ylabel('$J_z (\mathrm{\frac{A}{m}})$','interpreter','latex')
legend('Real Part', 'Imaginary Part');
grid on
cleanfigure();
% matlab2tikz('filename',sprintf('dielectric_e_current_RI.tex'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Electric Current Plot Absolute Value
figure(2)
% x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
H = plot(x, abs(Ie));
ax = gca;
H(1).Color = 'black';
H(1).LineWidth = 1.4;
% title(['Current on the wire of half-length $ h = .35\lambda$ at f =  ',int2str(f/1e6), ' MHz'],'Interpreter','latex')
set(gcf,'Color','white'); % Set background color to white
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% ax.XTick = [-.3498 -0.2625  -0.1750 -0.0875 0 0.0875 0.1750 0.2625 0.3498];
% ax.XTickLabel = { '-h','-.75h','-.5h','-.25h' , '0' ,'.25h', '5h', '.75h', 'h'};
% ax.YTick = [-2e-3   -1e-3 0 1e-3 2e-3];
% ax.YTickLabel = { '-.002','-.001','0' , '.001', '.002'};
% axis([ -.3498 .3498 -2.5e-3 2.5e-3]);
hold on
title(['Current Distribution on a Dielectric plate of length ',int2str(len), '$\lambda$ at $\phi_i = \pi/2$'],'Interpreter','latex')
xlabel('$\frac{x}{\lambda}$','interpreter','latex')
ylabel('$J_z (\mathrm{\frac{A}{m}})$','interpreter','latex')
grid on
cleanfigure();
% matlab2tikz('filename',sprintf('dielectric_e_current_abs.tex'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic Current Plot
figure(3)
% x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
H = plot(x, real(Im),x, imag(Im));
ax = gca;
H(1).Color = 'black';
H(1).LineWidth = 1.4;
H(2).Color = 'black';
H(2).LineWidth = 1.4;
H(2).LineStyle = '--';
% title(['Current on the wire of half-length $ h = .35\lambda$ at f =  ',int2str(f/1e6), ' MHz'],'Interpreter','latex')
set(gcf,'Color','white'); % Set background color to white
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% ax.XTick = [-.3498 -0.2625  -0.1750 -0.0875 0 0.0875 0.1750 0.2625 0.3498];
% ax.XTickLabel = { '-h','-.75h','-.5h','-.25h' , '0' ,'.25h', '5h', '.75h', 'h'};
% ax.YTick = [-2e-3   -1e-3 0 1e-3 2e-3];
% ax.YTickLabel = { '-.002','-.001','0' , '.001', '.002'};
% axis([ -.3498 .3498 -2.5e-3 2.5e-3]);
hold on
title(['Magnetic Current Distribution on a Dielectric plate of length ',int2str(len), '$\lambda$ at $\phi_i = \pi/2$'],'Interpreter','latex')
xlabel('$\frac{x}{\lambda}$','interpreter','latex')
ylabel('$M_x (\mathrm{\frac{V}{m}})$','interpreter','latex')
legend('Real Part', 'Imaginary Part');
grid on
cleanfigure();
% matlab2tikz('filename',sprintf('dielectric_m_current_RI.tex'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Magnetic Current Plot Absolute Value
figure(4)
% x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
H = plot(x, abs(Im));
ax = gca;
H(1).Color = 'black';
H(1).LineWidth = 1.4;
% title(['Current on the wire of half-length $ h = .35\lambda$ at f =  ',int2str(f/1e6), ' MHz'],'Interpreter','latex')
set(gcf,'Color','white'); % Set background color to white
set(gca,'FontName','times new roman','FontSize',11) % Set axes fonts to Times New Roman
% ax.XTick = [-.3498 -0.2625  -0.1750 -0.0875 0 0.0875 0.1750 0.2625 0.3498];
% ax.XTickLabel = { '-h','-.75h','-.5h','-.25h' , '0' ,'.25h', '5h', '.75h', 'h'};
% ax.YTick = [-2e-3   -1e-3 0 1e-3 2e-3];
% ax.YTickLabel = { '-.002','-.001','0' , '.001', '.002'};
% axis([ -.3498 .3498 -2.5e-3 2.5e-3]);
hold on
title(['Magnetic Current Distribution on a Dielectric plate of length ',int2str(len), '$\lambda$ at $\phi_i = \pi/2$'],'Interpreter','latex')
xlabel('$\frac{x}{\lambda}$','interpreter','latex')
ylabel('$M_x (\mathrm{\frac{V}{m}})$','interpreter','latex')
grid on
cleanfigure();
% matlab2tikz('filename',sprintf('dielectric_e_current_abs.tex'));