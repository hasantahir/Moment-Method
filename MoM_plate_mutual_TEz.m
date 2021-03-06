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
L1 = 6*w;
L2 = 6*w;
len1 = 3;
len2 = 1;
M1 = 400; % Sections on the strip
M2 = 400;
deltax1 = L1 / M1;
deltax2 = L2 / M1;
gamma = exp(vpa(eulergamma));
e = exp(1);
euler = vpa(eulergamma);
% x = linspace(0 + deltax , L - deltax, M);
x1 = linspace(-.5*L1 , .5*L1 , M1);
x2 = linspace(-.5*L2, .5*L2, M2);

% y-spacing
y1 = 1*lambda/2;
Z = zeros(1,M1 + M2);
%%



for j = 1 : 1 % Only one iteration is needed as impedance matrix can be constructed through its Toeplitz property
    for k = 1 : M1
        
        % Make a symbolic Function of x_prime
        % For a dielectric sheet
        fun1 = @(xp) (besselh(0,2,K0*abs(x1(j) - xp)) + ...
            besselh(0,2,K1*abs(x1(j) - xp)) + ...
            besselh(2,2,K0*abs(x1(j) - xp)) + ...
            besselh(2,2,K1*abs(x1(j) - xp)));

%         fun1 = @(xp1) besselh(0,2,K*abs(x1(j) - xp1)) + ...
%             besselh(2,2,K*abs(x1(j) - xp1));
        % This calculation is based on the reference cited in the code
        % introduction.
        % Mathematical Modeling of Eq. 4.63
        xp_upper = x1(k) + deltax1; % Upper limit of integration
        xp_lower = x1(k); % Lower limit of integration
        
        %% Numerical Integration Using Gauss_Quadratures
        if abs( j - k) >= 1
            
            int_part1 =  integral(fun1,xp_lower,xp_upper,'RelTol',1e-12,'AbsTol',1e-12);
            %
            Z(k) =  int_part1; % Pocklington Integral for j not equal to k
            
        elseif j == k
            YZ1 = deltax1*(1 - 1i/pi*(3 + 2 *log(1.781*K0*deltax1/(4*e)) + 16/(K0*deltax1)^2));
            YZ2 = deltax1*(1 - 1i/pi*(3 + 2 *log(1.781*K1*deltax1/(4*e)) + 16/(K1*deltax1)^2));
            Z(k) = 1*(YZ1 + YZ2);
        end
    end
    
    for k = 1 : M2
        
        fun1 = @(xp2) besselh(0,2,K0*sqrt((x1(j) - xp2).^2 + y1^2)) + ...
            besselh(2,2,K0*sqrt((x1(j) - xp2).^2 + y1^2));
        % Mathematical Modeling of Eq. 4.63
        xp_upper = x2(k) + deltax2; % Upper limit of integration
        xp_lower = x2(k); % Lower limit of integration
        
        %% Numerical Integration Using Gauss_Quadratures
        %         if abs( j - k) >= 1
%         if abs( j - k) >= 1
            
            int_part1 =  integral(fun1,xp_lower,xp_upper,'RelTol',1e-12,'AbsTol',1e-12);
            %
            Z(k + M1) =  int_part1; % Pocklington Integral for j not equal to k
            
       if j == k && y1 <= lambda / 4
            Z(j,k + M2) = deltax1*(1 - 1i/pi*(3 + 2 *log(1.781*K1*deltax1/(4*e)) + 16/(K1*deltax1)^2));
        end
        
    end
end
% Ignore omega/4 factor
Z_self = Z(1:M1);
Z_mutual  = Z(M1+1 : M1 + M2);

Z = toeplitz(real(Z(:))) + 1i*toeplitz(imag(Z(:))); % Make a Toeplitz matrix out of a row vector
Z_self = toeplitz(real(Z_self(:))) + 1i*toeplitz(imag(Z_self(:))); % Make a Toeplitz matrix out of a row vector
Z_mutual = toeplitz(real(Z_mutual(1,:))) + 1i*toeplitz(imag(Z_mutual(1,:))); % Make a Toeplitz matrix out of a row vector
% Z = toeplitz(Z);
%
phi = 1*pi/2;
X = [x1,zeros(1,length(x2))];
%     X = [x1,x2];
V1 = exp(1i*K0*x1'*cos(phi));
V2 = exp(1i*K0*sqrt(x2.^2 + y1^2)'*cos(phi));
%     V = eta0*sin(phi)*[V1',V2'];
V = exp(1i*K0*X'*cos(phi));
V = eta0*sin(phi)*[sin(phi)*ones(M1,1).*exp(1i*K0*x1'*cos(phi)); zeros(M2,1)];
% V = zeros(M1 + M2 -2,1); % Initialize Source (RHS)
% V = eta0*sin(phi)*exp(1i*K*X'*cos(phi));

% V_self = eta0*sin(phi)*exp(1i*K*x1'*cos(phi));
% V_mutual = eta0*sin(phi)*exp(1i*K*(x2)'*cos(phi));
% V_mutual(1) = 0;
% V_mutual(M2) = 0;

I = Z\V;
% I_self = Z_self\ V_self;
% I_mutual = Z_mutual\ V_mutual;

     I_self = I(1:M1);
     I_mutual = I(M1+1:M1+M2);

% Current Plot
figure(1)
% x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
H = plot(x2, real(I_mutual),x2, imag(I_mutual));
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
title(['Current Distribution Induced on a PEC plate of length ',int2str(len2), '$\lambda$ at $\phi_i = \pi/2$'],'Interpreter','latex')
xlabel('$\frac{x}{\lambda}$','interpreter','latex')
ylabel('$J_z \mathrm{A}$','interpreter','latex')
legend('Real Part', 'Imaginary Part');
grid on


% Current Plot
figure(2)
% x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
H = plot(x2, abs(I_mutual));
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
title(['Current Distribution Induced on a PEC plate of length ',int2str(len2), '$\lambda$ at $\phi_i = \pi/2$'],'Interpreter','latex')
xlabel('$\frac{x}{\lambda}$','interpreter','latex')
ylabel('$J_z \mathrm{A}$','interpreter','latex')
grid on



% Current Plot
figure(3)
% x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
H = plot(x1, real(I_self),x1, imag(I_self));
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
title(['Self Current Distribution on a Dielectric plate of length ',int2str(len1), '$\lambda$ at $\phi_i = \pi/2$'],'Interpreter','latex')
xlabel('$\frac{x}{\lambda}$','interpreter','latex')
ylabel('$J_z \mathrm{A}$','interpreter','latex')
legend('Real Part', 'Imaginary Part');
grid on


% Current Plot
figure(4)
% x = linspace(-.7*lambda/2, .7*lambda/2, M) ;
H = plot(x1, abs(I_self));
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
title(['Self Current Distribution on a Dielectric plate of length ',int2str(len1), '$\lambda$ at $\phi_i = \pi/2$'],'Interpreter','latex')
xlabel('$\frac{x}{\lambda}$','interpreter','latex')
ylabel('$J_z \mathrm{A}$','interpreter','latex')
grid on
