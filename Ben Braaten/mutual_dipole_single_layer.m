%**************************************************************************
% Author:       Ben Braaten
% Date:         11/23/08
% Description:  This code calculates the mutual impedance between two
%               printed dipoles on a single anisotropic dielectric layer.
%               The dipoles can have length L and width W.
%**************************************************************************
%                           FFFFF   IIIII   N     N      A       L
%                           F         I     N N   N     A A      L
%                           FFF       I     N  N  N    AAAAA     L
%                           F         I     N   N N   A     A    L
%                           F       IIIII   N     N  A       A   LLLLL

% DDDD    IIIII      SSS     SSS    EEEE    RRRR    TTTTT       A       TTTTT  IIIII    OOOO    N     N
% D   D     I       S       S       E       R   R     T        A A        T      I     O    O   N N   N
% D   D     I         S       S     EEE     RRRR      T       AAAAA       T      I     O    O   N  N  N
% D   D     I           S       S   E       R  R      T      A     A      T      I     O    O   N   N N
% DDDD    IIIII      SSS     SSS    EEEE    R   R     T     A       A     T    IIIII    OOOO    N     N

%                                CCC    OOOO    DDDD    EEEE
%                               C      O    O   D   D   E
%                               C      O    O   D   D   EEE
%                               C      O    O   D   D   E
%                                CCC    OOOO    DDDD    EEEE




clc
clear all


s=0;                                                                        %counter reset
f=.5e9;                                                                     %frequency
w=2*pi*f;                                                                   %omega
c=3e8;                                                                      %speed of light
lambda0=c/f;                                                                %free-space wavelength
a=.4e-3;                                                                    %dipole radius
W=.5e-3;                                                                    %dipole width


e11_matrix=[3.25 5.12 3.25 5.12 5.12];                                      %aniso. permittivity matrix
e12_matrix=[3.25 3.25 5.12 3.25 3.25];                                      %aniso. permittivity matrix
d1_matrix=[1.58e-3 1.58e-3 1.58e-3 5*1.58e-3 10*1.58e-3];                   %substrate thickness matrix

for count=1:5;
    
    
    for G=0e-3:10e-3:600e-3;                                                %steps through values of G
        S=10e-3;                                                            %separation distance (Echelon)
        M=14;                                                               %number of dipole segments
        L=15e-2+2*W;                                                        %length of dipole
        d1=d1_matrix(count);                                                %substrate thickness
        e0=8.854e-12;                                                       %constant
        e11=e11_matrix(count);                                              %aniso. permittivity
        e12=e12_matrix(count);                                              %aniso. permittivity
        u0=1.2566e-6;                                                       %constant
        k=w*sqrt(u0*e0*e11)                                                 %wavenumber in the material
        k0=w*sqrt(u0*e0)                                                    %free-space wavenumber
        
        s=s+1;                                                              %counter
        N=M;                                                                %number of source points
        Vs=1;                                                               %magnitude of delta-source
        delx=L/N;                                                           %x-dim. step size
        delz=W;                                                             %z-dim. step size
        y=[d1*ones(1,M-1) d1*ones(1,M-1)];                                  %both dipole locations
        xm=-L/2+delx:delx:L/2-delx;                                         %x-dim. match points
        xm=[xm xm+G];
        xn=xm;                                                              %x-dim. source points
        zm=[zeros(1,M-1),S*ones(1,M-1)];                                    %z-dim. match points
        zn=zm;                                                              %z-dim. source points
        w=2*pi*f;                                                           %omega
        
        %root finding routine
        %The following routine uses the secant method
        %to find the poles in the spectral domain
        %immittance fuctions.
        
        iteration=100;
        root1(1)=k;
        root1(2)=k-.01;
        root2(1)=0;
        root2(2)=.9*k0*j;
        
        p=2;
        while p<iteration;
            if e11==1;
                p=iteration+1;
                root1(p)=k0;
            else
                p=p+1;
                gamma01=sqrt(root1(p-1)^2-w^2*u0*e0);
                gamma02=sqrt(root1(p-2)^2-w^2*u0*e0);
                gammae1=sqrt((e12)/(e11)*(root1(p-1)^2-w^2*u0*e0*e11));
                gammae2=sqrt((e12)/(e11)*(root1(p-2)^2-w^2*u0*e0*e11));
                gammah1=sqrt(root1(p-1)^2-w^2*u0*e0*e12);
                gammah2=sqrt(root1(p-2)^2-w^2*u0*e0*e12);
                f11=gammae1+gamma01*e12*coth(gammae1*d1);
                f12=gammae2+gamma02*e12*coth(gammae2*d1);
                root1(p)=root1(p-1)-(f11)*(root1(p-1)-root1(p-2))/((f11)-(f12));
                
                
                gamma01=sqrt(root2(p-1)^2-w^2*u0*e0);
                gamma02=sqrt(root2(p-2)^2-w^2*u0*e0);
                gammae1=sqrt((e12)/(e11)*(root1(p-1)^2-w^2*u0*e0*e11));
                gammae2=sqrt((e12)/(e11)*(root1(p-2)^2-w^2*u0*e0*e11));
                gammah1=sqrt(root1(p-1)^2-w^2*u0*e0*e12);
                gammah2=sqrt(root1(p-2)^2-w^2*u0*e0*e12);
                f11=gamma01+gammah1*coth(gammah1*d1);
                f12=gamma02+gammah2*coth(gammah2*d1);
                root2(p)=root2(p-1)-(f11)*(root2(p-1)-root2(p-2))/((f11)-(f12));
                
                if abs(root1(p)-root1(p-1))<.001&p>1
                    p=iteration+1;
                else
                    % do nothing
                end
                
            end
        end
        
        root=real(root1(length(root1)));                                     %saves root
        
        %polar points
        %The following routines setup the polar integration points.
        mult=500+(L-.5*lambda0)*(-500/lambda0)+200;
        kmax_alpha=mult*k0;
        kmax_beta=mult*k0;
        Nalpha=kmax_alpha/k0;
        Nbeta=kmax_beta/k0;
        
        if root>=k0
            del_r_fine=k0/50;
            del_r_very_fine=(k-k0)*.01;
            r_fine_lower=del_r_fine:del_r_fine:k0-del_r_fine/2;
            r_very_fine=[fliplr([root-del_r_very_fine:-del_r_very_fine:k0+del_r_very_fine+del_r_fine/2]) ...
                root+del_r_very_fine:del_r_very_fine:k-del_r_very_fine];
            r_fine=[r_fine_lower r_very_fine];
            r=r_fine;
            del_r=[del_r_fine*ones(1,length(r_fine_lower)) del_r_very_fine*ones(1,length(r_very_fine))];
        else
            del_r_fine=k0/30;
            del_r_very_fine=root/30;
            r_very_fine=[[del_r_very_fine:del_r_very_fine:root-del_r_very_fine] ...
                [root+del_r_very_fine:del_r_very_fine:k-del_r_very_fine]];
            r_fine=[r_very_fine];
            r=r_fine;
            del_r=[del_r_very_fine*ones(1,length(r_very_fine))];
        end
        
        del_theta=pi/50;
        theta=[del_theta/2:del_theta:pi/2-del_theta/2 pi/2+del_theta/2:del_theta:pi-del_theta/2 ...
            pi+del_theta/2:del_theta:3*pi/2-del_theta/2 3*pi/2+del_theta/2:del_theta:2*pi-del_theta/2];
        
        Kxx_matrix=zeros(length(r),length(theta));
        Zxx=zeros(length(r),length(theta));
        
        
        %rectangular points
        %The following routines setup the rectangular integration points.
        del_alpha_lower=kmax_alpha/Nalpha;
        del_alpha_upper=del_alpha_lower;
        del_alpha_middle=k/30;
        
        alpha_lower=-fliplr([k+del_alpha_lower/2:del_alpha_lower:kmax_alpha]);
        alpha_middle=-k+del_alpha_middle/2:del_alpha_middle:k-del_alpha_middle/2;
        alpha_upper=k+del_alpha_upper/2:del_alpha_upper:kmax_alpha;
        alpha=[alpha_lower alpha_middle alpha_upper];
        del_alpha=[del_alpha_lower*ones(1,length(alpha_lower)) del_alpha_middle*ones(1,length(alpha_middle)) ...
            del_alpha_upper*ones(1,length(alpha_upper))];
        
        del_beta=kmax_beta/Nbeta;
        
        
        %Polar integration routine
        %The following routine integrates in the polar coordinates around and
        %near the poles between k0 and k.
        for m=1:2*(M-1);
            m
            for n=1:2*(N-1);
                Kxx_matrix_fine=0;
                Kxx_matrix=0;
                for g=1:length(r);
                    
                    ks=k;
                    
                    %  2-D PWS
                    rxn_PWS=(2./sin(ks*delx)).*(1./(ks-r(g).*sin(theta))+1./(ks+r(g).*sin(theta))).*exp(-j*r(g).* ...
                        sin(theta)*xn(n)).*sin((ks*delx+r(g).*sin(theta)*delx)/2).*sin((ks*delx-r(g).*sin(theta)*delx)/2);
                    rxm_PWS=(2./sin(ks*delx)).*(1./(ks--r(g).*sin(theta))+1./(ks+-r(g).*sin(theta))).*exp(j*r(g).* ...
                        sin(theta)*xm(m)).*sin((ks*delx+-r(g).*sin(theta)*delx)/2).*sin((ks*delx--r(g).*sin(theta)*delx)/2);
                    
                    rxn_pulse=exp(j*-r(g).*cos(theta)*zn(n)).*(delz*sin(r(g).*cos(theta)*delz/2)./(delz*r(g).* ...
                        cos(theta)/2)).*exp(j*-r(g).*cos(theta)*delz/2);
                    
                    rxm_pulse=exp(j*r(g).*cos(theta)*zm(m)).*(delz*sin(r(g).*cos(theta)*delz/2)./(delz*r(g).* ...
                        cos(theta)/2)).*exp(j*r(g).*cos(theta)*delz/2);
                    
                    rxn=rxn_PWS.*rxn_pulse;
                    rxm=rxm_PWS.*rxm_pulse;
                    
                    
                    gamma0=sqrt(r(g)^2-w^2*u0*e0);
                    gammae1=sqrt((e12/e11)*(r(g)^2-w^2*u0*e0*e11));
                    gammah1=sqrt(r(g)^2-w^2*u0*e0*e12);
                    
                    R=1;
                    
                    Zxx=(cos(theta).^2).*(-j*w*u0./(gamma0+gammah1.*coth(gammah1*d1)))+ ...
                        (sin(theta).^2).*(-gammae1.*gamma0./(j*w*e0*(gammae1+gamma0*e12.*coth(gammae1*d1))));
                    
                    Kxx_matrix=Kxx_matrix+sum(Zxx.*rxm.*rxn.*R.*del_theta.*del_r(g)*r(g));
                end
                matrix_polar(m,n)=(1/(2*pi)^2)*Kxx_matrix;
            end
        end
        
        
        
        clear Kxx_matrix Zxx
        
        %rectangular integration routine
        %The following routine integrates the rectangular region outside of the
        %k circle.  This routine avoids the poles associated with +/-k in the
        %basis functions
        
        for m=1:2*(M-1);
            m
            for n=1:2*(N-1);
                Kxx_matrix_fine=0;
                Kxx_matrix=0;
                for g=1:length(alpha);
                    
                    if abs(alpha(g))<k
                        beta_jump=sqrt(k^2-alpha(g)^2)+del_beta/2;
                        beta=beta_jump:del_beta:kmax_beta;
                        beta=[beta -beta];
                    else
                        beta=del_beta/2:del_beta:kmax_beta;
                        beta=[beta -beta];
                    end
                    
                    ks=k;
                    
                    %  2-D PWS
                    
                    rxn_PWS=(2./sin(ks*delx)).*(1./(ks-alpha(g))+1./(ks+alpha(g))).*exp(-j*alpha(g)*xn(n)).* ...
                        sin((ks*delx+alpha(g)*delx)/2).*sin((ks*delx-alpha(g)*delx)/2);
                    rxm_PWS=(2./sin(ks*delx)).*(1./(ks--alpha(g))+1./(ks+-alpha(g))).*exp(j*alpha(g)*xm(m)).* ...
                        sin((ks*delx+-alpha(g)*delx)/2).*sin((ks*delx--alpha(g)*delx)/2);
                    
                    rxn_pulse=exp(j*-beta.*zn(n)).*(delz*sin(beta.*delz/2)./(delz*beta./2)).*exp(j*-beta.*delz/2);
                    
                    rxm_pulse=exp(j*beta.*zm(m)).*(delz*sin(beta.*delz/2)./(delz*beta./2)).*exp(j*beta.*delz/2);
                    
                    
                    rxn=rxn_PWS.*rxn_pulse;
                    rxm=rxm_PWS.*rxm_pulse;
                    
                    R=1;
                    
                    gamma0=sqrt(alpha(g)^2+beta.^2-w^2*u0*e0);
                    gammae1=sqrt((e12/e11)*(alpha(g)^2+beta.^2-w^2*u0*e0*e11));
                    gammah1=sqrt(alpha(g)^2+beta.^2-w^2*u0*e0*e12);
                    
                    Zxx=(beta.^2./(beta.^2+alpha(g)^2)).*(-j*w*u0./(gamma0+gammah1.*coth(gammah1*d1)))+ ...
                        (alpha(g).^2./(beta.^2+alpha(g)^2)).*(-gammae1.*gamma0./(j*w*e0*(gammae1+gamma0*e12.*coth(gammae1*d1))));
                    
                    Kxx_matrix=Kxx_matrix+sum(Zxx.*rxm.*rxn.*R.*del_beta.*del_alpha(g));
                    
                end
                matrix_rect(m,n)=(1/(2*pi)^2)*Kxx_matrix;
            end
        end
        
        matrix=matrix_polar+matrix_rect;                                    %saves the rectangular
        %and polar integration results
        %in the impedance matrix.
        V=zeros(2*(M-1),1);                                                 %fills in voltage matrix
        ws=2*(1-cos(k*delx))/(k*sin(k*delx))*delz;
        V(ceil(M/2))=-Vs*ws/(2*delx);
        V(ceil(M/2)-1)=.5*V(ceil(M/2));
        V(ceil(M/2)+1)=.5*V(ceil(M/2));
        
        J=inv(matrix)*V                                                     %solves for surface current (A/m)
        I=J*delz;                                                           %solves for total current (A)
        Z_in(s)=Vs/I(ceil(M/2))                                             %calculates short-circuit impedance
        
        matrix_mutual=matrix(1:M-1,1:M-1);
        
        V_mutual=zeros(M-1,1);                                              %removes parasitic dipole
        ws=2*(1-cos(k*delx))/(k*sin(k*delx))*delz;
        V_mutual(ceil(M/2))=-Vs*ws/(2*delx);
        V_mutual(ceil(M/2)-1)=.5*V_mutual(ceil(M/2));
        V_mutual(ceil(M/2)+1)=.5*V_mutual(ceil(M/2));
        I_mutual=inv(matrix_mutual)*V_mutual*delz;
        
        I1p=I_mutual(ceil(M/2));
        Z_oc(s)=Vs/I1p                                                      %solves for open-circuit voltage
        Z_m(s)=sqrt(Z_oc(s)*(Z_oc(s)-Z_in(s)))                              %solves for mutual impedance
        
    end
    
end