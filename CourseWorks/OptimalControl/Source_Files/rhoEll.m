function [val, point] = rhoEll(psi,ell_conf)
%RHOELL Summary of this function goes here
%   this is support function of ELLIPSE implentation   
%Detailed explanation goes here
%   psi - vec arg (2xN); point - support vec ;
%   val - value of support function rho
%   ell_conf = [b,c,phi,q]
%   see more comments in latex attached file

%CONST:
    b = ell_conf{1}; %semimajor axis
    c = ell_conf{2}; %semiminor axis
    phi = ell_conf{3}; %anticlockwise rotation angle
    q = ell_conf{4}; %vector(column) of shift
    
    %2x2 matrix of ellipse configuration
        a_11 = b;   a_12 = 0;
        a_21 = 0;   a_22 = c;
        A = [a_11 a_12; a_21 a_22];
    %2x2 matrix of anticlockwise rotation with angle phi
        r_11 = cos(phi);    r_12 = -sin(phi);
        r_21 = sin(phi);    r_22 = cos(phi);
        R = [r_11 r_12; r_21 r_22];
    %result matrix W=RCB
        W = R*A;
        W_t = W';
        
%therefore Ell_q = W*Circ_1(0) + q
%sup func of Circ_1(0) = sqrt((psi_1)^2 +(psi_2)^2)
%value of support function
    W_t_psi = W_t*psi; %2xN
    val = sqrt(W_t_psi(1,:).^2 + W_t_psi(2,:).^2) + q'*psi; %1xN 
   
    normW_t_psi = sqrt(W_t_psi(1,:).^2 + W_t_psi(2,:).^2); %1xN
    x_0 =  W_t_psi./normW_t_psi; %2xN
    y_0 = W*x_0 + q; %2xN
    point = y_0; %2xN
end

