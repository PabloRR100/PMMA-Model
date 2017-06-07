function k = constantes(X, T, R, B, C)
    
%     kd  = k(1);
%     ki0 = k(2);
%     ki1 = k(3);
%     kp  = k(4);
%     ktd = k(5);
%     ktc = k(6);
%     kfM = k(7);

    %k(1)  = 2.5061 * 10^-5;   % Datos Mejico
    k(1)  = 1.6*10^10 * exp(-111380/(R*T));
    k(2)  = 2.8711 * 10^-11;  % Emilio
    k(3)  = 4.92 * 10^5 * exp(-18195.54/(R*T));      % Mahabadi

    k(4)  = 4.92 * 10^5 * exp(-18195.54/(R*T));      % Mahabadi
    %k(4)  = 4.92 * 10^7 * exp(-18153.74/(R*T));      % Matthew Justin --> SE DISPARA MUCHISIMO

    %k(5)  = 9.80 * 10^7 * exp(-2930.180/(R*T));      % Mahabadi
        ktd0 = 9.80 * 10^7 * exp(-2930.180/(R*T)); 
    k(5)  = ktd0 * exp( B*X + C*X^2);                 % Friis
    %k(5)   = ktd0 * (1 + B*X + C*X^2);                   % Pablo
    %ktd  = (kp * M)^2 / (2*kd*ef*I3*(2/(1-X))^2);        % Pablo2

    %k(6) = 0;                                          % Solo desproporción
    k(6) = 0.9 * k(5);                                 % Matthew Justin
    %k(6) = (X/(1-X)) * k(5);                          % Pablo 2
    %k(6) = 3.956*10^4 * exp(17096/(R*T));             % Bevington
    %k(6) = k(7);                                       % Hacerlas iguales

    %k(7)  = 0.05;
    %k(7) = k(4) * 9480 * exp(-58113/(R*T));
    k(7)  = 2.41 * 10^9 * exp(-2930.180/(R*T));      % Matthew Justin



end