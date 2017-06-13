function k = constantes(X, T, R, B, C)
    
%     kd  = k(1);
%     ki0 = k(2);
%     ki1 = k(3);
%     kp  = k(4);
%     ktd = k(5);
%     ktc = k(6);
%     kfM = k(7);

    %k(1)  = 1.4282 * 10^10 * exp(-99064.12/(R*T));   % Datos Mejico
    k(1)  = 1.43 * 10^10 * exp(-111380/(R*T));         % Ajuste
    k(2) = 0;
    k(3) = 4.92 * 10^5 * exp(-18195.54/(R*T));      % Mahabadi

    %k(4)  = 4.92 * 10^5 * exp(-18195.54/(R*T));      % Mahabadi
    k(4)  = 3.4 * 10^5 * exp(-18195.54/(R*T));      % Mahabadi

    %k(5)  = 9.80 * 10^7 * exp(-2930.180/(R*T));      % Mahabadi
        ktd0 = 9.80 * 10^7 * exp(-2930.180/(R*T)); 
    k(5)  = ktd0 * exp( B*X + C*X^2);                 % Friis

    %k(6) = 0.9 * k(5);                                 % Matthew Justin
    %k(6) = (X/(1-X)) * k(5);                          % Pablo I
    k(6) = (X/(1-X^2)) * k(5)/10;
    %k(6) = (X/(1-X)^2) * k(5);

    %k(7)  = 2.41 * 10^9 * exp(-2930.180/(R*T));      % Matthew Justin
    k(7) = 2.012 * 10^9 * exp(-70454.84/(R*T));     % Ajuste del trabajo

end