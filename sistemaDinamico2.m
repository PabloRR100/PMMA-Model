
% Sistema dinamico en derivadas parciales para la polimerizacion a 130ºC

function f = sistemaDinamico2 (t,x)

    % CONDICIONES INICIALES
    
    global Rjul
    global Mo                           % Monomero inicial
    %global I3o                         % Iniciador inicial
    global ef                           % efiencia iniciador
    global T2
%     global B
%     global C
    
%     B = -2;
%     C = -1;
    
    % COMPONENTES DEL SISTEMA
    
        I3   = x(1);        % x(1) = Iniciador trifuncional 
        M    = x(2);        % x(2) = Monomero   
        I2p2 = x(3);        % x(3) -> Diradical, 2 peroxidos sin descomponer

        R1   = x(4);        % x(4) -> Monoradicales
        R2   = x(5);        % x(5) -> Diradicales

        P    = x(6);        % x(6) -> Polimero 

        PeP  = x(7);        % x(7) -> Peroxidos totales en PMMA

    
        X = (Mo - M) / (Mo);
        
    % CONSTANTES CINÉTICAS
   
        kd   = 2.5061 * 10^-4;   % Datos Mejico
        ki0  = 2.8711 * 10^-11;  % Emilio
        ki1  = 4.92 * 10^5 * exp(-18195.54/(Rjul*T2)); % Mahabadi
        
        kp   = 4.92 * 10^5 * exp(-18195.54/(Rjul*T2)); % Mahabadi
        %kp   = 4.92 * 10^7 * exp(-18153.74/(Rjul*T2)); % Matthew Justin --> SE DISPARA MUCHISIMO
         
        ktd  = 9.80 * 10^7 * exp(-2930.180/(Rjul*T2)); % Mahabadi
%             ktd0 = 9.80 * 10^7 * exp(-2930.180/(Rjul*T2)); 
%         ktd  = ktd0 * ( exp( B*X + C*X^2));            % Friis
        %ktd  = (kp * M)^2 / (2*kd*ef*I3*(2/(1-X))^2);  % Pablo
        
        ktc = 0.9 * ktd;                              % Matthew Justin
        %ktc = 9.80 * 10^7 * exp(-2930.180/(Rjul*T1)); % Hacerlas iguales
        
        %kfM = 0.9;              % Emilio
        kfM  = 2.41 * 10^9 * exp(-2930.180/(Rjul*T2)); % Matthew Justin

            
    % BALANCES DE COMPONENTES 
    
                            % Evolucion de DEKTP (I3)                               - A.1
        f(1) = -3*ef*kd*I3; 
                            % Evolucion de monomero (M)                             - A.2
        f(2) = -kp*M*(R1+2*R2);    
                            % Evolucion de I2p2                                     - A.5
        f(3) = 3*ef*kd*I3 - 2*ki1*M*I2p2;
                            % Evolucion de R1                                       - A.11
        f(4) = 2*ki0*M^3 + 2*kfM*M*R2 - ktc*R1^2 - ktd*R1*(R1+2*R2) + 2*ktd*R2*(R1+2*R2) + 2*ef*kd*PeP;
                            % Evolucion de R2                                       - A.12
        f(5) = 2*(ki1*I2p2 - kfM*R2)*M - 2*ktc*R2*(R1+2*R2) - 2*ktd*R2*(R1+2*R2);          
                            % Evolucion de P                                        - A.16
        f(6) = kfM*M*R1 + ktc/2*R1^2 + 2*ktd*R1*(R1*R2) - kd*PeP;                  
                            % Evolucion de PeP                                      - A.23
        f(7) = ki1*I2p2*M - kd*PeP;

        f = f';

end