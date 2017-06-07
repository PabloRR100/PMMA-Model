
% Sistema dinamico en derivadas parciales para la polimerizacion a 120�C

function f = sistemaDinamico120(t,x)

    % CONDICIONES INICIALES
    
    global Rjul
    global Mo                           % Monomero inicial
    global ef                           % efiencia iniciador
    global T2
    global B
    global C

    
    % COMPONENTES DEL SISTEMA
    
        I3   = x(1);        % x(1) = Iniciador trifuncional 
        M    = x(2);        % x(2) = Monomero   
        I2p2 = x(3);        % x(3) -> Diradical, 2 peroxidos sin descomponer

        R1   = x(4);        % x(4) -> Monoradicales
        R2   = x(5);        % x(5) -> Diradicales

        %P    = x(6);        % x(6) -> Polimero 

        PeP  = x(7);        % x(7) -> Peroxidos totales en PMMA

        Rt = R1 + 2*R2;         % Radicales totales
        X  = (Mo - M) / (Mo);   % Conversi�n
        
    % CONSTANTES CIN�TICAS
        
        T = T2;
        R = Rjul;
        k = constantes(X, T, R, B, C);
        
        kd  = k(1);
        ki0 = k(2);
        ki1 = k(3);
        kp  = k(4);
        ktd = k(5);
        ktc = k(6);
        kfM = k(7);
   
            
    % BALANCES DE COMPONENTES 
    
                            % Evolucion de DEKTP (I3)                               - A.1
        f(1) = -3*ef*kd*I3; 
                            % Evolucion de monomero (M)                             - A.2
        f(2) = -kp*M*Rt;    
                            % Evolucion de I2p2                                     - A.5
        f(3) = 3*ef*kd*I3 - 2*ki1*M*I2p2;
                            % Evolucion de R1                                       - A.11
        f(4) = ki0*M + 4*kfM*M*R2 - 2*(ktc+ktd)*R1^2 + 8*ktd*R2^2 + 2*ef*kd*PeP;
                            % Evolucion de R2                                       - A.12
        f(5) = 2*(ki1*I2p2 - kfM*R2)*M - 2*ktc*R2*Rt - 2*ktd*R2*(R1+2*R2);
                            % Evolucion de P                                        - A.16
        %f(6) = kfM*M*R1 + ktc/2*R1^2 + ktd*R1*Rt - kd*PeP;
                            % Evolucion de PeP                                      - A.23
        f(7) = ki1*I2p2*M - kd*PeP;

        f = f';

end