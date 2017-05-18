
% Sistema dinamico en derivadas parciales para la polimerizacion

function f = sistemaDinamico1 (t,x)

    %Condiciones iniciales
    
    global Mo                           % Monomero inicial
    global I3o                          % Iniciador inicial
    global ef                           % Eficiencia iniciador
    global T1
    global T2
    
    % Variables
    
        I3   = x(1);        % x(1) = Iniciador trifuncional 
        M    = x(2);        % x(2) = Monomero   
        I2p2 = x(3);        % x(3) -> Diradical, 2 peroxidos sin descomponer

        R1   = x(4);        % x(4) -> Monoradicales
        R2   = x(5);        % x(5) -> Diradicales

        P    = x(6);        % x(6) -> Polimero 

        PeP  = x(7);        % x(7) -> Peroxidos totales en PMMA

    
        X = (Mo - M) / (Mo);
        
    % Parametros cineticos
   
        kd  = 0;
        ki0 = 0;
        ki1 = 0;
        kp  = 0;
        kfM = 0;
        ktc = 0;
        ktd = 0;

            
    % Balances ------------------------------------------------------------
    
                            % Evolucion de DEKTP (I3)                               - A.1
        f(1) = -3*efic*kd*I3; 
                            % Evolucion de monomero (M)                             - A.2
        f(2) = -kp*M*(R1+2*R2);    
                            % Evolucion de I2p2                                     - A.5
        f(3) = 3*efic*kd*I3 - 2*ki1*M*I2p2;
                            % Evolucion de R1                                       - A.11
        f(4) = 2*ki0*M^3 + 2*kfM*M*R2 - ktc*R1^2 - ktd*R1*(R1+2*R2) + 2*ktd*R2*(R1+2*R2) + 2*efic*kd*PeP;
                            % Evolucion de R2                                       - A.12
        f(5) = 2*(ki1*I2p2 - kfM*R2)*M - 2*ktc*R2*(R1+2*R2) - 2*ktd*R2*(R1+2*R2);          
                            % Evolucion de P                                        - A.16
        f(6) = kfM*M*R1 + ktc/2*R1^2 + 2*ktd*R1*(R1*R2) - kd*PeP;                  
                            % Evolucion de PeP                                      - A.23
        f(7) = ki1*I2p2*M - kd*PeP;

        f = f';

end