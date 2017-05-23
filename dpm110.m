
global tiempo
global T1
global Rjul

nmax = 12000;   % Máxima longitud de cadena
n = 1:nmax;     % Vector de 1 en 1 hasta n

denominador(1:length(tiempo)) = 0;    % Inicializar A con el tamaño del vector tiempo
alp(1:length(tiempo)) = 0;  % Inicializar alp con el tamaños del vector tiempo

% Vectores de las especies

    % Diradicales --> Inicializamos a 0
    
        R0 (1:nmax) = 0;
        R1 (1:nmax) = 0;
        R2 (1:nmax) = 0;
        R3 (1:nmax) = 0;
        R4 (1:nmax) = 0;
        R5 (1:nmax) = 0;
        R6 (1:nmax) = 0;
        R7 (1:nmax) = 0;
        R8 (1:nmax) = 0;
        R9 (1:nmax) = 0;
        R10(1:nmax) = 0;
        
    % Monoradicales --> Inicializamos a 0
    
        r0 (1:nmax) = 0;
        r1 (1:nmax) = 0;
        r2 (1:nmax) = 0;
        r3 (1:nmax) = 0;
        r4 (1:nmax) = 0;
        r5 (1:nmax) = 0;
        r6 (1:nmax) = 0;
        r7 (1:nmax) = 0;
        r8 (1:nmax) = 0;
        r9 (1:nmax) = 0;
        r10(1:nmax) = 0; 
 
% Cálculo
        
    for t = 1:length(tiempo)-1
        
        % Constantes cinéticas

            B = -4;
            C = -5;

            T = T1;
            R = Rjul;
            k = constantes(X, T, R, B, C);

            kd  = k(1);
            ki0 = k(2);
            ki1 = k(3);
            kp  = k(4);
            ktd = k(5);
            ktc = k(6);
            kfM = k(7);
        
      denominador(t)   = (kp+kfM) * M(t) + (ktc+ktd)*(R1(t) + 2*R2(t));
      alp(t) = denominador(t) / (kp*M(t));
      
      R0(1) = 2*ki1*Ip2p*M(t) / denominador(t);
      
      r0(1) = (ki0 - 4*kfM*R2(t))*M(t) ;
        
    end
    
    
