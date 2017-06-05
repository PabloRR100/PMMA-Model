function f = conversion()

    global tiempo
    global Rjul
    global T1
    global T2
    global T3

    global Mo
    global I3o

    to = 0;         % s
    tf = 5400;      % s 

    tiempo  = linspace(to, tf, 100000);     % Vector tiempo de reacción
    step = 100;                             % Máximo intervalo de tiempo que se dará a ODE

    inicial = zeros(1, 7);  % Vector iniciado a 0 de las 8 variables del sistema
    inicial(1) = I3o;
    inicial(2) = Mo;
    
    [t1,x1] = ode23s(@(t,x)sistemaDinamico110(t,x), tiempo, inicial, odeset('Maxstep', step)); % 110 �C
    [t2,x2] = ode23s(@(t,x)sistemaDinamico120(t,x), tiempo, inicial, odeset('Maxstep', step)); % 120 �C
    [t3,x3] = ode23s(@(t,x)sistemaDinamico130(t,x), tiempo, inicial, odeset('Maxstep', step)); % 130 �C

    
    %iniciador110  = x1(:,1);                     % Iniciador(t)
    
    k110 = constantes(0, T1, Rjul, 0, 0);
        kp110 = k110(4);
    monomero110   = x1(:,2);                                % Monomero(t)   a 110�C
    polimero110   = x1(:,6);                                % Polimero(t)   a 110�C
    PeP110        = x1(:,7);                                % PeP(t)        a 110�C
    conversion110 = (Mo - monomero110(:))/Mo;               % Conversion(t) a 110�C
    
    rad110        = x1(:,4);                                % R1(t)         a 110�C
    rads110       = x1(:,5);                                % R2(t)         a 110�C
    radicales110  = rad110 + 2.*rads110;                    % R(t) totales  a 110�C
    rp110         = kp110 .* monomero110 .* radicales110;   % Rp(t)         a 110�C
    
    % --- %
    
    %iniciador120  = x2(:,1);                     % Iniciador(t)
    
    k120 = constantes(0, T2, Rjul, 0, 0);
        kp120 = k120(4);
    monomero120   = x2(:,2);                                % Monomero(t)   a 120�C
    polimero120   = x2(:,6);                                % Polimero(t)   a 120�C
    PeP120        = x2(:,7);                                % PeP(t)        a 120�C
    conversion120 = (Mo - monomero120(:))/Mo;               % Conversion(t) a 120�C
    
    rad120        = x2(:,4);                                % R1(t)         a 120�C
    rads120       = x2(:,5);                                % R2(t)         a 120�C
    radicales120  = rad120 + 2.*rads120;                    % R(t) totales  a 120�C
    rp120         = kp120 .* monomero120 .* radicales120;   % Rp(t)         a 120�C
    
    % --- %
    
    %iniciador130  = x3(:,1);                     % Iniciador(t)
    
    k130 = constantes(0, T3, Rjul, 0, 0);
        kp130 = k130(4);
    monomero130   = x3(:,2);                                % Monomero(t)   a 130�C
    polimero130   = x3(:,6);                                % Polimero(t)   a 130�C
    PeP130        = x3(:,7);                                % PeP(t)        a 130�C
    conversion130 = (Mo - monomero130(:))/Mo;               % Conversion(t) a 130�C
    
    rad130        = x3(:,4);                                % R1(t)         a 130�C
    rads130       = x3(:,5);                                % R2(t)         a 130�C
    radicales130  = rad130 + 2.*rads130;                        % R(t) totales  a 130�C
    rp130         = kp130 .* monomero130 .* radicales130;   % Rp(t)         a 130�C
    
   
    % Creaci�n del gr�fico
    
    %pintarIniciador(t2, iniciador120, iniciador130)
       
    pintarConversiones(t3, conversion110, conversion120, conversion130, rp110, rp120, rp130)
    %pintarPolimero(t1, polimero110, polimero120, polimero130, PeP110, PeP120, PeP130)
        % Estamos pintando mal, esto viene despu�s del M�dulo de Dist.
      
        % Comentado porque lo estamos ejecutando desde main para ahorrar
        % hacer conversion en cada prueba
%     dpm110(x1)
%     dpm120(x2)
%     dpm130(x3)
   
end