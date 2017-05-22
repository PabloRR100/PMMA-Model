function conversion()

    global tiempo

    global Mo
    global I3o

    to = 0;         % s
    tf = 5400;      % s 

    tiempo  = linspace(to, tf, 100000);     % Vector tiempo de reacci√≥n
    step = 100;                             % M√°ximo intervalo de tiempo que se dar√° a ODE

    inicial = zeros(1, 7);  % Vector iniciado a 0 de las 8 variables del sistema
    inicial(1) = I3o;
    inicial(2) = Mo;
    
    [t1,x1] = ode23s(@(t,x)sistemaDinamico110(t,x), tiempo, inicial, odeset('Maxstep', step)); % 110 ∫C
    [t2,x2] = ode23s(@(t,x)sistemaDinamico120(t,x), tiempo, inicial, odeset('Maxstep', step)); % 120 ∫C
    [t3,x3] = ode23s(@(t,x)sistemaDinamico130(t,x), tiempo, inicial, odeset('Maxstep', step)); % 130 ∫C

    
    %iniciador110  = x1(:,1);                     % Iniciador(t)
    
    monomero110   = x1(:,2);                     % Monomero(t)   a 110∫C
    conversion110 = (Mo - monomero110(:))/Mo;    % Conversion(t) a 110∫C
    
    rad110        = x1(:,4);                     % R1(t)         a 110∫C
    rads110       = x1(:,5);                     % R2(t)         a 110∫C
    radicales110  = rad110 + rads110;            % R(t) totales  a 110∫C
    rp110         = monomero110 .* radicales110; % Rp(t)         a 110∫C
    
    % --- %
    
    %iniciador120  = x2(:,1);                     % Iniciador(t)
    
    monomero120   = x2(:,2);                     % Monomero(t)   a 120∫C
    conversion120 = (Mo - monomero120(:))/Mo;    % Conversion(t) a 120∫C
    
    rad120        = x2(:,4);                     % R1(t)         a 120∫C
    rads120       = x2(:,5);                     % R2(t)         a 120∫C
    radicales120  = rad120 + rads120;            % R(t) totales  a 120∫C
    rp120         = monomero120 .* radicales120; % Rp(t)         a 120∫C
    
    % --- %
    
    %iniciador130  = x3(:,1);                     % Iniciador(t)
    
    monomero130   = x3(:,2);                     % Monomero(t)   a 130∫C
    conversion130 = (Mo - monomero130(:))/Mo;    % Conversion(t) a 130∫C
    
    rad130        = x3(:,4);                     % R1(t)         a 130∫C
    rads130       = x3(:,5);                     % R2(t)         a 130∫C
    radicales130  = rad130 + rads130;            % R(t) totales  a 130∫C
    rp130         = monomero130 .* radicales130; % Rp(t)         a 130∫C
    
   
    % CreaciÛn del gr·fico
    
    %pintarIniciador(t2, iniciador120, iniciador130)
       
    pintarConversiones(t3, conversion110, conversion120, conversion130, rp110, rp120, rp130)
    
end