function conversion()

    global tiempo

    global Mo
    global I3o

    to = 0;         % s
    tf = 4000;      % s 

    tiempo  = linspace(to, tf, 100000);     % Vector tiempo de reacci√≥n
    step = 100;                             % M√°ximo intervalo de tiempo que se dar√° a ODE

    inicial = zeros(1, 7);  % Vector iniciado a 0 de las 8 variables del sistema
    inicial(1) = I3o;
    inicial(2) = Mo;
    
    [t1,x1] = ode23s(@(t,x)sistemaDinamico120(t,x), tiempo, inicial, odeset('Maxstep', step)); % 120 ∫C
    [t2,x2] = ode23s(@(t,x)sistemaDinamico130(t,x), tiempo, inicial, odeset('Maxstep', step)); % 130 ∫C
    
    iniciador120  = x1(:,1);                     % Iniciador(t)
    
    monomero120   = x1(:,2);                     % Monomero(t)   a 120∫C
    conversion120 = (Mo - monomero120(:))/Mo;    % Conversion(t) a 120∫C
    
    rad120        = x1(:,4);                     % R1(t)         a 120∫C
    rads120       = x1(:,5);                     % R2(t)         a 120∫C
    radicales120  = rad120 + rads120;            % R(t) totales  a 120∫C
    rp120         = monomero120 .* radicales120; % Rp(t)         a 120∫C
    
    % --- %
    
    iniciador130  = x2(:,1);                     % Iniciador(t)
    
    monomero130   = x2(:,2);                     % Monomero(t)   a 120∫C
    conversion130 = (Mo - monomero130(:))/Mo;    % Conversion(t) a 120∫C
    
    rad130        = x2(:,4);                     % R1(t)         a 120∫C
    rads130       = x2(:,5);                     % R2(t)         a 120∫C
    radicales130  = rad130 + rads130;            % R(t) totales  a 120∫C
    rp130         = monomero130 .* radicales130; % Rp(t)         a 130∫C
    
   
    % CreaciÛn del gr·fico
    
    pantalla = get(groot, 'ScreenSize');
    figure('Position',[0 0 pantalla(3)/2 pantalla(4)/2])
        pintarIniciador(t1, iniciador120, iniciador130)
        
    pantalla = get(groot, 'ScreenSize');
    figure('Position',[pantalla(3)/2 0 pantalla(3)/2 pantalla(4)])
        pintarConversiones(t2, conversion120, conversion130, rp120, rp130)
    
end