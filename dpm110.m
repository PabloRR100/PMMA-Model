function dpm110(x)

    global tiempo
    global T1
    global Rjul
    global B
    global C
    global Mo
    global Vo
    global nmax
    
    nmax = 1000;           % M�xima longitud de cadena
    long = length(tiempo);  % Longitud de los vectores de variables (t)
    V = Vo;
    
    denominador(1:long) = 0;      % Inicializar denominador con el tama�o del vector tiempo
    alp(1:long) = 0;              % Inicializar alp con el tama�os del vector tiempo

    % Vectores de las especies

        % Diradicales --> Inicializamos a 0

           %Rn = Diradical con N grupos per�xidos
            R0 = zeros(1, nmax);
            R1 = zeros(1, nmax);
            R2 = zeros(1, nmax);
            R3 = zeros(1, nmax);
            R4 = zeros(1, nmax);
            R5 = zeros(1, nmax);
            R6 = zeros(1, nmax);
            R7 = zeros(1, nmax);
            R8 = zeros(1, nmax);

        % Monoradicales --> Inicializamos a 0

           %rn = Monoradical con n grupos per�xidos
            r0 = zeros(1, nmax);
            r1 = zeros(1, nmax);
            r2 = zeros(1, nmax);
            r3 = zeros(1, nmax);
            r4 = zeros(1, nmax);
            r5 = zeros(1, nmax);
            r6 = zeros(1, nmax);
            r7 = zeros(1, nmax);
            r8 = zeros(1, nmax);

                % Generaci�n de monoradicales por desproporci�n
        
                   %termDrn = Monoradical con n grupos per�xidos por Desproporci�n 
                    termdr0 = zeros(1, nmax-1);
                    termdr1 = zeros(1, nmax-1);
                    termdr2 = zeros(1, nmax-1);
                    termdr3 = zeros(1, nmax-1);
                    termdr4 = zeros(1, nmax-1);
                    termdr5 = zeros(1, nmax-1);
                    termdr6 = zeros(1, nmax-1);
                    termdr7 = zeros(1, nmax-1);
                    termdr8 = zeros(1, nmax-1);

        % Pol�meros --> Inicializamos a 0

           %pn = Pol�mero con n grupos per�xidos
            p0 = zeros(long, nmax);
            p1 = zeros(long, nmax);
            p2 = zeros(long, nmax);
            p3 = zeros(long, nmax);
            p4 = zeros(long, nmax);
            p5 = zeros(long, nmax);
            p6 = zeros(long, nmax);
            p7 = zeros(long, nmax);
            p8 = zeros(long, nmax);
            
               %N grupos per�xidos sin descomponer
                NPS0 = zeros(long, nmax);
                NPS1 = zeros(long, nmax);
                NPS2 = zeros(long, nmax);
                NPS3 = zeros(long, nmax);
                NPS4 = zeros(long, nmax);
                NPS5 = zeros(long, nmax);
                NPS6 = zeros(long, nmax);
                NPS7 = zeros(long, nmax);
                NPS8 = zeros(long, nmax);
            

    % C�lculo
    
        M    = x(:,2);
        I2p2 = x(:,3);
        RT   = x(:,4) + 2.*x(:,5);
    
    for t = 1:long-1                    % BUCLE PARA AVANZAR EN EL TIEMPO

        % Constantes cin�ticas

            T = T1;
            R = Rjul;
            X    = (x(t,2) - Mo) ./ x(t,2);
            k = constantes(X, T, R, B, C);

            kd  = k(1);
            ki0 = k(2);
            ki1 = k(3);
            kp  = k(4);
            ktd = k(5);
            ktc = k(6);
            kfM = k(7);

            denominador(t) = (kp+kfM)*M(t) + (ktc+ktd)*(RT(t));   % Factor com�n al radical
            alp(t) = denominador(t) / (kp*M(t));                  % Para poner 1 sobre kp*M(t)

        % Dar el valor inicial antes del bucle;

            % Diradicales

                %R0(1) = 2*ki1*I2p0(t)*M(t) / denominador(t);
                %R1(1) = 2*ki1*I2p1(t)*M(t) / denominador(t);
                R2(1) = 2*ki1*I2p2(t)*M(t) / denominador(t);

            % Monoradicales

                %r0(1) = (ki1*Ip0(t) + ki0 + 2*kfM*R0(t)) * M(t) / denominador(t);
                %r1(1) = (ki1*Ip1(t) + ki0 + 2*kfM*R1(t)) * M(t) / denominador(t);
                r2(1) = (ki1*I2p2(t) + ki0 + 2*kfM*R2(t)) * M(t) / denominador(t);

            % Pol�meros

                p0(t,1) = kfM*r0(1)*M(t);
                p1(t,1) = kfM*r1(1)*M(t);
                p2(t,1) = kfM*r2(1)*M(t);
                p3(t,1) = kfM*r3(1)*M(t);
                p4(t,1) = kfM*r4(1)*M(t);
                p5(t,1) = kfM*r5(1)*M(t);
                p6(t,1) = kfM*r6(1)*M(t);
                p7(t,1) = kfM*r7(1)*M(t);
                p8(t,1) = kfM*r8(1)*M(t);
              
        % Generaci�n de diradicales por terminaci�n (Solo combinaci�n)
              
            termcR0 = ktc / (kp*M(t)) * (  conv(R0, R0));
            termcR1 = ktc / (kp*M(t)) * (2*conv(R0, R1));
            termcR2 = ktc / (kp*M(t)) * (2*conv(R2, R0) +   conv(R1, R1));
            termcR3 = ktc / (kp*M(t)) * (2*conv(R3, R0) + 2*conv(R2, R1));
            termcR4 = ktc / (kp*M(t)) * (2*conv(R4, R0) + 2*conv(R3, R1) +   conv(R2, R2));
            termcR5 = ktc / (kp*M(t)) * (2*conv(R5, R0) + 2*conv(R4, R1) + 2*conv(R3, R2));
            termcR6 = ktc / (kp*M(t)) * (2*conv(R6, R0) + 2*conv(R5, R1) + 2*conv(R4, R2) +   conv(R3, R3));
            termcR7 = ktc / (kp*M(t)) * (2*conv(R7, R0) + 2*conv(R6, R1) + 2*conv(R5, R2) + 2*conv(R4, R3));
            termcR8 = ktc / (kp*M(t)) * (2*conv(R8, R0) + 2*conv(R7, R1) + 2*conv(R6, R2) + 2*conv(R5, R3) + conv(R4, R4));
                
        % Generaci�n de monoradicales por terminaci�n por combinaci�n
              
            termcr0 = 2*ktc / (kp*M(t)) * (conv(R0, r0));
            termcr1 = 2*ktc / (kp*M(t)) * (conv(R1, r0) + conv(R0, r1));
            termcr2 = 2*ktc / (kp*M(t)) * (conv(R2, r0) + conv(R1, r1) + conv(R0, r2));
            termcr3 = 2*ktc / (kp*M(t)) * (conv(R3, r0) + conv(R2, r1) + conv(R1, r2) + conv(R0, r3));
            termcr4 = 2*ktc / (kp*M(t)) * (conv(R4, r0) + conv(R3, r1) + conv(R2, r3) + conv(R1, r3) + conv(R0, r4));
            termcr5 = 2*ktc / (kp*M(t)) * (conv(R5, r0) + conv(R4, r1) + conv(R3, r4) + conv(R2, r3) + conv(R1, r4) + conv(R0, r5));
            termcr6 = 2*ktc / (kp*M(t)) * (conv(R6, r0) + conv(R5, r1) + conv(R4, r5) + conv(R3, r3) + conv(R2, r4) + conv(R1, r5) + conv(R0, r6));
            termcr7 = 2*ktc / (kp*M(t)) * (conv(R7, r0) + conv(R6, r1) + conv(R5, r6) + conv(R4, r3) + conv(R3, r4) + conv(R2, r5) + conv(R1, r6) + conv(R0, r7));
            termcr8 = 2*ktc / (kp*M(t)) * (conv(R8, r0) + conv(R7, r1) + conv(R6, r7) + conv(R5, r3) + conv(R4, r4) + conv(R3, r5) + conv(R2, r6) + conv(R1, r7) + conv(R0, r8));
                
        % Generaci�n de monoradicales ��??
        
%             genr1 = zeros(nmax-1, 1:1);
%             genr2 = zeros(nmax-1, 1:2);
%             genr3 = zeros(nmax-1, 1:3);
%             genr4 = zeros(nmax-1, 1:4);
%             genr5 = zeros(nmax-1, 1:5);
%             genr6 = zeros(nmax-1, 1:6);
%             genr7 = zeros(nmax-1, 1:7);
%             genr8 = zeros(nmax-1, 1:8);
            
            
        for n = 2:nmax                      % Bucle para sacar el n a partir del n-1 de diradicales

            % Diradicales en cada t (no se almacenan)

                R0(n) = (R0(n-1) + termcR0(n-1)) / alp(t);
                R1(n) = (R1(n-1) + termcR1(n-1)) / alp(t);
                R2(n) = (R2(n-1) + termcR2(n-1)) / alp(t);
                R3(n) = (R3(n-1) + termcR3(n-1)) / alp(t);
                R4(n) = (R4(n-1) + termcR4(n-1)) / alp(t);
                R5(n) = (R5(n-1) + termcR5(n-1)) / alp(t);
                R6(n) = (R6(n-1) + termcR6(n-1)) / alp(t);
                R7(n) = (R7(n-1) + termcR7(n-1)) / alp(t);
                R8(n) = (R8(n-1) + termcR8(n-1)) / alp(t);

                % Generaci�n de monoradicales por terminaci�n por desproporci�n (no se almacenan)

                    termdr0(n) = 2*ktd / (kp*M(t)) * R0(n)*RT(t);
                    termdr1(n) = 2*ktd / (kp*M(t)) * R1(n)*RT(t);
                    termdr2(n) = 2*ktd / (kp*M(t)) * R2(n)*RT(t);
                    termdr3(n) = 2*ktd / (kp*M(t)) * R3(n)*RT(t);
                    termdr4(n) = 2*ktd / (kp*M(t)) * R4(n)*RT(t);
                    termdr5(n) = 2*ktd / (kp*M(t)) * R5(n)*RT(t);
                    termdr6(n) = 2*ktd / (kp*M(t)) * R6(n)*RT(t);
                    termdr7(n) = 2*ktd / (kp*M(t)) * R7(n)*RT(t);
                    termdr8(n) = 2*ktd / (kp*M(t)) * R8(n)*RT(t);

            % Monoradicales en cada t (no se almacenan)

                r0(n) = (r0(n-1) + termdr0(n) + 2*kfM/kp*R0(n) + termcr0(n-1)) / alp(t);
                r1(n) = (r1(n-1) + termdr1(n) + 2*kfM/kp*R1(n) + termcr1(n-1)) / alp(t);
                r2(n) = (r2(n-1) + termdr2(n) + 2*kfM/kp*R2(n) + termcr2(n-1)) / alp(t);
                r3(n) = (r3(n-1) + termdr3(n) + 2*kfM/kp*R3(n) + termcr3(n-1)) / alp(t);
                r4(n) = (r4(n-1) + termdr4(n) + 2*kfM/kp*R4(n) + termcr4(n-1)) / alp(t);
                r5(n) = (r5(n-1) + termdr5(n) + 2*kfM/kp*R5(n) + termcr5(n-1)) / alp(t);
                r6(n) = (r6(n-1) + termdr6(n) + 2*kfM/kp*R6(n) + termcr6(n-1)) / alp(t);
                r7(n) = (r7(n-1) + termdr7(n) + 2*kfM/kp*R7(n) + termcr7(n-1)) / alp(t);
                r8(n) = (r8(n-1) + termdr8(n) + 2*kfM/kp*R8(n) + termcr8(n-1)) / alp(t);           

        end
        
        % Generaci�n de cadena de pol�mero por terminaci�n por combinaci�n
            
            termp0 = 0.5*ktc * (  conv(r0, r0));
            termp1 = 0.5*ktc * (2*conv(r1, r0));
            termp2 = 0.5*ktc * (2*conv(r2, r0) + conv(r1,r1));
            termp3 = 0.5*ktc * (  conv(r3, r0) + conv(r2,r1));
            termp4 = 0.5*ktc * (  conv(r4, r0) + conv(r3,r1) + conv(r2,r2));
            termp5 = 0.5*ktc * (  conv(r5, r0) + conv(r4,r1) + conv(r3,r2));
            termp6 = 0.5*ktc * (  conv(r6, r0) + conv(r5,r1) + conv(r4,r2) + conv(r3,r3));
            termp7 = 0.5*ktc * (  conv(r7, r0) + conv(r6,r1) + conv(r5,r2) + conv(r4,r3));
            termp8 = 0.5*ktc * (  conv(r8, r0) + conv(r7,r1) + conv(r6,r2) + conv(r5,r3) + conv(r4,r4));
            
        
        for n = 2:nmax
                
             % GENERACI�N Pol�mero en cada t

                p0(t,n) = kfM*r0(n)*M(t) + termp0(n-1) + ktd*r0(n)*RT(t);
                p1(t,n) = kfM*r1(n)*M(t) + termp1(n-1) + ktd*r1(n)*RT(t);
                p2(t,n) = kfM*r2(n)*M(t) + termp2(n-1) + ktd*r2(n)*RT(t);
                p3(t,n) = kfM*r3(n)*M(t) + termp3(n-1) + ktd*r3(n)*RT(t);
                p4(t,n) = kfM*r4(n)*M(t) + termp4(n-1) + ktd*r4(n)*RT(t);
                p5(t,n) = kfM*r5(n)*M(t) + termp5(n-1) + ktd*r5(n)*RT(t);
                p6(t,n) = kfM*r6(n)*M(t) + termp6(n-1) + ktd*r6(n)*RT(t);
                p7(t,n) = kfM*r7(n)*M(t) + termp7(n-1) + ktd*r7(n)*RT(t);
                p8(t,n) = kfM*r8(n)*M(t) + termp8(n-1) + ktd*r8(n)*RT(t);

                % MOLES de pol�mero 
                NPS0(n+1,1) = NPS0(n,1) + p0(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS1(n+1,1) = NPS1(n,1) + p1(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS2(n+1,1) = NPS2(n,1) + p2(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS3(n+1,1) = NPS3(n,1) + p3(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS4(n+1,1) = NPS4(n,1) + p4(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS5(n+1,1) = NPS5(n,1) + p5(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS6(n+1,1) = NPS6(n,1) + p6(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS7(n+1,1) = NPS7(n,1) + p7(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS8(n+1,1) = NPS8(n,1) + p8(t,n)*V*(tiempo(t+1)-tiempo(t));

        end
            
        NPS0(NPS0<0)=0;
        NPS1(NPS1<0)=0;
        NPS2(NPS2<0)=0;
        NPS3(NPS3<0)=0;
        NPS4(NPS4<0)=0;
        NPS5(NPS5<0)=0;
        NPS6(NPS6<0)=0;
        NPS7(NPS7<0)=0;
        NPS8(NPS8<0)=0;

    end
    
    pesos(NPS0, NPS1, NPS2, NPS3, NPS4, NPS5, NPS6, NPS7, NPS8)
    
end
