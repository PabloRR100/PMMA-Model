function [Mn, Mw, X] = dpm110(tiempo, x, T)


    global Rjul
    global B
    global C
    global Mo
    global Vo
    global nmax

    long = length(tiempo);   % Longitud de los vectores de variables (t)
    V = Vo;

    denominador(1:long) = 0;      % Inicializar denominador con el tamaño del vector tiempo
    alp(1:long) = 0;              % Inicializar alp con el tamaños del vector tiempo

    % Vectores de las especies

        % Diradicales --> Inicializamos a 0

           %Rn = Diradical con N grupos peróxidos
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

           %rn = Monoradical con n grupos peróxidos
            r0 = zeros(1, nmax);
            r1 = zeros(1, nmax);
            r2 = zeros(1, nmax);
            r3 = zeros(1, nmax);
            r4 = zeros(1, nmax);
            r5 = zeros(1, nmax);
            r6 = zeros(1, nmax);
            r7 = zeros(1, nmax);
            r8 = zeros(1, nmax);
            
            %genrn = Generaci�n de monoradicales por ruptura de cadena
            genr1 = zeros(nmax, 1);
            genr2 = zeros(nmax, 2);
            genr3 = zeros(nmax, 3);
            genr4 = zeros(nmax, 4);
            genr5 = zeros(nmax, 5);
            genr6 = zeros(nmax, 6);
            genr7 = zeros(nmax, 7);
            genr8 = zeros(nmax, 8);

                % Generación de monoradicales por desproporción

                   %termdrn = Monoradical con n grupos peróxidos por Desproporción 
                    termdr0 = zeros(1, nmax-1);
                    termdr1 = zeros(1, nmax-1);
                    termdr2 = zeros(1, nmax-1);
                    termdr3 = zeros(1, nmax-1);
                    termdr4 = zeros(1, nmax-1);
                    termdr5 = zeros(1, nmax-1);
                    termdr6 = zeros(1, nmax-1);
                    termdr7 = zeros(1, nmax-1);
                    termdr8 = zeros(1, nmax-1);

        % Polímeros --> Inicializamos a 0

           %genPn = Término de --> Generación de Polímero con n grupos peróxidos
            genP0 = zeros(long, nmax);
            genP1 = zeros(long, nmax);
            genP2 = zeros(long, nmax);
            genP3 = zeros(long, nmax);
            genP4 = zeros(long, nmax);
            genP5 = zeros(long, nmax);
            genP6 = zeros(long, nmax);
            genP7 = zeros(long, nmax);
            genP8 = zeros(long, nmax);

               %NPSn --> Moles de Polímero con N grupos peróxidos sin descomponer
                NPS0 = zeros(long, nmax);
                NPS1 = zeros(long, nmax);
                NPS2 = zeros(long, nmax);
                NPS3 = zeros(long, nmax);
                NPS4 = zeros(long, nmax);
                NPS5 = zeros(long, nmax);
                NPS6 = zeros(long, nmax);
                NPS7 = zeros(long, nmax);
                NPS8 = zeros(long, nmax);


    % Cálculo

        M    = x(:,2);
        I2p2 = x(:,3);
        RT   = x(:,4) + 2.*x(:,5);
        X    = (Mo - x(:,2)) / Mo;

    for t = 1:long-1                    % BUCLE PARA AVANZAR EN EL TIEMPO

        % Constantes cinéticas

            R = Rjul;

            k = constantes(X(t), T, R, B, C);

            ki0 = k(2);
            ki1 = k(3);
            kp  = k(4);
            ktd = k(5);
            ktc = k(6);
            kfM = k(7);
            kdp = k(8);

            denominador(t) = (kp+kfM)*M(t) + (ktc+ktd)*RT(t);     % Factor común al radical
            alp(t) = denominador(t) / (kp*M(t));                  % Para poner 1 sobre kp*M(t)

        % Dar el valor inicial antes del bucle;

            % Diradicales

                %R0(1) = 2*ki1*I2p0(t)*M(t) / denominador(t); % Si hay iniciador monofuncional
                %R1(1) = 2*ki1*I2p1(t)*M(t) / denominador(t); % Si hay iniciador difuncional
                R2(1) = 2*ki1*I2p2(t)*M(t) / denominador(t);

            % Monoradicales

                %       r0(1) = (ki0 + 2*kfM*R0(t)) * M(t) / denominador(t);  % Si hay iniciador monofuncional
                r0(1) = (ki0 + kfM*RT(t)) * M(t) / denominador(t);
                %r1(1) = (ki1*Ip1(t) + ki0 + 2*kfM*R1(t)) * M(t) / denominador(t);  % Si hay iniciador difuncional
                r2(1) = (2*kfM*R2(t)) * M(t) / denominador(t);


            % Polímeros

                genP0(t,1) = kfM*r0(1)*M(t);
                genP1(t,1) = kfM*r1(1)*M(t);
                genP2(t,1) = kfM*r2(1)*M(t);
                genP3(t,1) = kfM*r3(1)*M(t);
                genP4(t,1) = kfM*r4(1)*M(t);
                genP5(t,1) = kfM*r5(1)*M(t);
                genP6(t,1) = kfM*r6(1)*M(t);
                genP7(t,1) = kfM*r7(1)*M(t);
                genP8(t,1) = kfM*r8(1)*M(t);

        % Generación de diradicales por terminación (Solo combinación)

            termcR0 = ktc / (kp*M(t)) * (  conv(R0, R0));
            termcR1 = ktc / (kp*M(t)) * (2*conv(R0, R1));
            termcR2 = ktc / (kp*M(t)) * (2*conv(R2, R0) +   conv(R1, R1));
            termcR3 = ktc / (kp*M(t)) * (2*conv(R3, R0) + 2*conv(R2, R1));
            termcR4 = ktc / (kp*M(t)) * (2*conv(R4, R0) + 2*conv(R3, R1) +   conv(R2, R2));
            termcR5 = ktc / (kp*M(t)) * (2*conv(R5, R0) + 2*conv(R4, R1) + 2*conv(R3, R2));
            termcR6 = ktc / (kp*M(t)) * (2*conv(R6, R0) + 2*conv(R5, R1) + 2*conv(R4, R2) +   conv(R3, R3));
            termcR7 = ktc / (kp*M(t)) * (2*conv(R7, R0) + 2*conv(R6, R1) + 2*conv(R5, R2) + 2*conv(R4, R3));
            termcR8 = ktc / (kp*M(t)) * (2*conv(R8, R0) + 2*conv(R7, R1) + 2*conv(R6, R2) + 2*conv(R5, R3) + conv(R4, R4));

        % Generación de monoradicales por terminación por combinación

            termcr0 = 2*ktc / (kp*M(t)) * (  conv(R0, r0));
            termcr1 = 2*ktc / (kp*M(t)) * (2*conv(R1, r0) + 2*conv(R0, r1));
            termcr2 = 2*ktc / (kp*M(t)) * (2*conv(R2, r0) +   conv(R1, r1) + 2*conv(R0, r2));
            termcr3 = 2*ktc / (kp*M(t)) * (2*conv(R3, r0) + 2*conv(R2, r1) + 2*conv(R1, r2) + 2*conv(R0, r3));
            termcr4 = 2*ktc / (kp*M(t)) * (2*conv(R4, r0) + 2*conv(R3, r1) +   conv(R2, r2) + 2*conv(R1, r3) + 2*conv(R0, r4));
            termcr5 = 2*ktc / (kp*M(t)) * (2*conv(R5, r0) + 2*conv(R4, r1) + 2*conv(R3, r2) + 2*conv(R2, r3) + 2*conv(R1, r4) + 2*conv(R0, r5));
            termcr6 = 2*ktc / (kp*M(t)) * (2*conv(R6, r0) + 2*conv(R5, r1) + 2*conv(R4, r2) +   conv(R3, r3) + 2*conv(R2, r4) + 2*conv(R1, r5) + 2*conv(R0, r6));
            termcr7 = 2*ktc / (kp*M(t)) * (2*conv(R7, r0) + 2*conv(R6, r1) + 2*conv(R5, r2) + 2*conv(R4, r3) + 2*conv(R3, r4) + 2*conv(R2, r5) + 2*conv(R1, r6) + 2*conv(R0, r7));
            termcr8 = 2*ktc / (kp*M(t)) * (2*conv(R8, r0) + 2*conv(R7, r1) + 2*conv(R6, r2) + 2*conv(R5, r3) +   conv(R4, r4) + 2*conv(R3, r5) + 2*conv(R2, r6) + 2*conv(R1, r7) + 2*conv(R0, r8));


        for n = 2:nmax                      % Bucle para sacar el n a partir del n-1 de radicales

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

                % Generación de monoradicales por terminación por desproporción (no se almacenan)

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

        % Generación de cadena de polímero por terminación por combinación

            termp0 = 0.5*ktc * (  conv(r0, r0));
            termp1 = 0.5*ktc * (2*conv(r1, r0));
            termp2 = 0.5*ktc * (2*conv(r2, r0) + conv(r1,r1));
            termp3 = 0.5*ktc * (2*conv(r3, r0) + 2*conv(r2,r1));
            termp4 = 0.5*ktc * (2*conv(r4, r0) + 2*conv(r3,r1) + conv(r2,r2));
            termp5 = 0.5*ktc * (2*conv(r5, r0) + 2*conv(r4,r1) + 2*conv(r3,r2));
            termp6 = 0.5*ktc * (2*conv(r6, r0) + 2*conv(r5,r1) + 2*conv(r4,r2) + conv(r3,r3));
            termp7 = 0.5*ktc * (2*conv(r7, r0) + 2*conv(r6,r1) + 2*conv(r5,r2) + 2*conv(r4,r3));
            termp8 = 0.5*ktc * (2*conv(r8, r0) + 2*conv(r7,r1) + 2*conv(r6,r2) + 2*conv(r5,r3) + conv(r4,r4));


        for n = 2:nmax

             % Término de GENERACIÓN de Polímero de longitud de cadena n en cada t

                genP0(t,n) = kfM*r0(n)*M(t) + termp0(n-1) + ktd*r0(n)*RT(t);
                genP1(t,n) = kfM*r1(n)*M(t) + termp1(n-1) + ktd*r1(n)*RT(t);
                genP2(t,n) = kfM*r2(n)*M(t) + termp2(n-1) + ktd*r2(n)*RT(t);
                genP3(t,n) = kfM*r3(n)*M(t) + termp3(n-1) + ktd*r3(n)*RT(t);
                genP4(t,n) = kfM*r4(n)*M(t) + termp4(n-1) + ktd*r4(n)*RT(t);
                genP5(t,n) = kfM*r5(n)*M(t) + termp5(n-1) + ktd*r5(n)*RT(t);
                genP6(t,n) = kfM*r6(n)*M(t) + termp6(n-1) + ktd*r6(n)*RT(t);
                genP7(t,n) = kfM*r7(n)*M(t) + termp7(n-1) + ktd*r7(n)*RT(t);
                genP8(t,n) = kfM*r8(n)*M(t) + termp8(n-1) + ktd*r8(n)*RT(t);

                % MOLES de polímero 
                NPS0(t+1,n) = NPS0(t,n) + genP0(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS1(t+1,n) = NPS1(t,n) + genP1(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS2(t+1,n) = NPS2(t,n) + genP2(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS3(t+1,n) = NPS3(t,n) + genP3(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS4(t+1,n) = NPS4(t,n) + genP4(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS5(t+1,n) = NPS5(t,n) + genP5(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS6(t+1,n) = NPS6(t,n) + genP6(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS7(t+1,n) = NPS7(t,n) + genP7(t,n)*V*(tiempo(t+1)-tiempo(t));
                NPS8(t+1,n) = NPS8(t,n) + genP8(t,n)*V*(tiempo(t+1)-tiempo(t));

        end
        
        for n = 2:nmax          % Cortes por ruptura secuencial
           
            y1 = corte(1,n);
            a = y1(1,1);
            b = y1(1,2) + 1;
            c = y1(2,1);
            d = y1(2,2) + 1;
            genr1(a,b) = genr1(a,b) + kdp*NPS1(t+1,n) / V;
            genr1(c,d) = genr1(c,d) + kdp*NPS1(t+1,n) / V;
            
            y2 = corte(2,n);
            a = y2(1,1);
            b = y2(1,2) + 1;
            c = y2(2,1);
            d = y2(2,2) + 1;
            genr2(a,b) = genr2(a,b) + kdp*NPS2(t+1,n) / V;
            genr2(c,d) = genr2(c,d) + kdp*NPS2(t+1,n) / V;
            
            y3 = corte(3,n);
            a = y3(1,1);
            b = y3(1,2) + 1;
            c = y3(2,1);
            d = y3(2,2) + 1;
            genr3(a,b) = genr3(a,b) + kdp*NPS3(t+1,n) / V;
            genr3(c,d) = genr3(c,d) + kdp*NPS3(t+1,n) / V;
            
            y4 = corte(4,n);
            a = y4(1,1);
            b = y4(1,2) + 1;
            c = y4(2,1);
            d = y4(2,2) + 1;
            genr4(a,b) = genr4(a,b) + kdp*NPS4(t+1,n) / V;
            genr4(c,d) = genr4(c,d) + kdp*NPS4(t+1,n) / V;
            
            y5 = corte(5,n);
            a = y5(1,1);
            b = y5(1,2) + 1;
            c = y5(2,1);
            d = y5(2,2) + 1;
            genr5(a,b) = genr5(a,b) + kdp*NPS5(t+1,n) / V;
            genr5(c,d) = genr5(c,d) + kdp*NPS5(t+1,n) / V;
            
            y6 = corte(6,n);
            a = y6(1,1);
            b = y6(1,2) + 1;
            c = y6(2,1);
            d = y6(2,2) + 1;
            genr6(a,b) = genr6(a,b) + kdp*NPS6(t+1,n) / V;
            genr6(c,d) = genr6(c,d) + kdp*NPS6(t+1,n) / V;
            
            y7 = corte(7,n);
            a = y7(1,1);
            b = y7(1,2) + 1;
            c = y7(2,1);
            d = y7(2,2) + 1;
            genr7(a,b) = genr7(a,b) + kdp*NPS7(t+1,n) / V;
            genr7(c,d) = genr7(c,d) + kdp*NPS7(t+1,n) / V;
            
            y8 = corte(8,n);
            a = y8(1,1);
            b = y8(1,2) + 1;
            c = y8(2,1);
            d = y8(2,2) + 1;
            genr8(a,b) = genr8(a,b) + kdp*NPS8(t+1,n) / V;
            genr8(c,d) = genr8(c,d) + kdp*NPS8(t+1,n) / V;
            
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
    
    [Mn, Mw] = pesos110(tiempo, NPS0, NPS1, NPS2, NPS3, NPS4, NPS5, NPS6, NPS7, NPS8);

end