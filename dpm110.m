function dpm110

    global tiempo
    global T1
    global Rjul
    global B
    global C

            % Quitar cuando se acople a conversion
             
                to = 0;         % s
                tf = 5400;      % s 
                tiempo  = linspace(to, tf, 100000);    

                B = -4;
                C = -5;
    
    nmax = 12000;   % M�xima longitud de cadena
    
    denominador(1:length(tiempo)) = 0;      % Inicializar denominador con el tama�o del vector tiempo
    alp(1:length(tiempo)) = 0;              % Inicializar alp con el tama�os del vector tiempo

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

           %rn = Diradical con n grupos per�xidos
            r0 = zeros(1, nmax);
            r1 = zeros(1, nmax);
            r2 = zeros(1, nmax);
            r3 = zeros(1, nmax);
            r4 = zeros(1, nmax);
            r5 = zeros(1, nmax);
            r6 = zeros(1, nmax);
            r7 = zeros(1, nmax);
            r8 = zeros(1, nmax);

        % Pol�meros --> Inicializamos a 0

           %pn = Pol�mero con n grupos per�xidos
            p0 = zeros(length(tiempo), nmax);
            p1 = zeros(length(tiempo), nmax);
            p2 = zeros(length(tiempo), nmax);
            p3 = zeros(length(tiempo), nmax);
            p4 = zeros(length(tiempo), nmax);
            p5 = zeros(length(tiempo), nmax);
            p6 = zeros(length(tiempo), nmax);
            p7 = zeros(length(tiempo), nmax);
            p8 = zeros(length(tiempo), nmax);
            
           %N grupos per�xidos sin descomponer
            NPS0 = zeros(length(tiempo), nmax);
            NPS1 = zeros(length(tiempo), nmax);
            NPS2 = zeros(length(tiempo), nmax);
            NPS3 = zeros(length(tiempo), nmax);
            NPS4 = zeros(length(tiempo), nmax);
            NPS5 = zeros(length(tiempo), nmax);
            NPS6 = zeros(length(tiempo), nmax);
            NPS7 = zeros(length(tiempo), nmax);
            NPS8 = zeros(length(tiempo), nmax);
            
        % Generaci�n de monoradicales por desproporci�n
        
            termdr0(n) = zeros(1, nmax-1);
            termdr1(n) = zeros(1, nmax-1);
            termdr2(n) = zeros(1, nmax-1);
            termdr3(n) = zeros(1, nmax-1);
            termdr4(n) = zeros(1, nmax-1);
            termdr5(n) = zeros(1, nmax-1);
            termdr6(n) = zeros(1, nmax-1);
            termdr7(n) = zeros(1, nmax-1);
            termdr8(n) = zeros(1, nmax-1);

    % C�lculo

    for t = 1:length(tiempo)-1

        % Constantes cin�ticas

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

            denominador(t) = (kp+kfM)*M(t) + (ktc+ktd)*(R1(t) + 2*R2(t));   % Factor com�n al radical
            alp(t) = denominador(t) / (kp*M(t));                            % Para poner 1 sobre kp*M(t)

        % Dar el valor inicial antes del bucle;

            % Diradicales

                R0(1) = 2*ki1*I2p0(t)*M(t) / denominador(t);
                R1(1) = 2*ki1*I2p1(t)*M(t) / denominador(t);
                R2(1) = 2*ki1*I2p2(t)*M(t) / denominador(t);

            % Monoradicales

                r0(1) = (ki1*Ip0(t) + ki0 + 2*kfM*R0(t)) * M(t) / denominador(t);
                r1(1) = (ki1*Ip1(t) + ki0 + 2*kfM*R1(t)) * M(t) / denominador(t);
                r2(1) = (ki1*Ip2(t) + ki0 + 2*kfM*R2(t)) * M(t) / denominador(t);

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
              
            termR0 = ktc / (kp*M(t)) * (  conv(R0, R0));
            termR1 = ktc / (kp*M(t)) * (2*conv(R0, R1));
            termR2 = ktc / (kp*M(t)) * (2*conv(R2, R0) +   conv(R1, R1));
            termR3 = ktc / (kp*M(t)) * (2*conv(R3, R0) + 2*conv(R2, R1));
            termR4 = ktc / (kp*M(t)) * (2*conv(R4, R0) + 2*conv(R3, R1) +   conv(R2, R2));
            termR5 = ktc / (kp*M(t)) * (2*conv(R5, R0) + 2*conv(R4, R1) + 2*conv(R3, R2));
            termR6 = ktc / (kp*M(t)) * (2*conv(R6, R0) + 2*conv(R5, R1) + 2*conv(R4, R2) +   conv(R3, R3));
            termR7 = ktc / (kp*M(t)) * (2*conv(R7, R0) + 2*conv(R6, R1) + 2*conv(R5, R2) + 2*conv(R4, R3));
            termR8 = ktc / (kp*M(t)) * (2*conv(R8, R0) + 2*conv(R7, R1) + 2*conv(R6, R2) + 2*conv(R5, R3) + conv(R4, R4));
                
        % Generaci�n de monoradicales por terminaci�n por combinaci�n
              
            termr0 = 2*ktc / (kp*M(t)) * (conv(R0, r0));
            termr1 = 2*ktc / (kp*M(t)) * (conv(R1, r0) + conv(R0, r1));
            termr2 = 2*ktc / (kp*M(t)) * (conv(R2, r0) + conv(R1, r1) + conv(R0, r2));
            termr3 = 2*ktc / (kp*M(t)) * (conv(R3, r0) + conv(R2, r1) + conv(R1, r2) + conv(R0, r3));
            termr4 = 2*ktc / (kp*M(t)) * (conv(R4, r0) + conv(R3, r1) + conv(R2, r3) + conv(R1, r3) + conv(R0, r4));
            termr5 = 2*ktc / (kp*M(t)) * (conv(R5, r0) + conv(R4, r1) + conv(R3, r4) + conv(R2, r3) + conv(R1, r4) + conv(R0, r5));
            termr6 = 2*ktc / (kp*M(t)) * (conv(R6, r0) + conv(R5, r1) + conv(R4, r5) + conv(R3, r3) + conv(R2, r4) + conv(R1, r5) + conv(R0, r6));
            termr7 = 2*ktc / (kp*M(t)) * (conv(R7, r0) + conv(R6, r1) + conv(R5, r6) + conv(R4, r3) + conv(R3, r4) + conv(R2, r5) + conv(R1, r6) + conv(R0, r7));
            termr8 = 2*ktc / (kp*M(t)) * (conv(R8, r0) + conv(R7, r1) + conv(R6, r7) + conv(R5, r3) + conv(R4, r4) + conv(R3, r5) + conv(R2, r6) + conv(R1, r7) + conv(R0, r8));
                
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
            
        % Generaci�n de monoradicales ��??
        
%             genr1 = zeros(nmax-1, 1:1);
%             genr2 = zeros(nmax-1, 1:2);
%             genr3 = zeros(nmax-1, 1:3);
%             genr4 = zeros(nmax-1, 1:4);
%             genr5 = zeros(nmax-1, 1:5);
%             genr6 = zeros(nmax-1, 1:6);
%             genr7 = zeros(nmax-1, 1:7);
%             genr8 = zeros(nmax-1, 1:8);
            
                
        % Bucle para sacar el n a partir del n-1 de diradicales
            
            for n = 2:nmax
                
                % Diradicales en cada t (no se almacenan)

                    R0(n) = (R0(n-1) + termR0(n-1)) / alp(n);
                    R1(n) = (R1(n-1) + termR0(n-1)) / alp(n);
                    R2(n) = (R2(n-1) + termR0(n-1)) / alp(n);
                    R3(n) = (R3(n-1) + termR0(n-1)) / alp(n);
                    R4(n) = (R4(n-1) + termR0(n-1)) / alp(n);
                    R5(n) = (R5(n-1) + termR0(n-1)) / alp(n);
                    R6(n) = (R6(n-1) + termR0(n-1)) / alp(n);
                    R7(n) = (R7(n-1) + termR0(n-1)) / alp(n);
                    R8(n) = (R8(n-1) + termR0(n-1)) / alp(n);

                    % Generaci�n de monoradicales por terminaci�n por desproporci�n (no se almacenan)

                        termdr0(n) = 2*ktd / (kp*M(t)) * R0(n)*R(t);
                        termdr1(n) = 2*ktd / (kp*M(t)) * R1(n)*R(t);
                        termdr2(n) = 2*ktd / (kp*M(t)) * R2(n)*R(t);
                        termdr3(n) = 2*ktd / (kp*M(t)) * R3(n)*R(t);
                        termdr4(n) = 2*ktd / (kp*M(t)) * R4(n)*R(t);
                        termdr5(n) = 2*ktd / (kp*M(t)) * R5(n)*R(t);
                        termdr6(n) = 2*ktd / (kp*M(t)) * R6(n)*R(t);
                        termdr7(n) = 2*ktd / (kp*M(t)) * R7(n)*R(t);
                        termdr8(n) = 2*ktd / (kp*M(t)) * R8(n)*R(t);
                    
                % Monoradicales en cada t (no se almacenan)

                    r0(n) = (r0(n-1) + termdr(n)/kp*M(t) + 2*kfM/kp*R0(n) + termr0(n-1)) / alp(t);
                    r1(n) = (r1(n-1) + termdr(n)/kp*M(t) + 2*kfM/kp*R1(n) + termr1(n-1)) / alp(t);
                    r2(n) = (r2(n-1) + termdr(n)/kp*M(t) + 2*kfM/kp*R2(n) + termr2(n-1)) / alp(t);
                    r3(n) = (r3(n-1) + termdr(n)/kp*M(t) + 2*kfM/kp*R3(n) + termr3(n-1)) / alp(t);
                    r4(n) = (r4(n-1) + termdr(n)/kp*M(t) + 2*kfM/kp*R4(n) + termr4(n-1)) / alp(t);
                    r5(n) = (r5(n-1) + termdr(n)/kp*M(t) + 2*kfM/kp*R5(n) + termr5(n-1)) / alp(t);
                    r6(n) = (r6(n-1) + termdr(n)/kp*M(t) + 2*kfM/kp*R6(n) + termr6(n-1)) / alp(t);
                    r7(n) = (r7(n-1) + termdr(n)/kp*M(t) + 2*kfM/kp*R7(n) + termr7(n-1)) / alp(t);
                    r8(n) = (r8(n-1) + termdr(n)/kp*M(t) + 2*kfM/kp*R8(n) + termr8(n-1)) / alp(t);           

                 % Pol�mero en cada t
                 
                    p0(t,n) = kfM*r0(n)*M(t) + termp0(t-1);
                    p1(t,n) = kfM*r1(n)*M(t) + termp1(t-1);
                    p2(t,n) = kfM*r2(n)*M(t) + termp2(t-1);
                    p3(t,n) = kfM*r3(n)*M(t) + termp3(t-1);
                    p4(t,n) = kfM*r4(n)*M(t) + termp4(t-1);
                    p5(t,n) = kfM*r5(n)*M(t) + termp5(t-1);
                    p6(t,n) = kfM*r6(n)*M(t) + termp6(t-1);
                    p7(t,n) = kfM*r7(n)*M(t) + termp7(t-1);
                    p8(t,n) = kfM*r8(n)*M(t) + termp8(t-1);
                 
                    NPS0(n+1,1) = NPS0(n,1) + p0(t+1)-t(t);
                    NPS1(n+1,1) = NPS1(n,1) + p1(t+1)-t(t);
                    NPS2(n+1,1) = NPS2(n,1) + p2(t+1)-t(t);
                    NPS3(n+1,1) = NPS3(n,1) + p3(t+1)-t(t);
                    NPS4(n+1,1) = NPS4(n,1) + p4(t+1)-t(t);
                    NPS5(n+1,1) = NPS5(n,1) + p5(t+1)-t(t);
                    NPS6(n+1,1) = NPS6(n,1) + p6(t+1)-t(t);
                    NPS7(n+1,1) = NPS7(n,1) + p7(t+1)-t(t);
                    NPS8(n+1,1) = NPS8(n,1) + p8(t+1)-t(t);
                
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

end
