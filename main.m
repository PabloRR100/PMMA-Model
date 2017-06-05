clear variables
clear globals
        
% DATOS EMPIRICOS

    global tiempos110
    global convers110

    global tiempos120
    global convers120
    
    global tiempos130
    global convers130
    
%     global tInic120
%     global iniciador120
%     
%     global tInic130
%     global iniciador130
    
    datosMejico = '/Users/pablorr10/OneDrive/UNIVERSIDAD/Ingeniería Química/MIQ/PROYECTO/Matlab/data/Resultados de Pablo.xlsx';
    
    	[tiempos110, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'B8:B17');
        [convers110, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'D8:D17');
        [tiempos120, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'B25:B34');
        [convers120, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'D25:D34');
        [tiempos130, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'B42:B51');
        [convers130, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'D42:D51');

%         [tInic120, ~, ~]     = xlsread(datosMejico, 'DescomposicionDEKTP', 'B32:B35'); 
%         [iniciador120, ~, ~] = xlsread(datosMejico, 'DescomposicionDEKTP', 'D32:D35'); 
%         [tInic130, ~, ~]     = xlsread(datosMejico, 'DescomposicionDEKTP', 'B21:B25'); 
%         [iniciador130, ~, ~] = xlsread(datosMejico, 'DescomposicionDEKTP', 'D21:D25'); 
    
        
        
% DEFINICION DE LOS PARAMETROS DEL SISTEMA

    global Rjul
    global Rcal
        Rjul = 8.314;       % J/molÂ·K
        Rcal = 1.987;       % cal/molÂ·K

    global ef
        ef = 0.75;
        
%     global B
%     global C


% DEFINICION DE LAS VARIABLES DEL SISTEMA
    
    global T1               % 110 ºC
    global T2               % 120 ºC
    global T3               % 130 ºC
        T1 = 110 + 273;     % K
        T2 = 120 + 273;     % K
        T3 = 130 + 273;     % K
     
    global PM
    global p

        PM = 100;   % g/mol
        p  = 940;  % g/L

    global Mo
    global Vo
    global I3o
    
        I3o = 0.01;     % mol/L
        %I3o = 0;     % Para las pruebas sin iniciador
        Mo  = p / PM;   % mol/L ( densidad / PM) 
        Vo  = 0.9;      % L

    global B
    global C
    
        B = -4;
        C = -5;
        
% FUNCION PARA CALCULAR

    conversion      
