
clear variables
clear globals
        
% DATOS EMPIRICOS

    global tiempos110
    global convers110
    global MN110
    global MW110

    global tiempos120
    global convers120
    global MN120
    global MW120
    
    global tiempos130
    global convers130
    global MN130
    global MW130
    
    
    datosMejico = '/Users/pablorr10/OneDrive/UNIVERSIDAD/Ingeniería Química/MIQ/PROYECTO/Matlab/data/Resultados de Pablo.xlsx';
    
    	[tiempos110, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'B8:B17');
        [convers110, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'D8:D17');
        [MN110, ~, ~]      = xlsread(datosMejico, 'PolimerizacionDEKTP', 'E8:E17');
            MN110 = MN110*10^5;                                                         % En el excel están por 10^-5
        [MW110, ~, ~]      = xlsread(datosMejico, 'PolimerizacionDEKTP', 'F8:F17');
            MW110 = MW110*10^5;
        
        
        [tiempos120, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'B25:B34');
        [convers120, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'D25:D34');
        [MN120, ~, ~]      = xlsread(datosMejico, 'PolimerizacionDEKTP', 'E25:E34');
            MN120 = MN120*10^5;
        [MW120, ~, ~]      = xlsread(datosMejico, 'PolimerizacionDEKTP', 'F25:F34');
            MW120 = MW120*10^5;
        
        [tiempos130, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'B42:B51');
        [convers130, ~, ~] = xlsread(datosMejico, 'PolimerizacionDEKTP', 'D42:D51');
        [MN130, ~, ~]      = xlsread(datosMejico, 'PolimerizacionDEKTP', 'E42:E51');
            MN130 = MN130*10^5;
        [MW130, ~, ~]      = xlsread(datosMejico, 'PolimerizacionDEKTP', 'F42:F51');
            MW130 = MW130*10^5;
        
        
% DEFINICION DE LOS PARAMETROS DEL SISTEMA

    global Rjul
    global Rcal
        Rjul = 8.314;       % J/molÂ·K
        Rcal = 1.987;       % cal/molÂ·K

    global ef
        ef = 0.75;


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
        