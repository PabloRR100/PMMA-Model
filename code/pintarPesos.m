function pintarPesos(t1, x1, Mn1, Mw1, t2, x2, Mn2, Mw2, t3, x3, Mn3, Mw3)

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
    
    pantalla = get(groot, 'ScreenSize');
    
    figure('Position',[0 pantalla(4)/2 pantalla(3) pantalla(4)/2])
    hold on
        plot(x1*100, Mn1, 'g')
        plot(x1*100, Mw1, 'g--')
        plot(convers110, MN110, 'gs')
        plot(convers110, MW110, 'g*')
        plot(x2*100, Mn2, 'b')
        plot(x2*100, Mw2, 'b--')
        plot(convers120, MN120, 'bs')
        plot(convers120, MW120, 'b*')
        plot(x3*100, Mn3, 'r')
        plot(x3*100, Mw3, 'r--')
        plot(convers130, MN130, 'rs')
        plot(convers130, MW130, 'r*')
            title('Dpm con la conversi�n')
            xlabel('Conversi�n (%)')
            ylabel('M (g/mol)')
            xlim([0, 100])
            ylim([0.5*10^5, 12*10^5])
            legend('Mn 110�C Modelo', 'Mw 110�C Modelo', 'Mn 110�C Experimental', 'Mw 110�C Experimental', 'Mn 120�C Modelo', 'Mw 120�C Modelo', 'Mn 120�C Experimental', 'Mw 120�C Experimental', 'Mn 130�C Modelo', 'Mw 130�C Modelo', 'Mn 130�C Experimental', 'Mw 130�C Experimental')
        
    hold off
    
    figure('Position',[0 0 pantalla(3) pantalla(4)/2])
    hold on
        plot(t1/60, Mn1, 'g')
        plot(t1/60, Mw1, 'g--')
        plot(tiempos110, MN110, 'gs')
        plot(tiempos110, MW110, 'g*')
        plot(t2/60, Mn2, 'b')
        plot(t2/60, Mw2, 'b--')
        plot(tiempos120, MN120, 'bs')
        plot(tiempos120, MW120, 'b*')
        plot(t3/60, Mn3, 'r')
        plot(t3/60, Mw3, 'r--')
        plot(tiempos130, MN130, 'rs')
        plot(tiempos130, MW130, 'r*')
            title('Dpm en el tiempo')
            xlabel('Tiempo (min)')
            ylabel('M (g/mol)')
            xlim([5, 90])
            ylim([0.5*10^5, 12*10^5])
            legend('Mn 110�C Modelo', 'Mw 110�C Modelo', 'Mn 110�C Experimental', 'Mw 110�C Experimental', 'Mn 120�C Modelo', 'Mw 120�C Modelo', 'Mn 120�C Experimental', 'Mw 120�C Experimental', 'Mn 130�C Modelo', 'Mw 130�C Modelo', 'Mn 130�C Experimental', 'Mw 130�C Experimental')
        
    hold off
    
end