function pintarConversiones(t1, t2, t3, x110, x120, x130, r110, r120, r130)

    global tiempos110
    global convers110
    global tiempos120
    global convers120
    global tiempos130
    global convers130
    
    pantalla = get(groot, 'ScreenSize');
    figure('Position',[0 0 pantalla(3)/2 pantalla(4)])
    hold on
        plot(tiempos110,convers110,'g*')
        plot(tiempos120,convers120,'b*')
        plot(tiempos130,convers130,'r*')
        plot(t1/60, x110*100, 'g')
        plot(t2/60, x120*100, 'b')
        plot(t3/60, x130*100, 'r')
            title('Evolución de la Conversión')
            xlabel('Tiempo en el reactor (min)')
            ylabel('Conversión (%)')
            ylim([0 100]);
            legend('Experimental 110ºC', 'Experimental 120ºC', 'Experimental 130ºC', 'Modelo 110ºC', 'Modelo 120ºC', 'Modelo 130ºC')
            legend('Location','east')
    
    hold off
    
    pantalla = get(groot, 'ScreenSize');
    figure('Position',[pantalla(3)/2 pantalla(4)/2 pantalla(3)/2 pantalla(4)/2])
    hold on
        plot(x110*100, r110, 'g')
        plot(x120*100, r120, 'b')
        plot(x130*100, r130, 'r')
    title('Evolución de la Rp')
    xlabel('Conversión (%)')
    ylabel('Rp')
    legend('110 ºC', '120 ºC', '130 ºC')
    legend('Location','northwest')
    
    hold off
    
    
    pantalla = get(groot, 'ScreenSize');
    figure('Position',[pantalla(3)/2 0 pantalla(3)/2 pantalla(4)/2])
    hold on
        plot(t/60, r110, 'g')
        plot(t/60, r120, 'b')
        plot(t/60, r130, 'r')
    title('Evolución de la Rp')
    xlabel('Tiempo en el reactor (min)')
    ylabel('Rp')
    legend('110 ºC', '120 ºC', '130 ºC')
    legend('Location','east')

end    