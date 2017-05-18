function pintarConversiones(t, x120, x130, r120, r130)

    global tiempos120
    global convers120
    global tiempos130
    global convers130
    
    subplot(4,2,1:6)
    hold on
        plot(tiempos120,convers120,'b--')
        plot(tiempos130,convers130,'r--')
        plot(t/60, x120*100, 'b')
        plot(t/60, x130*100, 'r')
    title('Evolución de la Conversión')
    xlabel('Tiempo en el reactor (min)')
    ylabel('Conversión (%)')
    ylim([0 100]);
    legend('Experimental 120ºC', 'Experimental 130ºC', 'Modelo 120ºC', 'Modelo 130ºC')
    
    hold off
    
    subplot(4,2,7:8)
    hold on
        plot(t/60, r120, 'b')
        plot(t/60, r130, 'r')
    title('Evolución de la Rp')
    xlabel('Tiempo en el reactor (min)')
    ylabel('Rp')
    legend('120 ºC', '130 ºC')


    
end    