function pintarIniciador(t, i1, i2)

    global tInic120
    global iniciador120
    global tInic130
    global iniciador130
    
    hold on 
    
    plot(tInic120, iniciador120, 'b--')
    plot(tInic130, iniciador130, 'r--')
    plot(t/60, i1, 'r')
    plot(t/60, i2, 'b')
    
    title('Evolución del Iniciador')
    xlabel('Tiempo en el reactor (min)')
    xlim([0 max(max(tInic120),max(tInic130))])
    ylabel('Concentración Iniciador (M)')
    ylim([0 iniciador120(1)]);
    legend('Experimental 120ºC', 'Experimental 130ºC', 'Modelo 120', 'Modelo130')


    
end    