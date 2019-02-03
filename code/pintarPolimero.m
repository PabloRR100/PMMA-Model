function pintarPolimero(t, p1, p2, p3, PeP1, PeP2, PeP3)

    figure
    hold on
        plot(t/60, p1, 'g')
        plot(t/60, p2, 'b')
        plot(t/60, p3, 'r')
        plot(t/60, PeP1, 'g--')
        plot(t/60, PeP2, 'b--')
        plot(t/60, PeP3, 'r--')
    title('Evoluci�n de la Concentraci�n')
    xlabel('Tiempo en el reactor (min)')
    ylabel('Concentraci�n (g/L)')
    legend('Pol�mero 110�C', 'Pol�mero 120�C', 'Pol�mero 130�C', 'Pol�mero con grupos per�xidos 110�C', 'Pol�mero con grupos per�xidos 120�C', 'Pol�mero con grupos per�xidos 130�C')
    legend('Location','east')

end