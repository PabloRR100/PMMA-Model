function pintarPolimero(t, p1, p2, p3)

    figure
    hold on
        plot(t/60, p1, 'g')
        plot(t/60, p2, 'b')
        plot(t/60, p3, 'r')
    title('Evoluci�n de la Concentraci�n')
    xlabel('Tiempo en el reactor (min)')
    ylabel('Concentraci�n (g/L)')
    legend('110�C', '120�C', '130�C')
    legend('Location','east')

end