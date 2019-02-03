function pintarPolimero(t, p1, p2, p3, PeP1, PeP2, PeP3)

    figure
    hold on
        plot(t/60, p1, 'g')
        plot(t/60, p2, 'b')
        plot(t/60, p3, 'r')
        plot(t/60, PeP1, 'g--')
        plot(t/60, PeP2, 'b--')
        plot(t/60, PeP3, 'r--')
    title('Evolución de la Concentración')
    xlabel('Tiempo en el reactor (min)')
    ylabel('Concentración (g/L)')
    legend('Polímero 110ºC', 'Polímero 120ºC', 'Polímero 130ºC', 'Polímero con grupos peróxidos 110ºC', 'Polímero con grupos peróxidos 120ºC', 'Polímero con grupos peróxidos 130ºC')
    legend('Location','east')

end