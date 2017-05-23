function pintarPolimero(t, p1, p2, p3)

    figure
    hold on
        plot(t/60, p1, 'g')
        plot(t/60, p2, 'b')
        plot(t/60, p3, 'r')
    title('Evolución de la Concentración')
    xlabel('Tiempo en el reactor (min)')
    ylabel('Concentración (g/L)')
    legend('110ºC', '120ºC', '130ºC')
    legend('Location','east')

end