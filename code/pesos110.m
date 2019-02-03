function [Mn, Mw] = pesos110(t, NPS0, NPS1, NPS2, NPS3, NPS4, NPS5, NPS6, NPS7, NPS8)

    global PM
    global nmax
        n = 1:nmax;

    GPS0(1:length(t),1:nmax)=0;
    GPS1(1:length(t),1:nmax)=0;
    GPS2(1:length(t),1:nmax)=0;
    GPS3(1:length(t),1:nmax)=0;
    GPS4(1:length(t),1:nmax)=0;
    GPS5(1:length(t),1:nmax)=0;
    GPS6(1:length(t),1:nmax)=0;
    GPS7(1:length(t),1:nmax)=0;
    GPS8(1:length(t),1:nmax)=0;

    NPS(1:length(t))=0;
    GPS(1:length(t))=0;
    Mn(1:length(t))=0;
    Mw(1:length(t))=0;
    D(1:length(t))=0;
    PeP2(1:length(t))=0;


    for i=2:length(t)

        for j=1:nmax

            GPS0(i,j)=j*PM*NPS0(i,j);
            GPS1(i,j)=j*PM*NPS1(i,j);
            GPS2(i,j)=j*PM*NPS2(i,j);
            GPS3(i,j)=j*PM*NPS3(i,j);
            GPS4(i,j)=j*PM*NPS4(i,j);
            GPS5(i,j)=j*PM*NPS5(i,j);
            GPS6(i,j)=j*PM*NPS6(i,j);
            GPS7(i,j)=j*PM*NPS7(i,j);
            GPS8(i,j)=j*PM*NPS8(i,j);

        end

        NPS(i)=sum(NPS0(i,:)+NPS1(i,:)+NPS2(i,:)+NPS3(i,:)+NPS4(i,:)+NPS5(i,:)+NPS6(i,:)+NPS7(i,:)+NPS8(i,:));
        GPS(i)=sum(GPS0(i,:)+GPS1(i,:)+GPS2(i,:)+GPS3(i,:)+GPS4(i,:)+GPS5(i,:)+GPS6(i,:)+GPS7(i,:)+GPS8(i,:));

        Mn(i)=GPS(i)./NPS(i);
        Mw(i)=PM*sum((GPS0(i,:)+GPS1(i,:)+GPS2(i,:)+GPS3(i,:)+GPS4(i,:)+GPS5(i,:)+GPS6(i,:)+GPS7(i,:)+GPS8(i,:)).*n)./GPS(i);

        D(i)=Mw(i)./Mn(i);
        PeP2(i)=sum(NPS1(i,:)+2*NPS2(i,:)+3*NPS3(i,:)+4*NPS4(i,:)+5*NPS5(i,:)+6*NPS6(i,:)+7*NPS7(i,:)+8*NPS8(i,:));

    end
    
    
end