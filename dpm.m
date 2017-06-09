
global T1
global T2
global T3

global nmax;
    nmax = 30000;

[Mn1, Mw1, X1] = dpm110(t1, x1, T1);
[Mn2, Mw2, X2] = dpm120(t2, x2, T2);
[Mn3, Mw3, X3] = dpm130(t3, x3, T3);

pintarPesos(t1, X1, Mn1, Mw1, t2, X2, Mn2, Mw2, t3, X3, Mn3, Mw3)
    
