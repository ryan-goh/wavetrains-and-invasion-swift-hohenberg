x = [-200:0.01:100];
nu = -2.8765 - 0.3882i;
ep = -0.2;
ep1 = -0.2;

unu = real((x+1).*exp(min(nu*x,200)));
    EXp   =  exp(min((x-0)/ep1,500));
    CHIp  =  1./(1+EXp);
    CHIpm = 1-CHIp;
    EX   =  exp(min(-(x-(-20))/ep,500));
    %  EX   =  exp(-(X-X0)/eps);
    CHI  =  1./(1+EX);

    figure(4)
    plot(x,[CHI',CHIp'])

    figure(1)
    plot(x,CHIp.*unu)