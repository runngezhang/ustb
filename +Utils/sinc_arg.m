function x = sinc_arg(y)
x = fminbnd(@(arg) abs(sinc(arg)-y),0,1);