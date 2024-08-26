
NPV = 0.27;
PPV = 0.96;
sen=0.64;

x1=(NPV)/(1-NPV);
x2=sen*(1-PPV)/PPV/(1-sen);
x3=NPV/(1-NPV);

sp = x1/(x2+x3);
disp(sp);