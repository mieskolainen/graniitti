clear; close all;

a = sym('a',[4 1]);
b = sym('b',[4 1]);

c = sym('c',[4 1]);
d = sym('d',[4 1]);

syms mA mB mC mD alpha1 alpha2 E

eq1 = alpha1*c(1) + alpha2*d(1) == E

eq3 = mC^2 == c(2)^2 + c(3)^2 + c(4)^2
eq4 = mD^2 == d(2)^2 + d(3)^2 + d(4)^2

