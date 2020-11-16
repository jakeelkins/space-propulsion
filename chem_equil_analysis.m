% solve sys of eq ofr chem equil

syms x y z w

eq1 = x + y;

eq2 = 2*x + y + 2*z + w;

eq3 = ((y*z^(0.5))/x)*2/(x+y+z+w);

eq4 = ((w^(2))/z)*2/(x+y+z+w);

[solx, soly, solz, solw] = vpasolve(eq1==2,eq2==6,eq3==0.65116,eq4==0.0463);
