p := 2;
q := 3;
r := 7;
mu1 := 1/2*(1-1/p-1/q+1/r);
mu2 := 1/2*(1-1/p+1/q-1/r);
mu3 := 1/2*(1+1/p-1/q-1/r);
mu4 := 1/2*(1+1/p+1/q+1/r);
m := 2*LCM([p,q,r]);
[-1+&+[mu*a-Floor(mu*a): mu in [mu1,mu2,mu3,mu4]]: a in [1..m]|GCD(a, m) eq 1];
