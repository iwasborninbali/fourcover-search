\\ Rank upper bound for y^2 = x^3 + 730033053750 with certified cubic-field data.
default(nbthreads,1);
E = ellinit([0,0,0,0,730033053750]); rk = ellrankinit(E);
{ for(i = 1, #rk[3], print("field ", rk[3][i].pol, " disc ", rk[3][i].disc, " class group ", rk[3][i].clgp.cyc, " bnfcertify=", bnfcertify(rk[3][i]))); }
print("rank bounds from certified structure: ", ellrank(rk, 0)[1..3]);
