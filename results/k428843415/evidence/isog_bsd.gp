\\ Regulator prediction for y^2 = x^3 + 428843415 and its 3-isogenous partner (L-values agree across an isogeny).
default(nbthreads, 1);
E = ellinit([0,0,0,0,428843415]);
IS = ellisomat(E); L = IS[1];
{
for(i = 1, #L,
  my(F = ellinit(L[i][1]), M = ellminimalmodel(F), Lr = ellanalyticrank(M), b = ellbsd(M), tam = elltamagawa(M), t = elltors(M)[1]);
  print("curve ", i, ": ainvs ", M[1..5], " conductor ", ellglobalred(M)[1], " torsion ", t,
        " analytic_rank ", Lr[1], " L^(r)(1)/r! ", Lr[2]/Lr[1]!, " tamagawa ", tam,
        " Reg*Sha ", Lr[2]/Lr[1]!/b));
}
print("isogeny matrix ", IS[2]);
