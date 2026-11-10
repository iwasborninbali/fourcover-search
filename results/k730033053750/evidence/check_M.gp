\\ Independent PARI check of the point found on y^2 = x^3 + 730033053750.
default(nbthreads,1); default(realprecision,60);
E = ellinit([0,0,0,0,730033053750]);
Q = [-127864540025905752836434190966519096512662851319608332803792635443060325/2732848049858431036813881997543283958563840470801202716717782207848409, -122065928537596706107228249110895040861755243563538764427408281909699946209705677476565196173561467921763824675/142864138774358023251040652571655140848277152021135716683738993203046179270987130401371701590136836137827]; t = -652261919404834031/162343171503721777;
print("Q on E: ", ellisoncurve(E, Q)); print("torsion: ", elltors(E)[1]);
R = (26565)*t^0 + (-50820)*t^1 + (-41580)*t^2 + (-21000)*t^3 + (4620)*t^4; N = (345512475)*t^0 + (455023800)*t^1 + (-112058100)*t^2 + (190159200)*t^3 + (59579100)*t^4; T = (14823325613250)*t^0 + (11389585284000)*t^1 + (30768352807500)*t^2 + (-135488430000)*t^3 + (-12825049545000)*t^4 + (-1959590556000)*t^5 + (-1064848302000)*t^6;
print("R(t) is a square: ", issquare(R));
v = sqrtint(numerator(R)) / sqrtint(denominator(R)); xq = N / v^2;
print("pi(t, +-v) == Q: ", [xq, T/(2*v^3)] == Q || [xq, -T/(2*v^3)] == Q);
hQ = ellheight(E, Q); print("h(Q) = ", hQ);
ar = ellanalyticrank(E); B = ar[2] / ar[1]! / ellbsd(E); print("analytic rank ", ar[1], "; BSD Reg*Sha = ", B); print("h(Q) / (Reg*Sha) = ", hQ / B);
red(Ep, S, p, k) = if(denominator(S[1]) % p == 0, 0, ellmul(Ep, [Mod(S[1], p), Mod(S[2], p)], k) != [0]);
{
forprime(q = 2, 50, ok = 0;
  forprime(p = 5, 50000, if(ellglobalred(E)[1] % p == 0 || ok, next); G = ellgroup(E, p);
    if(#G == 1 && G[1] % q == 0, Ep = ellinit(E, p); if(red(Ep, Q, p, G[1] / q), ok = p)));
  print("Q not in ", q, "E(Q): ", if(ok, Str("witness p=", ok), "undecided")));
}
