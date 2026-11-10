\\ Independent PARI check of the point found on y^2 = x^3 + 428843415.
default(nbthreads,1); default(realprecision,60);
E = ellinit([0,0,0,0,428843415]);
Q = [-21749823817314505417329398578203428726044389786520594015230853682787114924753276288107768790071/35589106081171735697954614495653836664456516348375009203368133815059732677049884911719472225, -3006983147937359986718680977260908084857042951090561334004195430921637930842654352957809370917734940721827109177778158814798527002308734122092/212312527003382662900836631130652382841261077438584845825565459225013474360190212744967783408671947123118390643672707622161170983175995375]; t = 149103631471577628622190/87526920682404705781569;
print("Q on E: ", ellisoncurve(E, Q)); print("torsion: ", elltors(E)[1]);
R = (14105)*t^0 + (-4108)*t^1 + (1872)*t^2 + (1144)*t^3 + (-104)*t^4; N = (-3346031)*t^0 + (-8708908)*t^1 + (2346396)*t^2 + (-35152)*t^3 + (114244)*t^4; T = (68292321344)*t^0 + (-39359685612)*t^1 + (-8330101260)*t^2 + (-12634243960)*t^3 + (3437601960)*t^4 + (-457889952)*t^5 + (-63519664)*t^6;
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
