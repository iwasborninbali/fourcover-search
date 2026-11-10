\\ Export the 2-cover basis of an elliptic curve E/Q, with exact maps to E, as JSON (input for search4.py).
\\ Usage:  echo 'A = [a1,a2,a3,a4,a6]; NAME = "myname";' | cat - pipeline/export_world.gp | gp -q -f > world.json
\\ Each cover is v^2 = R(t) (ascending coefficients); a point (t, v) maps to E by x = N(t)/v^2, y = (-(a1 x + a3) + T(t)/v^3)/2.
\\ The map formula is checked at random points before it is written.
default(nbthreads, 1);
E = ellinit(A); a1 = E.a1; a3 = E.a3;
print1("{\"world\": \"", NAME, "\", \"ainvs\": ", A, ", \"conductor\": ", ellglobalred(E)[1], ", \"known_point\": null, ");
print1("\"torsion\": ", elltors(E)[1], ", \"coefficient_order\": \"ascending\", \"map_to_original_E\": [\"N(t)/v^2\", \"(-(a1*x + a3) + T(t)/v^3)/2\"], \"covers\": [");
C = ell2cover(E);
{
for(i = 1, #C,
  R = C[i][1]; M = C[i][2];
  Nn = subst(M[1] * y^2, y, 1); Tt = subst(simplify((2 * M[2] + a1 * M[1] + a3) * y^3), y, 1);
  if(type(Tt) != "t_POL" || poldegree(Tt) > 6, error("T is not a polynomial"));
  for(j = 0, 10, my(t0 = random(100) - 50, v0 = 1 + random(9));
    my(xE = subst(Nn, x, t0) / v0^2, yE = (-(a1 * xE + a3) + subst(Tt, x, t0) / v0^3) / 2);
    if([xE, yE] != subst(subst(M, x, t0), y, v0), error("map formula mismatch")));
  print1(if(i > 1, ", ", ""), "{\"label\": \"c", i, "\", \"R\": ", Vecrev(R), ", \"N\": ", Vecrev(Nn), ", \"T\": ", Vecrev(Tt), "}"));
}
print("]}");
