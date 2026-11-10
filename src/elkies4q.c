/* Elkies-style lattice search for rational points on a four-cover Q2 = Q3 = 0 in P^3, second version.
   Same contract as elkies4.c (input, POINT/DENSE/COMPONENT/DONE output, exact int128 check, completeness NOT certified),
   with two changes that matter for large H:
   1. Osculating box. The real locus is cut into segments [p_k, p_{k+1}] on S^3. For each segment the frame at its
      midpoint m is (m, tangent t, n1 = direction of the second difference p_k + p_{k+1} - 2m, n2 = the rest).
      Along a segment of length h the components are ~h/2 (tangent), ~kappa h^2/8 (n1) and ~c h^3/48 (n2),
      so the box volume is ~H^4 h^6 and the number of segments grows like H^(2/3) (elkies4.c: H^(4/5), with both
      normals given the n1 width). The half-widths are the endpoint components times a safety factor: for
      f(s) = b s^2 + c s^3 the interior extremum is at most half the endpoint maximum, and the tangent component
      is monotone.
   2. Precision. Geometry (projection onto the curve, frames, the scaled basis W*T and LLL) is done in __float128.
      elkies4.c added H*1e-16 to the normal widths, which dominates the box once H >~ 1e10.
   Each segment: |<u,m>| <= H, |<u,t>| <= w1, |<u,n1>| <= w2, |<u,n2>| <= w3; the box lies in the radius-2 ball of the
   scaled form; LLL plus Fincke-Pohst enumerates the ball; every candidate is checked exactly.
   A component on S^3 is traced from a seed until it returns to +seed or -seed (u and -u are the same point).
   Usage: elkies4q H VOL SEEDS [F0 F1] < coeffs.txt   (optional: only the fraction [F0, F1) of each component's length,
   so that several cores can share one curve). Build: gcc -O2 -o elkies4q elkies4q.c -lquadmath -lm */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <quadmath.h>
typedef long double LD;
typedef __float128 Q;
typedef __int128 i128;
typedef long long ll;

static ll c2[10], c3[10];
static Q S2[4][4], S3[4][4];
static Q H, VOL;
static long n_arcs = 0, n_enum = 0, n_points = 0, n_dense = 0, n_retry = 0;
#define DENSE_CAP 2e5L

static void symm(const ll *c, Q S[4][4]) {
    Q m = 0; for (int k = 0; k < 10; k++) if (fabsq((Q)c[k]) > m) m = fabsq((Q)c[k]);
    int k = 0;
    for (int i = 0; i < 4; i++) for (int j = i; j < 4; j++, k++) {
        if (i == j) S[i][i] = (Q)c[k] / m; else S[i][j] = S[j][i] = (Q)c[k] / (2 * m);
    }
}
static Q quad(Q S[4][4], const Q *u) { Q r = 0; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) r += S[i][j] * u[i] * u[j]; return r; }
static void grad(Q S[4][4], const Q *u, Q *g) { for (int i = 0; i < 4; i++) { g[i] = 0; for (int j = 0; j < 4; j++) g[i] += 2 * S[i][j] * u[j]; } }
static Q dot(const Q *a, const Q *b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]; }
static void normalize(Q *a) { Q n = sqrtq(dot(a, a)); for (int i = 0; i < 4; i++) a[i] /= n; }
static Q dist2(const Q *a, const Q *b) { Q r = 0; for (int i = 0; i < 4; i++) r += (a[i]-b[i])*(a[i]-b[i]); return r; }
static Q det3(Q a[3][3]) {
    return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}
/* Unit vector orthogonal to three 4-vectors (generalized cross product). */
static void cross4(const Q *r0, const Q *r1, const Q *r2, Q *t) {
    const Q *R[3] = {r0, r1, r2};
    for (int k = 0; k < 4; k++) {
        Q M[3][3];
        for (int r = 0; r < 3; r++) for (int c = 0, cc = 0; c < 4; c++) if (c != k) M[r][cc++] = R[r][c];
        t[k] = ((k & 1) ? -1 : 1) * det3(M);
    }
    normalize(t);
}
/* Newton (minimum-norm steps) onto Q2 = Q3 = 0, |u| = 1. */
static int project(Q *u) {
    for (int it = 0; it < 60; it++) {
        Q F[3] = {quad(S2, u), quad(S3, u), (dot(u, u) - 1) / 2}, J[3][4];
        grad(S2, u, J[0]); grad(S3, u, J[1]); for (int i = 0; i < 4; i++) J[2][i] = u[i];
        if (fabsq(F[0]) + fabsq(F[1]) + fabsq(F[2]) < 1e-31Q) return 1;
        Q G[3][3]; for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) G[a][b] = dot(J[a], J[b]);
        Q d = det3(G); if (fabsq(d) < 1e-60Q) return 0;
        Q y[3];
        for (int c = 0; c < 3; c++) { Q M[3][3]; memcpy(M, G, sizeof M); for (int r = 0; r < 3; r++) M[r][c] = F[r]; y[c] = det3(M) / d; }
        for (int i = 0; i < 4; i++) u[i] -= J[0][i]*y[0] + J[1][i]*y[1] + J[2][i]*y[2];
    }
    return 0;
}
static void tangent(const Q *p, Q *t) { Q g2[4], g3[4]; grad(S2, p, g2); grad(S3, p, g3); cross4(p, g2, g3, t); }
/* Next curve point at chord distance ~h along t; its tangent oriented like t. */
static int step_to(const Q *p, const Q *t, Q h, Q *q, Q *tq) {
    for (int i = 0; i < 4; i++) q[i] = p[i] + h * t[i];
    normalize(q);
    if (!project(q)) return 0;
    tangent(q, tq); if (dot(tq, t) < 0) for (int i = 0; i < 4; i++) tq[i] = -tq[i];
    return 1;
}
/* p has just passed target (+seed or -seed) along t. */
static int passed(const Q *p, const Q *t, const Q *seed, Q h) {
    for (int sgn = 1; sgn >= -1; sgn -= 2) {
        Q d[4]; for (int i = 0; i < 4; i++) d[i] = sgn * seed[i] - p[i];
        if (dot(d, d) < 16 * h * h && dot(d, t) <= 0) return 1;
    }
    return 0;
}

/* ---- lattice part ---- */
static ll T[4][4];          /* current unimodular transform (warm start between segments) */
static Q W[4][4];           /* scaled frame: rows e_a / w_a */
static void apply(const ll *row, Q *out) { for (int a = 0; a < 4; a++) { out[a] = 0; for (int j = 0; j < 4; j++) out[a] += W[a][j] * (Q)row[j]; } }

static void gso(Q b[4][4], Q bs[4][4], Q mu[4][4], Q *B) {
    for (int i = 0; i < 4; i++) {
        memcpy(bs[i], b[i], sizeof(Q) * 4);
        for (int j = 0; j < i; j++) { mu[i][j] = dot(b[i], bs[j]) / B[j]; for (int k = 0; k < 4; k++) bs[i][k] -= mu[i][j] * bs[j][k]; }
        B[i] = dot(bs[i], bs[i]);
    }
}
static void lll(Q b[4][4]) {
    Q bs[4][4], mu[4][4], B[4];
    int k = 1, guard = 0;
    gso(b, bs, mu, B);
    while (k < 4 && guard++ < 20000) {
        for (int j = k - 1; j >= 0; j--) {
            Q q = roundq(mu[k][j]);
            if (q != 0) {
                ll qi = (ll)q;
                for (int a = 0; a < 4; a++) { b[k][a] -= q * b[j][a]; T[k][a] -= qi * T[j][a]; }
                gso(b, bs, mu, B);
            }
        }
        if (B[k] >= (0.99Q - mu[k][k-1] * mu[k][k-1]) * B[k-1]) k++;
        else {
            for (int a = 0; a < 4; a++) { Q s = b[k][a]; b[k][a] = b[k-1][a]; b[k-1][a] = s; ll t = T[k][a]; T[k][a] = T[k-1][a]; T[k-1][a] = t; }
            gso(b, bs, mu, B);
            if (k > 1) k--;
        }
    }
}

#define MAXP 4096
static ll found[MAXP][4]; static int nfound = 0;
static ll gcdll(ll a, ll b) { if (a < 0) a = -a; if (b < 0) b = -b; while (b) { ll t = a % b; a = b; b = t; } return a; }
static void report(ll *u) {
    ll g = gcdll(gcdll(u[0], u[1]), gcdll(u[2], u[3])); if (g == 0) return;
    ll v[4]; for (int i = 0; i < 4; i++) v[i] = u[i] / g;
    int s = 0; for (int i = 0; i < 4; i++) if (v[i]) { s = v[i] < 0; break; }
    if (s) for (int i = 0; i < 4; i++) v[i] = -v[i];
    for (int k = 0; k < nfound; k++) if (!memcmp(found[k], v, sizeof v)) return;
    if (nfound < MAXP) memcpy(found[nfound++], v, sizeof v);
    printf("POINT %lld %lld %lld %lld\n", v[0], v[1], v[2], v[3]); fflush(stdout); n_points++;
}
static int exact_on_curve(const ll *u) {
    i128 r2 = 0, r3 = 0; int k = 0;
    for (int i = 0; i < 4; i++) for (int j = i; j < 4; j++, k++) { i128 m = (i128)u[i] * u[j]; r2 += c2[k] * m; r3 += c3[k] * m; }
    return r2 == 0 && r3 == 0;
}

/* Enumerate integer x with |sum x_i b_i|^2 <= R2 (Fincke-Pohst on the reduced basis, long double). */
static void ldot_gso(LD b[4][4], LD bs[4][4], LD mu[4][4], LD *B) {
    for (int i = 0; i < 4; i++) {
        memcpy(bs[i], b[i], sizeof(LD) * 4);
        for (int j = 0; j < i; j++) {
            LD d = 0; for (int k = 0; k < 4; k++) d += b[i][k] * bs[j][k];
            mu[i][j] = d / B[j]; for (int k = 0; k < 4; k++) bs[i][k] -= mu[i][j] * bs[j][k];
        }
        B[i] = 0; for (int k = 0; k < 4; k++) B[i] += bs[i][k] * bs[i][k];
    }
}
static void enumerate(LD b[4][4], LD R2) {
    LD bs[4][4], mu[4][4], B[4]; ldot_gso(b, bs, mu, B);
    LD x[4], c[4], l[4]; int i = 3;
    c[3] = 0; l[3] = 0; x[3] = ceill(-sqrtl(R2 / B[3]));
    LD hi[4]; hi[3] = floorl(sqrtl(R2 / B[3]));
    for (;;) {
        if (x[i] > hi[i]) { if (++i > 3) break; x[i] += 1; continue; }
        LD d = x[i] - c[i]; LD li = l[i] + d * d * B[i];
        if (li > R2 * (1 + 1e-12L)) { x[i] += 1; continue; }
        if (i == 0) {
            n_enum++;
            ll u[4] = {0, 0, 0, 0}; int nz = 0;
            for (int a = 0; a < 4; a++) { ll xa = (ll)x[a]; if (xa) nz = 1; for (int j = 0; j < 4; j++) u[j] += xa * T[a][j]; }
            if (nz && exact_on_curve(u)) report(u);
            x[0] += 1; continue;
        }
        int j = i - 1; l[j] = li;
        if (j == 1) {   /* the remaining (x0, x1) slice: if dense, hand it to the exact plane solver */
            LD area = M_PI * (R2 * (1 + 1e-12L) - li) / sqrtl(B[0] * B[1]);
            if (area > DENSE_CAP) {
                ll w[4] = {0, 0, 0, 0};
                for (int a = 2; a < 4; a++) { ll xa = (ll)x[a]; for (int q = 0; q < 4; q++) w[q] += xa * T[a][q]; }
                printf("DENSE %lld %lld %lld %lld  %lld %lld %lld %lld  %lld %lld %lld %lld\n",
                       T[0][0], T[0][1], T[0][2], T[0][3], T[1][0], T[1][1], T[1][2], T[1][3], w[0], w[1], w[2], w[3]);
                n_dense++; x[i] += 1; continue;
            }
        }
        LD s = 0; for (int k = i; k < 4; k++) s += x[k] * mu[k][j];
        c[j] = -s;
        LD r = sqrtl((R2 * (1 + 1e-12L) - li) / B[j]);
        x[j] = ceill(c[j] - r); hi[j] = floorl(c[j] + r);
        i = j;
    }
}

static void lattice_segment(Q e[4][4], const Q *w) {
    for (int a = 0; a < 4; a++) for (int j = 0; j < 4; j++) W[a][j] = e[a][j] / w[a];
    Q b[4][4];
    for (int pass = 0; pass < 2; pass++) { for (int i = 0; i < 4; i++) apply(T[i], b[i]); lll(b); }
    for (int i = 0; i < 4; i++) apply(T[i], b[i]);   /* exact integer transform, quad frame */
    LD bl[4][4]; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) bl[i][j] = (LD)b[i][j];
    long e0 = n_enum;
    enumerate(bl, 4.04L);   /* box [-1,1]^4 lies in radius 2; 1% slack for rounding */
    n_arcs++;
    if (getenv("ELK_DEBUG") && n_arcs % 20000 == 1)
        fprintf(stderr, "seg %ld w1=%.3Le w2=%.3Le w3=%.3Le count=%ld\n", n_arcs, (LD)w[1], (LD)w[2], (LD)w[3], n_enum - e0);
}

static unsigned long long rng = 88172645463325252ULL;
static LD urand(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return (rng >> 11) * (1.0L / 9007199254740992.0L); }
static LD gauss(void) { LD a = urand() + 1e-30L, b = urand(); return sqrtl(-2 * logl(a)) * cosl(2 * M_PI * b); }

#define MAXS 200000
static Q samp[MAXS][4]; static int nsamp = 0;
static int near_traced(const Q *p, Q tol) {
    for (int k = 0; k < nsamp; k++) {
        Q d1 = 0, d2 = 0; for (int i = 0; i < 4; i++) { d1 += (p[i]-samp[k][i])*(p[i]-samp[k][i]); d2 += (p[i]+samp[k][i])*(p[i]+samp[k][i]); }
        if (d1 < tol*tol || d2 < tol*tol) return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 4 && argc != 6) { fprintf(stderr, "usage: elkies4q H VOL SEEDS [F0 F1] < coeffs\n"); return 2; }
    H = strtoflt128(argv[1], 0); VOL = strtoflt128(argv[2], 0); int seeds = atoi(argv[3]);
    Q F0 = argc == 6 ? strtoflt128(argv[4], 0) : 0, F1 = argc == 6 ? strtoflt128(argv[5], 0) : 1;
    for (int k = 0; k < 10; k++) if (scanf("%lld", &c2[k]) != 1) return 3;
    for (int k = 0; k < 10; k++) if (scanf("%lld", &c3[k]) != 1) return 3;
    for (int k = 0; k < 10; k++) if (llabs(c2[k]) >= (1LL << 40) || llabs(c3[k]) >= (1LL << 40)) return 4;
    /* exact check: candidates have |u_i| <= 2H (radius-2 ball), so |Q(u)| <= 10 * maxc * 4 H^2 must stay below 2^126 */
    ll maxc = 1; for (int k = 0; k < 10; k++) { if (llabs(c2[k]) > maxc) maxc = llabs(c2[k]); if (llabs(c3[k]) > maxc) maxc = llabs(c3[k]); }
    if (40 * (Q)maxc * H * H > 8.5e37Q || H > 1e15Q) { fprintf(stderr, "H too large for exact int128 check\n"); return 4; }
    symm(c2, S2); symm(c3, S3);
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) T[i][j] = (i == j);
    const Q HC = 2e-4Q, SAMPLE_STEP = 0.004Q, TOL = 0.01Q, DMAX = 1e-3Q, HMIN = 1e-13Q, MARGIN = 1e-3Q;
    const Q SAFE = 1.3Q, FLOOR = H * 1e-27Q + 1e-30Q;
    int comp = 0;
    for (int s = 0; s < seeds; s++) {
        Q seed[4]; for (int i = 0; i < 4; i++) seed[i] = (Q)gauss(); normalize(seed);
        if (!project(seed) || near_traced(seed, TOL)) continue;
        comp++;
        /* coarse pass: length to the first return to +-seed; samples for the seed test */
        Q p[4], t[4], q[4], tq[4], L = 0, next_sample = 0;
        memcpy(p, seed, sizeof p); tangent(p, t);
        int closed = 0;
        while (L < 50) {
            if (!step_to(p, t, HC, q, tq)) { fprintf(stderr, "projection failed (coarse) at %Lg\n", (LD)L); return 5; }
            L += sqrtq(dist2(p, q));
            if (L >= next_sample && nsamp < MAXS) { memcpy(samp[nsamp++], q, sizeof q); next_sample += SAMPLE_STEP; }
            memcpy(p, q, sizeof p); memcpy(t, tq, sizeof t);
            if (L > 0.05Q && passed(p, t, seed, HC)) { closed = 1; break; }
        }
        if (!closed) fprintf(stderr, "component %d did not close in the coarse pass (L=%Lg)\n", comp, (LD)L);
        Q lo = F0 > 0 ? F0 * L - MARGIN : -1, hi = F1 < 1 ? F1 * L + MARGIN : 1e30Q;
        /* fine pass */
        memcpy(p, seed, sizeof p); tangent(p, t);
        Q total = 0;
        while (total + HC < lo) {
            if (!step_to(p, t, HC, q, tq)) { fprintf(stderr, "projection failed (skip) at %Lg\n", (LD)total); return 5; }
            total += sqrtq(dist2(p, q)); memcpy(p, q, sizeof p); memcpy(t, tq, sizeof t);
        }
        long arcs0 = n_arcs;
        Q h = powq(48 * VOL / (H*H*H*H), 1 / 6.0Q); if (h > DMAX) h = DMAX;
        Q maxV = 0;
        for (;;) {
            Q m[4], tm[4], e[4][4], w[4];
            if (!step_to(p, t, h, q, tq)) { fprintf(stderr, "projection failed at %Lg\n", (LD)total); return 5; }
            for (int i = 0; i < 4; i++) m[i] = p[i] + q[i];
            normalize(m);
            if (!project(m)) { fprintf(stderr, "midpoint projection failed at %Lg\n", (LD)total); return 5; }
            tangent(m, tm); if (dot(tm, t) < 0) for (int i = 0; i < 4; i++) tm[i] = -tm[i];
            Q d[4]; for (int i = 0; i < 4; i++) d[i] = p[i] + q[i] - 2 * m[i];
            Q dm = dot(d, m), dt = dot(d, tm); for (int i = 0; i < 4; i++) d[i] -= dm * m[i] + dt * tm[i];
            if (dot(d, d) < 1e-60Q) {   /* no measurable curvature: any unit normal */
                for (int k = 0; k < 4; k++) {
                    for (int i = 0; i < 4; i++) d[i] = (i == k);
                    Q a = dot(d, m), b2 = dot(d, tm);
                    for (int i = 0; i < 4; i++) d[i] -= a * m[i] + b2 * tm[i];
                    if (dot(d, d) > 0.1Q) break;
                }
            }
            normalize(d);
            memcpy(e[0], m, sizeof m); memcpy(e[1], tm, sizeof tm); memcpy(e[2], d, sizeof d); cross4(m, tm, d, e[3]);
            w[0] = H;
            for (int a = 1; a < 4; a++) {
                Q c = fabsq(dot(p, e[a])), c2v = fabsq(dot(q, e[a])); if (c2v > c) c = c2v;
                w[a] = H * c * (a == 1 ? 1.02Q : SAFE) + FLOOR;
            }
            Q V = 16 * w[0] * w[1] * w[2] * w[3];
            if (V > 3 * VOL && h > HMIN) { Q f = powq(VOL / V, 1 / 6.0Q); if (f < 0.3Q) f = 0.3Q; if (f > 0.9Q) f = 0.9Q; h *= f; n_retry++; continue; }
            lattice_segment(e, w);
            if (V > maxV) maxV = V;
            total += sqrtq(dist2(p, q));
            Q f = powq(VOL / V, 1 / 6.0Q); if (f < 0.7Q) f = 0.7Q; if (f > 1.4Q) f = 1.4Q;
            h *= f; if (h > DMAX) h = DMAX;
            memcpy(p, q, sizeof p); memcpy(t, tq, sizeof t);
            if (total > hi) break;
            if (F1 >= 1 && total > 0.05Q && passed(p, t, seed, h)) break;
            if (total > L + 0.05Q) { fprintf(stderr, "fine pass overran the coarse length (total=%Lg L=%Lg)\n", (LD)total, (LD)L); break; }
            if ((n_arcs & ((1 << 21) - 1)) == 0) { fprintf(stderr, "progress comp=%d total=%Lg/%Lg segments=%ld enum=%ld points=%ld h=%.3Le\n", comp, (LD)total, (LD)L, n_arcs, n_enum, n_points, (LD)h); fflush(stderr); }
        }
        printf("COMPONENT %d length %.6Lf range %.4Lf-%.4Lf segments %ld max_box_volume %.3Lf\n", comp, (LD)L, (LD)F0, (LD)F1, n_arcs - arcs0, (LD)maxV); fflush(stdout);
    }
    printf("DONE H=%Lg VOL=%Lg components=%d arcs=%ld enumerated=%ld points=%ld dense_slices=%ld retries=%ld\n", (LD)H, (LD)VOL, comp, n_arcs, n_enum, n_points, n_dense, n_retry);
    return 0;
}
