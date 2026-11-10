/* Elkies-style lattice search for rational points on a four-cover Q2 = Q3 = 0 in P^3 (world 3, U3).
   Intended to find every primitive integer u with |u|_2 <= H on the curve; completeness is NOT certified (real components
   are found from random seeds, curvature is estimated numerically, H is not an enforced cutoff); emitted points are exact:
   the real locus is traced on S^3; for each short arc (half-length delta, centre p, unit tangent tau,
   normals n1, n2) every integer u in the cone over the arc with |u| <= H satisfies
       |<u,p>| <= H,  |<u,tau>| <= H delta,  |<u,n_i>| <= eps = H (kappa delta^2 / 2) + margin,
   a box of volume 4 kappa^2 H^4 delta^5. The box lies in the ball of radius 2 for the scaled form;
   LLL plus Fincke-Pohst enumerates all integer points of that ball, and each is checked exactly.
   Input: 20 integers (Q2 then Q3; coefficients of u_i u_j, i <= j, order 00 01 02 03 11 12 13 22 23 33).
   A dense two-dimensional slice {x0 v0 + x1 v1 + w} (near a rational point of small height) is not enumerated:
   it is printed as DENSE v0 v1 w and solved exactly by dense_solve.py (a plane meets the quartic curve in <= 4 points).
   Usage: elkies4 H VOL SEEDS < coeffs.txt     Output: POINT u0 u1 u2 u3 (primitive), DENSE ..., COMPONENT ..., DONE ... */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
typedef long double LD;
typedef __int128 i128;
typedef long long ll;

static ll c2[10], c3[10];
static LD S2[4][4], S3[4][4];
static LD H, VOL;
static long n_arcs = 0, n_enum = 0, n_points = 0, n_dense = 0;
#define DENSE_CAP 2e5L

static void symm(const ll *c, LD S[4][4]) {
    LD m = 0; for (int k = 0; k < 10; k++) if (fabsl((LD)c[k]) > m) m = fabsl((LD)c[k]);
    int k = 0;
    for (int i = 0; i < 4; i++) for (int j = i; j < 4; j++, k++) {
        if (i == j) S[i][i] = c[k] / m; else S[i][j] = S[j][i] = c[k] / (2 * m);
    }
}
static LD quad(LD S[4][4], const LD *u) { LD r = 0; for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) r += S[i][j] * u[i] * u[j]; return r; }
static void grad(LD S[4][4], const LD *u, LD *g) { for (int i = 0; i < 4; i++) { g[i] = 0; for (int j = 0; j < 4; j++) g[i] += 2 * S[i][j] * u[j]; } }
static LD dot(const LD *a, const LD *b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]; }
static void normalize(LD *a) { LD n = sqrtl(dot(a, a)); for (int i = 0; i < 4; i++) a[i] /= n; }
static LD det3(LD a[3][3]) {
    return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}

/* Newton (minimum-norm steps) onto Q2 = Q3 = 0, |u| = 1. */
static int project(LD *u) {
    for (int it = 0; it < 80; it++) {
        LD F[3] = {quad(S2, u), quad(S3, u), (dot(u, u) - 1) / 2}, J[3][4];
        grad(S2, u, J[0]); grad(S3, u, J[1]); for (int i = 0; i < 4; i++) J[2][i] = u[i];
        if (fabsl(F[0]) + fabsl(F[1]) + fabsl(F[2]) < 100 * LDBL_EPSILON) return 1;
        LD G[3][3]; for (int a = 0; a < 3; a++) for (int b = 0; b < 3; b++) G[a][b] = dot(J[a], J[b]);
        LD d = det3(G); if (fabsl(d) < 1e-30L) return 0;
        LD y[3];
        for (int c = 0; c < 3; c++) { LD M[3][3]; memcpy(M, G, sizeof M); for (int r = 0; r < 3; r++) M[r][c] = F[r]; y[c] = det3(M) / d; }
        for (int i = 0; i < 4; i++) u[i] -= J[0][i]*y[0] + J[1][i]*y[1] + J[2][i]*y[2];
    }
    return 0;
}
/* Unit tangent: generalized cross product of p, grad Q2, grad Q3. */
static void tangent(const LD *p, LD *t) {
    LD R[3][4]; memcpy(R[0], p, sizeof(LD) * 4); grad(S2, p, R[1]); grad(S3, p, R[2]);
    for (int k = 0; k < 4; k++) {
        LD M[3][3];
        for (int r = 0; r < 3; r++) for (int c = 0, cc = 0; c < 4; c++) if (c != k) M[r][cc++] = R[r][c];
        t[k] = ((k & 1) ? -1 : 1) * det3(M);
    }
    normalize(t);
}
static void normals(const LD *p, const LD *t, LD *n1, LD *n2) {
    LD cand[4][4]; LD best[4]; int idx[4] = {0, 1, 2, 3};
    for (int e = 0; e < 4; e++) {
        for (int i = 0; i < 4; i++) cand[e][i] = (i == e);
        LD a = dot(cand[e], p), b = dot(cand[e], t);
        for (int i = 0; i < 4; i++) cand[e][i] -= a * p[i] + b * t[i];
        best[e] = dot(cand[e], cand[e]);
    }
    for (int i = 0; i < 4; i++) for (int j = i + 1; j < 4; j++) if (best[idx[j]] > best[idx[i]]) { int s = idx[i]; idx[i] = idx[j]; idx[j] = s; }
    memcpy(n1, cand[idx[0]], sizeof(LD) * 4); normalize(n1);
    memcpy(n2, cand[idx[1]], sizeof(LD) * 4);
    LD a = dot(n2, n1); for (int i = 0; i < 4; i++) n2[i] -= a * n1[i];
    normalize(n2);
}

/* ---- lattice part ---- */
static ll T[4][4];          /* current unimodular transform (warm start between arcs) */
static LD W[4][4];          /* scaled frame: rows p/H, tau/(H delta), n1/eps, n2/eps */
static void apply(const ll *row, LD *out) { for (int a = 0; a < 4; a++) { out[a] = 0; for (int j = 0; j < 4; j++) out[a] += W[a][j] * (LD)row[j]; } }

static void gso(LD b[4][4], LD bs[4][4], LD mu[4][4], LD *B) {
    for (int i = 0; i < 4; i++) {
        memcpy(bs[i], b[i], sizeof(LD) * 4);
        for (int j = 0; j < i; j++) { mu[i][j] = dot(b[i], bs[j]) / B[j]; for (int k = 0; k < 4; k++) bs[i][k] -= mu[i][j] * bs[j][k]; }
        B[i] = dot(bs[i], bs[i]);
    }
}
static void lll(LD b[4][4]) {
    LD bs[4][4], mu[4][4], B[4];
    int k = 1, guard = 0;
    gso(b, bs, mu, B);
    while (k < 4 && guard++ < 10000) {
        for (int j = k - 1; j >= 0; j--) {
            LD q = roundl(mu[k][j]);
            if (q != 0) {
                ll qi = (ll)q;
                for (int a = 0; a < 4; a++) { b[k][a] -= q * b[j][a]; T[k][a] -= qi * T[j][a]; }
                gso(b, bs, mu, B);
            }
        }
        if (B[k] >= (0.99L - mu[k][k-1] * mu[k][k-1]) * B[k-1]) k++;
        else {
            for (int a = 0; a < 4; a++) { LD s = b[k][a]; b[k][a] = b[k-1][a]; b[k-1][a] = s; ll t = T[k][a]; T[k][a] = T[k-1][a]; T[k-1][a] = t; }
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

/* Enumerate integer x with |sum x_i b_i|^2 <= R2 (Fincke-Pohst on the reduced basis). */
static void enumerate(LD b[4][4], LD R2) {
    LD bs[4][4], mu[4][4], B[4]; gso(b, bs, mu, B);
    LD x[4], c[4], l[4]; int i = 3;
    c[3] = 0; l[3] = 0; x[3] = ceill(-sqrtl(R2 / B[3]));
    /* simple recursive-style enumeration with explicit stack */
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

static void lattice_arc(const LD *p, const LD *t, const LD *n1, const LD *n2, LD delta, LD eps) {
    for (int j = 0; j < 4; j++) { W[0][j] = p[j] / H; W[1][j] = t[j] / (H * delta); W[2][j] = n1[j] / eps; W[3][j] = n2[j] / eps; }
    LD b[4][4]; for (int i = 0; i < 4; i++) apply(T[i], b[i]);
    lll(b);
    for (int i = 0; i < 4; i++) apply(T[i], b[i]);   /* recompute from the exact integer transform */
    long e0 = n_enum;
    enumerate(b, 4.04L);   /* box [-1,1]^4 lies in radius 2; 1% slack for rounding */
    n_arcs++;
    if (getenv("ELK_DEBUG") && n_arcs % 50000 == 1) {
        LD bs[4][4], mu[4][4], B[4]; gso(b, bs, mu, B);
        fprintf(stderr, "arc %ld delta=%.3Le eps=%.3Le count=%ld B=%.2Le %.2Le %.2Le %.2Le T0=%lld %lld %lld %lld\n", n_arcs, delta, eps, n_enum - e0, B[0], B[1], B[2], B[3], T[0][0], T[0][1], T[0][2], T[0][3]);
    }
}

static unsigned long long rng = 88172645463325252ULL;
static LD urand(void) { rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return (rng >> 11) * (1.0L / 9007199254740992.0L); }
static LD gauss(void) { LD a = urand() + 1e-30L, b = urand(); return sqrtl(-2 * logl(a)) * cosl(2 * M_PI * b); }

#define MAXS 200000
static LD samp[MAXS][4]; static int nsamp = 0;
static int near_traced(const LD *p, LD tol) {
    for (int k = 0; k < nsamp; k++) {
        LD d1 = 0, d2 = 0; for (int i = 0; i < 4; i++) { d1 += (p[i]-samp[k][i])*(p[i]-samp[k][i]); d2 += (p[i]+samp[k][i])*(p[i]+samp[k][i]); }
        if (d1 < tol*tol || d2 < tol*tol) return 1;
    }
    return 0;
}

int main(int argc, char **argv) {
    if (argc != 4) { fprintf(stderr, "usage: elkies4 H VOL SEEDS < coeffs\n"); return 2; }
    H = strtold(argv[1], 0); VOL = strtold(argv[2], 0); int seeds = atoi(argv[3]);
    for (int k = 0; k < 10; k++) if (scanf("%lld", &c2[k]) != 1) return 3;
    for (int k = 0; k < 10; k++) if (scanf("%lld", &c3[k]) != 1) return 3;
    for (int k = 0; k < 10; k++) if (llabs(c2[k]) >= (1LL << 40) || llabs(c3[k]) >= (1LL << 40)) return 4;
    if (H > 1e12L) { fprintf(stderr, "H too large for exact int128 check\n"); return 4; }
    symm(c2, S2); symm(c3, S3);
    for (int i = 0; i < 4; i++) for (int j = 0; j < 4; j++) T[i][j] = (i == j);
    const LD SAMPLE_STEP = 0.004L, TOL = 0.01L, KFLOOR = 0.05L, DMAX = 1e-3L;
    int comp = 0;
    for (int s = 0; s < seeds; s++) {
        LD p[4]; for (int i = 0; i < 4; i++) p[i] = gauss(); normalize(p);
        if (!project(p) || near_traced(p, TOL)) continue;
        comp++;
        LD seed[4]; memcpy(seed, p, sizeof seed);
        LD t[4]; tangent(p, t);
        LD total = 0, next_sample = 0, kprev = KFLOOR, maxk = 0;
        long arcs0 = n_arcs;
        for (;;) {
            /* trial step to measure curvature ahead */
            LD dguess = powl(VOL / (4 * (2*kprev) * (2*kprev) * H*H*H*H), 0.2L); if (dguess > DMAX) dguess = DMAX;
            LD h = 1.8L * dguess, q[4], tq[4];
            for (int i = 0; i < 4; i++) q[i] = p[i] + h * t[i];
            normalize(q);
            if (!project(q)) { fprintf(stderr, "projection failed at total=%Lg\n", total); return 5; }
            tangent(q, tq); if (dot(tq, t) < 0) for (int i = 0; i < 4; i++) tq[i] = -tq[i];
            LD step = sqrtl((q[0]-p[0])*(q[0]-p[0]) + (q[1]-p[1])*(q[1]-p[1]) + (q[2]-p[2])*(q[2]-p[2]) + (q[3]-p[3])*(q[3]-p[3]));
            LD a[4]; for (int i = 0; i < 4; i++) a[i] = (tq[i] - t[i]) / step;
            LD ap = dot(a, p), at = dot(a, t); for (int i = 0; i < 4; i++) a[i] -= ap * p[i] + at * t[i];
            LD kap = sqrtl(dot(a, a)); if (kap > maxk) maxk = kap;
            LD keff = 2 * (kap > kprev ? kap : kprev); if (keff < KFLOOR) keff = KFLOOR;
            LD delta = powl(VOL / (4 * keff * keff * H*H*H*H), 0.2L); if (delta > DMAX) delta = DMAX;
            if (step > 1.8L * delta * 1.0001L) { kprev = keff; continue; }   /* curvature grew: retry with a shorter step */
            LD eps = H * keff * delta * delta / 2 + H * 1e-16L + 1e-9L;
            LD n1[4], n2[4]; normals(p, t, n1, n2);
            lattice_arc(p, t, n1, n2, delta, eps);
            total += step;
            if (total >= next_sample && nsamp < MAXS) { memcpy(samp[nsamp++], p, sizeof(LD) * 4); next_sample += SAMPLE_STEP; }
            memcpy(p, q, sizeof p); memcpy(t, tq, sizeof t); kprev = kap > KFLOOR ? kap : KFLOOR;
            LD ds = 0; for (int i = 0; i < 4; i++) ds += (p[i]-seed[i])*(p[i]-seed[i]);
            if (total > 0.05L && ds < 4 * h * h) break;
            if (total > 50) { fprintf(stderr, "component did not close (total=%Lg)\n", total); break; }
            if ((n_arcs & ((1 << 22) - 1)) == 0) { fprintf(stderr, "progress comp=%d total=%Lg arcs=%ld enum=%ld points=%ld\n", comp, total, n_arcs, n_enum, n_points); fflush(stderr); }
        }
        printf("COMPONENT %d length %.6Lf arcs %ld max_kappa %.4Lf\n", comp, total, n_arcs - arcs0, maxk); fflush(stdout);
    }
    printf("DONE H=%Lg VOL=%Lg components=%d arcs=%ld enumerated=%ld points=%ld dense_slices=%ld\n", H, VOL, comp, n_arcs, n_enum, n_points, n_dense);
    return 0;
}
