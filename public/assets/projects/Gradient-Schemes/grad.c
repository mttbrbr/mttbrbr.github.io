/**
 * Gradient schemes: Gauss linear vs least squares, on five kinds of mesh.
 *
 * Periodic quadrilateral meshes on [0,1) x [0,Ly), N x N cells. The field and
 * its exact gradient are known analytically. For every cell the gradient is
 * computed with
 *   Gauss linear : grad = 1/V sum_f phi_f S_f, phi_f linearly interpolated
 *                  with OpenFOAM's face weights
 *   leastSquares : minimise sum_f w_f (phi_N - phi_P - d_f . grad)^2,
 *                  w_f = 1/|d_f|^2 (OpenFOAM's leastSquares weights)
 * Meshes: uniform, skewed (rows leaning alternately), random (nodes jittered),
 * stretched (aspect ratio 50), stretched + curved (thin cells following a
 * wavy surface, like a boundary layer over a bump).
 *
 * Build: gcc -O3 -march=native grad.c -o grad -lm     Run: ./grad outdir
 */
#define _DEFAULT_SOURCE
#ifndef ARV
#define ARV 50
#endif
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static int N, MESH; static double Lx = 1, Ly = 1, AR = 1, AMP = 0;
static double *jx, *jy;       /* random jitter of the nodes */
enum { UNIFORM, SKEWED, RANDOM, STRETCHED, CURVED };

/* node (i, j), any integers (periodic images) */
static void node(int i, int j, double *x, double *y) {
    int I = ((i % N) + N) % N, J = ((j % N) + N) % N;
    double ox = (double)((i - I) / N) * Lx, oy = (double)((j - J) / N) * Ly;   /* integer number of periods */
    double hx = Lx / N, hy = Ly / N;
    double X = I * hx, Y = J * hy;
    switch (MESH) {
    case SKEWED: X += 0.3 * hx * ((J % 2) ? -1 : 1); break;
    case RANDOM: X += jx[J * N + I] * hx; Y += jy[J * N + I] * hy; break;
    case CURVED: Y += AMP * sin(2 * M_PI * X / Lx); break;
    default: break;
    }
    *x = X + ox; *y = Y + oy;
}
/* field: phi = sin(2 pi x) sin(2 pi eta), eta = (y - AMP sin(2 pi x)) / Ly  (follows the mesh layers) */
static void field(double x, double y, double *p, double *gx, double *gy) {
    double s = AMP * sin(2 * M_PI * x), ds = AMP * 2 * M_PI * cos(2 * M_PI * x);
    double eta = (y - s) / Ly;
    double A = sin(2 * M_PI * x), B = sin(2 * M_PI * eta), dA = 2 * M_PI * cos(2 * M_PI * x), dB = 2 * M_PI * cos(2 * M_PI * eta);
    *p = A * B;
    *gx = dA * B + A * dB * (-ds / Ly);
    *gy = A * dB / Ly;
}
typedef struct { double cx, cy, V; } Cell;
static void cell_geom(int i, int j, Cell *c) {
    double x[4], y[4];
    node(i, j, &x[0], &y[0]); node(i + 1, j, &x[1], &y[1]); node(i + 1, j + 1, &x[2], &y[2]); node(i, j + 1, &x[3], &y[3]);
    double A = 0, cx = 0, cy = 0;
    for (int k = 0; k < 4; k++) {
        int l = (k + 1) % 4; double cr = x[k] * y[l] - x[l] * y[k];
        A += cr; cx += (x[k] + x[l]) * cr; cy += (y[k] + y[l]) * cr;
    }
    A *= 0.5; c->V = A; c->cx = cx / (6 * A); c->cy = cy / (6 * A);
}

int main(int argc, char **argv) {
    const char *dir = argc > 1 ? argv[1] : ".";
    char path[512]; snprintf(path, sizeof path, "%s/gradients.csv", dir);
    FILE *f = fopen(path, "w"); fprintf(f, "mesh,N,gauss,leastsquares,leastsquares_cellLimited\n");
    const char *names[] = {"uniform", "skewed", "random", "stretched", "stretched+curved"};
    printf("relative L2 error of the gradient: Gauss linear | least squares | least squares + cellLimited\n");
    for (MESH = 0; MESH < 5; MESH++) {
        AR = (MESH == STRETCHED || MESH == CURVED) ? ARV : 1;
        Ly = 1.0 / AR;                      /* square domain in cells, cells of aspect ratio AR */
        printf("%-18s", names[MESH]);
        for (N = 16; N <= 256; N *= 2) {
            AMP = (MESH == CURVED) ? 0.05 : 0;     /* wavy layers: amplitude 0.05, much larger than Ly = 0.02 */
            jx = malloc(N * N * 8); jy = malloc(N * N * 8);
            unsigned long long r = 99;
            for (int k = 0; k < N * N; k++) {
                r = r * 6364136223846793005ULL + 1; jx[k] = 0.25 * (2.0 * (r >> 11) / 9007199254740992.0 - 1);
                r = r * 6364136223846793005ULL + 1; jy[k] = 0.25 * (2.0 * (r >> 11) / 9007199254740992.0 - 1);
            }
            Cell *C = malloc(N * N * sizeof(Cell)); double *p = malloc(N * N * 8);
            for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
                cell_geom(i, j, &C[j * N + i]);
                double gx, gy; field(C[j * N + i].cx, C[j * N + i].cy, &p[j * N + i], &gx, &gy);
            }
            double eg2 = 0, el2 = 0, ec2 = 0, n2 = 0;
            int di[4] = {1, 0, -1, 0}, dj[4] = {0, 1, 0, -1};
            for (int j = 0; j < N; j++) for (int i = 0; i < N; i++) {
                Cell P = C[j * N + i];
                double gx = 0, gy = 0, Gxx = 0, Gxy = 0, Gyy = 0, bx = 0, by = 0;
                for (int k = 0; k < 4; k++) {
                    int in = i + di[k], jn = j + dj[k];
                    Cell Nb; cell_geom(in, jn, &Nb);            /* neighbour with its periodic image coordinates */
                    double pn = p[((jn + N) % N) * N + (in + N) % N];
                    /* face nodes (counter-clockwise around P) */
                    double x1, y1, x2, y2;
                    if (k == 0) { node(i + 1, j, &x1, &y1); node(i + 1, j + 1, &x2, &y2); }
                    else if (k == 1) { node(i + 1, j + 1, &x1, &y1); node(i, j + 1, &x2, &y2); }
                    else if (k == 2) { node(i, j + 1, &x1, &y1); node(i, j, &x2, &y2); }
                    else { node(i, j, &x1, &y1); node(i + 1, j, &x2, &y2); }
                    double Sx = y2 - y1, Sy = -(x2 - x1), xf = 0.5 * (x1 + x2), yf = 0.5 * (y1 + y2);
                    double dx = Nb.cx - P.cx, dy = Nb.cy - P.cy;
                    /* OpenFOAM linear weight: w = (S . (C_N - x_f)) / (S . (C_N - C_P)) */
                    double w = (Sx * (Nb.cx - xf) + Sy * (Nb.cy - yf)) / (Sx * dx + Sy * dy);
                    double pf = w * p[j * N + i] + (1 - w) * pn;
                    gx += pf * Sx; gy += pf * Sy;
                    double ww = 1 / (dx * dx + dy * dy);
                    Gxx += ww * dx * dx; Gxy += ww * dx * dy; Gyy += ww * dy * dy;
                    bx += ww * dx * (pn - p[j * N + i]); by += ww * dy * (pn - p[j * N + i]);
                }
                gx /= P.V; gy /= P.V;
                double det = Gxx * Gyy - Gxy * Gxy, lx = (Gyy * bx - Gxy * by) / det, ly = (Gxx * by - Gxy * bx) / det;
                /* cellLimited 1 (applied to the least-squares gradient): scale the gradient so that
                   extrapolated face values stay within the min/max of the cell and its neighbours */
                double pmin = p[j * N + i], pmax = pmin, alpha = 1;
                for (int k = 0; k < 4; k++) { double pn = p[((j + dj[k] + N) % N) * N + (i + di[k] + N) % N]; if (pn < pmin) pmin = pn; if (pn > pmax) pmax = pn; }
                for (int k = 0; k < 4; k++) {
                    double x1, y1, x2, y2;
                    if (k == 0) { node(i + 1, j, &x1, &y1); node(i + 1, j + 1, &x2, &y2); }
                    else if (k == 1) { node(i + 1, j + 1, &x1, &y1); node(i, j + 1, &x2, &y2); }
                    else if (k == 2) { node(i, j + 1, &x1, &y1); node(i, j, &x2, &y2); }
                    else { node(i, j, &x1, &y1); node(i + 1, j, &x2, &y2); }
                    double dpf = lx * (0.5 * (x1 + x2) - P.cx) + ly * (0.5 * (y1 + y2) - P.cy), lim = 1;
                    if (dpf > 1e-300) lim = (pmax - p[j * N + i]) / dpf; else if (dpf < -1e-300) lim = (pmin - p[j * N + i]) / dpf;
                    if (lim < alpha) alpha = lim;
                }
                if (alpha < 0) alpha = 0;
                double cxg = alpha * lx, cyg = alpha * ly;
                double pe, ex, ey; field(P.cx, P.cy, &pe, &ex, &ey);
                ec2 += ((cxg - ex) * (cxg - ex) + (cyg - ey) * (cyg - ey)) * P.V;
                eg2 += ((gx - ex) * (gx - ex) + (gy - ey) * (gy - ey)) * P.V;
                el2 += ((lx - ex) * (lx - ex) + (ly - ey) * (ly - ey)) * P.V;
                n2 += (ex * ex + ey * ey) * P.V;
            }
            double eg = sqrt(eg2 / n2), el = sqrt(el2 / n2), ec = sqrt(ec2 / n2);
            fprintf(f, "%s,%d,%.6e,%.6e,%.6e\n", names[MESH], N, eg, el, ec);
            printf("  %4d: %.1e | %.1e | %.1e", N, eg, el, ec);
            free(jx); free(jy); free(C); free(p);
        }
        printf("\n");
    }
    fclose(f);
    return 0;
}
