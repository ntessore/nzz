// nzz: compute two-dimensional per-bin redshift distribution
// ---
// author:  Nicolas Tessore <nicolas.tessore@manchester.ac.uk>
// date:    28 May 2019

#define _XOPEN_SOURCE 600

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <signal.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static const int RW = 5;

int mapsort(const void* a, const void* b)
{
    const double* x = a;
    const double* y = b;
    
    if(x[4] < y[4])
        return -1;
    if(x[4] > y[4])
        return +1;
    
    if(x[2] < y[2])
        return -1;
    if(x[2] > y[2])
        return +1;
    if(x[1] < y[1])
        return -1;
    if(x[1] > y[1])
        return +1;
    if(x[0] < y[0])
        return -1;
    if(x[0] > y[0])
        return +1;
    
    return 0;
}

static inline int index(double x, double y, double z, double s, int w, int h)
{
    return (int)(z/s)*(w*h) + (int)(y/s)*w + (int)(x/s);
}

static inline int query(int q, const int ma[], int gx, int gy, int gz,
                                                    int gr, int* c, int v[])
{
    int i, il, ih, j, jl, jh, k, kl, kh, l, m, n, p;
    
    i = q/(gx*gy);
    j = (q/gx)%gy;
    k = q%gx;
    
    il = i > gr ? i-gr : 0;
    ih = i+gr < gz ? i+gr+1 : gz;
    
    jl = j > gr ? j-gr : 0;
    jh = j+gr < gy ? j+gr+1 : gy;
    
    kl = k > gr ? k-gr : 0;
    kh = k+gr < gx ? k+gr+1 : gx;
    
    n = 0;
    p = -1;
    
    for(i = il; i < ih; ++i)
    {
        for(j = jl; j < jh; ++j)
        {
            k = (i*gy + j)*gx;
            l = ma[k + kl];
            m = ma[k + kh];
            
            if(l == p)
                p = (v[2*n-1] = m);
            else
                p = (v[2*n+0] = l, v[2*n+1] = m), ++n;
        }
    }
    
    *c = n;
    
    return q;
}

#include "io.c" // yes, really

static const char* ANIM[] = {
    "\033[34m\xe2\xa0\xb7\033[0m", "\033[34m\xe2\xa0\xaf\033[0m",
    "\033[34m\xe2\xa0\x9f\033[0m", "\033[34m\xe2\xa0\xbb\033[0m",
    "\033[34m\xe2\xa0\xbd\033[0m", "\033[34m\xe2\xa0\xbe\033[0m"
};
static const int NANIM = sizeof(ANIM)/sizeof(*ANIM);

volatile sig_atomic_t AL;
volatile sig_atomic_t QQ;

void handler(int s)
{
    AL = (s == SIGALRM);
    QQ = (s == SIGQUIT);
    signal(s, handler);
}

int main(int argc, char* argv[])
{
    char* cfgfile;
    struct config cfg;
    
    bool ls, sc, tc;
    int nc, nd, nr;
    double dl, dh, d0, dm, Dl, Dh, rl, rh, rm;
    double ui, uo;
    
    int rc[2];
    double* rv[2];
    int cc, xc;
    double* cv;
    double* xv;
    
    double gs;
    int gr;
    double xl, xh, yl, yh, zl, zh;
    int gx, gy, gz, ng;
    int* ma;
    
    double* Z;
    double* N;
    double s;
    
    time_t st;
    int dt;
    
    int i, j;
    
    char* bf, *nf, *sv, *sx;
    
    if(isatty(fileno(stdout)))
    {
        bf = "\033[1m";
        nf = "\033[0m";
        sv = "\033[32m\xe2\x9c\x94\033[0m";
        sx = "\033[31m\xe2\x9c\x98\033[0m";
    }
    else
    {
        bf = nf = "";
        sv = sx = ">";
    }
    
    cfgfile = NULL;
    memset(&cfg, 0, sizeof(cfg));
    
    if(argc > 5)
        goto err_usage;
    
    if(argc > 1 && strcmp(argv[1], "--") != 0)
        cfgfile = strdup(argv[1]);
    if(argc > 2 && strcmp(argv[2], "--") != 0)
        cfg.output = strdup(argv[2]);
    for(i = 3; i < argc; ++i)
        if(strcmp(argv[i], "--") != 0)
            cfg.catv[cfg.catc++] = strdup(argv[i]);
    
    if(!cfgfile)
        cfgfile = strdup("nzz.cfg");
    
    readcfg(cfgfile, &cfg);
    
    printf("\n");
    printf("%sconfiguration file %s%s\n", bf, cfgfile, nf);
    printf("\n");
    printcfg(&cfg);
    printf("\n");
    
    sc = cfg.coords >= COORDS_LONLAT;
    ls = cfg.spacing == SPACING_LOG;
    
    ui = UCONV[cfg.units];
    uo = UCONV[cfg.thunit];
    
    nd = cfg.nth;
    dl = cfg.thmin*uo;
    dh = cfg.thmax*uo;
    
    nr = cfg.nz;
    rl = cfg.zmin;
    rh = cfg.zmax;
    
#ifdef _OPENMP
    if(cfg.num_threads)
        omp_set_num_threads(cfg.num_threads);
    tc = cfg.thread_data == TDATA_COPY;
#else
    tc = false;
#endif
    
    if(sc)
    {
        dl = 2*sin(0.5*dl);
        dh = 2*sin(0.5*dh);
    }
    
    if(ls)
    {
        d0 = log(dl);
        dm = nd/(log(dh) - d0);
    }
    else
    {
        d0 = dl;
        dm = nd/(dh - d0);
    }
    
    Dl = dl*dl;
    Dh = dh*dh;
    
    rm = nr/(rh - rl);
    
    for(nc = 0; nc < cfg.catc; ++nc)
    {
        printf("%sread catalog %d%s\n", bf, nc, nf);
        fflush(stdout);
        
        rc[nc] = 0;
        rv[nc] = NULL;
        readc(cfg.catv[nc], cfg.coords, ui, &rc[nc], &rv[nc]);
        
        printf("%s done with %d points\n", sv, rc[nc]);
        printf("\n");
    }
    
    printf("%sbuild index%s\n", bf, nf);
    fflush(stdout);
    
    gs = 0.25*dh;
    gr = ceil(dh/gs);
    
    xl = xh = rv[0][0];
    yl = yh = rv[0][1];
    zl = zh = rv[0][2];
    for(j = 0; j < nc; ++j)
    {
        for(i = 1; i < rc[j]; ++i)
        {
            if(rv[j][i*RW+0] < xl) xl = rv[j][i*RW+0];
            if(rv[j][i*RW+0] > xh) xh = rv[j][i*RW+0];
            if(rv[j][i*RW+1] < yl) yl = rv[j][i*RW+1];
            if(rv[j][i*RW+1] > yh) yh = rv[j][i*RW+1];
            if(rv[j][i*RW+2] < zl) zl = rv[j][i*RW+2];
            if(rv[j][i*RW+2] > zh) zh = rv[j][i*RW+2];
        }
    }
    
    gx = floor((xh - xl)/gs) + 1;
    gy = floor((yh - yl)/gs) + 1;
    gz = floor((zh - zl)/gs) + 1;
    
    ng = gx*gy*gz;
    
    for(j = 0; j < nc; ++j)
    {
        for(i = 0; i < rc[j]; ++i)
            rv[j][i*RW+4] = index(rv[j][i*RW+0]-xl, rv[j][i*RW+1]-yl,
                                                rv[j][i*RW+2]-zl, gs, gx, gy);
        qsort(rv[j], rc[j], RW*sizeof(double), mapsort);
    } 
    
    ma = malloc((ng+1)*sizeof(int));
    if(!ma)
        goto err_alloc;
    
    cc = rc[0];
    cv = rv[0];
    xc = rc[nc-1];
    xv = rv[nc-1];
    
    for(i = 0, j = 0; i < ng; ++i)
    {
        while(j < xc && xv[j*RW+4] < i)
            j += 1;
        ma[i] = j;
    }
    ma[ng] = xc;
    
    printf("%s done with %d x %d x %d grid cells\n", sv, gx, gy, gz);
    printf("\n");
    
    Z = calloc(nd, sizeof(double));
    N = calloc(nd*nr*nr, sizeof(double));
    if(!Z || !N)
        goto err_alloc;
    
    s = 0;
    
    signal(SIGALRM, handler);
    signal(SIGQUIT, handler);
    AL = QQ = 0;
    
    printf("%sworking%s\n", bf, nf);
    fflush(stdout);
    
    st = time(NULL);
    dt = 0;
    
    #pragma omp parallel default(none) shared(Z, N, s, AL, QQ, st, dt) \
        firstprivate(ls, sc, tc, nd, nc, d0, dm, Dl, Dh, nr, rl, rm, \
            gr, gx, gy, gz, ng, cc, cv, xc, xv, ma, ANIM, NANIM, stdout)
    {
        int q, qc, nq;
        int* qr;
        
        double* cv_;
        double* xv_;
        int* ma_;
        
        double* Z_;
        double* N_;
        double s_;
        
        bool fb;
        
        int i, j, jh;
        
        nq = 0;
        qr = malloc((2*gr+1)*(2*gr+1)*2*sizeof(int));
        if(!qr)
            perror(NULL), abort();
        
        if(tc)
        {
            cv_ = malloc(cc*RW*sizeof(double));
            if(cv != xv)
                xv_ = malloc(xc*RW*sizeof(double));
            else
                xv_ = cv_;
            ma_ = malloc((ng+1)*sizeof(int));
            if(!cv_ || !xv_ || !ma_)
                perror(NULL), abort();
            
            memcpy(cv_, cv, cc*RW*sizeof(double));
            if(cv != xv)
                memcpy(xv_, xv, xc*RW*sizeof(double));
            memcpy(ma_, ma, (ng+1)*sizeof(int));
        }
        else
        {
            cv_ = cv;
            xv_ = xv;
            ma_ = ma;
        }
        
        Z_ = calloc(nd, sizeof(double));
        N_ = calloc(nd*nr*nr, sizeof(double));
        if(!Z_ || !N_)
            perror(NULL), abort();
        
        s_ = 0;
        
        fb = false;
        
        #pragma omp master
        if(isatty(fileno(stdout)))
        {
            fb = true;
            AL = false;
            alarm(1);
            
#ifdef _OPENMP
            printf("\r%s %d thread(s) ", ANIM[0], omp_get_num_threads());
            fflush(stdout);
#endif
        }
        
        qc = -1;
        
        #pragma omp for schedule(dynamic, 1) nowait
        for(i = 0; i < cc; ++i)
        {
            const double xi = cv_[i*RW+0];
            const double yi = cv_[i*RW+1];
            const double zi = cv_[i*RW+2];
            const double ri = cv_[i*RW+3];
            const int    qi = cv_[i*RW+4];
            
            const int ni = rm*(ri - rl);
            
            if(QQ)
                continue;
            
            if(AL && fb)
            {
                dt = difftime(time(NULL), st);
                
                printf("\r%s %.2f%%", ANIM[dt%NANIM], 100.*i/cc);
                printf(" in %02d:%02d:%02d ", dt/3600, (dt/60)%60, dt%60);
                fflush(stdout);
                
                AL = false;
                alarm(1);
            }
            
            if(ni < 0 || ni >= nr)
                continue;
            
            if(qi != qc)
                qc = query(qi, ma, gx, gy, gz, gr, &nq, qr);
            
            for(q = 0; q < nq; ++q)
            {
                for(j = qr[2*q+0], jh = qr[2*q+1]; j < jh; ++j)
                {
                    const double xj = xv_[j*RW+0];
                    const double yj = xv_[j*RW+1];
                    const double zj = xv_[j*RW+2];
                    const double rj = xv_[j*RW+3];
                    
                    const int nj = rm*(rj - rl);
                    
                    const double dx = xi - xj;
                    const double dy = yi - yj;
                    const double dz = zi - zj;
                    
                    const double D = dx*dx + dy*dy + dz*dz;
                    
                    if(nj >= 0 && nj < nr && D >= Dl && D < Dh)
                    {
                        const int k = dm*((ls ? 0.5*log(D) : sqrt(D)) - d0);
                        const int l = k*nr*nr + ni*nr + nj;
                        
                        Z_[k] += 1;
                        N_[l] += 1;
                        s_ += 1;
                    }
                }
            }
        }
        
        #pragma omp critical
        {
            for(j = 0; j < nd; ++j)
            {
                const int k = j*nr*nr;
                Z[j] += Z_[j];
                for(i = 0; i < nr*nr; ++i)
                    N[k+i] += (Z_[j]/Z[j])*(N_[k+i]/Z_[j] - N[k+i]);
            }
            s += s_;
        }
        
        free(qr);
        free(Z_);
        free(N_);
        
        if(tc)
        {
            free(cv_);
            if(cv_ != xv_)
                free(xv_);
            free(ma_);
        }
    }
    
    dt = difftime(time(NULL), st);
    
    if(isatty(fileno(stdin)))
        printf("\r");
    printf("%s done with %.0f pairs", sv, s);
    printf(" in %02d:%02d:%02d  \n", dt/3600, (dt/60)%60, dt%60);
    printf("\n");
    
    output(cfg.output, nd, nr, N);
    
    free(Z);
    free(N);
    free(ma);
    for(j = 0; j < nc; ++j)
        free(rv[j]);
    
    free(cfgfile);
    freecfg(&cfg);
    
    return EXIT_SUCCESS;
    
err_usage:
    fprintf(stderr, "usage: nzz [config] [output] [cat ...]\n");
    return EXIT_FAILURE;
    
err_alloc:
    perror(NULL);
    return EXIT_FAILURE;
}
