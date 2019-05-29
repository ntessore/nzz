#include <string.h>

enum { SPACING_LIN, SPACING_LOG, NUM_SPACING };
char* CFG_SPACING[] = { "lin", "log" };
char* PRN_SPACING[] = { "linear", "logarithmic" };

enum { COORDS_FLAT, COORDS_LONLAT, NUM_COORDS };
char* CFG_COORDS[] = { "flat", "lonlat" };
char* PRN_COORDS[] = { "flat", "lon, lat" };

enum { TDATA_SHARE, TDATA_COPY, NUM_TDATA };
char* CFG_TDATA[] = { "share", "copy" };
char* PRN_TDATA[] = { "share", "copy" };

enum { UNIT_RAD, UNIT_DEG, UNIT_ARCMIN, UNIT_ARCSEC, NUM_UNIT };
char* CFG_UNIT[] = { "rad", "deg", "arcmin", "arcsec" };
char* PRN_UNIT[] = { "radian", "degree", "arcmin", "arcsec" };

const double UCONV[] = {
    1.0,
    0.017453292519943295769,
    0.00029088820866572159615,
    0.0000048481368110953599359
};

#ifndef LINELEN
#define LINELEN 1024
#endif

struct config {
    int catc;
    char* catv[2];
    int units;
    int coords;
    char* output;
    int nth;
    double thmin;
    double thmax;
    int thunit;
    int spacing;
    int nz;
    double zmin;
    double zmax;
    int num_threads;
    int thread_data;
};

static int findstr(const char* s, int n, char* v[])
{
    int i;
    for(i = 0; i < n; ++i)
        if(strcmp(s, v[i]) == 0)
            break;
    return i;
}

void readcfg(const char* f, struct config* cfg)
{
    FILE* fp;
    char buf[LINELEN];
    size_t l;
    char* key;
    char* val;
    
    fp = fopen(f, "r");
    if(!fp)
        goto err_fopen;
    
    for(l = 1; fgets(buf, sizeof(buf), fp); ++l)
    {
        key = strtok(buf, " \t\r\n");
        val = strtok(NULL, " \t\r\n");
        
        if(!key || *key == '#')
            continue;
        if(!val || *val == '#')
            goto err_no_value;
        
        if(strcmp(key, "cat") == 0)
        {
            if(cfg->catc == 2)
                goto err_bad_key;
            cfg->catv[cfg->catc++] = strdup(val);
        }
        else if(strcmp(key, "units") == 0)
        {
            cfg->units = findstr(val, NUM_UNIT, CFG_UNIT);
            if(cfg->units == NUM_UNIT)
                goto err_bad_value;
        }
        else if(strcmp(key, "coords") == 0)
        {
            cfg->coords = findstr(val, NUM_COORDS, CFG_COORDS);
            if(cfg->coords == NUM_COORDS)
                goto err_bad_value;
        }
        else if(strcmp(key, "output") == 0)
            cfg->output = strdup(val);
        else if(strcmp(key, "nth") == 0)
        {
            cfg->nth = atoi(val);
            if(cfg->nth <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "thmin") == 0)
        {
            cfg->thmin = atof(val);
            if(cfg->thmin <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "thmax") == 0)
        {
            cfg->thmax = atof(val);
            if(cfg->thmax <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "thunit") == 0)
        {
            cfg->thunit = findstr(val, NUM_UNIT, CFG_UNIT);
            if(cfg->thunit == NUM_UNIT)
                goto err_bad_value;
        }
        else if(strcmp(key, "spacing") == 0)
        {
            cfg->spacing = findstr(val, NUM_SPACING, CFG_SPACING);
            if(cfg->spacing == NUM_SPACING)
                goto err_bad_value;
        }
        else if(strcmp(key, "nz") == 0)
        {
            cfg->nz = atoi(val);
            if(cfg->nz <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "zmin") == 0)
        {
            cfg->zmin = atof(val);
            if(cfg->zmin <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "zmax") == 0)
        {
            cfg->zmax = atof(val);
            if(cfg->zmax <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "num_threads") == 0)
        {
            cfg->num_threads = atoi(val);
            if(cfg->num_threads <= 0)
                goto err_bad_value;
        }
        else if(strcmp(key, "thread_data") == 0)
        {
            cfg->thread_data = findstr(val, NUM_TDATA, CFG_TDATA);
            if(cfg->thread_data == NUM_TDATA)
                goto err_bad_value;
        }
        else
            goto err_bad_key;
    }
    
    fclose(fp);
    
    if(cfg->catc == 0)
        { key = "cat"; goto err_missing_key; }
    if(!cfg->output)
        { key = "output"; goto err_missing_key; }
    if(!cfg->nth)
        { key = "nth"; goto err_missing_key; }
    if(!cfg->thmin)
        { key = "thmin"; goto err_missing_key; }
    if(!cfg->thmax)
        { key = "thmax"; goto err_missing_key; }
    if(!cfg->nz)
        { key = "nz"; goto err_missing_key; }
    if(!cfg->zmin)
        { key = "zmin"; goto err_missing_key; }
    if(!cfg->zmax)
        { key = "zmax"; goto err_missing_key; }
    
    return;
    
err_fopen:
    perror(f);
    exit(EXIT_FAILURE);
    
err_bad_key:
    fprintf(stderr, "error: %s:%zu: invalid key `%s`\n", f, l, key);
    exit(EXIT_FAILURE);
    
err_missing_key:
    fprintf(stderr, "error: %s: missing required `%s` key\n", f, key);
    exit(EXIT_FAILURE);
    
err_no_value:
    fprintf(stderr, "error: %s:%zu: missing value\n", f, l);
    exit(EXIT_FAILURE);
    
err_bad_value:
    fprintf(stderr, "error: %s:%zu: invalid value `%s`\n", f, l, val);
    exit(EXIT_FAILURE);
}

void freecfg(struct config* cfg)
{
    for(int i = 0; i < cfg->catc; ++i)
        free(cfg->catv[i]);
    free(cfg->output);
}

void printcfg(const struct config* cfg)
{
    printf("catalog ......... %s%s\n", cfg->catc > 1 ? "0: " : "", *cfg->catv);
    for(int i = 1; i < cfg->catc; ++i)
        printf("                  %d: %s\n", i, cfg->catv[i]);
    printf("units ........... %s\n", PRN_UNIT[cfg->units]);
    printf("coordinates ..... %s\n", PRN_COORDS[cfg->coords]);
    printf("\n");
    printf("output file ..... %s\n", cfg->output);
    printf("num. points ..... %u\n", cfg->nth);
    printf("point range ..... %lg to %lg %s\n", cfg->thmin, cfg->thmax,
                                                    PRN_UNIT[cfg->thunit]);
    printf("point spacing ... %s\n", PRN_SPACING[cfg->spacing]);
    printf("num. redshifts .. %u\n", cfg->nz);
    printf("redshift range .. %lg to %lg\n", cfg->zmin, cfg->zmax);
#ifdef _OPENMP
    printf("\n");
    printf("num. threads .... %d\n", omp_get_max_threads());
    printf("thread data ..... %s\n", PRN_TDATA[cfg->thread_data]);
#endif
}

int readc(const char* f, int co, double ui, int* n, double** c)
{
    FILE* fp;
    int l;
    char buf[LINELEN];
    int i, a;
    double* d;
    char* sx;
    char* sy;
    char* sr;
    char* sw;
    double x, y, r, w;
    
    if(strcmp(f, "-") == 0)
        fp = stdin;
    else
    {
        fp = fopen(f, "r");
        if(!fp)
            goto err_fopen;
    }
    
    i = *n;
    a = 2*i + 1;
    d = realloc(*c, a*RW*sizeof(double));
    if(!d)
        goto err_alloc;
    
    for(l = 1; fgets(buf, sizeof buf, fp); ++l)
    {
        sx = strtok(buf, " \t\r\n");
        sy = strtok(NULL, " \t\r\n");
        sr = strtok(NULL, " \t\r\n");
        sw = strtok(NULL, " \t\r\n");
        
        if(!sx || *sx == '#')
            continue;
        
        if(!sy || *sy == '#')
        {
            fprintf(stderr, "error: %s:%d: missing `y` value\n", f, l);
            exit(EXIT_FAILURE);
        }
        if(!sr || *sr == '#')
        {
            fprintf(stderr, "error: %s:%d: missing `r` value\n", f, l);
            exit(EXIT_FAILURE);
        }
        if(!sw || *sw == '#')
        {
            fprintf(stderr, "error: %s:%d: missing `w` value\n", f, l);
            exit(EXIT_FAILURE);
        }
        
        x = atof(sx)*ui;
        y = atof(sy)*ui;
        r = atof(sr);
        w = atof(sw);
        
        switch(co)
        {
        case COORDS_FLAT:
            d[i*RW+0] = 1;
            d[i*RW+1] = x;
            d[i*RW+2] = y;
            break;
        case COORDS_LONLAT:
            d[i*RW+0] = cos(x)*cos(y);
            d[i*RW+1] = sin(x)*cos(y);
            d[i*RW+2] = sin(y);
            break;
        }
        
        d[i*RW+3] = r;
        d[i*RW+4] = w;
        
        i += 1;
        
        if(i == a)
        {
            a *= 2;
            d = realloc(d, a*RW*sizeof(double));
            if(!d)
                goto err_alloc;
        }
    }
    
    if(fp != stdin)
        fclose(fp);
    
    d = realloc(d, i*RW*sizeof(double));
    if(!d)
        goto err_alloc;
    
    *n = i;
    *c = d;
    
    return 0;
    
err_fopen:
    perror(f);
    exit(EXIT_FAILURE);
    
err_alloc:
    perror(NULL);
    exit(EXIT_FAILURE);
}

void output(const char* f, int n, int m, const double* N)
{
    FILE* fp;
    int i, j, k;
    
    if(strcmp(f, "-") == 0)
        fp = stdout;
    else
    {
        fp = fopen(f, "w");
        if(!fp)
            goto err_fopen;
    }

    for(k = 0; k < n; ++k)
    {
        for(j = 0; j < m; ++j)
        {
            for(i = 0; i < m; ++i)
                fprintf(fp, " % .18e", N[k*m*m+j*m+i]);
            fprintf(fp, "\n");
        }
        fprintf(fp, "\n");
    }
    
    if(fp != stdout)
        fclose(fp);
    
    return;
    
err_fopen:
    perror(f);
    exit(EXIT_FAILURE);
}
