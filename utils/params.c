#include "cado.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>              /* isdigit isspace */
#include <limits.h>             /* INT_MIN INT_MAX */
#include <errno.h>
#include <stdarg.h>
#include <inttypes.h>

#include "params.h"
#include "macros.h"

void param_list_init(param_list pl)
{
    pl->size = 0;
    pl->alloc = 16;
    pl->p = (parameter *) malloc(pl->alloc * sizeof(parameter));
    pl->consolidated = 1;
    pl->aliases = NULL;
    pl->naliases = 0;
    pl->naliases_alloc = 0;
    pl->knobs = NULL;
    pl->nknobs = 0;
    pl->nknobs_alloc = 0;
}

void param_list_clear(param_list pl)
{
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        if (pl->p[i]->key) free(pl->p[i]->key);
        free(pl->p[i]->value);
    }
    free(pl->p);
    free(pl->aliases);
    for(int i = 0 ; i < pl->nknobs ; i++) {
        free(pl->knobs[i]->knob);
    }
    free(pl->knobs);
    memset(pl, 0, sizeof(pl));
}

static void make_room(param_list pl, unsigned int more)
{
    if (pl->size + more <= pl->alloc) {
        return;
    }

    for( ; pl->size + more > pl->alloc ; ) {
      pl->alloc += 8 + (pl->alloc >> 1);
    }

    pl->p = (parameter *) realloc(pl->p, pl->alloc * sizeof(parameter));
}

static int param_list_add_key_nostrdup(param_list pl,
        char * key, char * value, enum parameter_origin o)
{
    make_room(pl, 1);
    int r = pl->size;
    pl->p[pl->size]->key = key;
    pl->p[pl->size]->value = value;
    pl->p[pl->size]->origin = o;
    // knobs always count as parsed, of course. Hence the (value==NULL) thing
    pl->p[pl->size]->parsed = (value == NULL);
    // values above 1 are built within the sorting step.
    pl->p[pl->size]->seen = 1;
    pl->size++;
    pl->consolidated = 0;
    return r;
}

int param_list_add_key(param_list pl,
        const char * key, const char * value, enum parameter_origin o)
{
    int r = param_list_add_key_nostrdup(pl,
            key ? strdup(key) : NULL, value ? strdup(value) : NULL, o);
    return r;
}

struct sorting_data {
    parameter_srcptr s;
    int p;
};

typedef int (*sortfunc_t) (const void *, const void *);

int strcmp_or_null(const char * a, const char * b)
{
    if (a == NULL) { return b ? -1 : 0; }
    if (b == NULL) { return a ? 1 : 0; }
    return strcmp(a, b);
}

int paramcmp(const struct sorting_data * a, const struct sorting_data * b)
{
    int r;
    r = strcmp_or_null(a->s->key, b->s->key);
    if (r) return r;
    r = a->s->origin - b->s->origin;
    if (r) return r;
    /* Compare by pointer position so that the sort is stable */
    return a->p - b->p;
}

/* Sort (for searching), and remove duplicates. */
void param_list_consolidate(param_list pl)
{
    if (pl->consolidated) {
        return;
    }
    struct sorting_data * intermediate;
    intermediate = malloc(pl->size * sizeof(struct sorting_data));
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        intermediate[i].p = i;
        intermediate[i].s = pl->p[i];
    }
    qsort(intermediate, pl->size, sizeof(struct sorting_data),
            (sortfunc_t) &paramcmp);

    parameter * np = (parameter *) malloc(pl->alloc * sizeof(parameter));
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        memcpy(np[i], pl->p[intermediate[i].p], sizeof(parameter));
    }
    free(pl->p);
    pl->p = np;
    free(intermediate);

    // now remove duplicates. The sorting has priorities right.
    unsigned int j = 0;
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        if (pl->p[i]->key != NULL && i + 1 < pl->size && strcmp(pl->p[i]->key, pl->p[i+1]->key) == 0) {
            /* The latest pair in the list is the one having highest
             * priority. Do we don't do the copy at this moment.
             */
            free(pl->p[i]->key);
            free(pl->p[i]->value);
            // this value is useful for knobs
            pl->p[i+1]->seen += pl->p[i]->seen;
        } else {
            // I can't see why there could conceivably be a problem if i == j,
            // but valgrind complains...
            if (i != j) {
                memcpy(pl->p[j], pl->p[i], sizeof(parameter));
            }
            j++;
        }
    }
    pl->size = j;

    pl->consolidated = 1;
}

void param_list_remove_key(param_list pl, const char * key)
{
    unsigned int j = 0;
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        if (strcmp_or_null(pl->p[i]->key, key) == 0) {
            if (pl->p[i]->key) free(pl->p[i]->key);
            free(pl->p[i]->value);
        } else {
            if (i != j) {
                memcpy(pl->p[j], pl->p[i], sizeof(parameter));
            }
            j++;
        }
    }
    pl->size = j;

    pl->consolidated = 1;
}


int param_list_read_stream(param_list pl, FILE *f)
{
    int all_ok=1;
    const int linelen = 512;
    char line[linelen];
    char * newkey;
    char * newvalue;
    while (!feof(f)) {
        if (fgets(line, linelen, f) == NULL)
            break;
        if (line[0] == '#')
            continue;
        // remove possible comment at end of line.
        char * hash;
        if ((hash = strchr(line, '#')) != NULL) {
            *hash = '\0';
        }

        char * p = line;

        // trailing space
        int l = strlen(p);
        for( ; l && isspace(p[l-1]) ; l--);
        p[l] = '\0';

        // leading space.
        for( ; *p && isspace(*p) ; p++, l--);

        // empty ps are ignored.
        if (l == 0)
            continue;

        // look for a left-hand-side.
        l = 0;
        if (!(isalpha(p[l]) || p[l] == '_' || p[l] == '-')) {
            param_list_add_key(pl, NULL, line, PARAMETER_FROM_FILE);
            continue;
        }
        for( ; p[l] && (isalnum(p[l]) || p[l] == '_' || p[l] == '-') ; l++);

        int lhs_length = l;

        if (lhs_length == 0) {
            fprintf(stderr, "Parse error, no usable key for config line:\n%s\n",
                    line);
            all_ok=0;
            continue;
        }

        /* Now we can match (whitespace*)(separator)(whitespace*)(data)
         */
        char * q = p + lhs_length;
        for( ; *q && isspace(*q) ; q++);

        /* match separator, which is one of : = := */
        if (*q == '=') {
            q++;
        } else if (*q == ':') {
            q++;
            if (*q == '=')
                q++;
        } else if (q == p + lhs_length) {
            fprintf(stderr, "Parse error, no separator for config line:\n%s\n",
                    line);
            all_ok=0;
            continue;
        }
        for( ; *q && isspace(*q) ; q++);

        newkey = malloc(lhs_length + 1);
        memcpy(newkey, p, lhs_length);
        newkey[lhs_length]='\0';

        newvalue = strdup(q);

        param_list_add_key_nostrdup(pl, newkey, newvalue, PARAMETER_FROM_FILE);
    }
    param_list_consolidate(pl);

    return all_ok;
}

int param_list_read_file(param_list pl, const char * name)
{
    FILE * f;
    f = fopen(name, "r");
    if (f == NULL) {
        fprintf(stderr, "Cannot read %s\n", name);
        exit(1);
    }
    int r = param_list_read_stream(pl, f);
    fclose(f);
    return r;
}

int param_list_configure_alias(param_list pl, const char * key, const char * alias)
{
    size_t len = strlen(alias);

    ASSERT_ALWAYS(alias != NULL);
    ASSERT_ALWAYS(key != NULL);
    /* A knob may be aliases, but only as another knob !!! */
    ASSERT_ALWAYS(key[0] != '-' || (alias[0] == '-' && alias[len-1] != '='));

    if (alias[0] != '-' && alias[len-1] != '=') {
        /* Then, accept both the --xxx and xxx= forms */
        char * tmp;
        tmp = malloc(len + 4);
        snprintf(tmp, len+4, "--%s", alias);
        param_list_configure_alias(pl, key, strdup(tmp));
        snprintf(tmp, len+4, "%s=", alias);
        param_list_configure_alias(pl, key, strdup(tmp));
        free(tmp);
        return 0;
    }

    if (pl->naliases == pl->naliases_alloc) {
        pl->naliases_alloc += 1;
        pl->naliases_alloc <<= 1;
        pl->aliases = realloc(pl->aliases, pl->naliases_alloc * sizeof(param_list_alias));
    }
    pl->aliases[pl->naliases]->alias = alias;
    pl->aliases[pl->naliases]->key = key;
    pl->naliases++;
    return 0;
}

int param_list_configure_knob(param_list pl, const char * knob, int * ptr)
{
    ASSERT_ALWAYS(knob != NULL);
    if ((pl->nknobs + 1) >= pl->nknobs_alloc) {
        pl->nknobs_alloc += 2;
        pl->nknobs_alloc <<= 1;
        pl->knobs = realloc(pl->knobs, pl->nknobs_alloc * sizeof(param_list_knob));
    }
    char * tmp = (char *) malloc(strlen(knob)+2);
    tmp[0]='-';
    if (knob[1] == '-') { // have -- in the knob
        strncpy(tmp, knob, strlen(knob) + 1);
    } else {
        strncpy(tmp+1, knob, strlen(knob) + 1);
    }
    // put the -- version
    pl->knobs[pl->nknobs]->knob = strdup(tmp);
    pl->knobs[pl->nknobs]->ptr = ptr;
    if (ptr) *ptr = 0;
    pl->nknobs++;
    // put the - version
    pl->knobs[pl->nknobs]->knob = strdup(tmp+1);
    pl->knobs[pl->nknobs]->ptr = ptr;
    if (ptr) *ptr = 0;
    pl->nknobs++;
    free(tmp);

    return 0;
}


static int param_list_update_cmdline_alias(param_list pl,
        param_list_alias al,
        int * p_argc, char const *** p_argv)
{
    const char * a = (*p_argv[0]);
    if (al->alias[strlen(al->alias)-1] == '=') {
        // since knobs are aliased only by knobs, we know we have a plain
        // option here.
        if (strncmp(a, al->alias, strlen(al->alias)) != 0)
            return 0;
        a += strlen(al->alias);
        if (a[1] == '\0')
            return 0;
        // no +1 , since the alias contains the = sign already.
        param_list_add_key(pl, al->key, a, PARAMETER_FROM_CMDLINE);
        (*p_argv)+=1;
        (*p_argc)-=1;
        return 1;
    }
    if (strcmp(a, al->alias) == 0) {
        if (al->key[0] == '-') {
            /* This is a knob ; we have to treat it accordingly. The
             * difficult part is to properly land on
             * param_list_update_cmdline_knob at the proper time. This
             * means in particular not necessarily there. It's
             * considerably easier to simply change the value in the
             * command line.
             */
            (*p_argv)[0] = al->key;
            /* leave argv and argc unchanged. */
            return 0;
        }
        (*p_argv)+=1;
        (*p_argc)-=1;
        if (*p_argc == 0) {
            fprintf(stderr, "Option %s requires an argument\n", a);
            exit(1);
        }
        param_list_add_key(pl, al->key, (*p_argv[0]), PARAMETER_FROM_CMDLINE);
        (*p_argv)+=1;
        (*p_argc)-=1;
        return 1;
    }
    return 0;
}

static int param_list_update_cmdline_knob(param_list pl,
        param_list_knob knob,
        int * p_argc, char const *** p_argv)
{
    const char * a = (*p_argv[0]);
    if (strcmp(a, knob->knob) == 0) {
        int r;
        r = param_list_add_key(pl, knob->knob, NULL, PARAMETER_FROM_CMDLINE);
        (*p_argv)+=1;
        (*p_argc)-=1;
        if (knob->ptr) (*(knob->ptr))++;
        return 1;
    }
    return 0;
}

int param_list_update_cmdline(param_list pl,
        int * p_argc, char *** pp_argv)
{
    char const *** p_argv = (char const ***) pp_argv;
    if (*p_argc == 0)
        return 0;

    int i;

    /* We rely on having alias scanning first, because this incurs a
     * command line changed (could get along without except in the case
     * of knobs where it's particularly handy).
     */
    for(i = 0 ; i < pl->naliases ; i++) {
        if (param_list_update_cmdline_alias(pl, pl->aliases[i], p_argc, p_argv))
            return 1;
    }

    for(i = 0 ; i < pl->nknobs ; i++) {
        if (param_list_update_cmdline_knob(pl, pl->knobs[i], p_argc, p_argv))
            return 1;
    }

    const char * a = (*p_argv[0]);
    if (*p_argc >= 2 && a[0] == '-') {
        a++;
        a+= *a == '-';
        int x=0;
        /* parameters _must_ begin by alphabetic characters, or
         * otherwise we suffer to distinguish immediate numerical
         * entries.
         */
        if (!isalpha(a[x]))
            return 0;
        for( ; a[x] && (isalnum(a[x]) || a[x] == '_' || a[x] == '-') ; x++);
        if (a[x] == '\0') {
            param_list_add_key(pl, a, (*p_argv)[1], PARAMETER_FROM_CMDLINE);
            (*p_argv)+=2;
            (*p_argc)-=2;
            return 1;
        }
    } else {
        /* Check for <key>=<value> syntax */
        int x=0;
        for( ; a[x] && (isalnum(a[x]) || a[x] == '_' || a[x] == '-') ; x++);
        if (a[x] == '=' && a[x+1]) {
            char * newkey = malloc(x+1);
            memcpy(newkey, a, x);
            newkey[x]='\0';
            char * newvalue = strdup(a+x+1);
            param_list_add_key_nostrdup(pl,
                    newkey, newvalue, PARAMETER_FROM_CMDLINE);
            (*p_argv)+=1;
            (*p_argc)-=1;
            return 1;
        }
    }
    return 0;
}

int param_strcmp(const char * a, parameter_srcptr b)
{
    return strcmp_or_null(a, b->key);
}

static int assoc(param_list pl, const char * key)
{
    void * found;

    param_list_consolidate(pl);
    found = bsearch(key, pl->p, pl->size,
            sizeof(parameter), (sortfunc_t) param_strcmp);
    if (found == NULL)
        return -1;
    parameter * c = (parameter *) found;
    return c-pl->p;
}

int param_list_parse_long(param_list pl, const char * key, long * r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    long res;
    res = strtol(value, &end, 0);
    if (*end != '\0') {
        fprintf(stderr, "Parse error: parameter for key %s is not a long: %s\n",
                key, value);
        exit(1);
    }
    if (r)
        *r = res;
    return pl->p[v]->seen;
}

int param_list_parse_int(param_list pl, const char * key, int * r)
{
    long res;
    int rc;
    rc = param_list_parse_long(pl, key, &res);
    if (rc == 0)
        return 0;
    if (res > INT_MAX || res < INT_MIN) {
        fprintf(stderr, "Parse error:"
                " parameter for key %s does not fit within an int: %ld\n",
                key, res);
        exit(1);
    }
    if (r)
        *r  = res;
    return rc;
}

int param_list_parse_int_and_int(param_list pl, const char * key, int * r, const char * sep)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    long res[2];
    res[0] = strtol(value, &end, 0);
    if (strncmp(end, sep, strlen(sep)) != 0) {
        fprintf(stderr, "Parse error: parameter for key %s"
                " must match %%d%s%%d; got %s\n",
                key, sep, pl->p[v]->value);
        exit(1);
    }
    value = end + strlen(sep);
    res[1] = strtol(value, &end, 0);
    if (*end != '\0') {
        fprintf(stderr, "Parse error: parameter for key %s"
                " must match %%d%s%%d; got %s\n",
                key, sep, pl->p[v]->value);
        exit(1);
    }
    if (r) {
        r[0] = res[0];
        r[1] = res[1];
    }
    return pl->p[v]->seen;
}

int param_list_parse_intxint(param_list pl, const char * key, int * r)
{
    return param_list_parse_int_and_int(pl, key, r, "x");
}


int param_list_parse_ulong(param_list pl, const char * key, unsigned long * r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    unsigned long res;
    res = strtoul(value, &end, 0);
    if (*end != '\0') {
        fprintf(stderr, "Parse error:"
                " parameter for key %s is not an ulong: %s\n",
                key, value);
        exit(1);
    }
    if (r)
        *r = res;
    return pl->p[v]->seen;
}

int param_list_parse_size_t(param_list pl, const char * key, size_t * r)
{
    unsigned long t;
    int res;
    res = param_list_parse_ulong(pl, key, &t);
    if (res && r) { *r = t; }
    return res;
}


int param_list_parse_uint64(param_list pl, const char * key, uint64_t * r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    uint64_t res;
    res = strtoumax(value, &end, 0);
    if (*end != '\0') {
        fprintf(stderr, "Parse error:"
                " parameter for key %s is not an ulong: %s\n",
                key, value);
        exit(1);
    }
    if (r)
        *r = res;
    return pl->p[v]->seen;
}

int param_list_parse_uint(param_list pl, const char * key, unsigned int * r)
{
    unsigned long res;
    int rc = param_list_parse_ulong(pl, key, &res);
    if (rc == 0)
        return 0;
    if (res > UINT_MAX) {
        fprintf(stderr, "Parse error:"
                " parameter for key %s does not fit within an uinsigned int: %ld\n",
                key, res);
        exit(1);
    }
    if (r)
        *r  = res;
    return rc;
}


int param_list_parse_double(param_list pl, const char * key, double * r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    double res;
    res = strtod(value, &end);
    if (*end != '\0') {
        fprintf(stderr, "Parse error: parameter for key %s is not an int: %s\n",
                key, value);
        exit(1);
    }
    if (r)
        *r = res;
    return pl->p[v]->seen;
}

int param_list_parse_string(param_list pl, const char * key, char * r, size_t n)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    pl->p[v]->parsed=1;
    char * value = pl->p[v]->value;
    if (r && strlen(value) > n-1) {
        fprintf(stderr, "Parse error:"
                " parameter for key %s does not fit within string buffer"
                " of length %lu\n", key, (unsigned long) n);
        exit(1);
    }
    if (r)
        strncpy(r, value, n);
    return pl->p[v]->seen;
}

int param_list_parse_int_list(param_list pl, const char * key, int * r, size_t n, const char * sep)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    char * value = pl->p[v]->value;
    pl->p[v]->parsed=1;
    char * end;
    int * res = malloc(n * sizeof(int));
    memset(res, 0, n * sizeof(int));
    size_t parsed = 0;
    for( ;; ) {
        res[parsed] = strtol(value, &end, 0);
        if (parsed++ == n)
            break;
        if (parsed && *end == '\0')
            break;
        if (strncmp(end, sep, strlen(sep)) != 0) {
            fprintf(stderr, "Parse error: parameter for key %s"
                    " must match %%d(%s%%d)*; got %s\n",
                    key, sep, pl->p[v]->value);
            exit(1);
        }
        value = end + strlen(sep);
    }
    if (*end != '\0') {
        fprintf(stderr, "Parse error: parameter for key %s"
                " must match %%d(%s%%d)*; got %s\n",
                key, sep, pl->p[v]->value);
        exit(1);
    }
    if (r) {
        memcpy(r, res, n * sizeof(int));
    }
    free(res);
    return parsed;
}

const char * param_list_lookup_string(param_list pl, const char * key)
{
    int v = assoc(pl, key);
    if (v < 0)
        return NULL;
    pl->p[v]->parsed=1;
    const char * value = pl->p[v]->value;
    return value;
}

int param_list_parse_mpz(param_list pl, const char * key, mpz_ptr r)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    pl->p[v]->parsed=1;
    unsigned int nread;
    int rc;
    char * value = pl->p[v]->value;
    if (r) {
        rc = gmp_sscanf(value, "%Zd%n", r, &nread);
    } else {
        /* scan even when the result is not wanted */
        rc = gmp_sscanf(value, "%*Zd%n", &nread);
    }
    if (rc != 1 || value[nread] != '\0') {
        fprintf(stderr, "Parse error: parameter for key %s is not an mpz: %s\n",
                key, value);
        exit(1);
    }
    return pl->p[v]->seen;
}

int param_list_parse_knob(param_list pl, const char * key)
{
    int v = assoc(pl, key);
    if (v < 0)
        return 0;
    pl->p[v]->parsed=1;
    if (pl->p[v]->value != NULL) {
        fprintf(stderr, "Parse error: option %s accepts no argument\n", key);
        exit(1);
    }
    return pl->p[v]->seen;
}

int param_list_all_consumed(param_list pl, char ** extraneous)
{
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        if (!pl->p[i]->parsed) {
            pl->p[i]->parsed = 1;
            if (extraneous) {
                *extraneous = pl->p[i]->key;
            }
            return 0;
        }
    }
    return 1;
}

int param_list_warn_unused(param_list pl)
{
    int u = 0;
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        if (pl->p[i]->origin != PARAMETER_FROM_FILE && !pl->p[i]->parsed) {
            fprintf(stderr, "Warning: unused command-line parameter %s\n",
                    pl->p[i]->key);
            u++;
        }
    }
    return u;
}

void param_list_display(param_list pl, FILE *f)
{
    param_list_consolidate(pl);
    for(unsigned int i = 0 ; i < pl->size ; i++) {
        fprintf(f,"%s=%s\n", pl->p[i]->key, pl->p[i]->value);
    }
}

void param_list_save(param_list pl, const char * filename)
{
    FILE * f = fopen(filename, "w");
    if (f == NULL) {
        fprintf(stderr, "fopen(%s): %s\n", filename, strerror(errno));
        exit(1);
    }

    param_list_display(pl, f);
    fclose(f);
}

int param_list_save_parameter(param_list pl, enum parameter_origin o, 
        const char * key, const char * format, ...)
{
    va_list ap;
    va_start(ap, format);

    char * tmp;
    int rc;
    rc = vasprintf(&tmp, format, ap);
    param_list_add_key(pl, key, tmp, o);
    free(tmp);

    return rc;
}
