#include "cado.h" // IWYU pragma: keep

#include <cstdio>
#include <cstdlib>

#include "fm.hpp"
#include "macros.h"
#include "utils_cxx.hpp"

fm_t * fm_create()
{
    fm_t * t = (fm_t *)malloc(sizeof(*t));
    ASSERT_ALWAYS(t != NULL);
    t->len_method = 4;
    t->len_proba = 1;
    t->len_time = 1;
    t->method = (unsigned long *)calloc(t->len_method, sizeof(unsigned long));
    ASSERT_ALWAYS(t->method != NULL);
    t->proba = (double *)calloc(t->len_proba, sizeof(double));
    ASSERT_ALWAYS(t->proba != NULL);
    t->time = (double *)calloc(t->len_time, sizeof(double));
    ASSERT_ALWAYS(t->time != NULL);
    t->len_p_min = 0;
    return t;
}

void fm_free(fm_t * t)
{
    if (t != NULL) {
        free(t->time);
        free(t->proba);
        free(t->method);
        free(t);
    }
}

unsigned long const * fm_get_method(fm_t const * t)
{
    return t->method;
}

double const * fm_get_proba(fm_t const * t)
{
    return t->proba;
}

double const * fm_get_time(fm_t const * t)
{
    return t->time;
}

unsigned int fm_get_len_method(fm_t const * t)
{
    return t->len_method;
}

unsigned int fm_get_len_proba(fm_t const * t)
{
    return t->len_proba;
}

unsigned int fm_get_len_time(fm_t const * t)
{
    return t->len_time;
}

unsigned int fm_get_len_p_min(fm_t const * t)
{
    return t->len_p_min;
}

void fm_set_method(fm_t * t, unsigned long const * value, unsigned int len)
{
    if (len != t->len_method) { // realloc
        checked_realloc(t->method, len);
        t->len_method = len;
    }

    for (unsigned int i = 0; i < t->len_method; i++)
        t->method[i] = value[i];
}

void fm_set_proba(fm_t * t, double const * value, unsigned int len,
                  unsigned int len_p_min)
{
    t->len_p_min = len_p_min;
    if (len != t->len_proba) { // realloc
        checked_realloc(t->proba, len);
        t->len_proba = len;
    }

    for (unsigned int i = 0; i < t->len_proba; i++)
        t->proba[i] = value[i];
}

void fm_set_time(fm_t * t, double const * value, unsigned int len)
{
    if (len == 0)
        return;

    if (len != t->len_time) { // realloc
        checked_realloc(t->time, len);
        t->len_time = len;
    }

    for (unsigned int i = 0; i < t->len_time; i++)
        t->time[i] = value[i];
}

fm_t * fm_copy(fm_t const * t)
{
    fm_t * cop = fm_create();
    fm_set_method(cop, t->method, t->len_method);
    fm_set_proba(cop, t->proba, t->len_proba, t->len_p_min);
    fm_set_time(cop, t->time, t->len_time);
    return cop;
}

void fm_put_zero(fm_t * t)
{
    t->method[2] = 0; // B1
    t->method[3] = 0; // B2
    for (unsigned int i = 0; i < t->len_proba; i++)
        t->proba[i] = 0;
    for (unsigned int i = 0; i < t->len_time; i++)
        t->time[i] = 0;
}

bool fm_is_zero(fm_t const * t)
{
    return (t->method[2] == 0 && t->method[3] == 0);
}

int fm_is_equal(fm_t const * c1, fm_t const * c2)
{
    const unsigned int len = c1->len_method;
    for (unsigned int i = 0; i < len; i++)
        if (c1->method[i] != c2->method[i])
            return false;
    return true;
}

int fm_print(fm_t const * t)
{
    return fm_fprint(stdout, t);
}

int fm_fprint(FILE * file, fm_t const * elem)
{
    if (file == nullptr)
        return -1;

    unsigned long const * method = fm_get_method(elem);
    const unsigned int len_method = fm_get_len_method(elem);
    for (unsigned int i = 0; i < len_method; i++)
        fprintf(file, "%lu ", method[i]);
    fputs("| ", file);

    fprintf(file, "%d ", elem->len_p_min);
    double const * proba = fm_get_proba(elem);
    const unsigned int len_proba = fm_get_len_proba(elem);
    for (unsigned int i = 0; i < len_proba; i++)
        fprintf(file, "%lf ", proba[i]);
    fputs("| ", file);

    double const * time = fm_get_time(elem);
    const unsigned int len_time = fm_get_len_time(elem);
    for (unsigned int i = 0; i < len_time; i++)
        fprintf(file, "%lf ", time[i]);

    fputs("|\n", file);
    return 0;
}
