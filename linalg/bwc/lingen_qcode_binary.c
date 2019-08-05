#include "cado.h"
#include "lingen_qcode_binary.h"
#include "utils.h"

void lingen_qcode_init(lingen_qcode_data_ptr qq, unsigned int m, unsigned int b, unsigned int length, unsigned int outlength)
{
    qq->m = m;
    qq->b = b;
    qq->length = length;
    qq->outlength = outlength;
    qq->iptrs = malloc(m * b * sizeof(unsigned long *));
    qq->optrs = malloc(b * b * sizeof(unsigned long *));
    memset(qq->iptrs, 0, m * b * sizeof(unsigned long *));
    memset(qq->optrs, 0, b * b * sizeof(unsigned long *));

    qq->local_delta = malloc(b * sizeof(unsigned int));
    memset(qq->local_delta, 0, b * sizeof(unsigned int));
}
void lingen_qcode_clear(lingen_qcode_data_ptr qq)
{
    free(qq->iptrs);
    free(qq->optrs);
    free(qq->local_delta);
    memset(qq, 0, sizeof(*qq));
}

#if 0
unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq)
{
    /* It's not technically necessary, but the typical use case is really
     * like this, so I'm taking the liberty to shorten the code a bit
     * using this fact */
    ASSERT_ALWAYS(qq->length <= ULONG_BITS);

    /* read the complete input */
    {
        unsigned long tmask = 1;
        for(unsigned int k = 0 ; k < qq->length ; k++, tmask <<= 1) {
            unsigned long jmask = 0;
            for(unsigned int j = 0 ; j < qq->b ; j++, jmask <<= 1) {
                unsigned long * Ac = qq->A[k] + j / ULONG_BITS;
                if (!jmask) jmask = 1;
                for(unsigned int i = 0 ; i < qq->m ; i++) {
                    if (qq->iptrs[i * qq->b + j][0] & tmask) {
                        *Ac |= jmask;
                    }
                    Ac += iceildiv(qq->b, ULONG_BITS);
                }
            }
        }
    }

    /* please implement me ! */
    abort();


    /* Now, store the complete output */
    {
        unsigned long tmask = 1;
        for(unsigned int k = 0 ; k < qq->outlength ; k++, tmask <<= 1) {
            unsigned long jmask = 0;
            for(unsigned int j = 0 ; j < qq->b ; j++, jmask <<= 1) {
                unsigned long * Xc = qq->X[k] + j / ULONG_BITS;
                if (!jmask) jmask = 1;
                for(unsigned int i = 0 ; i < qq->b ; i++) {
                    if (*Xc & jmask) {
                        qq->optrs[i * qq->b + j][0] |= tmask;
                    }
                    Xc += iceildiv(qq->b, ULONG_BITS);
                }
            }
        }
    }
    return qq->length;
}
#endif

unsigned int lingen_qcode_output_column_length(lingen_qcode_data_srcptr qq, unsigned int j)
{
    return qq->local_delta[j];
}

