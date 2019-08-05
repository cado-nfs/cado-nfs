#ifndef LINGEN_QCODE_BINARY_H_
#define LINGEN_QCODE_BINARY_H_

#ifdef __cplusplus
extern "C" {
#endif

struct lingen_qcode_data_s {
    /* we have a matrix of size m*b, and one of size b*b. The second
     * dimension is not called n in order to avoid confusion with the n
     * in BW ; in fact we have b = m + n */
    unsigned int m, b;
    unsigned int t;
    unsigned long length, outlength;

    unsigned int * local_delta;

    /* where we grab our input and store our output */
    /* Note that we don't own the data corresponding to the innermost
     * level (while we do own the outermost level for iptrs and optrs).
     */
    unsigned int * delta;
    unsigned int * ch;

    const unsigned long ** iptrs;
    unsigned long ** optrs;
};


typedef struct lingen_qcode_data_s lingen_qcode_data[1];
typedef struct lingen_qcode_data_s * lingen_qcode_data_ptr;
typedef const struct lingen_qcode_data_s * lingen_qcode_data_srcptr;

void lingen_qcode_init(lingen_qcode_data_ptr qq, unsigned int m, unsigned int b, unsigned int length, unsigned int outlength);
void lingen_qcode_clear(lingen_qcode_data_ptr qq);

unsigned int lingen_qcode_output_column_length(lingen_qcode_data_srcptr qq, unsigned int j);

unsigned int lingen_qcode_do(lingen_qcode_data_ptr qq);

static inline void lingen_qcode_hook_delta(lingen_qcode_data_ptr qq, unsigned int * delta)
{
    qq->delta = delta;
}

static inline void lingen_qcode_hook_chance_list(lingen_qcode_data_ptr qq, unsigned int * ch)
{
    qq->ch = ch;
}

static inline void lingen_qcode_hook_input(lingen_qcode_data_ptr qq, unsigned int i, unsigned int j, unsigned long * poly)
{
    qq->iptrs[i * qq->b + j] = poly;
}

static inline void lingen_qcode_hook_output(lingen_qcode_data_ptr qq, unsigned int i, unsigned int j, unsigned long * poly)
{
    qq->optrs[i * qq->b + j] = poly;
}

#ifdef __cplusplus
}
#endif

#endif	/* LINGEN_QCODE_BINARY_H_ */
