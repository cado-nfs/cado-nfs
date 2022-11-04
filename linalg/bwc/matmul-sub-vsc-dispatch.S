#if 0
static inline void matmul_sub_vsc_dispatch(abobj_ptr x, abt * dst, abt const * src, const uint16_t * q, unsigned int count)
{
    for( ; count-- ; ) {
        abcopy(x, dst, src + aboffset(x, *q++), 1);
        dst += aboffset(x, 1);
    }
}
#endif
.section .note.GNU-stack,"",%progbits
	.text
        .p2align 4,,15
.globl matmul_sub_vsc_dispatch
        .type   matmul_sub_vsc_dispatch, @function
matmul_sub_vsc_dispatch:
        testl   %r8d, %r8d
        movl    %r8d, %eax
        je      .L134
        subl    $1, %eax
        leaq    2(%rcx,%rax,2), %r8
        .p2align 4,,10
        .p2align 3
.L133:
        movzwl  (%rcx), %eax
        addq    $2, %rcx
        mov     %eax, %eax
        movq    (%rdx,%rax,8), %rax
        movq    %rax, (%rsi)
        addq    $8, %rsi
        cmpq    %r8, %rcx
        jne     .L133
.L134:
        rep
        ret
        .size   matmul_sub_vsc_dispatch, .-matmul_sub_vsc_dispatch

