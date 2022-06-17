#if 0
static inline void matmul_sub_vsc_dispatch_asm(abobj_ptr x, abt * dst, abt const * src, const uint16_t * q, unsigned int count)
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
.globl matmul_sub_vsc_dispatch_asm
        .type   matmul_sub_vsc_dispatch_asm, @function
matmul_sub_vsc_dispatch_asm:
        // arg order: %rdi, %rsi, %rdx, %rcx, %r8, and %r9
        // %rdi dst
        // %rsi src
	// %rdx q
	// %rcx ncoeffs
        testl   %ecx, %ecx
        je      .L134
        subl    $1, %ecx
        leaq    2(%rdx,%rcx,2), %r8
        .p2align 4,,10
        .p2align 3
.L133:
        movzwl  (%rdx), %ecx
        addq    $2, %rdx
        movq    (%rsi,%rcx,8), %rax
        movq    %rax, (%rdi)
        addq    $8, %rdi
        cmpq    %r8, %rdx
        jne     .L133
.L134:
        rep
        ret
        .size   matmul_sub_vsc_dispatch_asm, .-matmul_sub_vsc_dispatch_asm

