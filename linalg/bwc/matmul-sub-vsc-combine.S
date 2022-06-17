#if 0
static inline void matmul_sub_vsc_combine_asm(abobj_ptr x, abt * dst, const abt * * mptrs, const uint8_t * q, unsigned int count)
{
    for( ; count-- ; ) {
        uint8_t c = *q++;
        abadd(x, dst, mptrs[c]);
        mptrs[c] += aboffset(x, c != 0);
        dst += aboffset(x, c == 0);
    }
}
#endif
.section .note.GNU-stack,"",%progbits
	.text
	.p2align 4,,15
.globl matmul_sub_vsc_combine_asm
	.type	matmul_sub_vsc_combine_asm, @function
matmul_sub_vsc_combine_asm:
        pushq %rbp
        pushq %rbx
        pushq %r15
        pushq %r14
        pushq %r13
        pushq %r12
        // arg order: %rdi, %rsi, %rdx, %rcx, %r8, and %r9
        // %rdi dst
        // %rsi src
	// %rdx matrix  --> rbp
        movq    %rdx, %rbp
	// %rcx ncoeffs  --> r8 and r9
        movq    %rcx, %r8
        movq    %rcx, %r9
        test    %r9,%r9
	je	.Lfinished
	.p2align 4,,10
	.p2align 3

#define one_simple_loop							\
	movzbl	(%rbp), %eax                                          ; \
	leaq	(%rsi,%rax,8), %rdx                                   ; \
	movq	(%rdx), %rcx                                          ; \
        movq    (%rdi), %rbx					      ; \
	xorq	(%rcx), %rbx                                          ; \
	movq	%rbx, (%rdi)                                          ; \
        xorq    %rbx, %rbx					      ; \
        testb   %al, %al					      ; \
        setne   %bl						      ; \
        leaq    (%rcx, %rbx, 8), %rcx				      ; \
        movq    %rcx, (%rdx)					      ; \
        xorb    $1, %bl						      ; \
        leaq    (%rdi, %rbx, 8), %rdi

#define ALIGN_MASK      7

.Lunaligned:
        movq    %rbp, %rax
        andq    $ALIGN_MASK, %rax
        je .Laligned
        one_simple_loop
        addq $1, %rbp
        subq $1, %r9
        je .Lfinished
        jmp .Lunaligned

.Laligned:
	leaq	(%rbp,%r9), %r15
        cmpq    $8, %r9
        jb      .Ltail
        subq    $8, %r9
	leaq	(%rbp,%r9), %r15

        movq    (%rbp), %rax
        addq    $8, %rbp
        cmpq    %rbp, %r15
        jb .Lfixup

.Lmain:
        // rax contains eight bytes ahead.

#define inner_ops							\
	leaq	(%rsi,%rcx,8), %rdx                                   ; \
	movq	(%rdx), %r8                                           ; \
        movq    (%rdi), %r9					      ; \
	xorq	(%r8), %r9                                            ; \
	movq	%r9, (%rdi)                                           ; \
        xorq    %rbx, %rbx					      ; \
        testb   %cl, %cl					      ; \
        setne   %bl						      ; \
        leaq    (%r8, %rbx, 8), %r8				      ; \
        movq    %r8, (%rdx)					      ; \
        xorb    $1, %bl						      ; \
        leaq    (%rdi, %rbx, 8), %rdi

        // do this 8 times.
        movzbq  %al, %rcx
        shrq    $8, %rax
        inner_ops
        movzbq  %al, %rcx
        shrq    $8, %rax
        inner_ops
        movzbq  %al, %rcx
        shrq    $8, %rax
        inner_ops
        movzbq  %al, %rcx
        shrq    $8, %rax
        inner_ops

        movzbq  %al, %rcx
        shrq    $8, %rax
        inner_ops
        movzbq  %al, %rcx
        shrq    $8, %rax
        inner_ops

        movzbq  %al, %rcx
        shrq    $8, %rax
        inner_ops
        movzbq  %al, %rcx
        inner_ops

	movq    (%rbp), %rax

        addq    $8, %rbp
        cmpq    %rbp, %r15
        ja      .Lmain
.Lfixup:
        movq    $8, %r9
        leaq    (%r15,%r9,1), %r15
        subq    $8, %rbp

.Ltail:
        one_simple_loop
        addq    $1, %rbp
        cmpq    %rbp, %r15
        jne .Ltail

.Lfinished:
        popq %r12
        popq %r13
        popq %r14
        popq %r15
        popq %rbx
        popq %rbp
	ret
	.size	matmul_sub_vsc_combine_asm, .-matmul_sub_vsc_combine_asm

