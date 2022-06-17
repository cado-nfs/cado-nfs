.section .note.GNU-stack,"",%progbits
	.text
	.p2align 4,,15
.globl matmul_sub_small1_asm
	.type	matmul_sub_small1_asm, @function
matmul_sub_small1_asm:
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

        // at .Lfinished, rbp is shipped as a return value.
	je	.Lfinished

        // we want our reads aligned, we want to unroll our main loop,
        // and leave at least some space in the tail because the tight
        // loop does some readahead.

#define ALIGN_MASK      63

#define one_simple_xor(bufreg1, bufreg2, bufreg3)               \
	movzwq	(%rbp), % ## bufreg1                            ; \
	movzwq	2(%rbp), % ## bufreg2                           ; \
	movq	(%rsi,% ## bufreg1,8), %bufreg3                 ; \
	xorq	%bufreg3, (%rdi,% ## bufreg2,8)                 ; \
        addq    $4, %rbp

.Lunaligned_loop:
        movq    %rbp, %rax
        andq     $ALIGN_MASK, %rax
        je .Laligned
        one_simple_xor(r10, r11, rdx)
        subq    $1, %r9
	je	.Lfinished
        jmp .Lunaligned_loop

.Laligned:
        // now rbp is aligned.

        // four bytes of readahead seem to provide timings a lot more
        // stable than with just two.

#define READAHEAD_COEFFS        4
        // we're reading READAHEAD_COEFFS coeffs ahead of time, so the
        // tail must contain at least that many coefficients.

        leaq    (%rbp,%r9,4),%r15

        cmpq $READAHEAD_COEFFS, %r9

        // the tail loop is just simple and stupid. rbp points to a
        // location where exactly r9 coeffs are to be read. Note that
        // some coefficients are already in registers at this point. But
        // this does not matter much (we could avoid re-reading, but
        // really we don't care).
        jb .Ltail_loop

        // we're going to advance rbp one loop iteration ahead. This
        // eases the handling of the tail.

        // LOOP_LENGTH is the number of bytes read in rbp per turn. 16 or
        // 64 make no difference here.
#define LOOP_LENGTH     16

        subq $READAHEAD_COEFFS, %r9
        // r15+4*READAHEAD_COEFFS : end data pointer.
        leaq    (%rbp,%r9,4),%r15

        movl 0(%rbp), %eax
        movl 4(%rbp), %ebx
        movl 8(%rbp), %ecx
        movl 12(%rbp), %edx
        addq    $LOOP_LENGTH, %rbp

	cmpq	%rbp, %r15
	jb	.Lfixup_before_tail

.Lmain_loop:

#define one_xor(idxreg, bufreg1, bufreg2, offset)               \
        movzwq % ## idxreg, % ## bufreg1                        ; \
        shrq $16, %r ## idxreg                                  ; \
	movq	(%rsi,% ## bufreg1, 8), % ## bufreg2            ; \
	xorq	% ## bufreg2, (%rdi,%r ## idxreg, 8)            ; \
        movl offset(%rbp), %e ## idxreg
        
/*      This expands as (for example) :

        movzwq %ax, %r14
        shrq $16, %rax
	movq	(%rsi,%r14, 8), %r8
	xorq	%r8, (%rdi,%rax,8)
        movl 0(%rbp), %eax
        */

        /*
        one_xor(ax, r14, r8, -48)
        one_xor(bx, r12, r9, -44)
        one_xor(cx, r14, r10, -40)
        one_xor(dx, r12, r11, -36)
        one_xor(ax, r14, r8, -32)
        one_xor(bx, r12, r9, -28)
        one_xor(cx, r14, r10, -24)
        one_xor(dx, r12, r11, -20)
        one_xor(ax, r14, r8, -16)
        one_xor(bx, r12, r9, -12)
        one_xor(cx, r14, r10, -8)
        one_xor(dx, r12, r11, -4)
        */
        one_xor(ax, r14, r8, 0)
        one_xor(bx, r12, r9, 4)
        one_xor(cx, r14, r10, 8)
        one_xor(dx, r12, r11, 12)

        // note that we need 4*READAHEAD_COEFFS bytes of readahead at the
        // end of the array.

	addq	$LOOP_LENGTH, %rbp
	cmpq	%rbp, %r15
	ja	.Lmain_loop

.Lfixup_before_tail:
        // fix r15 for real value
        movq    $READAHEAD_COEFFS, %r9
        leaq    (%r15,%r9,4), %r15
        // also fix rbp
        subq    $LOOP_LENGTH, %rbp

.Ltail_loop:
        cmpq    %rbp, %r15
        je      .Lfinished
        one_simple_xor(r10, r11, rdx)
        jmp .Ltail_loop

.Lfinished:
	movq	%rbp, %rax

        popq %r12
        popq %r13
        popq %r14
        popq %r15
        popq %rbx
        popq %rbp
	ret
	.size	matmul_sub_small1_asm, .-matmul_sub_small1_asm
