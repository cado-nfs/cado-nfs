.section .note.GNU-stack,"",%progbits
	.text
	.p2align 4,,15
.globl matmul_sub_small2_asm
	.type	matmul_sub_small2_asm, @function
matmul_sub_small2_asm:
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

// first arg must be one of ax bx cx dx ; second: r[8-15] ; third: any
#define one_simple_xor(bufreg1, bufreg2, bufreg3)               \
	movzwl	(%rbp), %e ## bufreg1                           ; \
        movl    %e ## bufreg1, % ## bufreg2 ## d                ; \
        andl    $4095, % ## bufreg2 ## d                        ; \
        shrq    $12, %r ## bufreg1                              ; \
        leaq    (%rsi, %r ## bufreg1, 8), %rsi                  ; \
	movq	(%rsi), % ## bufreg3                            ; \
	xorq	% ## bufreg3, (%rdi,% ## bufreg2,8)             ; \
        addq    $2, %rbp

.Lunaligned_loop:
        movq    %rbp, %rax
        andq     $ALIGN_MASK, %rax
        je .Laligned
        one_simple_xor(dx, r10, r11)
        subq    $1, %r9
	je	.Lfinished
        jmp .Lunaligned_loop

.Laligned:
        // now rbp is aligned.

        // four bytes of readahead seem to provide timings a lot more
        // stable than with just two.

#define READAHEAD_COEFFS        8
        // we're reading READAHEAD_COEFFS coeffs ahead of time, so the
        // tail must contain at least that many coefficients.

        leaq    (%rbp,%r9,2),%r15

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
        leaq    (%rbp,%r9,2),%r15

        movl 0(%rbp), %eax
        movl 4(%rbp), %ebx
        movl 8(%rbp), %ecx
        movl 12(%rbp), %edx
        addq    $LOOP_LENGTH, %rbp

	cmpq	%rbp, %r15
	jb	.Lfixup_before_tail

.Lmain_loop:

        // We have two consecutive entries stored in eax.
#define two_xors(idxreg, breg1, breg2, breg3, offset)		\
        movzwl  % ## idxreg, % ## breg1 ## d			; \
        movl    % ## breg1 ## d, % ## breg2 ## d		; \
        shrl    $12, % ## breg1 ## d		        	; \
        leaq    (%rsi, % ## breg1, 8), %rsi			; \
        andl    $0x0fff, % ## breg2 ## d			; \
        shrl    $16, %e ## idxreg		        	; \
        movq    (%rsi), % ## breg3		        	; \
        xorq    % ## breg3, (%rdi, % ## breg2, 8)		; \
        movzwl  % ## idxreg, % ## breg1 ## d			; \
        shrl    $12, % ## breg1 ## d			        ; \
        leaq    (%rsi, % ## breg1, 8), %rsi			; \
        andl    $0x0fff, %e ## idxreg   			; \
        movq    (%rsi), % ## breg3		        	; \
        xorq    % ## breg3, (%rdi, %r ## idxreg, 8)		; \
        movl    offset(%rbp), %e ## idxreg


/* expansion for two_xors(ax,r8,r9,r10,0):
        movzwl  % ## idxreg, %r8d
        movl    %r8d, %r9d
        shrl    $12, $r8d
        leal    (%rsi, %r8, 8), %rsi
        andl    $0x0fff, %r9d
        shll    $3, $r9d
        shrl    $16, %eax
        movq    (%rsi), %r10
        xorq    %r10, (%rdi, %r9, 8)
        movzwl  %ax, %r8d
        shrl    $12, $r8d
        leal    (%rsi, %r8, 8), %rsi
        andl    $0x0fff, %eax
        shll    $3, %eax
        movq    (%rsi), %r10
        xorq    %r10, (%rdi, %rax, 8)
        movl    0(%rbp), %eax
 */

        two_xors(ax, r8,  r9 , r10, 0)
        two_xors(bx, r11, r12, r13, 4)
        two_xors(cx, r8,  r9,  r10, 8)
        two_xors(dx, r11, r12, r13, 12)

        // note that we need 4*READAHEAD_COEFFS bytes of readahead at the
        // end of the array.

	addq	$LOOP_LENGTH, %rbp
	cmpq	%rbp, %r15
	ja	.Lmain_loop

.Lfixup_before_tail:
        // fix r15 for real value
        movq    $READAHEAD_COEFFS, %r9
        leaq    (%r15,%r9,2), %r15
        // also fix rbp
        subq    $LOOP_LENGTH, %rbp

.Ltail_loop:
        cmpq    %rbp, %r15
        je      .Lfinished
        one_simple_xor(dx, r10, r11)
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
	.size	matmul_sub_small2_asm, .-matmul_sub_small2_asm
