#if 0
static void
matmul_sub_large_fbd_asm(abt ** sb, const abt * z, const uint8_t * q, unsigned int n)
{
    /* Dispatch data found in z[0]...z[n] such that z[i] is in array
     * pointed to by sb[q[i]]. Exactly n coefficients are expected. All
     * the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        *sb[q[c]]++ = *z++;
    }
}
#endif
.section .note.GNU-stack,"",%progbits
	.text
	.p2align 4,,15
.globl matmul_sub_large_fbd_asm
	.type	matmul_sub_large_fbd_asm, @function
matmul_sub_large_fbd_asm:
        pushq %rbp
        pushq %rbx
        pushq %r15
        pushq %r14
        pushq %r13
        pushq %r12
        // arg order: %rdi, %rsi, %rdx, %rcx, %r8, and %r9
        // %rdi sbp
        // %rsi z
	// %rdx q  --> rbp
        movq    %rdx, %rbp
	// %rcx q  --> r8 and r9
        movq    %rcx, %r8
        movq    %rcx, %r9

        test    %r9, %r9

        // at .Lfinished, rbp is shipped as a return value.
	je	.Lfinished

        // we want our reads aligned, we want to unroll our main loop,
        // and leave at least some space in the tail because the tight
        // loop does some readahead.

#define ALIGN_MASK      7

#define one_simple_move(reg2,reg3,reg4)   	\
        movq    (%rsi), % ## reg3			; \
        movzbl  (%rbp),  % ## reg2 ## d		        ; \
	movq	(%rdi,% ## reg2,8), % ## reg4		; \
	movq	% ## reg3, (% ## reg4)		       	; \
	addq	$8, (%rdi,% ## reg2,8)                  ; \
        addq    $1, %rbp                                ; \
        addq    $8, %rsi

.Lunaligned_loop:
        movq    %rbp, %rax
        andq     $ALIGN_MASK, %rax
        je .Laligned
        one_simple_move(r8, r10, r11)
        subq    $1, %r9
	je	.Lfinished
        jmp .Lunaligned_loop

.Laligned:
        // now rbp is aligned.

        // four bytes of readahead seem to provide timings a lot more
        // stable than with just two.

#define READAHEAD_COEFFS        32
        // we're reading READAHEAD_COEFFS coeffs ahead of time, so the
        // tail must contain at least that many coefficients.

        leaq    (%rbp,%r9,1),%r15

        cmpq $READAHEAD_COEFFS, %r9

        // the tail loop is just simple and stupid. rbp points to a
        // location where exactly r9 coeffs are to be read. Note that
        // some coefficients are already in registers at this point. But
        // this does not matter much (we could avoid re-reading, but
        // really we don't care).
        jb .Ltail_loop

        // we're going to advance rbp one loop iteration ahead. This
        // eases the handling of the tail.

#define LOOP_LENGTH     32

        subq $READAHEAD_COEFFS, %r9
        // r15+1*READAHEAD_COEFFS : end data pointer.
        leaq    (%rbp,%r9,1),%r15

        movq 0(%rbp), %rax
        movq 8(%rbp), %rbx
        movq 16(%rbp), %rcx
        movq 24(%rbp), %rdx

        addq    $LOOP_LENGTH, %rbp

	cmpq	%rbp, %r15
	jb	.Lfixup_before_tail

.Lmain_loop:

// reg1 : one of a..d, containing 8 coeffs, one per byte (only the bucket
// number).
// reg2 scratch buffer for storing bucket number
// reg3 scratch buffer for the coefficient block in the source vector
// reg4 scratch buffer for the address in the bucket
#define use_one_coefficient(reg1,reg2,reg3,reg4,soffset)   	\
        movq    soffset(%rsi), % ## reg3			; \
        movzbl  % ## reg1 ## l, % ## reg2 ## d		        ; \
	movq	(%rdi,% ## reg2,8), % ## reg4			; \
	movq	% ## reg3, (% ## reg4)		        	; \
	addq	$8, % ## reg4                                   ; \
        movq    % ## reg4, (%rdi,% ## reg2,8)

#define eight_moves_and_reload(reg1,reg2,reg3,reg4,offset,soffset)       \
        use_one_coefficient(reg1,reg2,reg3,reg4, soffset + 0)   ; \
        shrq    $8, %r ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4, soffset + 8)   ; \
        shrq    $8, %r ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4, soffset + 16)  ; \
        shrq    $8, %r ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4, soffset + 24)  ; \
        shrq    $8, %r ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4, soffset + 32)  ; \
        shrl    $8, %e ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4, soffset + 40)  ; \
        shrl    $8, %e ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4, soffset + 48)  ; \
        shrl    $8, %e ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4, soffset + 56)  ; \
	movq	offset (%rbp), %r ## reg1 ## x

        prefetchnta     0x100(%rsi)
        eight_moves_and_reload(a,r8,r9,r10,,0)
        addq    $0x40, %rsi
        eight_moves_and_reload(b,r11,r12,r13,8,0)
        addq    $0x40, %rsi
        eight_moves_and_reload(c,r8,r9,r10,16,0)
        addq    $0x40, %rsi
        eight_moves_and_reload(d,r11,r12,r13,24,0)
        addq    $0x40, %rsi

        // note that we need 2*READAHEAD_COEFFS bytes of readahead at the
        // end of the array.

	addq	$LOOP_LENGTH, %rbp
	cmpq	%rbp, %r15
	ja	.Lmain_loop

.Lfixup_before_tail:
        // fix r15 for real value
        movq    $READAHEAD_COEFFS, %r9
        leaq    (%r15,%r9,1), %r15
        // also fix rbp
        subq    $LOOP_LENGTH, %rbp

.Ltail_loop:
        cmpq    %rbp, %r15
        je      .Lfinished
        one_simple_move(r8, r10, r11)
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
	.size	matmul_sub_large_fbd_asm, .-matmul_sub_large_fbd_asm
