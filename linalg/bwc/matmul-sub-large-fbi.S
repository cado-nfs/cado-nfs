#if 0
static inline void matmul_sub_large_fbi_asm(abobj_ptr x, abt ** sb, const abt * z, const uint8_t * q, unsigned int n)
{
    /* Dispatch data found in z[0]...z[f(n-1)] such that z[f(i)] is in
     * array pointed to by sb[q[2*i+1]]. The function f(i) is given by
     * the sum q[0]+q[2]+...+q[2*(i-1)]. Exactly 2n coefficients are
     * expected in q[] All the sb[] pointers are increased */
    for(unsigned int c = 0 ; c < n ; c++) {
        z += aboffset(x, *q);
        q++;
        abcopy(x, sb[*q], z, 1);
        sb[*q]+= aboffset(x, 1);
        q++;
    }
}

static inline void matmul_sub_large_fbi_asm(ulong ** sb, const ulong * z, const uint8_t * q, unsigned int n)
    for(unsigned int c = 0 ; c < n ; c++) {
        z += *q++;
        *sb[*q++]++ = *z;
    }
}
#endif
.section .note.GNU-stack,"",%progbits
	.text
	.p2align 4,,15
.globl matmul_sub_large_fbi_asm
	.type	matmul_sub_large_fbi_asm, @function
matmul_sub_large_fbi_asm:
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
        test    %r9,%r9

        // at .Lfinished, rbp is shipped as a return value.
	je	.Lfinished

        // we want our reads aligned, we want to unroll our main loop,
        // and leave at least some space in the tail because the tight
        // loop does some readahead.

#define ALIGN_MASK      15

#define one_simple_move(reg2,reg3,reg4)   	\
        movzbl  (%rbp),  % ## reg2 ## d	        	        ; \
        movq    (%rsi,% ## reg2,8), % ## reg3			; \
	leaq	(%rsi,% ## reg2,8), %rsi			; \
        movzbl  1(%rbp),  % ## reg2 ## d		        ; \
	movq	(%rdi,% ## reg2,8), % ## reg4			; \
	movq	% ## reg3, (% ## reg4)		        	; \
	addq	$8, (%rdi,% ## reg2,8)                          ; \
        addq    $2, %rbp

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

#define READAHEAD_COEFFS        16
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

#define LOOP_LENGTH     16

        subq $READAHEAD_COEFFS, %r9
        // r15+2*READAHEAD_COEFFS : end data pointer.
        leaq    (%rbp,%r9,2),%r15

// do the readahead.

        // movq 0(%rbp), %rax
        movdqa (%rbp), %xmm0
        // movq 8(%rbp), %rbx
        // movq 16(%rbp), %rcx
        // movq 24(%rbp), %rdx

        addq    $LOOP_LENGTH, %rbp

	cmpq	%rbp, %r15
	jb	.Lfixup_before_tail

        movq    %rsp, %rcx
        subq    $16, %rsp
        andq    $-16, %rsp
        movq    %rcx, -8(%rsp)

.Lmain_loop:

// reg1 : one of a..d, containing 4 coeffs in two byte pairs (col index
// difference, bucket number).
// reg2 scratch buffer for storing index difference and (afterwards)
// bucket number
// reg3 scratch buffer for the coefficient block in the source vector
// reg4 scratch buffer for the address in the bucket
#define use_one_coefficient(reg1,reg2,reg3,reg4)   	        \
        movzbl  % ## reg1 ## l, % ## reg2 ## d		        ; \
        movq    (%rsi,% ## reg2,8), % ## reg3			; \
	leaq	(%rsi,% ## reg2,8), %rsi			; \
        shrq    $8, %r ## reg1 ## x		               	; \
        movzbl  % ## reg1 ## l, % ## reg2 ## d		        ; \
	movq	(%rdi,% ## reg2,8), % ## reg4			; \
	movq	% ## reg3, (% ## reg4)		        	; \
	addq	$8, % ## reg4                                   ; \
        movq    % ## reg4, (%rdi,% ## reg2,8)

#define four_moves_and_reload(reg1,reg2,reg3,reg4,offset)       \
        use_one_coefficient(reg1,reg2,reg3,reg4)                ; \
        shrq    $8, %r ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4)                ; \
        shrq    $8, %r ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4)                ; \
        shrl    $8, %e ## reg1 ## x                             ; \
        use_one_coefficient(reg1,reg2,reg3,reg4)                ; \
	movq	offset (%rbp), %r ## reg1 ## x

        // on core2 U9400, using vstrip 6 of c135b, 0x200 and 0x400 are
        // good choices.
        // prefetchnta     0x200(%rsi)
        // four_moves_and_reload(a,r8,r9,r10,)
        // four_moves_and_reload(b,r11,r12,r13,8)
        // four_moves_and_reload(c,r8,r9,r10,16)
        // four_moves_and_reload(d,r11,r12,r13,24)

//////////////////////////
    
        // prefetcht1     0x40(%rsi)
        // prefetcht1     0x80(%rsi)
        // prefetcht1     0xc0(%rsi)
        // prefetcht1     0x100(%rsi)

        movdqa  %xmm0, (%rsp)
        movq    (%rsp), %rdx
        movq    8(%rsp), %rax
        movdqa  0(%rbp), %xmm0
        // movq    %rax, %rdx
        // movq    0(%rbp), %rax

#define read_two_bytes_and_load(dreg,r8,r9d,mm0)              \
        movzbl  % ## dreg ## l, % ## r8 ## d       	; \
        leaq    (%rsi,% ## r8,8), %rsi  		; \
        movq    (%rsi), % ## mm0        		; \
        shrq    $8,%r ## dreg ## x                 	; \
        movzbl  % ## dreg ## l, % ## r9d       		; \
        shrq    $8,%r ## dreg ## x

#define read_two_bytes_and_load2(dreg,bx,r9d,mm0)              \
        movzbl  % ## dreg ## l, %e ## bx       		; \
        leaq    (%rsi,%r ## bx,8), %rsi  		; \
        movq    (%rsi), % ## mm0        		; \
        shrq    $8,%r ## dreg ## x                 	; \
        movzbl  % ## dreg ## l, % ## r9d       		; \
        shrq    $8,%r ## dreg ## x

        read_two_bytes_and_load(d,r8,r8d,mm0)
        read_two_bytes_and_load(d,r10,r10d,mm1)
        read_two_bytes_and_load(d,r12,r12d,mm2)
        read_two_bytes_and_load(d,r14,r14d,mm3)

        read_two_bytes_and_load(a,r9,r9d,mm4)
        read_two_bytes_and_load(a,r11,r11d,mm5)
        read_two_bytes_and_load(a,r13,r13d,mm6)
        read_two_bytes_and_load2(a,bx,ebx,mm7)

        // now do stores

#define do_one_store(r9,mm0)                    \
	movq	(%rdi,% ## r9,8), %rcx          ; \
	movq	% ## mm0, (%rcx)                ; \
	addq	$8, %rcx                        ; \
        movq    %rcx, (%rdi,% ## r9,8)

        do_one_store(r8,mm0)
        do_one_store(r10,mm1)
        do_one_store(r12,mm2)
        do_one_store(r14,mm3)

        do_one_store(r9,mm4)
        do_one_store(r11,mm5)
        do_one_store(r13,mm6)
        do_one_store(rbx,mm7)


//////////////////////////

        // note that we need 2*READAHEAD_COEFFS bytes of readahead at the
        // end of the array.

	addq	$LOOP_LENGTH, %rbp
	cmpq	%rbp, %r15
	ja	.Lmain_loop

        emms

        movq    -8(%rsp), %rsp

.Lfixup_before_tail:
        // fix r15 for real value
        movq    $READAHEAD_COEFFS, %r9
        leaq    (%r15,%r9,2), %r15
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
	.size	matmul_sub_large_fbi_asm, .-matmul_sub_large_fbi_asm
