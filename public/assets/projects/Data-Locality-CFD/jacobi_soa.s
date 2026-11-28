	.file	"jacobi_kernel.c"
	.text
	.p2align 4
	.globl	jacobi_row
	.type	jacobi_row, @function
jacobi_row:
.LFB0:
	.cfi_startproc
	pushq	%rbp
	.cfi_def_cfa_offset 16
	.cfi_offset 6, -16
	leaq	-1(%r8), %rax
	movq	%rsp, %rbp
	.cfi_def_cfa_register 6
	pushq	%r15
	pushq	%r14
	pushq	%r13
	pushq	%r12
	pushq	%rbx
	andq	$-32, %rsp
	.cfi_offset 15, -24
	.cfi_offset 14, -32
	.cfi_offset 13, -40
	.cfi_offset 12, -48
	.cfi_offset 3, -56
	movq	%rdi, -16(%rsp)
	movq	%rsi, -24(%rsp)
	movq	%rcx, -8(%rsp)
	cmpq	$1, %rax
	jbe	.L11
	leaq	-3(%r8), %rax
	vmovapd	%xmm0, %xmm18
	vbroadcastsd	%xmm0, %ymm15
	movq	%rdi, %rbx
	imulq	%r8, %r9
	vbroadcastsd	%xmm1, %ymm14
	vbroadcastsd	%xmm2, %ymm17
	movq	%rdx, %rdi
	vbroadcastsd	%xmm3, %ymm16
	vbroadcastsd	%xmm4, %ymm7
	vbroadcastsd	%xmm5, %ymm6
	cmpq	$2, %rax
	jbe	.L6
	movq	%r9, %rax
	leaq	1(%r9,%r8), %rdx
	leaq	16(,%r9,8), %rcx
	movq	%rsi, %r11
	subq	%r8, %rax
	salq	$3, %rdx
	leaq	-2(%r8), %r10
	leaq	(%rdi,%rcx), %r13
	leaq	8(,%rax,8), %rsi
	leaq	0(,%r9,8), %rax
	addq	%r11, %rcx
	movq	%r10, -40(%rsp)
	leaq	(%rdi,%rdx), %r15
	leaq	(%rdi,%rsi), %r14
	addq	%rbx, %rdx
	addq	%rbx, %rsi
	leaq	(%r11,%rax), %rbx
	movq	-8(%rsp), %r11
	shrq	$2, %r10
	leaq	(%rdi,%rax), %r12
	movq	%r10, -32(%rsp)
	vmovapd	%ymm16, %ymm12
	vmovapd	%ymm15, %ymm11
	vmovapd	%ymm14, %ymm10
	leaq	8(%r11,%rax), %r11
	vmovapd	%ymm7, %ymm9
	vmovapd	%ymm6, %ymm8
	salq	$5, %r10
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L4:
	vmovupd	(%rdx,%rax), %ymm20
	vmovupd	(%rcx,%rax), %ymm23
	vmovupd	0(%r13,%rax), %ymm0
	vmovupd	(%r15,%rax), %ymm19
	vsubpd	(%rsi,%rax), %ymm20, %ymm21
	vsubpd	(%rbx,%rax), %ymm23, %ymm24
	vaddpd	(%r12,%rax), %ymm0, %ymm0
	vaddpd	(%r14,%rax), %ymm19, %ymm19
	vmulpd	%ymm11, %ymm21, %ymm22
	vmulpd	%ymm10, %ymm24, %ymm25
	vmulpd	%ymm12, %ymm0, %ymm0
	vmulpd	%ymm17, %ymm19, %ymm19
	vaddpd	%ymm25, %ymm22, %ymm26
	vaddpd	%ymm19, %ymm0, %ymm0
	vmulpd	%ymm9, %ymm26, %ymm27
	vsubpd	%ymm27, %ymm0, %ymm0
	vdivpd	%ymm8, %ymm0, %ymm0
	vmovupd	%ymm0, (%r11,%rax)
	addq	$32, %rax
	cmpq	%r10, %rax
	jne	.L4
	movq	-32(%rsp), %rax
	salq	$2, %rax
	leaq	1(%rax), %rdx
	cmpq	-40(%rsp), %rax
	je	.L11
.L3:
	movq	%r8, %rcx
	subq	%rax, %rcx
	cmpq	$3, %rcx
	je	.L5
	movq	%r9, %r10
	leaq	1(%r9,%r8), %rsi
	movq	-16(%rsp), %rbx
	subq	%r8, %r10
	addq	%rax, %rsi
	leaq	1(%rax,%r10), %r10
	addq	%r9, %rax
	vmovupd	(%rdi,%rax,8), %xmm8
	vmovupd	(%rdi,%r10,8), %xmm0
	vaddpd	16(%rdi,%rax,8), %xmm8, %xmm8
	vaddpd	(%rdi,%rsi,8), %xmm0, %xmm0
	vmulpd	%xmm16, %xmm8, %xmm16
	vmovupd	(%rbx,%rsi,8), %xmm8
	movq	-24(%rsp), %rsi
	vmulpd	%xmm17, %xmm0, %xmm0
	vmovupd	16(%rsi,%rax,8), %xmm9
	vsubpd	(%rbx,%r10,8), %xmm8, %xmm8
	vsubpd	(%rsi,%rax,8), %xmm9, %xmm9
	movq	-8(%rsp), %rsi
	vmulpd	%xmm15, %xmm8, %xmm8
	vaddpd	%xmm16, %xmm0, %xmm0
	vmulpd	%xmm14, %xmm9, %xmm14
	vaddpd	%xmm14, %xmm8, %xmm8
	vmulpd	%xmm7, %xmm8, %xmm7
	vsubpd	%xmm7, %xmm0, %xmm0
	vdivpd	%xmm6, %xmm0, %xmm6
	vmovupd	%xmm6, 8(%rsi,%rax,8)
	leaq	-2(%rcx), %rax
	andq	$-2, %rax
	andl	$1, %ecx
	je	.L11
	addq	%rax, %rdx
.L5:
	addq	%r9, %rdx
	leaq	(%r8,%rdx), %rcx
	movq	%rdx, %rax
	leaq	1(%rdx), %rsi
	vmovsd	(%rdi,%rcx,8), %xmm0
	subq	%r8, %rax
	vaddsd	(%rdi,%rax,8), %xmm0, %xmm0
	vmulsd	%xmm2, %xmm0, %xmm0
	vmovsd	(%rdi,%rsi,8), %xmm2
	vaddsd	-8(%rdi,%rdx,8), %xmm2, %xmm2
	movq	-16(%rsp), %rdi
	vmulsd	%xmm3, %xmm2, %xmm2
	vaddsd	%xmm2, %xmm0, %xmm0
	vmovsd	(%rdi,%rcx,8), %xmm2
	vsubsd	(%rdi,%rax,8), %xmm2, %xmm2
	movq	-24(%rsp), %rax
	vmovsd	(%rax,%rsi,8), %xmm3
	vmulsd	%xmm18, %xmm2, %xmm2
	vsubsd	-8(%rax,%rdx,8), %xmm3, %xmm3
	movq	-8(%rsp), %rax
	vmulsd	%xmm1, %xmm3, %xmm3
	vaddsd	%xmm3, %xmm2, %xmm2
	vmulsd	%xmm4, %xmm2, %xmm2
	vsubsd	%xmm2, %xmm0, %xmm0
	vdivsd	%xmm5, %xmm0, %xmm0
	vmovsd	%xmm0, (%rax,%rdx,8)
.L11:
	vzeroupper
	leaq	-40(%rbp), %rsp
	popq	%rbx
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	popq	%rbp
	.cfi_remember_state
	.cfi_def_cfa 7, 8
	ret
.L6:
	.cfi_restore_state
	xorl	%eax, %eax
	movl	$1, %edx
	jmp	.L3
	.cfi_endproc
.LFE0:
	.size	jacobi_row, .-jacobi_row
	.ident	"GCC: (GNU) 16.2.1 20261122"
	.section	.note.GNU-stack,"",@progbits
