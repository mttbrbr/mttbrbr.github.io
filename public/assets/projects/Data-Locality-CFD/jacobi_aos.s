	.file	"jacobi_kernel.c"
	.text
	.p2align 4
	.globl	jacobi_row
	.type	jacobi_row, @function
jacobi_row:
.LFB0:
	.cfi_startproc
	leaq	-1(%rsi), %rax
	cmpq	$1, %rax
	jbe	.L9
	imulq	%rsi, %rdx
	subq	$24, %rsp
	.cfi_def_cfa_offset 32
	vmovapd	%xmm0, %xmm6
	movq	%rdi, %r11
	movq	%rbx, (%rsp)
	.cfi_offset 3, -32
	movq	%rsi, %rbx
	movq	%rbp, 8(%rsp)
	.cfi_offset 6, -24
	movq	%rdx, %rbp
	cmpq	$3, %rsi
	je	.L5
	leaq	(%rdx,%rsi), %rax
	movq	%r12, 16(%rsp)
	.cfi_offset 12, -16
	leaq	-2(%rsi), %r12
	vmovddup	%xmm2, %xmm14
	leaq	(%rax,%rax,2), %rax
	movq	%r12, %r9
	vmovddup	%xmm3, %xmm13
	movq	%r12, %r10
	salq	$4, %rax
	shrq	%r9
	andq	$-2, %r10
	vmovddup	%xmm0, %xmm12
	leaq	64(%rdi,%rax), %rcx
	movq	%rdx, %rax
	vmovddup	%xmm1, %xmm11
	addq	%r10, %r9
	subq	%rsi, %rax
	vmovddup	%xmm4, %xmm10
	vmovddup	%xmm5, %xmm9
	salq	$5, %r9
	leaq	(%rax,%rax,2), %rax
	vmovapd	%xmm0, %xmm16
	salq	$4, %rax
	leaq	64(%rdi,%rax), %rdx
	leaq	0(%rbp,%rbp,2), %rax
	salq	$4, %rax
	leaq	120(%rdi,%rax), %r8
	leaq	88(%r11,%rax), %rsi
	leaq	24(%rdi,%rax), %rdi
	xorl	%eax, %eax
	.p2align 4,,10
	.p2align 3
.L4:
	vmovupd	16(%rdx,%rax), %xmm7
	vmovupd	(%rdi,%rax), %xmm15
	vmovhpd	64(%rdx,%rax), %xmm7, %xmm0
	vmovupd	16(%rcx,%rax), %xmm7
	vmovhpd	48(%rdi,%rax), %xmm15, %xmm8
	vmovhpd	64(%rcx,%rax), %xmm7, %xmm6
	vmovupd	48(%r8,%rax), %xmm7
	vaddpd	%xmm6, %xmm0, %xmm0
	vmovlpd	8(%r8,%rax), %xmm7, %xmm6
	vmovupd	48(%rdi,%rax), %xmm7
	vmulpd	%xmm14, %xmm0, %xmm0
	vmovlpd	8(%rdi,%rax), %xmm7, %xmm7
	vaddpd	%xmm7, %xmm6, %xmm6
	vmovupd	(%rcx,%rax), %xmm7
	vmulpd	%xmm13, %xmm6, %xmm6
	vaddpd	%xmm6, %xmm0, %xmm0
	vmovhpd	48(%rcx,%rax), %xmm7, %xmm6
	vmovupd	(%rdx,%rax), %xmm7
	vmovhpd	48(%rdx,%rax), %xmm7, %xmm7
	vsubpd	%xmm7, %xmm6, %xmm6
	vmovupd	(%r8,%rax), %xmm7
	vmovhpd	48(%r8,%rax), %xmm7, %xmm7
	vsubpd	%xmm8, %xmm7, %xmm7
	vmulpd	%xmm12, %xmm6, %xmm6
	vmulpd	%xmm11, %xmm7, %xmm7
	vaddpd	%xmm7, %xmm6, %xmm6
	vmulpd	%xmm10, %xmm6, %xmm6
	vsubpd	%xmm6, %xmm0, %xmm0
	vdivpd	%xmm9, %xmm0, %xmm0
	vmovlpd	%xmm0, (%rsi,%rax)
	vmovhpd	%xmm0, 48(%rsi,%rax)
	addq	$96, %rax
	cmpq	%r9, %rax
	jne	.L4
	vmovapd	%xmm16, %xmm6
	cmpq	%r12, %r10
	je	.L12
	movq	16(%rsp), %r12
	.cfi_restore 12
	incq	%r10
.L3:
	addq	%rbp, %r10
	movq	8(%rsp), %rbp
	leaq	-3(%r10,%r10,2), %rcx
	leaq	(%rbx,%r10), %rax
	salq	$4, %rcx
	leaq	3(%r10,%r10,2), %rsi
	leaq	(%rax,%rax,2), %rdx
	movq	%r10, %rax
	addq	%r11, %rcx
	salq	$4, %rsi
	subq	%rbx, %rax
	leaq	(%r10,%r10,2), %rdi
	vmovsd	32(%rcx), %xmm0
	addq	%r11, %rsi
	leaq	(%rax,%rax,2), %rax
	salq	$4, %rdx
	salq	$4, %rax
	addq	%r11, %rdx
	salq	$4, %rdi
	movq	(%rsp), %rbx
	vaddsd	32(%rsi), %xmm0, %xmm0
	addq	%r11, %rax
	vmulsd	%xmm3, %xmm0, %xmm0
	vmovsd	32(%rax), %xmm3
	vaddsd	32(%rdx), %xmm3, %xmm3
	vmulsd	%xmm2, %xmm3, %xmm3
	vmovsd	24(%rsi), %xmm2
	vsubsd	24(%rcx), %xmm2, %xmm2
	vmulsd	%xmm1, %xmm2, %xmm2
	vmovsd	16(%rdx), %xmm1
	vaddsd	%xmm3, %xmm0, %xmm0
	vsubsd	16(%rax), %xmm1, %xmm1
	vmulsd	%xmm6, %xmm1, %xmm1
	vaddsd	%xmm1, %xmm2, %xmm2
	vmulsd	%xmm4, %xmm2, %xmm2
	vsubsd	%xmm2, %xmm0, %xmm0
	vdivsd	%xmm5, %xmm0, %xmm0
	vmovsd	%xmm0, 40(%r11,%rdi)
	addq	$24, %rsp
	.cfi_restore 3
	.cfi_restore 6
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L12:
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	.cfi_offset 6, -24
	.cfi_offset 12, -16
	movq	16(%rsp), %r12
	.cfi_restore 12
	movq	(%rsp), %rbx
	movq	8(%rsp), %rbp
	addq	$24, %rsp
	.cfi_restore 3
	.cfi_restore 6
	.cfi_def_cfa_offset 8
	ret
	.p2align 4,,10
	.p2align 3
.L9:
	ret
.L5:
	.cfi_def_cfa_offset 32
	.cfi_offset 3, -32
	.cfi_offset 6, -24
	movl	$1, %r10d
	jmp	.L3
	.cfi_endproc
.LFE0:
	.size	jacobi_row, .-jacobi_row
	.ident	"GCC: (GNU) 16.2.1 20261122"
	.section	.note.GNU-stack,"",@progbits
