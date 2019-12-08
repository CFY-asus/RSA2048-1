#include<stdio.h>
#include "rsa.h"

//use pk to encrypt msg, output is cipher
void rsa_encrypt(DTYPE *cipher, DTYPE cipher_len, DTYPE *msg, DTYPE msg_len, rsa_pk *pk)
{
	BN_MonExp(cipher, msg, msg_len, pk->e, MAX_PRIME_LENGTH, pk->n, MAX_MODULUS_LENGTH, pk->r_mod, pk->n_inv);
}

//use sk to decrype cipher, result is saved in output
void rsa_decrypt(DTYPE *output, DTYPE output_len, DTYPE *cipher, DTYPE cipher_len, rsa_sk *sk)
{
	int i;
	DTYPE M1[MAX_PRIME_LENGTH]={0};
	DTYPE M2[MAX_PRIME_LENGTH]={0};
	DTYPE M_tmp[MAX_PRIME_LENGTH]={0};
	DTYPE r1_tmp[MAX_PRIME_LENGTH]={0};
	DTYPE r2_tmp[MAX_MODULUS_LENGTH]={0};
	DTYPE out[MAX_MODULUS_LENGTH+1] = {0};

	BN_MonExp(M1, cipher, cipher_len, sk->d1, MAX_PRIME_LENGTH, sk->p, MAX_PRIME_LENGTH, sk->p_mod, sk->p0_inv);
	BN_MonExp(M2, cipher, cipher_len, sk->d2, MAX_PRIME_LENGTH, sk->q, MAX_PRIME_LENGTH, sk->q_mod, sk->q0_inv);
	if(BN_cmp(M2, MAX_PRIME_LENGTH, M1, MAX_PRIME_LENGTH)>=0) //M2>=M1
	{
		BN_sub(M_tmp, MAX_PRIME_LENGTH, M2, MAX_PRIME_LENGTH, M1, MAX_PRIME_LENGTH);
	}
	else //M1>M2
	{
		BN_sub(M_tmp, MAX_PRIME_LENGTH, M1, MAX_PRIME_LENGTH, M2, MAX_PRIME_LENGTH);
		BN_mod(M_tmp, M_tmp, MAX_PRIME_LENGTH, sk->q, MAX_PRIME_LENGTH);
		BN_sub(M_tmp, MAX_PRIME_LENGTH, sk->q, MAX_PRIME_LENGTH, M_tmp, MAX_PRIME_LENGTH);
	}
	// BN_mod_mul2(r1_tmp, M_tmp, MAX_PRIME_LENGTH, sk->p_inv, MAX_PRIME_LENGTH, sk->q, MAX_PRIME_LENGTH);
	BN_mod_mul(r1_tmp, M_tmp, MAX_PRIME_LENGTH, sk->p_inv, MAX_PRIME_LENGTH, sk->q, MAX_PRIME_LENGTH);
	BN_mul(r2_tmp, r1_tmp, MAX_PRIME_LENGTH, sk->p, MAX_PRIME_LENGTH);
	BN_add(out, MAX_MODULUS_LENGTH+1, r2_tmp, MAX_MODULUS_LENGTH, M1, MAX_PRIME_LENGTH);
	
	for(i=0; i<MAX_MODULUS_LENGTH; i++)
	{
		output[i] = out[i+1];
	}
}
