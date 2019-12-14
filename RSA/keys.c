#include<memory.h>
#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include "keys.h"

#define S_LENGTH 54
#define TRUE 1
#define FALSE 0

const int SP = 14;
DTYPE S[S_LENGTH] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251};

void seive(DTYPE *p, DTYPE *S, DTYPE *w)
{
	int i, j, flag=0;

	//first generate an odd number
	for(i=MAX_PRIME_LENGTH-1; i>=0; i--)
	{
		srand(time(NULL) * 1000 + (clock() % 1000) + i);
		p[i] = (DTYPE)rand(); 
		// p[i] = (DTYPE)rand() % MAX_VAL;
		if(i==MAX_PRIME_LENGTH-1) p[i] = (p[i] | 0x0001); //p is odd
		if(i==0) p[i] = (p[i] | 0x8000); //the most significant bit of p is 1
	}
	for(i=0; i<S_LENGTH; i++)
	{
		BN_mod(&w[i], p, MAX_PRIME_LENGTH, &S[i], 1);
	}
	
	while(flag==0)
	{
		for(i=0; i<S_LENGTH; i++)
		{
			if(w[i]==0)
			{
				for(j=0; j<S_LENGTH; j++)
				{
					w[j]=(w[j]+2)%S[j];
				}
				//p=p+2
				BN_inc(p, MAX_PRIME_LENGTH);
				BN_inc(p, MAX_PRIME_LENGTH);
				break;
			}
			else continue;
		}
		if(i==S_LENGTH) flag = 1;
	}
}

//running Miller-Rabin Probable Primality Test, length of n is MAX_PRIME_lENGTH
int probable_prime_test(DTYPE *n, DTYPE SecurityParam)
{
	int i, j, s=0, flag=0;
	DTYPE n_sub_1[MAX_PRIME_LENGTH]={0}, n_sub[MAX_PRIME_LENGTH]={0}, y[MAX_PRIME_LENGTH]={0};
	DTYPE a[2*MAX_PRIME_LENGTH]={0};
	DTYPE r[MAX_PRIME_LENGTH+1]={0}, r_mod[MAX_PRIME_LENGTH]={0};
	DTYPE inv, n0;

	BN_assign(n_sub_1, n, MAX_PRIME_LENGTH);
	BN_dec(n_sub_1, MAX_PRIME_LENGTH);
	BN_assign(n_sub, n_sub_1, MAX_PRIME_LENGTH); //n_sub = n-1

	for(i=2*MAX_PRIME_LENGTH-1, j=MAX_PRIME_LENGTH-1; i>=0 && j>=0; i--, j--)
	{
		a[i] = n_sub[j];
	} //a=n-1

	r[0]=1;
	BN_mod(r_mod, r, MAX_PRIME_LENGTH+1, n, MAX_PRIME_LENGTH);

	n0=n[MAX_PRIME_LENGTH-1];
	n0 = MAX_VAL-n0;
	n0 += 1;
	BN_MonMod_inv(&inv, &n0, 1);

	while(BN_LSB(n_sub_1[MAX_PRIME_LENGTH-1])!=1)
	{
		BN_right_shift(n_sub_1, MAX_PRIME_LENGTH, 1);
		s++;
	}

	for(i=1; i<=SecurityParam; i++)
	{
		BN_dec(a, 2*MAX_PRIME_LENGTH);
		BN_MonExp(y, a, 2*MAX_PRIME_LENGTH, n_sub_1, MAX_PRIME_LENGTH, n, MAX_PRIME_LENGTH, r_mod, inv);
		if(BN_isEqualOne(y, MAX_PRIME_LENGTH) == 0 && BN_cmp(y, MAX_PRIME_LENGTH, n_sub, MAX_PRIME_LENGTH) != 0)
		{
			for(j=0; j<s; j++)
			{
				if(BN_cmp(y, MAX_PRIME_LENGTH, n_sub, MAX_PRIME_LENGTH) != 0)
				{
					// BN_mod_mul2(y, y, MAX_PRIME_LENGTH, y, MAX_PRIME_LENGTH, n, MAX_PRIME_LENGTH);
					BN_mod_mul(y, y, MAX_PRIME_LENGTH, y, MAX_PRIME_LENGTH, n, MAX_PRIME_LENGTH);
					if(BN_isEqualOne(y, MAX_PRIME_LENGTH) == 1 || BN_cmp(y, MAX_PRIME_LENGTH, n_sub, MAX_PRIME_LENGTH) != 0) return FALSE;
				} //end of if
			} //end of for
		} //end of if
	} //end of for
	return TRUE;
}

//generate prime with length of MAX_PRIME_LENGTH, order=0->sk->p, order=1->sk->q
void prime_generate(rsa_sk *sk, int order)
{
	DTYPE p[MAX_PRIME_LENGTH], w[S_LENGTH]={0};

	do{
		seive(p, S, w);
	}while(probable_prime_test(p, SP)==FALSE); //SP, number of iterations in Miller Rabin Test

	if(order==0) BN_assign(sk->p, p, MAX_PRIME_LENGTH);
	else BN_assign(sk->q, p, MAX_PRIME_LENGTH);
}

//generate keys in RSA
void rsa_key_generation(rsa_pk *pk, rsa_sk *sk)
{
	BN_init(pk->e, MAX_PRIME_LENGTH);
	BN_init(pk->n, MAX_MODULUS_LENGTH);
	BN_init(sk->n, MAX_MODULUS_LENGTH);
	BN_init(sk->p, MAX_PRIME_LENGTH);
	BN_init(sk->q, MAX_PRIME_LENGTH);
	BN_init(sk->phi_n, MAX_MODULUS_LENGTH);
	BN_init(sk->d, MAX_MODULUS_LENGTH);
	BN_init(sk->d1, MAX_PRIME_LENGTH);
	BN_init(sk->d2, MAX_PRIME_LENGTH);
	BN_init(sk->p_inv, MAX_PRIME_LENGTH);
	BN_init(sk->r_mod, MAX_MODULUS_LENGTH);
	BN_init(sk->p_mod, MAX_PRIME_LENGTH);
	BN_init(sk->q_mod, MAX_PRIME_LENGTH);

	//e=2^16+1
	pk->e[MAX_PRIME_LENGTH - 2] = 1;
	pk->e[MAX_PRIME_LENGTH - 1] = 1;

	prime_generate(sk, 0); //sk->p
	prime_generate(sk, 1); //sk->q
	
	if(BN_cmp(sk->q, MAX_PRIME_LENGTH, sk->p, MAX_PRIME_LENGTH)==1)
	{
		BN_exchange(sk->q, sk->p, MAX_PRIME_LENGTH);
	} //if q is bigger than p, exchange them
	
	BN_mul(pk->n, sk->p, MAX_PRIME_LENGTH, sk->q, MAX_PRIME_LENGTH); //pk->n
	BN_assign(sk->n, pk->n, MAX_MODULUS_LENGTH); //sk->n

	DTYPE n0=pk->n[MAX_MODULUS_LENGTH-1];
	n0 = MAX_VAL-n0;
	n0 += 1;
	BN_MonMod_inv(&pk->n_inv, &n0, 1); //pk->n_inv
	sk->n_inv = pk->n_inv; //sk->n_inv

	DTYPE r[MAX_MODULUS_LENGTH+1] = {0};
	int bits = BN_valid_bits(pk->n, MAX_MODULUS_LENGTH);
	if(bits==MAX_MODULUS_BITS) r[0]=1;
	else
	{
		r[MAX_MODULUS_LENGTH] = 1;
		BN_left_shift(r, MAX_MODULUS_LENGTH + 1, bits);
	}
	BN_mod(pk->r_mod, r, MAX_MODULUS_LENGTH+1, pk->n, MAX_MODULUS_LENGTH); //pk->r_mod
	BN_assign(sk->r_mod, pk->r_mod, MAX_MODULUS_LENGTH); //sk->r_mod

	DTYPE p_dec[MAX_PRIME_LENGTH], q_dec[MAX_PRIME_LENGTH];
	BN_assign(p_dec, sk->p, MAX_PRIME_LENGTH);
	BN_assign(q_dec, sk->q, MAX_PRIME_LENGTH);
	BN_dec(p_dec, MAX_PRIME_LENGTH);
	BN_dec(q_dec, MAX_PRIME_LENGTH);
	BN_mul(sk->phi_n, p_dec, MAX_PRIME_LENGTH, q_dec, MAX_PRIME_LENGTH); //sk->phi_n
	BN_mod_inv(sk->d, pk->e, MAX_PRIME_LENGTH, sk->phi_n, MAX_MODULUS_LENGTH); //sk->d
	BN_mod(sk->d1, sk->d, MAX_MODULUS_LENGTH, p_dec, MAX_PRIME_LENGTH); //sk->d1
	BN_mod(sk->d2, sk->d, MAX_MODULUS_LENGTH, q_dec, MAX_PRIME_LENGTH); //sk->d2

	DTYPE p_mod_q[MAX_PRIME_LENGTH]={0};
	BN_sub(p_mod_q, MAX_PRIME_LENGTH, sk->p, MAX_PRIME_LENGTH, sk->q, MAX_PRIME_LENGTH);
	BN_mod_inv(sk->p_inv, p_mod_q, MAX_PRIME_LENGTH, sk->q, MAX_PRIME_LENGTH); //sk->p_inv

	DTYPE pp[MAX_PRIME_LENGTH+1] = {0};
	pp[0]=1;
	BN_mod(sk->p_mod, pp, MAX_PRIME_LENGTH+1, sk->p, MAX_PRIME_LENGTH); //sk->p_mod

	DTYPE qq[MAX_PRIME_LENGTH + 1] = { 0 };
	qq[0] = 1;
	BN_mod(sk->q_mod, qq, MAX_PRIME_LENGTH+1, sk->q, MAX_PRIME_LENGTH); //sk->q_mod

	DTYPE p0=sk->p[MAX_PRIME_LENGTH-1];
	p0 = MAX_VAL-p0;
	p0 += 1;
	BN_MonMod_inv(&sk->p0_inv, &p0, 1); //sk->p0_inv

	DTYPE q0=sk->q[MAX_PRIME_LENGTH-1];
	q0 = MAX_VAL-q0;
	q0 += 1;
	BN_MonMod_inv(&sk->q0_inv, &q0, 1); //sk->q0_inv
}
