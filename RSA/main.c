#include<stdio.h>
#include<time.h>
#include<string.h>
#include "rsa.h"

int main(void)
{
	rsa_pk pk;
	rsa_sk sk;

	clock_t key_s, key_e, rsa_s, rsa_e;
	double key_t, rsa_t;

	DTYPE msg[MAX_MODULUS_LENGTH] = { 0 };
	DTYPE cipher[MAX_MODULUS_LENGTH] = { 0 };
	DTYPE plaintext[MAX_MODULUS_LENGTH] = { 0 };

    char message[MAX_MODULUS_LENGTH * 2];
    char ciphertext[MAX_MODULUS_LENGTH * 4 + 1];
	char plainmsg[MAX_MODULUS_LENGTH * 2 + 1];

	char *m = "../../message.txt";
	char *c = "../../ciphertext.txt";
	char *p = "../../plaintext.txt";
	char *pkey = "../../publicKey.txt";
	char *skey = "../../privateKey.txt";

	FILE *fp;

	printf("Key generation starts...\n");
    key_s = clock();
	rsa_key_generation(&pk, &sk);
    key_e = clock();
	printf("Key generation done...\n\n");
    key_t = (key_e - key_s) / CLOCKS_PER_SEC;

    //write public keys into file
    fp = fopen(pkey, "w");
    fclose(fp);
    fp = fopen(pkey, "a");
    BN_printToFile(pk.n, MAX_MODULUS_LENGTH, fp);
    fputc('\n', fp);
    BN_printToFile(pk.e, MAX_PRIME_LENGTH, fp);
    fclose(fp);

    //write private keys into file
    fp = fopen(skey, "w");
    fclose(fp);
    fp = fopen(skey, "a");
    BN_printToFile(sk.p, MAX_PRIME_LENGTH, fp);
    fputc('\n', fp);
    BN_printToFile(sk.q, MAX_PRIME_LENGTH, fp);
    fputc('\n', fp);
    BN_printToFile(sk.d, MAX_MODULUS_LENGTH, fp);
    fclose(fp);

    //read message from file
    fp = fopen(m , "r");
    if(fp==NULL){
        printf("Cannot open file %s\n", m);
        return 0;
    }
    fgets(message, MAX_MODULUS_LENGTH * 2, fp);
    fclose(fp);
    string_to_hex(msg, message);
	
	printf("Encryption starts...\n");
    rsa_s = clock();
	rsa_encrypt(cipher, MAX_MODULUS_LENGTH, msg, MAX_MODULUS_LENGTH, &pk);
    rsa_e = clock();
    rsa_t = (double)(rsa_e-rsa_s) / CLOCKS_PER_SEC;
    printf("Encryption done...\n\n");

    //write ciphertext into file
    fp = fopen(c, "w");
    BN_printToFile(cipher, MAX_MODULUS_LENGTH, fp);
    fclose(fp);

    printf("Decryption starts...\n");
    rsa_s = clock();
    rsa_decrypt(plaintext, MAX_MODULUS_LENGTH, cipher, MAX_MODULUS_LENGTH, &sk);
    rsa_e = clock();
    printf("Decryption done...\n\n");
	rsa_t += (double)(rsa_e-rsa_s) / CLOCKS_PER_SEC;
    
    //write plaintext into file
    hex_to_string(plainmsg, plaintext);
    fp = fopen(p, "w");
    fputs(plainmsg, fp);
    fclose(fp);

	if(BN_cmp(msg, MAX_MODULUS_LENGTH, plaintext, MAX_MODULUS_LENGTH)==0)
    {
		printf("\nAfter decryption, plaintext equal to message.\n");
	}
	else
    {
		printf("\nAfter decryption, wrong answer.\n");
	}

	printf("\n%f seconds to generate keys, %f seconds to encrypt and decrypt\n", key_t, rsa_t);

	return 0;
}