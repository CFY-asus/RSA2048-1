#include<stdint.h>
#include<string.h>

#define WORD_SIZE 16
#define MAX_MODULUS_BITS 2048
#define MAX_MODULUS_LENGTH (MAX_MODULUS_BITS/WORD_SIZE)
#define MAX_PRIME_BITS (MAX_MODULUS_BITS/2)
#define MAX_PRIME_LENGTH (MAX_PRIME_BITS/WORD_SIZE)
#define DTYPE uint16_t
#define T_DTYPE uint32_t
#define MAX_VAL (uint16_t)0xffff
#define BN_MSB(x) (uint16_t)((x & 0x8000) >> (WORD_SIZE-1)) //most significant bit of x
#define BN_LSB(x) (uint16_t)(x & 0x0001) //least significant bit of x
#define STD_FORMAT_STR "%04x"

typedef struct{
    int flag;
    DTYPE length;
    DTYPE *bn;
}sm; //signed multiprecision integer

//print multiprecision integer in hex format
void BN_print_hex(DTYPE *data, DTYPE data_length);

void BN_printToFile(DTYPE *data, DTYPE data_length, FILE *fp);

//transform string to hex format
void string_to_hex(DTYPE *hex, char *str);

//transform hex to string
void hex_to_string(char *str, DTYPE *hex);

//init data to zero
void BN_init(DTYPE *data, DTYPE data_length);

//data>>=digit
void BN_right_shift(DTYPE *data, DTYPE data_length, int digit);

//data<<=digit
void BN_left_shift(DTYPE *data, DTYPE data_length, int digit);

//if data=1, return 1; else return 0
int BN_isEqualOne(DTYPE *data, DTYPE data_length);

//if data=0, return 1; else return 0
int BN_isEqualZero(DTYPE *data, DTYPE data_length);

//set data to 1
void BN_setOne(DTYPE *data, DTYPE data_length);

void BN_inc(DTYPE *data, DTYPE data_length); //increasa 1

void BN_dec(DTYPE *data, DTYPE data_length); //decrease 1

//return loc of the first non-zero number
int BN_valid_pos(DTYPE *BN, DTYPE BN_length);

int BN_valid_bits(DTYPE *m, DTYPE length);

//a>b=>1, a<b=>-1, a=b=>0
int BN_cmp(DTYPE *a, DTYPE a_length, DTYPE *b, DTYPE b_length);

//assign b to a
void BN_assign(DTYPE *a, DTYPE *b, DTYPE b_length);

//exchange values of a and b
void BN_exchange(DTYPE *a, DTYPE *b, DTYPE length);

//r=a+b a_length>=b_length by default, r_length=a_length+1
void BN_add(DTYPE *r, DTYPE r_length, DTYPE *a, DTYPE a_length, DTYPE *b, DTYPE b_length);

//r=a-b a_length>=b_length by default, r_length=a_length
void BN_sub(DTYPE *r, DTYPE r_length, DTYPE *a, DTYPE a_length, DTYPE *b, DTYPE b_length);

//r=a*a r_length=a_length*2
void BN_square(DTYPE *r, DTYPE *a, DTYPE a_length);

//r=a*b r_length=a_leng+b_length
void BN_mul(DTYPE *r, DTYPE *a, DTYPE a_length, DTYPE *b, DTYPE b_length);

//r=a(mod n) r_length=n_length, a_length>=n_length
void BN_mod(DTYPE *r, DTYPE *a, DTYPE a_length, DTYPE *n, DTYPE n_length);

//r=a*b(mod n) default r_length=n_length
//method - multiply then reduce
void BN_mod_mul(DTYPE *r, DTYPE *a, DTYPE a_length, DTYPE *b, DTYPE b_length, DTYPE *n, DTYPE n_length);

//r=a*b(mod n) using Blakley's Method
void BN_mod_mul2(DTYPE *r, DTYPE *a, DTYPE a_length, DTYPE *b, DTYPE b_length, DTYPE *n, DTYPE n_length);

//compute a_inv, a_inv*a=1(mod n)
void BN_mod_inv(DTYPE *a_inv, DTYPE *a, DTYPE a_length, DTYPE *n, DTYPE n_length);

//compute a_inv(a DTYPE), a_inv*a=1(mod 2^WORD_SIZE)
void BN_MonMod_inv(DTYPE *a_inv, DTYPE *a, DTYPE a_length);

//compute Montgomery Product of a and b with modulus n
void BN_MonPro(DTYPE *r, DTYPE *a, DTYPE a_length, DTYPE *b, DTYPE b_length, DTYPE *n, DTYPE n_length, DTYPE n_inv);

//r=a^e(mod n) using Montgomery Method
void BN_MonExp(DTYPE *result, DTYPE *a, DTYPE a_length, DTYPE *e, DTYPE e_length, DTYPE *n, DTYPE n_length, DTYPE *r, DTYPE inv);

//a_length=e_length=n_length=MAX_PRIME_LENGTH
void BN_pow_mod(DTYPE *result, DTYPE *a, DTYPE *e, DTYPE *n);

void signBN_init(sm *r, DTYPE len);

void BN_signAddShift(sm *a, DTYPE *n, DTYPE *tmp);

void signBN_init_assign(sm *r, DTYPE *a, DTYPE a_length);

//r->length=a->length=c->length
void BN_sign_sub(sm *r, sm *a, sm *c);