#include "keys.h"

void rsa_encrypt(DTYPE *cipher, DTYPE cipher_len, DTYPE *msg, DTYPE msg_len, rsa_pk *public_key);

void rsa_decrypt(DTYPE *output, DTYPE output_len, DTYPE *cipher, DTYPE cipher_len, rsa_sk *secret_key);