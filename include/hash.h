#ifndef HASH_H
#define HASH_H

typedef struct {
    unsigned long state[5];
    unsigned long count[2];
    unsigned char buffer[64];
} SHA1_CTX;

struct hash_ctx_s {
    SHA1_CTX context;
};
typedef struct hash_ctx_s hash_ctx_t[1];
typedef struct hash_ctx_s *hash_ctx_ptr;

void hash_init(hash_ctx_t context);
void hash_update(hash_ctx_t context, unsigned char *msg, unsigned int len);
void hash_final(unsigned char *digest, hash_ctx_t context);

enum {
    hash_length = 20,
};
#endif //HASH_H
