// Useful for tests.

// Read pairing from file specified as first argument, or from standard input
// if there is no first argument.
static inline void demo_get_pairing(pairing_t pairing, int argc, char **argv) {
  char s[1024];
  FILE *fp = stdin;

  if (argc > 1) {
    fp = fopen(argv[1], "r");
    if (!fp) pbc_die("error opening %s", argv[1]);
  }
  size_t count = fread(s, 1, 1024, fp);
  if (!count) pbc_die("error reading pairing");
  fclose(fp);

  pairing_init_set_buf(pairing, s, count);
}

#define EXPECT(cond) \
  if (cond); else pbc_err_count++, fprintf(stderr, "\n*** FAIL ***\n  %s:%d: %s\n\n", __FILE__, __LINE__, #cond)

int pbc_err_count;
