// Read pairing from file specified as first argument, or from standard input
// if there is no first argument.
static inline void demo_get_pairing(pairing_t pairing, int argc, char **argv) {
  char s[1024];
  FILE *fp = stdin;

  if (argc > 1) {
    fp = fopen(argv[1], "r");
    if (!fp) pbc_die("error opening %s", argv[1]);
  }
  if (!fread(s, 1, 1024, fp)) pbc_die("error reading pairing");
  fclose(fp);

  pairing_init_set_str(pairing, s);
}
