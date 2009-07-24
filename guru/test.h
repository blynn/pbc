#define EXPECT(cond) \
  if (cond); else errcount++, fprintf(stderr, "\n*** FAIL *** %s:%d: %s\n\n", __FILE__, __LINE__, #cond)

int errcount;
