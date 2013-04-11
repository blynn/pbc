#include <stdint.h> // for intptr_t
#include "pbc.h"
#include "pbc_test.h"

/* I've heard that sometimes automatic garbage collection can outperform
 * manual collection, so I briefly tried using the Boehm-Demers-Weiser GC
 * library. Both GMP and PBC support custom memory allocation routines so
 * incorporating the GC library is trivial.
 *
 * Automatic garbage collection appears to slow this program down a little,
 * even if only PBC collects automatically. (The case where PBC collects
 * manually but GMP collects automatically cannot be achieved with the GC
 * library because PBC objects point at GMP objects.)
 *
 * Perhaps specially-tailored memory allocation routines could shave off
 * some time, but one would have to thoroughly analyze PBC and GMP memory usage
 * patterns.
 *
 * Below is the commented-out code that collects garbage for PBC. Of course,
 * if you want to use it you must also tell the build system where to find
 * gc.h and to link with the GC library.
 *
 * Also, you may wish to write similar code for GMP (which I unfortunately
 * deleted before thinking that it might be useful for others).
 * Note GC_MALLOC_ATOMIC may be used for GMP since the mpz_t type does not
 * store pointers in the memory it allocates.
 *
 * The malloc and realloc functions should exit on failure but I didn't
 * bother since I was only seeing if GC could speed up this program.

#include <gc.h>
#include <pbc_utils.h>

void *gc_alloc(size_t size) {
  return GC_MALLOC(size);
}

void *gc_realloc(void *ptr, size_t size) {
  return GC_REALLOC(ptr, size);
}

void gc_free(void *ptr) {
  UNUSED_VAR(ptr);
}

 * The following should be the first two statements in main()

GC_INIT();
pbc_set_memory_functions(gc_alloc, gc_realloc, gc_free);

 */

int main(int argc, char **argv) {
  pairing_t pairing;
  element_t x, y, r, r2;
  int i, n;
  double t0, t1, ttotal, ttotalpp;
  pairing_pp_t pp;

  // Cheat for slightly faster times:
  // pbc_set_memory_functions(malloc, realloc, free);

  pbc_demo_pairing_init(pairing, argc, argv);

  element_init_G1(x, pairing);
  element_init_G2(y, pairing);
  element_init_GT(r, pairing);
  element_init_GT(r2, pairing);

  n = 10;
  ttotal = 0.0;
  ttotalpp = 0.0;
  for (i=0; i<n; i++) {
    element_random(x);
    element_random(y);

    pairing_pp_init(pp, x, pairing);
    t0 = pbc_get_time();
    pairing_pp_apply(r, y, pp);
    t1 = pbc_get_time();
    ttotalpp += t1 - t0;
    pairing_pp_clear(pp);

    t0 = pbc_get_time();

    element_pairing(r2, x, y);
    t1 = pbc_get_time();
    ttotal += t1 - t0;

    element_printf("x = %B\n", x);
    element_printf("y = %B\n", y);
    element_printf("e(x,y) = %B\n", r);
    if (element_cmp(r, r2)) {
      printf("BUG!\n");
      exit(1);
    }
  }
  printf("average pairing time = %f\n", ttotal / n);
  printf("average pairing time (preprocessed) = %f\n", ttotalpp / n);

  element_clear(x);
  element_clear(y);
  element_clear(r);
  element_clear(r2);

  pairing_clear(pairing);

  return 0;
}
