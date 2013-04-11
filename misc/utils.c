#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h> // for intptr_t
#include <gmp.h>

#include "pbc_utils.h"
#include "pbc_field.h"

static int pbc_msg_to_stderr = 1;

int pbc_set_msg_to_stderr(int i) {
  return pbc_msg_to_stderr = i;
}

static int out(const char *format, ...) {
  if (!pbc_msg_to_stderr) return 0;
  va_list params;

  va_start(params, format);
  int res = vfprintf(stderr, format, params);
  va_end(params);
  return res;
}

static void print_warning(void) {
  static int first = 1;
  if (first) {
    out("*** PBC asserts enabled: potential performance penalties ***\n");
    first = 0;
  }
}

void pbc_assert(int expr, char *msg, const char *func) {
  print_warning();
  if (!expr) {
    out("PBC assert failed: %s(): %s\n", func, msg);
    abort();
  }
}

void pbc_assert_match2(element_ptr a, element_ptr b, const char *func) {
  print_warning();
  if (a->field != b->field) {
    out("PBC assert failed: %s(): field mismatch\n", func);
    abort();
  }
}

void pbc_assert_match3(element_ptr a, element_ptr b, element_ptr c,
    const char *func) {
  print_warning();
  if (a->field != b->field) {
    out("PBC assert failed: %s(): first two args field mismatch\n", func);
    abort();
  }
  if (b->field != c->field) {
    out("PBC assert failed: %s(): last two args field mismatch\n", func);
    abort();
  }
}

// Print at most the first 1024 bytes of an error message.
static void report(const char *prefix, const char *err, va_list params) {
  char msg[1024];
  element_vsnprintf(msg, sizeof(msg), err, params);
  out("%s%s\n", prefix, msg);
}

void pbc_die(const char *err, ...) {
  va_list params;

  va_start(params, err);
  report("fatal: ", err, params);
  va_end(params);
  exit(128);
}

void pbc_info(const char *err, ...) {
  va_list params;

  va_start(params, err);
  report("", err, params);
  va_end(params);
}

void pbc_warn(const char *err, ...) {
  va_list params;

  va_start(params, err);
  report("warning: ", err, params);
  va_end(params);
}

void pbc_error(const char *err, ...) {
  va_list params;

  va_start(params, err);
  report("error: ", err, params);
  va_end(params);
}
