#ifndef __PBC_UTILS_H__
#define __PBC_UTILS_H__

#ifdef PBC_DEBUG

#define PBC_ASSERT(expr, msg) (pbc_assert(expr, msg, __func__))
#define PBC_ASSERT_MATCH2(a, b) (pbc_assert_match2(a, b, __func__))
#define PBC_ASSERT_MATCH3(a, b, c) (pbc_assert_match3(a, b, c, __func__))

#else

#define PBC_ASSERT(expr, msg) ((void) (0))
#define PBC_ASSERT_MATCH2(a, b) ((void) (0))
#define PBC_ASSERT_MATCH3(a, b, c) ((void) (0))

#endif

// die, warn and info based on Git code.

// Print error message to standard error and exit with code 128.
void pbc_die(const char *err, ...)
    __attribute__((__noreturn__))
    __attribute__((format (printf, 1, 2)));

// Print info message to standard error.
void pbc_info(const char *err, ...)
    __attribute__((format (printf, 1, 2)));

// Print warning message to standard error.
void pbc_warn(const char *err, ...)
    __attribute__((format (printf, 1, 2)));

// Print error message to standard error.
void pbc_error(const char *err, ...)
    __attribute__((format (printf, 1, 2)));

#ifndef UNUSED_VAR
#if defined(__GNUC__)
// We could use __attribute__((unused)) instead.
#define UNUSED_VAR(a) (void) a
#else
// From the ACE project: http://www.cs.wustl.edu/~schmidt/ACE.html
// silences warnings, and generates no code for many compilers
// See ACE_wrappers/ace/ace/config-macros.h:391
//
// Not anymore: gcc no longer likes it -blynn
#define UNUSED_VAR(a) do { /* nothing */ } while (&a == 0)
#endif
#endif

// For storing small integers in void *
// see http://www.gelato.unsw.edu.au/archives/linux-ia64/0008/0344.html
static inline void *int_to_voidp(int i) {
  // TODO: This won't work on some platforms.
  return (void *) (long) i;
}

#endif //__PBC_UTILS_H__
