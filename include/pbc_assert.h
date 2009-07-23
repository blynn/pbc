#ifndef __PBC_ASSERT_H__
#define __PBC_ASSERT_H__

#ifdef PBC_DEBUG

#define PBC_ASSERT(expr, msg) (pbc_assert(expr, msg, __func__))
#define PBC_ASSERT_MATCH2(a, b) (pbc_assert_match2(a, b, __func__))
#define PBC_ASSERT_MATCH3(a, b, c) (pbc_assert_match3(a, b, c, __func__))

#else

#define PBC_ASSERT(expr, msg) ((void) (0))
#define PBC_ASSERT_MATCH2(a, b) ((void) (0))
#define PBC_ASSERT_MATCH3(a, b, c) ((void) (0))

#endif

// Based on Git code.

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

#endif //__PBC_ASSERT_H__
