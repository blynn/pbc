#ifndef __PBC_UTILS_H__
#define __PBC_UTILS_H__

#ifdef PBC_DEBUG

/*@manual debug
Macro: if `expr` evaluates to 0, print `msg` and exit.
*/
#define PBC_ASSERT(expr, msg) \
    (pbc_assert(expr, msg, __func__))

/*@manual debug
Macro: if elements `a` and `b` are from different fields then exit.
*/
#define PBC_ASSERT_MATCH2(a, b) \
    (pbc_assert_match2(a, b, __func__))

/*@manual debug
Macro: if elements `a`, `b` and `c` are from different fields then exit.
*/
#define PBC_ASSERT_MATCH3(a, b, c) \
    (pbc_assert_match3(a, b, c, __func__))

#else

#define PBC_ASSERT(expr, msg) ((void) (0))
#define PBC_ASSERT_MATCH2(a, b) ((void) (0))
#define PBC_ASSERT_MATCH3(a, b, c) ((void) (0))

#endif

// die, warn and info based on Git code.

/*@manual log
By default error messages are printed to standard error.
Call `pbc_set_msg_to_stderr(0)` to suppress messages.
*/
int pbc_set_msg_to_stderr(int i);

/*@manual log
Reports error message and exits with code 128.
*/
void pbc_die(const char *err, ...)
    __attribute__((__noreturn__))
    __attribute__((format (printf, 1, 2)));

/*@manual log
Reports informational message.
*/
void pbc_info(const char *err, ...)
    __attribute__((format (printf, 1, 2)));

/*@manual log
Reports warning message.
*/
void pbc_warn(const char *err, ...)
    __attribute__((format (printf, 1, 2)));

/*@manual log
Reports error message.
*/
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
// C99 standard introduced the intptr_t and uintptr_t types,
// guaranteed to be able to hold pointers
static inline void *int_to_voidp(intptr_t i) {
  return (void *) i;
}

#endif //__PBC_UTILS_H__
