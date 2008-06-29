#ifndef __PBC_UTIL_H__
#define __PBC_UTIL_H__

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
static inline void *int_to_voidp(int i)
{
    //TODO: this won't work on some platforms 
    //assert(sizeof(long) == sizeof(void *));
    return (void *) (long) i;
}

#endif //__PBC_UTIL_H__
