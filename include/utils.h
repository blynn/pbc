#ifndef __UTIL_H__
#define __UTIL_H__

// from the ACE project: http://www.cs.wustl.edu/~schmidt/ACE.html
// silences warnings, and generates no code for many compilers
// See ACE_wrappers/ace/ace/config-macros.h:391
#define UNUSED_VAR(a) do { /* nothing */ } while (&a == 0)

// I was using this before (and without macros):
//#define UNUSED_VAR(a) (void) a

#endif
