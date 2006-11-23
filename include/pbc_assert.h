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

#endif //__PBC_ASSERT_H__
