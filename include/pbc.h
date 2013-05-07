#ifndef __PBC_H__
#define __PBC_H__

#include <stdarg.h>
#include <stdio.h>
#include <stdint.h> // for intptr_t
#include <stdlib.h>
#include <gmp.h>

#if defined (__cplusplus)
extern "C" {
#endif

#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_param.h"
#include "pbc_pairing.h"
#include "pbc_curve.h"
#include "pbc_mnt.h"
#include "pbc_a1_param.h"
#include "pbc_a_param.h"
#include "pbc_d_param.h"
#include "pbc_e_param.h"
#include "pbc_f_param.h"
#include "pbc_g_param.h"
#include "pbc_i_param.h"
#include "pbc_random.h"
#include "pbc_memory.h"

#if defined (__cplusplus)
}  // extern "C"
#endif

#endif //__PBC_H__
