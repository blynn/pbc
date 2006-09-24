/*
 * Extends printf with new conversion specifiers:
 * * 'B' for element_t,
 * * 'Y' for mpz_t (prefer 'Z' but this used in Linux libc5? [man 3 printf])
 *
 * element_printf, element_fprintf call printf after changing
 * occurrences of %Z to %Y, so that %Z can be used for mpz_t as in GMP.
 *
 * Cons of current approach:
 * * needs GNU C extension
 * * printf behaves differently
 * * possible spurious warnings if printf() is used
 * * gcc does not check format string if element_printf() is used.
 * * doesn't behave the same as gmp_printf on GMP types
 */
//TODO (eventually): write element_printf() explicitly, i.e. without using GNU C
//features and so that printf() is not affected,
//and I can directly use %Z for mpz_t

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <printf.h>
#include <gmp.h>
#include "utils.h"
#include "field.h"
#include "strclone.h"

static int print_element(FILE *stream, const struct printf_info *info,
	const void *const *args)
{
    UNUSED_VAR(info);
    return element_out_str(stream, 0, *((const element_ptr *) args[0]));
}

static int pointer_arginfo(const struct printf_info *info, size_t n,
	int *argtypes)
{
    UNUSED_VAR(info);
    if (n > 0) argtypes[0] = PA_POINTER;
    return 1;
}

static int print_mpz(FILE *stream, const struct printf_info *info,
	const void *const *args)
{
    UNUSED_VAR(info);
    return mpz_out_str(stream, 0, *((const mpz_ptr *) args[0]));
}

static int element_printf_init(void)
{
    if (register_printf_function('B', print_element, pointer_arginfo)) {
	return -1;
    }
    return register_printf_function('Y', print_mpz, pointer_arginfo);
}

static int printf_init_flag = 0;

int element_vfprintf(FILE *stream, const char *format, va_list ap)
{
    int status;
    char *copy, *c;

    copy = strclone(format);
    c = copy;
    while (*c) {
	if (*c == '%') {
	    c++;
	    if (!*c) break;
	    if (*c == 'Z') *c = 'Y';
	}
	c++;
    }

    if (!printf_init_flag) element_printf_init();
    status = vfprintf(stream, copy, ap);

    free(copy);
    return status;
}

int element_fprintf(FILE *stream, const char *format, ...)
{
    int status;
    va_list ap;

    va_start(ap, format);
    status = element_vfprintf(stream, format, ap);
    va_end(ap);
    return status;
}

int element_printf(const char *format, ...)
{
    int status;
    va_list ap;

    va_start(ap, format);
    status = element_vfprintf(stdout, format, ap);
    va_end(ap);
    return status;
}
