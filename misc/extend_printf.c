#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <printf.h>
#include "utils.h"
#include "field.h"

static int print_element(FILE *stream, const struct printf_info *info,
	const void *const *args)
{
    UNUSED_VAR(info);
    return element_out_str(stream, 0, *((const element_ptr *) args[0]));
}

static int print_element_arginfo(const struct printf_info *info, size_t n,
	int *argtypes)
{
    UNUSED_VAR(info);
    if (n > 0) argtypes[0] = PA_POINTER;
    return 1;
}

static int element_printf_init(void)
{
    return register_printf_function('B', print_element, print_element_arginfo);
}

static int printf_init_flag = 0;

int element_fprintf(FILE *stream, const char *format, ...)
{
    int status;
    va_list ap;

    if (!printf_init_flag) element_printf_init();
    va_start(ap, format);
    status = vfprintf(stream, format, ap);
    va_end(ap);
    return status;
}

int element_printf(const char *format, ...)
{
    int status;
    va_list ap;

    if (!printf_init_flag) element_printf_init();
    va_start(ap, format);
    status = vprintf(format, ap);
    va_end(ap);
    return status;
}
