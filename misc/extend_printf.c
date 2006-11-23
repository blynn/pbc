/*
 * Behaves as gmp_printf with new conversion specifier %B for element_t types
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <gmp.h>
#include "pbc_utils.h"
#include "pbc_field.h"
#include "pbc_memory.h"

//TODO: remove repeated code for error handling 
int element_vfprintf(FILE *stream, const char *format, va_list ap)
{
    int count = 0, status;
    char ch;
    char *copy, *c, *start, *next;
    element_ptr e;
    int found;

    copy = strdup(format);
    start = next = copy;

    for(;;) {
	for(;;) {
	    c = strchr(next, '%');
	    if (!c) {
		status = fprintf(stream, start);
		if (status < 0) {
		    count = -1;
		} else count += status;
		goto done;
	    }
	    if (!*(c + 1)) goto done;
	    if (!*(c + 1) != '%') break;
	    next = c + 2;
	}
	*c = 0;
	status = fprintf(stream, start);
	if (status < 0) {
	    count = -1;
	    goto done;
	} else count += status;
	*c = '%';
	start = c;
	found = 0;
	while(!found) {
	    c++;
	    switch (*c) {
		case '\0':
		    goto done;
		case 'B':
		    e = va_arg(ap, element_ptr);
		    status = element_out_str(stream, 0, e);
		    if (status < 0) {
			count = -1;
			goto done;
		    } else count += status;
		    found = 1;
		    break;
		/*
		case 'Z':
		    z = va_arg(ap, mpz_ptr);
		    ch = *(c+1);
		    *(c+1) = '\0';
		    printf("format string '%s'\n",start);
		    status = gmp_fprintf(stream, start, z);
		    if (status < 0) {
			count = -1;
			goto done;
		    } else count += status;
		    *(c+1) = ch;
		    found = 1;
		    break;
		    */
		default:
		    if (strchr("diouxXeEfFgGaAcspnm", *c)) {
			void *ptr = va_arg(ap, void *);
			ch = *(c+1);
			*(c+1) = '\0';
			status = gmp_fprintf(stream, start, ptr);
			if (status < 0) {
			    count = -1;
			    goto done;
			} else count += status;
			*(c+1) = ch;
			found = 1;
		    }
		    break;
	    }
	}
	next = start = c + 1;
    }

done:
    pbc_free(copy);

    return count;
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
