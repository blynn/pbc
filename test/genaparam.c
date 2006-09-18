#include "pbc.h"
#include "a_param.h"

int main(void)
{
    a_param_t ap;

    a_param_init(ap);
    a_param_gen(ap, 160, 512);

    a_param_out_str(stdout, ap);
    a_param_clear(ap);

    return 0;
}
