#include "pbc.h"
#include "e_param.h"

int main(void)
{
    e_param_t ep;

    e_param_init(ep);
    e_param_gen(ep, 160, 1024);
    e_param_out_str(stdout, ep);

    return 0;
}
