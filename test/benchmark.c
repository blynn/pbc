#include "pbc.h"
#include "get_time.h"

int main(void)
{
    pairing_t pairing;
    element_t x, y, r;
    int i, n;
    double t0, t1, ttotal;

    pairing_init_inp_str(pairing, stdin);
    element_init(x, pairing->G1);
    element_init(y, pairing->G2);
    element_init(r, pairing->GT);

    n = 10;
    ttotal = 0.0;
    for (i=0; i<n; i++) {
	element_random(x);
	element_random(y);

	element_printf("x = %B\n", x);
	element_printf("y = %B\n", y);
	t0 = get_time();
	bilinear_map(r, x, y, pairing);
	t1 = get_time();
	ttotal += t1 - t0;
	element_printf("e(x,y) = %B\n", r);
    }
    printf("average pairing time = %f\n", ttotal / n);
    element_clear(x);
    element_clear(y);
    element_clear(r);

    return 0;
}
