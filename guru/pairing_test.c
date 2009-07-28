#include <string.h>
#include "pbc.h"
#include "pbc_test.h"

// TODO: Probably better off as a test for the calculator.
static void a_test(void) {
  char buf[1024];
  char *gstr =
"[2382389466570123849673299401984867521337122094157231907755149435707124249269394670242462497382963719723036281844079382411446883273020125104982896098602669, 2152768906589770702756591740710760107878949212304343787392475836859241438597588807103470081101790991563152395601123682809718038151417122294066319979967168]";
  char *hstr =
"[5832612417453786541700129157230442590988122495898645678468800815872828277169950107203266157735206975228912899931278160262081308603240860553459187732968543, 5825590786822892934138376868455818413990615826926356662470129700411774690868351658310187202553513693344017463065909279569624651155563430675084173630054336]";
  char *astr = "171583727262251826931173602797951212789946235851";
  char *bstr = "233634857565210859330459959563397971304462340857";
  element_t g, h;
  element_t a, b;
  element_t ga, hb;
  element_t x1, x2;
  pairing_t pairing;
  pairing_init_set_str(pairing,
"type a\n"
"q 8780710799663312522437781984754049815806883199414208211028653399266475630880222957078625179422662221423155858769582317459277713367317481324925129998224791\n"
"h 12016012264891146079388821366740534204802954401251311822919615131047207289359704531102844802183906537786776\n"
"r 730750818665451621361119245571504901405976559617\n"
"exp2 159\n"
"exp1 107\n"
"sign1 1\n"
"sign0 1\n"
  );
  element_init_G1(g, pairing);
  element_init_G1(ga, pairing);
  element_init_G2(h, pairing);
  element_init_G2(hb, pairing);
  element_init_Zr(a, pairing);
  element_init_Zr(b, pairing);
  element_init_GT(x1, pairing);
  element_init_GT(x2, pairing);

  element_set_str(g, gstr, 0);
  element_snprint(buf, 1024, g);
  EXPECT(!strcmp(gstr, buf));

  element_set_str(h, hstr, 0);
  element_snprint(buf, 1024, h);
  EXPECT(!strcmp(hstr, buf));

  element_set_str(a, astr, 0);
  element_snprint(buf, 1024, a);
  EXPECT(!strcmp(astr, buf));

  element_set_str(b, bstr, 0);
  element_snprint(buf, 1024, b);
  EXPECT(!strcmp(bstr, buf));

  pairing_apply(x1, g, h, pairing);
  element_snprint(buf, 1024, x1);
  EXPECT(!strcmp(buf,
"[1352478452661998164151215014828915385601138645645403926287105573769451214277485326392786454433874957123922454604362337349978217917242114505658729401276644, 2809858014072341042857607405424304552357466023841122154308055820747972163307396014445308786731013691659356362568425895483877936945589613445697089590886519]"
  ));

  element_pow_zn(ga, g, a);
  element_snprint(buf, 1024, ga);
  EXPECT(!strcmp(buf,
"[3727290142167731134589933003026410141163353118002821914170365887139605219852868537686435214464927363733592858325260588072422405672197113236445369761687270, 8313413520789037477320458888316489483781506373846006723006557775349684878102042826049292521482530556981023752851151672326421296204733037418468523296005577]"
  ));

  element_pow_zn(hb, h, b);
  element_snprint(buf, 1024, hb);
  EXPECT(!strcmp(buf,
"[302169045606583472168811217560382970305157511680176350745436990853463473855962841196184541109617397027480204774682450915021848512168573082843355648090809, 7428193877404140917518137438384425427600294220905786853638038223349096573857683866658575603565175187399696035468569929483731011292133989973187846752806084]"
  ));

  pairing_apply(x2, ga, hb, pairing);
  element_snprint(buf, 1024, x2);
  EXPECT(!strcmp(buf,
"[5401677742232403160612802517983583823254857216272776607059355607024091426935935872461700304196658606704085604766577186374528948004140797833341187234647180, 4255900207739859478558185000995524505026245539159946661271849714832846423204570340979120001638894488614502770175520505048836617405342161594891740961421000]"
  ));

  element_pow_zn(x1, x1, a);
  element_pow_zn(x1, x1, b);
  element_snprint(buf, 1024, x1);
  EXPECT(!strcmp(buf,
"[5401677742232403160612802517983583823254857216272776607059355607024091426935935872461700304196658606704085604766577186374528948004140797833341187234647180, 4255900207739859478558185000995524505026245539159946661271849714832846423204570340979120001638894488614502770175520505048836617405342161594891740961421000]"
  ));

  EXPECT(!element_cmp(x1, x2));

  element_clear(g);
  element_clear(h);
  element_clear(a);
  element_clear(b);
  element_clear(ga);
  element_clear(hb);
  element_clear(x1);
  element_clear(x2);
  pairing_clear(pairing);
}

void oldtest(void) {
  element_t g, h;
  element_t x1, x2;
  element_t zg, zh, z;
  pairing_t pairing;
  pbc_param_t param;

  pbc_param_init_a_gen(param, 160, 512);
  pairing_init_pbc_param(pairing, param);
  pbc_param_clear(param);

  element_init_G1(g, pairing);
  element_init_G1(zg, pairing);
  element_init_G2(h, pairing);
  element_init_G2(zh, pairing);
  element_init_GT(x1, pairing);
  element_init_GT(x2, pairing);
  element_init_Zr(z, pairing);
  element_random(g);
  element_random(h);
  element_printf("g = %B\n", g);
  element_printf("h = %B\n", h);
  pairing_apply(x1, g, h, pairing);
  element_printf("f(g, h) = %B\n", x1);

  element_random(z);
  element_printf("z = %B\n", z);

  element_pow_zn(x1, x1, z);
  element_printf("f(g, h)^z = %B\n", x1);

  element_pow_zn(zg, g, z);
  element_printf("g^z = %B\n", zg);
  pairing_apply(x2, zg, h, pairing);
  element_printf("f(g^z, h) = %B\n", x2);

  EXPECT(!element_cmp(x1, x2));

  element_pow_zn(zh, h, z);
  element_printf("h^z = %B\n", zh);
  pairing_apply(x2, g, zh, pairing);
  element_printf("f(g, h^z) = %B\n", x2);

  EXPECT(!element_cmp(x1, x2));

  {
    int i;
    int len = element_length_in_bytes(h);
    printf("length_in_bytes(h) = %d\n", len);
    unsigned char *data = pbc_malloc(len);
    element_to_bytes(data, h);
    for (i=0; i<len; i++) {
      printf(" %02X", data[i]);
      if (15 == (i % 16)) printf("\n");
    }
    printf("\n");
    element_from_bytes(h, data);
    element_printf("from_bytes h = %B\n", h);
    free(data);
  }

  element_clear(g);
  element_clear(h);
  element_clear(x1);
  element_clear(x2);
  element_clear(zg);
  element_clear(zh);
  element_clear(z);
  pairing_clear(pairing);
}

int main(void) {
  a_test();
  oldtest();
  return pbc_err_count;
}
