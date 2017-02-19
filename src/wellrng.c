
#include <stdio.h>
#include <stdlib.h>

#include <rjmcmc/wellrng.h>

#define W 32
#define R 1391
#define P 15
#define MASKU (0xffffffffE >> (W - P))
#define MASKL (~MASKU)

#define M1 23
#define M2 481
#define M3 229

#define TEMPERB 0x93dd1400U
#define TEMPERC 0xfa118000U

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define MAT1(v) v
#define MAT2(a,v) ((v & 1U)?((v>>1)^a):(v>>1))
#define MAT3POS(t,v) (v>>t)
#define MAT3NEG(t,v) (v<<(-(t)))
#define MAT4POS(t,b,v) (v ^ ((v>>  t ) & b))
#define MAT4NEG(t,b,v) (v ^ ((v<<(-(t))) & b))
#define MAT5(r,a,ds,dt,v) ((v & dt)?((((v<<r)^(v>>(W-r)))&ds)^a):(((v<<r)^(v>>(W-r)))&ds))
#define MAT7(v) 0

#define FACT 2.32830643653869628906e-10


struct wellrng {
  unsigned int s[R];
  int i;
  int state;
};

static void kiss_init(unsigned int seed);
static unsigned int kiss_sample(void);
static struct {
  unsigned int x;
  unsigned int y;
  unsigned int z;
  unsigned int c;
} kiss_state = {
  123456789,
  362436000,
  521288629,
  7654321
};

wellrng_t *
wellrng_init_simple(unsigned int seed)
{
  int i;
  wellrng_t *r;

  kiss_init(seed);

  r = (wellrng_t*)malloc(sizeof(struct wellrng));
  if (r == NULL) {
    return NULL;
  }
  
  for (i = 0; i < R; i ++) {
    r->s[i] = kiss_sample();
  }

  r->state = 1;
  r->i = 0;
  return r;
}

unsigned int
wellrng_seed_size(void) 
{
  return R;
}

wellrng_t *
wellrng_init_direct(unsigned int *seed)
{
  int i;
  wellrng_t *r;

  r = (wellrng_t*)malloc(sizeof(struct wellrng));
  if (r == NULL) {
    return NULL;
  }
  
  for (i = 0; i < R; i ++) {
    r->s[i] = seed[i];
  }

  r->state = 1;
  r->i = 0;
  return r;
}

double 
wellrng_sample(wellrng_t *w)
{
#define V0            w->s[w->i]
#define VM1Over       w->s[w->i+M1-R]
#define VM1           w->s[w->i+M1]
#define VM2Over       w->s[w->i+M2-R]
#define VM2           w->s[w->i+M2]
#define VM3Over       w->s[w->i+M3-R]
#define VM3           w->s[w->i+M3]
#define Vrm1          w->s[w->i-1]
#define Vrm1Under     w->s[w->i+R-1]
#define Vrm2          w->s[w->i-2]
#define Vrm2Under     w->s[w->i+R-2]

#define newV0         w->s[w->i-1]
#define newV0Under    w->s[w->i-1+R]
#define newV1         w->s[w->i]
#define newVRm1       w->s[w->i-2]
#define newVRm1Under  w->s[w->i-2+R]

  unsigned int z0;
  unsigned int z1;
  unsigned int z2;
  unsigned int y;

  switch (w->state) {

  case 1:
    z0 = (Vrm1Under & MASKL) | (Vrm2Under & MASKU);
    z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
    z2 = MAT0NEG(-10,VM2) ^ MAT3NEG(-26,VM3);
    newV1  = z1 ^ z2;
    newV0Under = 
      MAT1(z0) ^ 
      MAT0POS(20,z1) ^  
      MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ 
      MAT1(newV1);
    w->i = R-1;
    w->state = 3;
    y = w->s[w->i] ^ ((w->s[w->i] << 7) & TEMPERB);
    y = y ^ ((y << 15) & TEMPERC);
    break;

  case 2:
    z0 = (Vrm1 & MASKL) | (Vrm2Under & MASKU);
    z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
    z2 = MAT0NEG(-10,VM2) ^ MAT3NEG(-26,VM3);
    newV1 = z1 ^ z2;
    newV0 =  MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
    w->i = 0;
    w->state = 1;
    y = w->s[w->i] ^ ((w->s[w->i] << 7) & TEMPERB);
    y = y ^ ((y << 15) & TEMPERC);
    break;

  case 3:
    z0 = (Vrm1 & MASKL) | (Vrm2 & MASKU);
    z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1Over);
    z2 = MAT0NEG(-10,VM2Over) ^ MAT3NEG(-26,VM3Over);
    newV1 = z1 ^ z2;
    newV0 = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
    w->i --;
    if((w->i + M1) < R) {
      w->state = 4;
    }

    y = w->s[w->i] ^ ((w->s[w->i] << 7) & TEMPERB);
    y = y ^ ((y << 15) & TEMPERC);
    break;

  case 4:
    z0 = (Vrm1 & MASKL) | (Vrm2 & MASKU);
    z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
    z2 = MAT0NEG(-10,VM2Over) ^ MAT3NEG(-26,VM3Over);
    newV1 = z1 ^ z2;
    newV0 = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
    w->i--;
    if ((w->i + M3) < R) {
      w->state = 5;
    }

    y = w->s[w->i] ^ ((w->s[w->i] << 7) & TEMPERB);
    y = y ^ ((y << 15) & TEMPERC);
    break;

  case 5:
    z0 = (Vrm1 & MASKL) | (Vrm2 & MASKU);
    z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
    z2 = MAT0NEG(-10,VM2Over) ^ MAT3NEG(-26,VM3);
    newV1 = z1 ^ z2;
    newV0 = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
    w->i--;
    if((w->i + M2) < R) {
      w->state = 6;
    }
    y = w->s[w->i] ^ ((w->s[w->i] << 7) & TEMPERB);
    y = y ^ ((y << 15) & TEMPERC);
    break;

  case 6:
    z0 = (Vrm1 & MASKL) | (Vrm2 & MASKU);
    z1 = MAT0NEG(-24,V0) ^ MAT0POS(30,VM1);
    z2 = MAT0NEG(-10,VM2) ^ MAT3NEG(-26,VM3);
    newV1 = z1 ^ z2;
    newV0 = MAT1(z0) ^ MAT0POS(20,z1) ^ MAT5(9,0xb729fcecU,0xfbffffffU,0x00020000U,z2) ^ MAT1(newV1);
    w->i--;
    if(w->i == 1) {
      w->state = 2;
    }
    y = w->s[w->i] ^ ((w->s[w->i] << 7) & TEMPERB);
    y = y ^ ((y << 15) & TEMPERC);
    break;

  default:
    fprintf(stderr, "wellrng_sample: invalid state %d\n", w->state);
    return -1.0;
  }

  return ((double) y * FACT);    
}

void
wellrng_destroy(wellrng_t *w)
{
  free(w);
}


static void kiss_init(unsigned int seed)
{
  kiss_state.x = seed;
  kiss_state.y = (seed + 1) * 8191;
  kiss_state.z = ~seed;
  kiss_state.c = (seed + 1) * 127;
}

static unsigned int kiss_sample(void)
{
  unsigned long long t, a = 698769069ULL;

  kiss_state.x = 69069*kiss_state.x + 12345;

  kiss_state.y ^= (kiss_state.y << 13);
  kiss_state.y ^= (kiss_state.y >> 17);
  kiss_state.y ^= (kiss_state.y << 5);

  t = a*kiss_state.y + kiss_state.c;
  kiss_state.c = t >> 32;

  kiss_state.z = t;
  return kiss_state.x + kiss_state.y + kiss_state.z;
}
