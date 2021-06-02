#include <stdint.h>
#include <stdio.h>
#define QINV -3327 // q^-1 mod 2^16
#define N 256
#define KYBER_Q 3329

const int16_t zetas[128] = {
  -1044,  -758,  -359, -1517,  1493,  1422,   287,   202,
   -171,   622,  1577,   182,   962, -1202, -1474,  1468,
    573, -1325,   264,   383,  -829,  1458, -1602,  -130,
   -681,  1017,   732,   608, -1542,   411,  -205, -1571,
   1223,   652,  -552,  1015, -1293,  1491,  -282, -1544,
    516,    -8,  -320,  -666, -1618, -1162,   126,  1469,
   -853,   -90,  -271,   830,   107, -1421,  -247,  -951,
   -398,   961, -1508,  -725,   448, -1065,   677, -1275,
  -1103,   430,   555,   843, -1251,   871,  1550,   105,
    422,   587,   177,  -235,  -291,  -460,  1574,  1653,
   -246,   778,  1159,  -147,  -777,  1483,  -602,  1119,
  -1590,   644,  -872,   349,   418,   329,  -156,   -75,
    817,  1097,   603,   610,  1322, -1285, -1465,   384,
  -1215,  -136,  1218, -1335,  -874,   220, -1187, -1659,
  -1185, -1530, -1278,   794, -1510,  -854,  -870,   478,
   -108,  -308,   996,   991,   958, -1460,  1522,  1628
};

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              centered representative congruent to a mod q in {-(q-1)/2,...,(q-1)/2}
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {-(q-1)/2,...,(q-1)/2} congruent to a modulo q.
**************************************************/
int16_t barrett_reduce(int16_t a) {
  int16_t t;
  const int16_t v = ((1<<26) + KYBER_Q/2)/KYBER_Q;

  t  = ((int32_t)v*a + (1<<25)) >> 26;
  t *= KYBER_Q;
  return a - t;
}


/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              16-bit integer congruent to a * R^-1 mod q, where R=2^16
*
* Arguments:   - int32_t a: input integer to be reduced;
*                           has to be in {-q2^15,...,q2^15-1}
*
* Returns:     integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
**************************************************/
int32_t montgomery_reduce_64(int64_t a)
{
  int32_t t;

  t = (int64_t)(int32_t)a*QINV;
  t = (a - (int64_t)t*KYBER_Q) >> 16;
  return t;
}

int16_t montgomery_reduce(int32_t a)
{
  int16_t t;

  t = (int16_t)a*QINV;
  t = (a - (int32_t)t*KYBER_Q) >> 16;
  return t;
}

static int16_t fqmul(int16_t a, int16_t b) {
  return montgomery_reduce((int32_t)a*b);
}

/*************************************************
* Name:        ntt
*
* Description: Inplace number-theoretic transform (NTT) in Rq.
*              input is in standard order, output is in bitreversed order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements of Zq
**************************************************/
void ntt(int16_t r[256]) {
  unsigned int len, start, j, k;
  int16_t t, zeta;

  k = 1;
  for(len = 128; len >= 2; len >>= 1) {
    for(start = 0; start < 256; start = j + len) {
      zeta = zetas[k++];
      for(j = start; j < start + len; j++) {
        t = fqmul(zeta, r[j + len]);
        r[j + len] = r[j] - t;
        r[j] = r[j] + t;
      }
    }
  }
}

void ntt_dit(int32_t r[256]) {
  unsigned int len, start, j, k;
  int32_t t, zeta;

  k = 1;
  for(len = 128; len >= 2; len >>= 1) {
    for(start = 0; start < 256; start = j + len) {
      zeta = zetas[k++];
      for(j = start; j < start + len; j++) {
        t = montgomery_reduce((int32_t)zeta* r[j + len]);
        r[j + len] = r[j] - t;
        r[j] = r[j] + t;
      }
    }
  }
}

/*************************************************
* Name:        invntt_tomont
*
* Description: Inplace inverse number-theoretic transform in Rq and
*              multiplication by Montgomery factor 2^16.
*              Input is in bitreversed order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements of Zq
**************************************************/
void invntt(int16_t r[256]) {
  unsigned int start, len, j, k;
  int16_t t, zeta;
  const int16_t f = 1441; // mont^2/128

  k = 127;
  for(len = 2; len <= 128; len <<= 1) {
    for(start = 0; start < 256; start = j + len) {
      zeta = zetas[k--];
      for(j = start; j < start + len; j++) {
        t = r[j];
        r[j] = barrett_reduce(t + r[j + len]);
        r[j + len] = r[j + len] - t;
        r[j + len] = fqmul(zeta, r[j + len]);
      }
    }
  }

  for(j = 0; j < 256; j++)
    r[j] = fqmul(r[j], f);
}

/*************************************************
* Name:        basemul
*
* Description: Multiplication of polynomials in Zq[X]/(X^2-zeta)
*              used for multiplication of elements in Rq in NTT domain
*
* Arguments:   - int16_t r[2]: pointer to the output polynomial
*              - const int16_t a[2]: pointer to the first factor
*              - const int16_t b[2]: pointer to the second factor
*              - int16_t zeta: integer defining the reduction polynomial
**************************************************/
void basemul(int16_t r[2], const int16_t a[2], const int16_t b[2], int16_t zeta)
{
  r[0]  = fqmul(a[1], b[1]);
  r[0]  = fqmul(r[0], zeta);
  r[0] += fqmul(a[0], b[0]);
  r[1]  = fqmul(a[0], b[1]);
  r[1] += fqmul(a[1], b[0]);
}

void invntt_tomont(int32_t r[N]) {
  unsigned int start, len, j, k;
  int32_t t, zeta;
  const int32_t f = 1441; // mont^2/128

  k = 128;
  for(len = 2; len <= 128; len <<= 1) {
    for(start = 0; start < 256; start = j + len) {
      zeta = zetas[--k];
      for(j = start; j < start + len; j++) {
        t = r[j];
        r[j] = barrett_reduce(t + r[j + len]);
        r[j + len] = r[j + len] - t;
        r[j + len] = montgomery_reduce((int32_t)zeta*r[j + len]);
      }
    }
  }

  for(j = 0; j < 256; j++)
    r[j] = montgomery_reduce((int32_t)f*r[j]);
}

int main()
{
  int i;
  int16_t a[N]={1290,-1477,-1468,1471,604,1454,-300,722,-1306,1086,-1171,-1105,-187,-1230,1064,-1482,1657,-640,-444,1452,-40,1642,738,438,1164,159,470,1058,1565,1351,265,-716,-276,-1330,1576,1074,1162,1473,-840,-889,-486,1291,1141,-1590,947,-1201,-461,86,-823,-1094,-565,765,-1053,9,-855,1534,767,565,-1255,-112,-372,560,-776,1487,1288,180,-282,-653,-298,-1384,1544,925,1311,391,630,-245,1402,-1555,549,1319,840,886,-683,1018,-547,482,-604,-1100,582,525,1471,1065,-1162,13,-1504,-1302,440,1163,785,-47,1637,1206,223,627,-1328,-27,833,-339,1340,45,-1555,-105,175,18,-489,-808,-488,-311,1539,1153,1110,756,252,-860,-1455,-1384,902,1259,1193,-1652,-184,293,-421,955,-1316,1184,273,-400,-625,-1488,1015,593,-749,-500,-1525,-1087,-831,-666,73,-245,1530,-1414,-241,-1605,20,-466,-1564,-1552,1227,-702,-1054,-126,1163,-1159,724,947,1095,1271,256,-1021,248,-1520,-733,1618,975,1236,-11,-884,-1244,-1309,400,1414,728,1060,-621,1093,-1616,1448,596,-235,1439,-361,1345,-1448,-1243,533,-1343,1449,1293,1186,-359,1236,-1205,-1415,-433,509,139,165,-1202,726,417,-418,1609,-222,-1030,-1393,-1191,-726,-1142,1388,-1223,1050,-1571,1590,704,107,570,-302,-1230,634,881,-1241,-1143,-908,488,1344,-1643,124,-790,-1142,-1514,-757,836,861,-863,1004,1406,-1015,-500,-1001,1139,-1308,59,-1348,-1417,-1491};
  int32_t a_1[N]={1290,-1477,-1468,1471,604,1454,-300,722,-1306,1086,-1171,-1105,-187,-1230,1064,-1482,1657,-640,-444,1452,-40,1642,738,438,1164,159,470,1058,1565,1351,265,-716,-276,-1330,1576,1074,1162,1473,-840,-889,-486,1291,1141,-1590,947,-1201,-461,86,-823,-1094,-565,765,-1053,9,-855,1534,767,565,-1255,-112,-372,560,-776,1487,1288,180,-282,-653,-298,-1384,1544,925,1311,391,630,-245,1402,-1555,549,1319,840,886,-683,1018,-547,482,-604,-1100,582,525,1471,1065,-1162,13,-1504,-1302,440,1163,785,-47,1637,1206,223,627,-1328,-27,833,-339,1340,45,-1555,-105,175,18,-489,-808,-488,-311,1539,1153,1110,756,252,-860,-1455,-1384,902,1259,1193,-1652,-184,293,-421,955,-1316,1184,273,-400,-625,-1488,1015,593,-749,-500,-1525,-1087,-831,-666,73,-245,1530,-1414,-241,-1605,20,-466,-1564,-1552,1227,-702,-1054,-126,1163,-1159,724,947,1095,1271,256,-1021,248,-1520,-733,1618,975,1236,-11,-884,-1244,-1309,400,1414,728,1060,-621,1093,-1616,1448,596,-235,1439,-361,1345,-1448,-1243,533,-1343,1449,1293,1186,-359,1236,-1205,-1415,-433,509,139,165,-1202,726,417,-418,1609,-222,-1030,-1393,-1191,-726,-1142,1388,-1223,1050,-1571,1590,704,107,570,-302,-1230,634,881,-1241,-1143,-908,488,1344,-1643,124,-790,-1142,-1514,-757,836,861,-863,1004,1406,-1015,-500,-1001,1139,-1308,59,-1348,-1417,-1491};
  ntt(a);
  ntt_dit(a_1);
  invntt(a);
  invntt_tomont(a_1);
  printf("INTT: ");
  for(i=0;i<N;i++)
  {
    printf("%d ",a[i]);
  }
  printf("\n\n");
  printf("INTT: ");
  for(i=0;i<N;i++)
  {
    printf("%d ",a_1[i]);
  }
  printf("\n\n");
}