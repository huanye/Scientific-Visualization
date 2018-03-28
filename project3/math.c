#include "rnd.h"
#include "rndPrivate.h"

/*
** rnd_m4_affine_inv: invert a given 4x4 matrix "mat", and store inverse in
** "inv", assuming that the bottom row of "mat" is [0 0 0 1].  This means the
** matrix represents an affine transform, and its inverse is simpler to
** compute than for a general 4x4 matrix.
**
** To do this inverse, you may use M34_UPPER, M3_INVERSE, V3_SET,
** MV3_MUL, M4_SET, and your awareness of how 3x3 and 4x4 matrices
** are represented as a 1D array, as described in rndMath.h.
**
** You may *not* use ell_4m_inv_d or any other function or macro that
** accomplishes 4x4 matrix inversion with a single line of code
*/
void
rnd_m4_affine_inv(real inv[16], const real mat[16]) {
  // v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v.v  begin student code
real left_upper[9];
M34_UPPER(left_upper,mat);
real lu_inv[9];
real tmp;
M3_INVERSE(lu_inv, left_upper, tmp);
real tmp_col[3];
real V3_col[3];
V3_SET(tmp_col,mat[3],mat[7],mat[11]);
MV3_MUL(V3_col,lu_inv,tmp_col);
M4_SET(inv,lu_inv[0],lu_inv[1],lu_inv[2],-V3_col[0],
           lu_inv[3],lu_inv[4],lu_inv[5],-V3_col[1],
           lu_inv[6],lu_inv[7],lu_inv[8],-V3_col[2],
           0,0,0,1);
  // ( 10 lines in reference implementation )
  // ^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^'^  end student code
  return;
}

void
rnd_m4_print(FILE *file, const real mm[16]) {
  uint nlen = rndSprintLenMax(0, mm, 16);
  char buff[256];
  sprintf(buff, "%% %u." RCD "g  %% %u." RCD "g  %% %u." RCD "g  %% %u." RCD "g\n",
          nlen, nlen, nlen, nlen);
  fprintf(file, buff, mm[ 0], mm[ 1], mm[ 2], mm[ 3]);
  fprintf(file, buff, mm[ 4], mm[ 5], mm[ 6], mm[ 7]);
  fprintf(file, buff, mm[ 8], mm[ 9], mm[10], mm[11]);
  fprintf(file, buff, mm[12], mm[13], mm[14], mm[15]);
  return;
}
