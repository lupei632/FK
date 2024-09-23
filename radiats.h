/****************************************************************
   Head file for radiats.c

   rad[4][3] is the horizontal radiation rad[DD,DS,SS,EX][Z,R,T]
*****************************************************************/

#ifndef __RADIATS_HEAD__
  #define __RADIATS_HEAD__

#include<stdio.h>
#include<math.h>

#define DEG2RAD  1.745329252e-2  /*degree to rad*/

// mapping the 3x3 stress/strain tensor to 1x6 array (11, 22, 33, 23, 13, 12)
static int indx_i[6] = {0,1,2,1,0,0};
static int indx_j[6] = {0,1,2,2,2,1};
static int indxij[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};

void    dc_radiat(float stk,float dip,float rak,float rad[4][3]);
void	sf_radiat(float, float, float rad[4][3]);
void	mt_radiat(float, float mt[3][3], float rad[4][3]);
float	nmtensor(float iso, float clvd, float stk, float dip, float rak, float stiff[6][6], float mt[3][3]);
void	potency2moment(float stiff[6][6], float ptt[6], float mtt[6]);
void	ThomsenTTI(float tti[7], float stiff[6][6]);
void	RotateStiff(float alpha, float phi, float stiff[6][6]);
void	rotationMatrix(float alpha, float phi, float rot[3][3]);
float	radpmt(float mt[3][3], float alpha, float az, int type);
void	fk2mtg(float az, int npt, float *fk[12], float *mt[18]);

#endif
