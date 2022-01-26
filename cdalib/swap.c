/***************************************************************/
/*	swap_ys.c:                                                 */
/*		00.11.21:                                              */
/*			ByteSwap                                           */
/*		00.12.15:                                              */
/*			Convert IBM->IEEE32                                */
/*	swap_c2(void *pt);		swap 2 byte data                   */
/*	swap_c4(void *pt);		swap 4 byte data                   */
/*	swap_c8(void *pt);		swap 8 byte data                   */
/***************************************************************/

#include  <stdio.h>
#include  <stdlib.h>
#include  <memory.h>
#include  <math.h>

/*Endian Big->Little <2byte>*/
void swap_c2(void *pt)
{
	char cbit[2],cswp;
	memcpy(cbit,pt,2);
	cswp =cbit[0];cbit[0]=cbit[1];cbit[1]=cswp;
	memcpy(pt,cbit,2);
}

/*Endian Big->Little <4byte>*/
void swap_c4(void *pt)
{
	char cbit[4],cswp;
	memcpy(cbit,pt,4);
	cswp =cbit[0];cbit[0]=cbit[3];cbit[3]=cswp;
	cswp =cbit[1];cbit[1]=cbit[2];cbit[2]=cswp;
	memcpy(pt,cbit,4);
}

/*Endian Big->Little <8byte>*/
void swap_c8(void *pt)
{
	char cbit[8],cswp;
	memcpy(cbit,pt,8);
	cswp =cbit[0];cbit[0]=cbit[7];cbit[7]=cswp;
	cswp =cbit[1];cbit[1]=cbit[6];cbit[6]=cswp;
	cswp =cbit[2];cbit[2]=cbit[5];cbit[5]=cswp;
	cswp =cbit[3];cbit[3]=cbit[4];cbit[4]=cswp;
	memcpy(pt,cbit,8);
}

void
swp2swp4(void *pt){
	char	*pc;
	char	c4[4];
	pc=(char *)pt;
	memcpy(&c4[0],&pc[2],2);
	memcpy(&c4[2],&pc[0],2);
	memcpy(pc,c4,4);
}


/********
 * FORTRAN INTERFACE
 ********/
void 
#if   defined C_FORTRAN_LOWER_
swap_c2_l_
#elif defined C_FORTRAN_LOWER
swap_c2_l
#elif defined C_FORTRAN_LOWER__
swap_c2_l__
#else
SWAP_C2_L
#endif
(void *pt)
{
 swap_c2(pt);
}

void 
#if   defined C_FORTRAN_LOWER_
swap_c4_l_
#elif defined C_FORTRAN_LOWER
swap_c4_l
#elif defined C_FORTRAN_LOWER__
swap_c4_l__
#else
SWAP_C4_L
#endif
(void *pt)
{
 swap_c4(pt);
}

void 
#if   defined C_FORTRAN_LOWER_
swap_c8_l_
#elif defined C_FORTRAN_LOWER
swap_c8_l
#elif defined C_FORTRAN_LOWER__
swap_c8_l__
#else
SWAP_C8_L
#endif
(void *pt)
{
 swap_c8(pt);
}

void
#if   defined C_FORTRAN_LOWER_
swp2swp4_l_
#elif defined C_FORTRAN_LOWER
swp2swp4_l
#elif defined C_FORTRAN_LOWER__
swp2swp4_l__
#else
SWP2SWP4_L
#endif
(void *pt)
{
  swp2swp4(pt);
}
