/***************************************************************
 * BIT UTILITY for C and FORTRAN    
 * 2003.10.24 recieve from SATO-san, SATELLITE-GROUP, NWP, JMA.
 ***************************************************************/

#include	<stdio.h>
#include	<stdlib.h>
#include	<memory.h>
/***************************************************************
* get bit template 
***************************************************************/
unsigned char
get_bit_tmpl(int ib){
	static unsigned char 
		bit_tmpl[]={
			0x80,0x40,0x20,0x10,
			0x08,0x04,0x02,0x01
		};
	if(ib>=0&&ib<8) return bit_tmpl[ib];
	else            return 0xff;
}

/***************************************************************
* get bit mask template
***************************************************************/
unsigned char
get_msk_tmpl(int ib){
	static unsigned char 
		msk_tmpl[]={
			0x7f,0xbf,0xdf,0xef,
			0xf7,0xfb,0xfd,0xfe
		};
	if(ib>=0&&ib<8) return msk_tmpl[ib];
	else            return 0xff;
}

/***************************************************************
* read bit data
***************************************************************/
int
read_bit(void *pc, int ib){
	unsigned char	tmpl,*pcpt;
	int				ibyt,ibit;

	ibyt=ib/8; ibit=ib%8;
	tmpl=get_bit_tmpl(ibit);
	pcpt=(unsigned char *)pc+ibyt;
	if(((*pcpt)&tmpl)==tmpl) return 1;
	else                     return 0;
}

/***************************************************************
* set 1 in bit data
***************************************************************/
void
set_bit1(void *pc, int ib){
	unsigned char	tmpl,*pcpt;
	int				ibyt,ibit;

	ibyt=ib/8; ibit=ib%8;
	tmpl=get_bit_tmpl(ibit);
	pcpt=(unsigned char *)pc+ibyt;
	(*pcpt)|=tmpl;
}

/***************************************************************
* set 0 in bit data
***************************************************************/
void
set_bit0(void *pc, int ib){
	unsigned char	tmpl,*pcpt;
	int				ibyt,ibit;

	ibyt=ib/8; ibit=ib%8;
	tmpl=get_msk_tmpl(ibit);
	pcpt=(unsigned char *)pc+ibyt;
	(*pcpt)&=tmpl;
}

/***************************************************************
* MOVEB for C
***************************************************************/
void
moveb_c(void* dst, int idst, void *src, int isrc, int leng){
	int	ic;
	for(ic=0;ic<leng;ic++){
		if(read_bit(src,isrc+ic)==1) set_bit1(dst,idst+ic);
		else                         set_bit0(dst,idst+ic);
	}
}

/***************************************************************
* MOVEB for Fortran
***************************************************************/
#ifndef HAVE_MOVEB
void
#  if   defined C_FORTRAN_LOWER_
moveb_
#  elif defined C_FORTRAN_LOWER
moveb
#  elif defined C_FORTRAN_LOWER__
moveb_
#  else
MOVEB
#  endif
(void* dst, int *idst, void *src, int *isrc, int *leng){
	moveb_c(dst,*idst-1,src,*isrc-1,*leng);
}
#endif

/***************************************************************
* MOVEC for C
***************************************************************/
void
movec_c(void* dst, int idst, void *src, int isrc, int leng){
	if(leng>0){
		unsigned char	*cdst,*csrc;
		cdst=(unsigned char *)dst+idst;
		csrc=(unsigned char *)src+isrc;
		memcpy(cdst,csrc,leng);
	}
}

/***************************************************************
* MOVEC for Fortran
***************************************************************/
#ifndef HAVE_MOVEC
void
#  if   defined C_FORTRAN_LOWER_
movec_
#  elif defined C_FORTRAN_LOWER
movec
#  elif defined C_FORTRAN_LOWER__
movec_
#  else
/* MOVEC */
movec_
#  endif
(void* dst, int *idst, void *src, int *isrc, int *leng)
{
	movec_c(dst, *idst-1, src, *isrc-1, *leng);
}
#endif

/***************************************************************
* show bit data for test
***************************************************************/
void show_bit(void *pc,int leng){
	int ic;
	for(ic=0;ic<leng;ic++){
		printf("%i",read_bit(pc,ic));
		if(ic%4==3) printf(" ");
	}
	printf("\n");
}
/***************************************************************
* main routien foe test
***************************************************************/
/*
void
main(void){
	int	src,dst,ic;
	src=65535;
	dst=0;
	show_bit(&src,16);
	show_bit(&dst,16);
	moveb((void *)&dst,4,(void *)&src,4,8);
	show_bit(&dst,16);
}
*/
