************************************************************************
*     -----------------
*     INCLUDE "guref.h"
*     -----------------
*
*     PMSREF ; ratio of Sei hansha on PMT
*     PMRREF ; ratio of Ran hansha on PMT
*              (not necessary, removed by Byeongsu Yang, 2008.08.24)
*     BSTREF ; ratio of reflection on Black Sheet
*     COREPMT; PMT correction efficiency
*-- Algorithm of reducing the number of photons is improved
*   and can be applied also LE mode. Byeongsu Yang 2008/04/23
*     maxfct ; Maximum of QE table values
*     REDUCFCT; photon-reducing factor
*
*     (Creation Date and Author)
*       1992.08.27 ; First version by Y.Koshio
*
************************************************************************
      REAL PMSREF, BSTREF, COREPMT, maxfct, REDUCFCT
      CHARACTER*3 cqetab
      COMMON /REFRTO/ PMSREF,BSTREF,COREPMT,maxfct,REDUCFCT,cqetab
