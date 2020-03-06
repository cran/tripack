
#include <R.h>
#include <Rinternals.h>

#include "tripack.h"

#include <R_ext/Rdynload.h>

/* Fortran interface descriptions: */

static R_NativePrimitiveArgType addcst_t[11] = {
  INTSXP,  /* NCC */
  INTSXP,  /* LCC */
  INTSXP,  /* N */
  REALSXP, /* X */
  REALSXP, /* Y */
  INTSXP,  /* LWK */
  INTSXP,  /* IWK */
  INTSXP,  /* LIST */
  INTSXP,  /* LPTR */
  INTSXP,  /* LEND */
  INTSXP   /* IER */
};

//not yet finished 
//static R_NativePrimitiveArgType awvoronoi_t[16] = {
//  INTSXP,  /* NCC */
//  INTSXP,  /* LCC */
//  INTSXP,  /* N */
//  REALSXP, /* X */
//  REALSXP, /* Y */
//  REALSXP, /* WT */
//  INTSXP,  /* LIST */
//  INTSXP,  /* LPTR */
//  INTSXP,  /* LEND */
//  INTSXP,  /* NT */
//  REALSXP, /* LCCC */
//  INTSXP,  /* ICCC */
//  INTSXP,  /* LCT */
//  INTSXP,  /* LTRI */
//  INTSXP,  /* MAXIT */
//  INTSXP   /* IER */
//};

static R_NativePrimitiveArgType troutq_t[17] = {
  INTSXP,  /* NCC */
  INTSXP,  /* LCC */
  INTSXP,  /* N */
  REALSXP, /* X */
  REALSXP, /* Y */
  INTSXP,  /* LIST */
  INTSXP,  /* LPTR */
  INTSXP,  /* LEND */
  INTSXP,  /* LOUT */
  INTSXP,  /* NNABS */
  INTSXP,  /* NPTR */
  INTSXP,  /* NPTR1 */
  INTSXP,  /* NABOR */
  INTSXP,  /* NBNOS */
  INTSXP,  /* NA */
  INTSXP,  /* NB */
  INTSXP   /* NT */
};

static R_NativePrimitiveArgType bnodes_t[8] = {
  INTSXP, /* N */
  INTSXP, /* LIST */
  INTSXP, /* LPTR */
  INTSXP, /* LEND */
  INTSXP, /* NODES */
  INTSXP, /* NB */
  INTSXP, /* NA */
  INTSXP  /* NT */
};

static R_NativePrimitiveArgType inhull_t[10] = {
  REALSXP, /* XP */
  REALSXP, /* YP */
  INTSXP,  /* NP */
  REALSXP, /* X */
  REALSXP, /* Y */
  INTSXP,  /* N */
  INTSXP,  /* LIST */
  INTSXP,  /* LPTR */
  INTSXP,  /* LEND */
  LGLSXP   /* INH */
};

static R_NativePrimitiveArgType onhull_t[11] = {
  REALSXP, /* XP */
  REALSXP, /* YP */
  INTSXP,  /* NP */
  REALSXP, /* X */
  REALSXP, /* Y */
  INTSXP,  /* N */
  INTSXP,  /* LIST */
  INTSXP,  /* LPTR */
  INTSXP,  /* LEND */
  LGLSXP,  /* ONH */
  REALSXP  /* EPS */
};

static R_NativePrimitiveArgType trlist_t[11] = {
  INTSXP, /* NCC */
  INTSXP, /* LCC */
  INTSXP, /* N */
  INTSXP, /* LIST */
  INTSXP, /* LPTR */
  INTSXP, /* LEND */
  INTSXP, /* NROW */
  INTSXP, /* NT */
  INTSXP, /* LTRI */
  INTSXP, /* LCT */
  INTSXP  /* IER */
};

static R_NativePrimitiveArgType trfind_t[12] = {
  INTSXP,  /* NST */
  REALSXP, /* PX */
  REALSXP, /* PY */
  INTSXP,  /* N */
  REALSXP, /* X */
  REALSXP, /* Y */
  INTSXP,  /* LIST */
  INTSXP,  /* LPTR */
  INTSXP,  /* LEND */
  INTSXP,  /* I1 */
  INTSXP,  /* I2 */
  INTSXP   /* I3 */
};

static R_NativePrimitiveArgType trmesh_t[11] = {
  INTSXP,  /* N */
  REALSXP, /* X */
  REALSXP, /* Y */
  INTSXP,  /* LIST */
  INTSXP,  /* LPTR */
  INTSXP,  /* LEND */
  INTSXP,  /* LNEW */
  INTSXP,  /* NEAR */
  INTSXP,  /* NEXT */
  REALSXP, /* DIST */
  INTSXP   /* IER */
};

static R_NativePrimitiveArgType voronoi_t[14] = {
  INTSXP,  /* NCC */
  INTSXP,  /* LCC */
  INTSXP,  /* N */
  REALSXP, /* X */
  REALSXP, /* Y */
  INTSXP,  /* LIST */
  INTSXP,  /* LPTR */
  INTSXP,  /* LEND */
  INTSXP,  /* NT */
  REALSXP, /* LCCC */
  INTSXP,  /* ICCC */
  INTSXP,  /* LCT */
  INTSXP,  /* LTRI */
  INTSXP   /* IER */
};

static R_FortranMethodDef fortranMethods[] = {
  {"addcst", (DL_FUNC) &F77_SUB(addcst), 11, addcst_t}, 
    /* add.constraint */
// not yet finished
//{"awvoronoi", (DL_FUNC) &F77_SUB(awvoronoi), 16, awvoronoi_t}, 
//  /* aw.voronoi.mosaic */
  {"troutq", (DL_FUNC) &F77_SUB(troutq), 17, troutq_t}, 
    /* cells, neighbours, plot.tri, print.tri, summary.tri */
  {"bnodes", (DL_FUNC) &F77_SUB(bnodes), 8, bnodes_t}, 
    /* convex.hull */
  {"inhull", (DL_FUNC) &F77_SUB(inhull), 10, inhull_t}, 
    /* in.convex.hull */
  {"onhull", (DL_FUNC) &F77_SUB(onhull), 11, onhull_t}, 
    /* on.convex.hull */
  {"trlist", (DL_FUNC) &F77_SUB(trlist), 11, trlist_t}, 
    /* triangles */
  {"trfind", (DL_FUNC) &F77_SUB(trfind), 12, trfind_t}, 
    /* tri.find */
  {"trmesh", (DL_FUNC) &F77_SUB(trmesh), 11, trmesh_t}, 
    /* tri.mesh */
  {"voronoi", (DL_FUNC) &F77_SUB(voronoi), 14, voronoi_t},
    /* voronoi.mosaic */
  {NULL, NULL, 0}
};

void
R_init_tripack(DllInfo *info)
{
  R_registerRoutines(info, 
		     NULL /*cMethods*/, NULL /*callMethods*/, 
		     fortranMethods, NULL/*externalMethods*/);
  R_useDynamicSymbols(info, TRUE);
  /* does not work with Fortran code:
   * R_useDynamicSymbols(info, FALSE);
   * R_forceSymbols(info, TRUE);
   */
}
