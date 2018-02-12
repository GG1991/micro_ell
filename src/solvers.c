#include "solvers.h"


#ifdef PETSC
int solvers_print_petsc_ksp_info(MPI_Comm COMM, KSP ksp) {

  int kspits, reason;
  double kspnorm;
  char *reason_s;

  KSPGetIterationNumber(ksp, &kspits);
  KSPGetConvergedReason(ksp, &reason);
  KSPGetResidualNorm(ksp, &kspnorm);
  switch(reason) {

    case KSP_CONVERGED_RTOL:
      reason_s = strdup("RTOL");
      break;
    case KSP_CONVERGED_ATOL:
      reason_s = strdup("ATOL");
      break;
    case KSP_DIVERGED_ITS:
      reason_s = strdup("DIV ITS");
      break;
    default :
      reason_s = strdup("DIV ITS");
      sprintf(reason_s, "%d", reason);
      break;
  }
  PetscPrintf(COMM,"kspits %D kspnorm %e kspreason %s", kspits, kspnorm, reason_s);

  return 0;
}
#endif
