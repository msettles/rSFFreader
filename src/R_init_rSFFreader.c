#include "rSFFreader.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

#define REGISTER_CCALLABLE(fun) \
        R_RegisterCCallable("Roche454Reads.h", #fun, (DL_FUNC) &fun)


static const R_CallMethodDef callMethods[] = {
/* RocheSFF-io.c */
        CALLMETHOD_DEF(read_roche_sff, 5),
		CALLMETHOD_DEF(read_roche_sff_header, 2),
		CALLMETHOD_DEF(sff_geometry, 1),

		{NULL, NULL, 0}
};



