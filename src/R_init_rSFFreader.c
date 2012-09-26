#include "rSFFreader.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

#define REGISTER_CCALLABLE(fun) \
        R_RegisterCCallable("rSFFreader.h", #fun, (DL_FUNC) &fun)


static const R_CallMethodDef callMethods[] = {
/* RocheSFF-io.c */
    CALLMETHOD_DEF(read_sff, 5),
		CALLMETHOD_DEF(read_sff_header, 2),
		CALLMETHOD_DEF(sff_geometry, 1),
		CALLMETHOD_DEF(write_phred_quality, 1),
		{NULL, NULL, 0}
};



