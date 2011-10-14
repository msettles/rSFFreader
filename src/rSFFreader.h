#ifndef _RSFFREADER_H_
#define _RSFFREADER_H_

#ifdef __cplusplus
extern "C" {
#endif

// R includes
//#include <R.h>
#include <Rdefines.h>
//#include <Rinternals.h> // Rprintf, SEXP
#include <R_ext/Rdynload.h>

// System includes
#include <zlib.h>
#include <stdint.h>		// uint64_t, uint32_t, uint16_t


SEXP read_sff(
		SEXP files,
		SEXP use_names,
		SEXP lkup_seq,
		SEXP lkup_qual,
		SEXP verbose
);

SEXP read_sff_header(
		SEXP files,
		SEXP verbose
);

SEXP sff_geometry(
		SEXP files
);

SEXP write_phred_quality(
		SEXP id,
		SEXP quality, 
		SEXP fname,
		SEXP fmode,
		SEXP max_width
);

#ifdef __cplusplus
}
#endif

#endif /* _RSFFHEADER_H_ */

