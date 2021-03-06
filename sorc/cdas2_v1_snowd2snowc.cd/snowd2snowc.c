/*
 * change snow depth to snow covert
 * 11/2019 W. Ebisuzaki
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>

#include "gribwlib.h"
#include "pds4.h"
#include "ncepopn.h"

int main(int argc, char **argv) {

    unsigned char *pds, *gds;
    FILE *input, *output;
    float *data;
    int ndata, i, count, cnvt, nblank;

    long int pos, len_grib;

    /* preliminaries .. open up all files */

    if (argc != 3) {
	fprintf(stderr, "%s [in SNOWD gribfile] [out SNOWC]\n", argv[0]);
	exit(8);
    }
    if ((input = fopen(argv[1],"rb")) == NULL) {
        fprintf(stderr,"could not open file: %s\n", argv[1]);
        exit(7);
    }
    if ((output = fopen(argv[2],"wb+")) == NULL) {
        fprintf(stderr,"could not open file: %s\n", argv[2]);
        exit(7);
    }

    for(nblank = count = cnvt = pos = 0;;) {

	len_grib = rd_grib_rec(input, pos, &pds, &gds, &data, &ndata);
	if (len_grib <= 0) break;
	count++;

	if (PDS_PARAM(pds) == SNOD) {
	    PDS_PARAM(pds) = SNOWC;
	    cnvt++;
	    for (i = 0; i < ndata; i++) {
		if (!UNDEFINED_VAL(data[i])) {
		    data[i] = (data[i] > 0) ? 1.0 : 0.0;
		}
	    }
	    pds = PDStool(pds, P_dec_scale(0), P_end);
	    wrt_grib_rec(pds, gds, data, ndata, output);
	}
	else {
	    wrt_grib_rec(pds, gds, data, ndata, output);
	}
	pos += len_grib;
    }
    fclose(output);
    fclose(input);
    printf("%d SNOWC records converted out of %d, %d blank\n", cnvt, count,
         nblank);
    return 0;
}
