/* Author: Matt Settles
   with original contributtion from:Kyu-Chul Cho (cho@vandals.uidaho.edu)
 * Last Modified: June 9, 2010
 */

#ifdef __cplusplus
extern "C" {
#endif

// Base rSFFreader includes
#include "rSFFreader.h"

// System Includes
#include <stdio.h>
#include <string.h>
#include <stdint.h>        // uint64_t, uint32_t, uint16_t
#include <stdlib.h>        // malloc(), free()
#include <arpa/inet.h>  // htons(), htonl()

//#include "encode.h"
#include "IRanges_interface.h"
#include "Biostrings_interface.h"


static int debug = 0;
static int flowgram_bytes_per_flow = 2;

/***************************** SFF file structures *****************************/
typedef struct CommonHeader {
    uint32_t magic_number;
    char version[4];
    uint64_t index_offset;
    uint32_t index_length;
    uint32_t number_of_reads;
    uint16_t commonHeader_length;
    uint16_t key_length;
    uint16_t number_of_flows_per_read;
    uint8_t flowgram_format_code;
    char *flow_chars;// flow_chars omitted
    char *key_sequence;// key_sequence omitted
    // eight_byte_padding is not necessary
} COMMONheader;

typedef struct ReadHeader {
    uint16_t read_header_length;
    uint16_t name_length;
    uint32_t number_of_bases;
    uint16_t clip_qual_left;
    uint16_t clip_qual_right;
    uint16_t clip_adapter_left;
    uint16_t clip_adapter_right;
    // name
    // eight_byte_padding
} READheader;

typedef struct ReadData {
    float* flows;
    int*   index;
    char*  name;
    char*  quality;
    char*  bases;
} READdata;

/***************************** Auxilary Functions *****************************/

int
isMagicNumberCorrect(uint32_t magic_number)
{
        if(0x2E736666==magic_number) return 1;
        else return 0;
}

// only version 1 is supported
int
isVersionInfoCorrect(char* version)
{
        if(version[0]==0 && version[1]==0 && version[2]==0 && version[3]==1)
                return 1;
        else return 0;
}


static void
debug_sffreader(int read, COMMONheader commonHeader,READheader header, READdata data)
{
    int i;
    Rprintf("read(%d) \n", read);
    Rprintf("read header length: %u\n", header.read_header_length);
    Rprintf("name length: %u\n", header.name_length);
    Rprintf("# of bases: %u\n", header.number_of_bases);
    Rprintf("clip qual left: %u\n", header.clip_qual_left);
    Rprintf("clip qual right: %u\n", header.clip_qual_right);
    Rprintf("clip adapter left: %u\n", header.clip_adapter_left);
    Rprintf("clip adapter right: %u\n", header.clip_adapter_right);
    Rprintf("Name:%s",data.name);
    Rprintf("\n");
    Rprintf("Flowgram: ");
    for(i=0; i<commonHeader.number_of_flows_per_read; i++) {
        Rprintf("%.2f\t", data.flows[i]);
    }
    Rprintf("\n");
    Rprintf("Flow Indexes: ");
    for(i=0; i<header.number_of_bases; i++) {
        Rprintf("%hu\t", data.index[i]);
    }
    Rprintf("\n");
    Rprintf("Bases:\t\t");
    for(i=0; i<header.number_of_bases; i++) {
        Rprintf("%c", data.bases[i]);
    }
    Rprintf("\n");
    Rprintf("Quality Scores: ");
    for(i=0; i<header.number_of_bases; i++) {
        Rprintf("%c", data.quality[i]);
    }
    Rprintf("\n");
}

static void
debug_headerreader(COMMONheader commonHeader)
{
    int i;
    Rprintf("magic number c:%d\n",
            isMagicNumberCorrect(commonHeader.magic_number));
    Rprintf("version: %c%c%c%c\n",
            commonHeader.version[0]+'0', commonHeader.version[1]+'0',
            commonHeader.version[2]+'0', commonHeader.version[3]+'0');
    Rprintf("correct version c:%d\n", isVersionInfoCorrect(commonHeader.version));
    Rprintf("index offset: %llu\n", commonHeader.index_offset);
    Rprintf("index length: %u\n", commonHeader.index_length);
    Rprintf("number of reads: %u\n", commonHeader.number_of_reads);
    Rprintf("commonHeader length: %u\n", commonHeader.commonHeader_length);
    Rprintf("key length: %u\n", commonHeader.key_length);
    Rprintf("number of flows per read: %u\n",
                commonHeader.number_of_flows_per_read);
    Rprintf("flowgram format code: %u\n",commonHeader.flowgram_format_code);
    Rprintf("flow chars:\n");
    for(i=0; i<commonHeader.number_of_flows_per_read; i++) {
        Rprintf("%c", commonHeader.flow_chars[i]);
    }
    Rprintf("\n");
    Rprintf("key:");
    for(i=0; i<commonHeader.key_length; i++) {
        Rprintf("%c", commonHeader.key_sequence[i]);
    }
    Rprintf("\n");
}

int
isBigEndian()
{
    uint64_t integerForTesting = UINT64_C(0x0123456789abcdef);
    unsigned char *ch = (unsigned char *) &integerForTesting;
    int returnValue;
    if (*ch == 0xef) // little endian
    {
        returnValue = 0;
    }
    else // big endian
    {
        returnValue = 1;
    }
    return returnValue;
}

/* convert 64 bit Big Endian integer to Native Endian(means current machine) */
// As far as I know, there is no standard function for 64 bit conversion
uint64_t
BE64toNA(uint64_t bigEndian)
{
    // if current machine uses Big Endian, then return the input value
    if(isBigEndian()) return bigEndian;
    else {    // else the machine must use little Endian. Convert to little Endian
        uint64_t littleEndian =
        ((bigEndian &   (0x00000000000000FF)) << 56) |
        ((bigEndian &   (0x000000000000FF00)) << 40) |
        ((bigEndian &   (0x0000000000FF0000)) << 24) |
        ((bigEndian &   (0x00000000FF000000)) << 8)  |
        ((bigEndian &   (0x000000FF00000000)) >> 8)  |
        ((bigEndian &   (0x0000FF0000000000)) >> 24) |
        ((bigEndian &   (0x00FF000000000000)) >> 40) |
        ((bigEndian &   (0xFF00000000000000)) >> 56);
        return littleEndian;
    }
}

int
commonHeaderPaddingSize(int number_of_flows_per_read, int key_length)
{
    int remaining = (7 + number_of_flows_per_read + key_length)%8;
    int padding_size;
    if(remaining==0) padding_size = 0;
    else padding_size = 8-remaining;
    return padding_size;
}

int
readHeaderPaddingSize(int read_header_len, int name_length)
{
    return read_header_len - name_length - 16;
}

int
readDataPaddingSize(int number_of_flows_per_read, int number_of_bases)
{
//     Rprintf("%i", (number_of_flows_per_read*flowgram_bytes_per_flow + (number_of_bases*3)));
    int remaining = ((number_of_flows_per_read*flowgram_bytes_per_flow) + (number_of_bases*3))%8;
    int padding_size;
    if(remaining==0) padding_size = 0;
    else padding_size = 8-remaining;
    return padding_size;
}

void
checkPadding (uint8_t *byte_padding, int padding_size)
{
    int i;
    for(i=0; i<padding_size; i++)
    {
        if(byte_padding[i] !=0 )
            Rf_error("The header does not have proper byte_padding");
    }
}

char
phredToChar(uint8_t score)
{
    return (char)((int)score + 33);
}

int
CharToPhred(char score)
{
    return ((int)score - 33);
}

/**************************** END Auxilary Functions **************************/

/**************************** SFF Header Functions ****************************/

/*
 * side effect: this function uses malloc.
 * ==> MUST call freeCommonHeader() function after calling readCommonHeader
 */
COMMONheader
readCommonHeader(const char *fname){
      int i, padding_size, fres;
    char ch;
    uint16_t uint16;
    uint8_t uint8;
    COMMONheader commonHeader;

    FILE* file = fopen (fname, "r");
    if(file==NULL)
        Rf_error("cannot open specified file: '%s'", fname);

    /** Read Common Header Section **/
    int size = fread( &commonHeader.magic_number, sizeof(uint32_t), 1, file);
    commonHeader.magic_number = htonl(commonHeader.magic_number);
    fres = fread( commonHeader.version, sizeof(char), 4, file);

    fres = fread( &commonHeader.index_offset, sizeof(uint64_t), 1, file);
    commonHeader.index_offset = BE64toNA(commonHeader.index_offset);

    fres = fread( &commonHeader.index_length, sizeof(uint32_t), 1, file);
    commonHeader.index_length = htonl(commonHeader.index_length);

    fres = fread( &commonHeader.number_of_reads, sizeof(uint32_t), 1, file);
    commonHeader.number_of_reads = htonl(commonHeader.number_of_reads);

    fres = fread( &commonHeader.commonHeader_length, sizeof(uint16_t), 1, file);
    commonHeader.commonHeader_length = htons(commonHeader.commonHeader_length);

    fres = fread( &commonHeader.key_length, sizeof(uint16_t), 1, file);
    commonHeader.key_length = htons(commonHeader.key_length);

    fres = fread( &commonHeader.number_of_flows_per_read, sizeof(uint16_t), 1, file);
    commonHeader.number_of_flows_per_read = htons(commonHeader.number_of_flows_per_read);

    fres - fread( &commonHeader.flowgram_format_code, sizeof(uint8_t),1, file);
    // check flowgram_format_code
    if (commonHeader.flowgram_format_code != 1)
        Rf_error("Unknown flowgram format code value:%u",commonHeader.flowgram_format_code);
    // flow chars
    // malloc must be freed before the function ends
    commonHeader.flow_chars =
        (char*) malloc(sizeof(char)*(commonHeader.number_of_flows_per_read+1));
    if (commonHeader.flow_chars==NULL) Rf_error("cannot allocate memory");

    for(i=0; i<commonHeader.number_of_flows_per_read; i++) {
        fres = fread(&ch, sizeof(char), 1, file);
        commonHeader.flow_chars[i] = ch;
    }
    commonHeader.flow_chars[i] = '\0';

    // key sequence
    // malloc must be freed before the function ends
    commonHeader.key_sequence =
        (char*) malloc(sizeof(char)*(commonHeader.key_length+1));
    if (commonHeader.key_sequence==NULL) Rf_error("cannot allocate memory");

    for(i=0; i<commonHeader.key_length; i++) {
        fres = fread(&ch, sizeof(char), 1, file);
        commonHeader.key_sequence[i] = ch;
    }
    commonHeader.key_sequence[i] = '\0';

    fclose(file);

    // if debug, print out header info
    if (debug == 1){
        debug_headerreader(commonHeader);
    }
    return commonHeader;
}


/* MUST called after using readCommonHeader(...) function */
void freeCommonHeader(COMMONheader commonHeader) {
    free(commonHeader.flow_chars);
    free(commonHeader.key_sequence);
}

void freeReadData(READdata readData) {
    free(readData.flows);
    free(readData.index);
    free(readData.name);
    free(readData.quality);
    free(readData.bases);
}

// Reads a single SFF file, converts the output to a R vector and returns the vector
SEXP
readOneSFFHeader (SEXP file)
{
    // we want to call readCommonHeader and convert to SEXP
    COMMONheader commonHeader = readCommonHeader(CHAR(file));
    SEXP headerList = allocVector(VECSXP, 12);
    SEXP eltnm;
    static const char *eltnames[] = {
        "filename", "magic_number", "version", "index_offset",
        "index_length", "number_of_reads", "header_length",
        "key_length", "number_of_flows_per_read", "flowgram_format_code", "flow_chars",
        "key_sequence"
    };
    PROTECT(headerList);
    SET_VECTOR_ELT(headerList, 0, mkString(CHAR(file)));
    SET_VECTOR_ELT(headerList, 1, Rf_ScalarInteger(commonHeader.magic_number));
    SET_VECTOR_ELT(headerList, 2, mkString(commonHeader.version));
    SET_VECTOR_ELT(headerList, 3, Rf_ScalarInteger(commonHeader.index_offset));
    SET_VECTOR_ELT(headerList, 4, Rf_ScalarInteger(commonHeader.index_length));
    SET_VECTOR_ELT(headerList, 5, Rf_ScalarInteger(commonHeader.number_of_reads));
    SET_VECTOR_ELT(headerList, 6, Rf_ScalarInteger(commonHeader.commonHeader_length));
    SET_VECTOR_ELT(headerList, 7, Rf_ScalarInteger(commonHeader.key_length));
    SET_VECTOR_ELT(headerList, 8,
                   Rf_ScalarInteger(commonHeader.number_of_flows_per_read));
    SET_VECTOR_ELT(headerList, 9, Rf_ScalarInteger(commonHeader.flowgram_format_code));
    SET_VECTOR_ELT(headerList, 10, mkString(commonHeader.flow_chars));
    SET_VECTOR_ELT(headerList, 11, mkString(commonHeader.key_sequence));

    PROTECT( eltnm = allocVector( STRSXP, 12 ) );
    for( int i = 0; i < 12; i++ )
        SET_STRING_ELT( eltnm, i, mkChar( eltnames[i] ) );
    namesgets( headerList, eltnm );
    UNPROTECT(1);
    freeCommonHeader(commonHeader);
    UNPROTECT(1);
    return headerList;
}

// Reads multiple SFF files, by calling readOneSFFHeader
SEXP
read_sff_header(SEXP files, SEXP verbose)
{
    int i, nfiles= 0,set_verbose;
    SEXP fname;
    if (!IS_CHARACTER(files))
        Rf_error("'%s' must be '%s'", "files", "character");

    set_verbose = LOGICAL(verbose)[0];

    nfiles = LENGTH(files);

    SEXP ans;
    PROTECT(ans = allocVector(VECSXP,nfiles));
    for (i = 0; i < nfiles; ++i) {
        R_CheckUserInterrupt();
        fname = STRING_ELT(files, i);
        if(set_verbose) {
            Rprintf("reading header for sff file:%s\n", CHAR(fname));
        }
        SET_VECTOR_ELT(ans,i,readOneSFFHeader(fname));
    }
    UNPROTECT(1);
    return ans;
}

int
count_reads_sum(SEXP files)
{
    COMMONheader header;
    int i, nfile = LENGTH(files);
    int nreads = 0;

    for(i=0; i<nfile; i++) {
        header = readCommonHeader(CHAR(STRING_ELT(files,i)));
        nreads += header.number_of_reads;
        freeCommonHeader(header);
    }
    return nreads;
}
/*************************** END SFF Header Functions *************************/

/****************************************************************************
 * Reading SFF files.
 */

typedef struct irange_values{
    int *width;
    int *start;
	int length;
} IRANGE_VALUES;

typedef struct sff_loader {
    void (*load_seqid)(struct sff_loader *loader,
        const cachedCharSeq *dataline);
    void (*load_seq)(struct sff_loader *loader,
        const cachedCharSeq *dataline);
    void (*load_qual)(struct sff_loader *loader,
        const cachedCharSeq *dataline);
    void (*load_qclip)(struct sff_loader *loader,
        int start, int width);
    void (*load_aclip)(struct sff_loader *loader,
        int start, int width);
    int nrec;
    void *ext;  /* loader extension (optional) */
} SFFloader;


typedef struct sff_loader_ext {
    CharAEAE ans_names_buf;
    cachedXVectorList cached_seq;
    cachedXVectorList cached_qual;
    IRANGE_VALUES cached_qual_clip;
    IRANGE_VALUES cached_adapt_clip;
    const int *lkup_seq;
    int lkup_length_seq;
    const int *lkup_qual;
    int lkup_length_qual;
} SFF_loaderExt;

static void SFF_load_seqid(SFFloader *loader,
        const cachedCharSeq *dataline)
{
    SFF_loaderExt *loader_ext;
    CharAEAE *ans_names_buf;

    loader_ext = loader->ext;
    ans_names_buf = &(loader_ext->ans_names_buf);
    // This works only because dataline->seq is nul-terminated!
    append_string_to_CharAEAE(ans_names_buf, dataline->seq);
    return;
}

static void SFF_load_seq(SFFloader *loader, const cachedCharSeq *dataline)
{
    SFF_loaderExt *loader_ext;
    cachedCharSeq cached_ans_elt;

    loader_ext = loader->ext;
    cached_ans_elt = get_cachedXRawList_elt(&(loader_ext->cached_seq),
                        loader->nrec);

    /* cached_ans_elt.seq is a const char * so we need to cast it to
       char * before we can write to it */
    Ocopy_bytes_to_i1i2_with_lkup(0, cached_ans_elt.length - 1,
        (char *) cached_ans_elt.seq, cached_ans_elt.length,
        dataline->seq, dataline->length,
        loader_ext->lkup_seq, loader_ext->lkup_length_seq);
    return;
}

static void SFF_load_qual(SFFloader *loader, const cachedCharSeq *dataline)
{
    SFF_loaderExt *loader_ext;
    cachedCharSeq cached_ans_elt;

    loader_ext = loader->ext;
    cached_ans_elt = get_cachedXRawList_elt(&(loader_ext->cached_qual),
                        loader->nrec);
    /* cached_ans_elt.seq is a const char * so we need to cast it to
       char * before we can write to it */
    Ocopy_bytes_to_i1i2_with_lkup(0, cached_ans_elt.length - 1,
        (char *) cached_ans_elt.seq, cached_ans_elt.length,
        dataline->seq, dataline->length,
        loader_ext->lkup_qual, loader_ext->lkup_length_qual);
    return;
}

static void SFF_load_qclip(SFFloader *loader, int start, int width)
{
    SFF_loaderExt *loader_ext;
    loader_ext = loader->ext;
    
	IRANGE_VALUES clip_result;
	clip_result = loader_ext->cached_qual_clip;

	clip_result.start[(int)loader->nrec] = start;
	clip_result.width[(int)loader->nrec] = width;

//	Rprintf("record: %i\tclip points qual:%i %i\t",
//            loader->nrec,clip_result.start[loader->nrec],clip_result.start[loader->nrec]);

    return;
}

static void SFF_load_aclip(SFFloader *loader, int start, int width)
{
    SFF_loaderExt *loader_ext;
    loader_ext = loader->ext;

	IRANGE_VALUES clip_result;
	clip_result = loader_ext->cached_adapt_clip;

	clip_result.start[loader->nrec] = start;
	clip_result.width[loader->nrec] = width;

//	Rprintf("clip points adapt:%i %i\n",
//            clip_result.start[loader->nrec],clip_result.start[loader->nrec]);
    return;
}

void freeLoader(SFF_loaderExt loader) {
    free(loader.cached_qual_clip.width);
    free(loader.cached_qual_clip.start);
    free(loader.cached_adapt_clip.width);
    free(loader.cached_adapt_clip.start);
}

static SFF_loaderExt new_SFF_loaderExt(SEXP seq, SEXP qual, SEXP lkup_seq, SEXP lkup_qual)
{
    SFF_loaderExt loader_ext;

    loader_ext.ans_names_buf =
        new_CharAEAE(get_XVectorList_length(seq), 0);
    loader_ext.cached_seq = cache_XVectorList(seq);
    loader_ext.cached_qual = cache_XVectorList(qual);

    loader_ext.cached_qual_clip.width = (int *) R_alloc((long) get_XVectorList_length(seq), sizeof(int));
    loader_ext.cached_qual_clip.start = (int *) R_alloc((long) get_XVectorList_length(seq), sizeof(int));

    loader_ext.cached_adapt_clip.width = (int *) R_alloc((long) get_XVectorList_length(seq), sizeof(int));
    loader_ext.cached_adapt_clip.start = (int *) R_alloc((long) get_XVectorList_length(seq), sizeof(int));

    if (lkup_seq == R_NilValue) {
        loader_ext.lkup_seq = NULL;
    } else {
        loader_ext.lkup_seq = INTEGER(lkup_seq);
        loader_ext.lkup_length_seq = LENGTH(lkup_seq);
    }
    if (lkup_qual == R_NilValue) {
        loader_ext.lkup_qual = NULL;
    } else {
        loader_ext.lkup_qual = INTEGER(lkup_qual);
        loader_ext.lkup_length_qual = LENGTH(lkup_qual);
    }
    return loader_ext;
}

static SFFloader new_SFF_loader(int load_seqids,
        SFF_loaderExt *loader_ext)
{
    SFFloader loader;

    loader.load_seqid = load_seqids ? &SFF_load_seqid : NULL;
    loader.load_seq = SFF_load_seq;
    loader.load_qual = SFF_load_qual;
    loader.load_qclip = SFF_load_qclip;
    loader.load_aclip = SFF_load_aclip;
    loader.nrec = 0;
    loader.ext = loader_ext;
    return loader;
}


/******************************* Main Function ********************************/
/* Read a sff file and print a commonHeader information
 * Argument:(0)    SEXP string - file path
 * Return:        SEXP string - file path
 */
static void
readSFF(SEXP string, int *recno, SFFloader *loader)
{
    // C Structures
    COMMONheader commonHeader;
    READheader header;

    cachedCharSeq dataline;

    // C declarations
    int i, padding_size, fres, load_record;
    char ch;
    uint16_t uint16;
    uint8_t uint8;
//    static uint8_t byte_padding[8];

    const char *string2 = CHAR(string);

    FILE* file = fopen (string2, "r");
    if(file==NULL)
        Rf_error("cannot open specified file: '%s'", string2);

    commonHeader = readCommonHeader(string2);
    // add fseek statement to bypass the commonHeader
    fseek(file, commonHeader.commonHeader_length , SEEK_SET);

    READdata data;

    // for every read,
    for(int read=0; read<commonHeader.number_of_reads; read++) {
        // Read Header Section
        fres = fread( &header.read_header_length, sizeof(uint16_t), 1, file);
        header.read_header_length = htons(header.read_header_length);
        fres = fread( &header.name_length, sizeof(uint16_t), 1, file);
        header.name_length = htons(header.name_length);
        fres = fread( &header.number_of_bases, sizeof(uint32_t), 1, file);
        header.number_of_bases = htonl(header.number_of_bases);
        fres = fread( &header.clip_qual_left, sizeof(uint16_t), 1, file);
        header.clip_qual_left = htons(header.clip_qual_left);
        fres = fread( &header.clip_qual_right, sizeof(uint16_t), 1, file);
        header.clip_qual_right = htons(header.clip_qual_right);
        fres = fread( &header.clip_adapter_left, sizeof(uint16_t), 1, file);
        header.clip_adapter_left = htons(header.clip_adapter_left);
        fres = fread( &header.clip_adapter_right, sizeof(uint16_t), 1, file);
        header.clip_adapter_right = htons(header.clip_adapter_right);

        // Determine whether or not to load the record
        load_record = loader != NULL;

        // save clip points
        if (load_record && loader->load_seq != NULL) {
            loader->load_qclip(loader, (int)header.clip_qual_left,(int)header.clip_qual_right-(int)header.clip_qual_left+1);
            loader->load_aclip(loader, (int)header.clip_adapter_left,(int)header.clip_adapter_right-(int)header.clip_adapter_left+1);
        }


        data.name = malloc(sizeof(char)*(header.name_length+1));
        if(data.name==NULL) Rf_error("id: cannot allocate memory");

        fres = fread( data.name, sizeof(char),header.name_length, file);
        data.name[header.name_length] = '\0';


        dataline.length = strlen(data.name);
        dataline.seq = data.name;
        if (load_record && loader->load_seqid != NULL) {
            loader->load_seqid(loader, &dataline);
        }

        padding_size =
        readHeaderPaddingSize(header.read_header_length,
                        header.name_length);
        fseek(file, padding_size , SEEK_CUR);
        //Read Data Section
        // Flows
        data.flows = (float*) malloc(sizeof(float)*(commonHeader.number_of_flows_per_read));
        for(i=0; i<commonHeader.number_of_flows_per_read; i++) {
            fres = fread( &uint16, sizeof(uint16_t), 1, file);
            uint16 = htons(uint16);
            data.flows[i] = uint16/100.0;
        }
        // indexes
        data.index = (int*) malloc(sizeof(int)*(header.number_of_bases));
        int cindex = 0;
        for(i=0; i<header.number_of_bases; i++) {
            fres = fread( &uint8, sizeof(uint8_t), 1, file);
            data.index[i] = (int)uint8 + cindex;
            cindex = data.index[i];
        }
        // bases
        data.bases = (char*) malloc(sizeof(char)*(header.number_of_bases+1));
        if (data.bases==NULL) Rf_error("cannot allocate memory");
        for(i=0; i<header.number_of_bases; i++) {
            fres = fread( &ch, sizeof(char), 1, file);
            data.bases[i] = ch;
            if(ch!='A' && ch!='T' && ch!='G' && ch!='C' && ch!='N') Rf_error("base error at %i:%c ",i, ch);
        }
        data.bases[header.number_of_bases]= '\0';

//TODO: SAMS function to adapter find/clip and assign MID

        dataline.length = strlen(data.bases);
        dataline.seq = data.bases;
        if (load_record && loader->load_seq != NULL)
            loader->load_seq(loader, &dataline);

        //quality
        data.quality = (char*) malloc(sizeof(char)*(header.number_of_bases+1));
        if (data.quality==NULL) Rf_error("cannot allocate memory");
        for(i=0; i<header.number_of_bases; i++) {
            fres = fread( &uint8, sizeof(uint8_t), 1, file);
            data.quality[i] = phredToChar(uint8);
        }
        data.quality[header.number_of_bases] = '\0';

        dataline.length = strlen(data.quality);
        dataline.seq = data.quality;
        if (load_record && loader->load_qual != NULL)
            loader->load_qual(loader, &dataline);

        if (load_record)
            loader->nrec++;

        if (debug == 1){
            debug_sffreader(read,commonHeader, header, data);
        }

        padding_size =
        readDataPaddingSize(commonHeader.number_of_flows_per_read,
                        header.number_of_bases);
        fseek(file, padding_size , SEEK_CUR);
        freeReadData(data);
        (*recno)++;
    } // end of for(every read)

    freeCommonHeader(commonHeader);
    fclose(file);
    return;
}
SEXP
sff_geometry(SEXP files)
{
    // C Structures
    COMMONheader commonHeader;
    READheader header;

    // C declarations
    int nrec, nfiles, i, recno, skip, padding_size, fres;
    const char *fname;
    nrec = recno = 0;
    static const char *names[] = {"number_of_reads","read_lengths"};
    // R declarations
    SEXP ans, ans_width, eltnm;

    if (!IS_CHARACTER(files))
        Rf_error("'%s' must be '%s'", "files", "character");

    nfiles = LENGTH(files);


    //  Retrieve number of records from the header(s)
    nrec = count_reads_sum(files);
    PROTECT(ans_width = NEW_INTEGER(nrec));

    for (i = 0; i < nfiles; ++i) {
        R_CheckUserInterrupt();
        fname = CHAR(STRING_ELT(files, i));

        FILE* file = fopen (fname, "r");
        if(file==NULL)
            Rf_error("cannot open specified file: '%s'", fname);

        commonHeader = readCommonHeader(fname);
        // add fseek statement to bypass the commonHeader
        fseek(file, commonHeader.commonHeader_length , SEEK_SET);

        // for every read,
        for(int read=0; read<commonHeader.number_of_reads; read++) {
            // Read Header Section
            fres = fread( &header.read_header_length, sizeof(uint16_t), 1, file);
            header.read_header_length = htons(header.read_header_length);
            fres = fread( &header.name_length, sizeof(uint16_t), 1, file);
            header.name_length = htons(header.name_length);
            fres = fread( &header.number_of_bases, sizeof(uint32_t), 1, file);
            header.number_of_bases = htonl(header.number_of_bases);
            fres = fread( &header.clip_qual_left, sizeof(uint16_t), 1, file);
            header.clip_qual_left = htons(header.clip_qual_left);
            fres = fread( &header.clip_qual_right, sizeof(uint16_t), 1, file);
            header.clip_qual_right = htons(header.clip_qual_right);
            fres = fread( &header.clip_adapter_left, sizeof(uint16_t), 1, file);
            header.clip_adapter_left = htons(header.clip_adapter_left);
            fres = fread( &header.clip_adapter_right, sizeof(uint16_t), 1, file);
            header.clip_adapter_right = htons(header.clip_adapter_right);

            skip = ((commonHeader.number_of_flows_per_read*flowgram_bytes_per_flow) +
                    (header.number_of_bases*3))+header.read_header_length-16;
            fseek(file,skip, SEEK_CUR);

            padding_size =
            readDataPaddingSize(commonHeader.number_of_flows_per_read,
                            header.number_of_bases);
            fseek(file, padding_size , SEEK_CUR);
            INTEGER(ans_width)[recno] = (int)(header.number_of_bases);
            recno++;

        } // end of for(every read)

        freeCommonHeader(commonHeader);
        fclose(file);
    }
    PROTECT(ans = allocVector(VECSXP,2));
    PROTECT(eltnm = allocVector( STRSXP, 2 ) );

    SET_VECTOR_ELT(ans, 0, ScalarInteger(nrec));
    SET_VECTOR_ELT(ans, 1, ans_width);
    for( int i = 0; i < 2; i++ )
        SET_STRING_ELT( eltnm, i, mkChar( names[i] ) );
    namesgets( ans, eltnm );
    UNPROTECT(3);
    return(ans);
}


SEXP
read_sff(SEXP files, SEXP use_names, SEXP lkup_seq, SEXP lkup_qual, SEXP verbose)
{
    int i, nfiles, recno,load_seqids,set_verbose, ans_length;
    SEXP fname, ans_geom, ans_names, header, nms;
    SEXP  ans = R_NilValue, reads = R_NilValue, quals = R_NilValue,
            qual_clip = R_NilValue, adapt_clip = R_NilValue;
	SEXP qclip_start, qclip_width, aclip_start, aclip_width;
    SFF_loaderExt loader_ext;
    SFFloader loader;

    if (!IS_CHARACTER(files))
        Rf_error("'%s' must be '%s'", "files", "character");

    nfiles = LENGTH(files);

    load_seqids = LOGICAL(use_names)[0];
    set_verbose = LOGICAL(verbose)[0];
    //  Retrieve SFF(s) Geometry
    PROTECT(ans_geom = sff_geometry(files));
    ans_length = INTEGER(VECTOR_ELT(ans_geom,0))[0];  //number of READS
    if(set_verbose) Rprintf("Total number of reads to be read: %d\n", ans_length);

    PROTECT(header = read_sff_header(files,verbose));
    PROTECT(reads = alloc_XRawList("DNAStringSet","DNAString",VECTOR_ELT(ans_geom,1)));
    PROTECT(quals = alloc_XRawList("BStringSet","BString",VECTOR_ELT(ans_geom,1)));
 
    loader_ext = new_SFF_loaderExt(reads, quals, lkup_seq,lkup_qual); //Biostrings/XStringSet_io.c
    loader = new_SFF_loader(load_seqids, &loader_ext); //Biostrings/XStringSet_io.c

    recno = 0;

    for (i = 0; i < nfiles; i++) {
        R_CheckUserInterrupt();
        fname = STRING_ELT(files, i);
        if(set_verbose) {
            Rprintf("reading file:%s\n", CHAR(fname));
        }
        readSFF(fname,&recno, &loader);
    }

    // load in the seq_ids
    if (load_seqids) {
        PROTECT(ans_names =
            new_CHARACTER_from_CharAEAE(&(loader_ext.ans_names_buf))); //IRanges
        set_XVectorList_names(reads, ans_names); //IRanges
//        set_XVectorList_names(quals, ans_names); //IRanges
        UNPROTECT(1);
    }

	PROTECT(qclip_start = NEW_INTEGER(ans_length));
	PROTECT(qclip_width = NEW_INTEGER(ans_length));
	memcpy(INTEGER(qclip_start), loader_ext.cached_qual_clip.start, sizeof(int) * ans_length);
	memcpy(INTEGER(qclip_width), loader_ext.cached_qual_clip.width, sizeof(int) * ans_length);
	PROTECT(qual_clip = new_IRanges("IRanges", qclip_start, qclip_width, R_NilValue));

	PROTECT(aclip_start = NEW_INTEGER(ans_length));
	PROTECT(aclip_width = NEW_INTEGER(ans_length));
	memcpy(INTEGER(aclip_start), loader_ext.cached_adapt_clip.start, sizeof(int) * ans_length);
	memcpy(INTEGER(aclip_width), loader_ext.cached_adapt_clip.width, sizeof(int) * ans_length);
	PROTECT(adapt_clip = new_IRanges("IRanges", aclip_start, aclip_width, R_NilValue));

    PROTECT(ans = NEW_LIST(5));
    SET_VECTOR_ELT(ans, 0, header);
    SET_VECTOR_ELT(ans, 1, reads); /* read */
    SET_VECTOR_ELT(ans, 2, quals); /* quality */
    SET_VECTOR_ELT(ans, 3, qual_clip); /* quality based clip points */
    SET_VECTOR_ELT(ans, 4, adapt_clip); /* adapter based clip points */
    UNPROTECT(11);

    PROTECT(nms = NEW_CHARACTER(5));
    SET_STRING_ELT(nms, 0, mkChar("header"));
	SET_STRING_ELT(nms, 1, mkChar("sread"));
    SET_STRING_ELT(nms, 2, mkChar("quality"));
    SET_STRING_ELT(nms, 3, mkChar("qualityClip"));
    SET_STRING_ELT(nms, 4, mkChar("adapterClip"));
    setAttrib(ans, R_NamesSymbol, nms);
	UNPROTECT(1);

    return ans;
}

/******************************* END Main Functions ********************************/

/******** Borrowing write_fastq from ShortRead *************/

/******************************* Write out phred qualities *************************/

char *
_cache_to_char(cachedXStringSet *cache, const int i,
               char *buf, const int width)
{
    cachedCharSeq roSeq = get_cachedXStringSet_elt(cache, i);
    if (roSeq.length > width)
        return NULL;
    strncpy(buf, roSeq.seq, roSeq.length);
    buf[roSeq.length] = '\0';
    return buf;
}

SEXP
write_phred_quality(SEXP id, SEXP quality, 
            SEXP fname, SEXP fmode, SEXP max_width)
{
    if (!(IS_S4_OBJECT(id) && 
          strcmp(get_classname(id), "BStringSet") == 0))
        Rf_error("'%s' must be '%s'", "id", "BStringSet");
    if (!(IS_S4_OBJECT(quality) && 
          strcmp(get_classname(quality), "BStringSet") == 0))
        Rf_error("'%s' must be '%s'", "quality", "BStringSet");
    const int len = get_XStringSet_length(id);
    if ((len != get_XStringSet_length(quality)))
        Rf_error("length() of %s must all be equal",
                 "'id', 'quality'");
    if (!(IS_CHARACTER(fname) && LENGTH(fname) == 1)) /* FIXME: nzchar */
        Rf_error("'%s' must be '%s'", "file", "character(1)");
    if (!(IS_CHARACTER(fmode) && LENGTH(fmode) == 1)) /* FIXME nchar()<3 */
        Rf_error("'%s' must be '%s'", "mode", "character(1)");
    if (!(IS_INTEGER(max_width) && LENGTH(max_width) == 1 &&
          INTEGER(max_width)[0] >= 0))
        Rf_error("'%s' must be %s", "max_width", "'integer(1)', >=0");
    const int width = INTEGER(max_width)[0];

    cachedXStringSet xid = cache_XStringSet(id),
                     xquality = cache_XStringSet(quality);

    FILE *fout = fopen(CHAR(STRING_ELT(fname, 0)), 
                       CHAR(STRING_ELT(fmode, 0)));
    if (fout == NULL)
        Rf_error("failed to open file '%s'", 
                 CHAR(STRING_ELT(fname, 0)));

    char *idbuf = (char *) R_alloc(sizeof(char), width + 1),
        *qualbuf = (char *) R_alloc(sizeof(char), width + 1);
    int i, j, phredval;
    for (i = 0; i < len; ++i) {
        idbuf = _cache_to_char(&xid, i, idbuf, width);
        if (idbuf == NULL){
			fclose(fout);
			Rf_error("failed to write record %d", i + 1);
		}
		fprintf(fout, ">%s\n",idbuf);
        
        qualbuf = _cache_to_char(&xquality, i, qualbuf, width);
        if (qualbuf == NULL){
			fclose(fout);
			Rf_error("failed to write record %d", i + 1);
		}
		j = 0;
        while(qualbuf[j] != '\0') {
			if (j != 0) fprintf(fout," ");
            phredval = CharToPhred(qualbuf[j]);
	        fprintf(fout, "%i", phredval);
			j++;
        }
		fprintf(fout,"\n");
    }
    fclose(fout);
    return R_NilValue;
}


#ifdef __cplusplus
}
#endif

