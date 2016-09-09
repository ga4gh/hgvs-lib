/**
 * Copyright 2016 Golden Hellix Inc.
 *
 * HGVS Library for taging genomic variants, projecting them to a
 * transcript alignment to a reference sequence, and producing HGVS
 * descriptions of those variants along with other complimentary
 * attributes.
 *
 */

#ifndef HGVS_H
#define HGVS_H

// C99 expected
#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>
#include <limits.h>
#include <float.h>
#include <math.h>

// Some C byob (bring your own buffer) size bounds
#define HGVS_MAX_ERR_MSG 1024
#define HGVS_MAX_OUTPUT_LEN 1024

/**
 * This structured must be constructed from some external source. Some
 * fields are optional. This structure resembles the UCSC table for gene
 * records and is capable of represent GFF/GTF gene records.
 */
typedef struct hgvs_tx_t {
  // Coordinates are 0-based, half-open
  const char* chr;
  int g_start;
  int g_stop;

  const char* tx_id;       // Used by HGVS c.
  const char* protein_id;  // Used by HGVS p., if NULL, use tx_name
  const char* gene_name;   // Optional, can be NULL

  char strand; //'-' if reverse, else forward

  // Coding sequence start/stop, in [g_start,g_stop]. 0 if non-coding.
  int cds_start;
  int cds_stop;

  // Pairs of exon start/stops, in [g_start,g_stop].
  // exon_count >= 1, and defines length of array of start/stops.
  int* exon_starts;
  int* exon_stops;
  int exon_count;

  // Following fields computed from above with hgvs_tx_validate
  bool validated; // Once true, hgvs_tx_validate is NO-OP
  bool is_valid;  // Failed constraints on record fields
  bool is_coding; // Has CDS start/stop set with valid ranges
  bool is_reverse_strand; // (strand == '-')

  // Following fields computed from above with hgvs_tx_encode
  bool encoded; // Once true, hgvs_tx_encode is NO-OP
  bool is_complete; // Has start codon and single stop codon where expected

  char* genomic_seq; // reference sequence from g_start to g_stop
  char* coding_seq; // Spliced coding transcript
  char* aa_seq; // Encoded amino acid sequence (translated coding_seq)

  struct hgvs_exon_t* exons; // Owned, constructed from above data
} hgvs_tx_t;

/*
        ==========[ Notes on coordinate systems ]==========

These classes use several different coordinate systems at various times. To avoid confusion, variable
and method names are often suffixed with the coordinate system they use. The systems are:

1. Chromosome - abbreviated as 'g'
These are half-open coordinates relative to the start of the reference sequence for a chromosome.

2. Transcript - abbreviated as 'ts'
These are half-open coordinates relative to a transcription start site and in the direction of
transcription (5' to 3' on the sense strand). They can also be thought of as coordinates within
the non-spliced primary RNA transcript.

3. Spliced - abbreviated as 'sp'
These are half-open coordinates inside of the exonic sequence for a transcript, and correspond to
the post-spliced RNA product. Spliced coordinates start at the 5' end of the sense strand and
stop at the 3' end.

4. Coding - abbreviated as 'cd'
These are half-open coordinates inside of the coding exonic region for a transcript. Coordinates start
at the beginning of the initiation codon. Dividing by 3 yields a coordinate in the amino acid sequence. 
*/

typedef struct hgvs_exon_t {
  hgvs_tx_t* tx; // Containing transcript
  int index; // 0-based. Human friendly ordinal is index+1

  int g_start;
  int g_stop;

  int ts_start;
  int ts_stop;

  int sp_start;
  int sp_stop;

  int cd_start;
  int cd_stop;

  bool is_coding;
} hgvs_exon_t;

typedef struct hgvs_var_t {
  const char* chr;
  int g_start;
  int g_stop;

  // Always set. Insertions will have g_ref = "", deletions will have g_alt = ""
  const char* g_ref;
  const char* g_alt;

  int ref_length; // Length of g_ref
  int alt_length; // Length of g_alt

  bool is_snp; // Single letter replacement, otherwise MNP or length polymorphism

} hgvs_var_t;

typedef struct hgvs_vartx_t {
  bool is_valid;

  hgvs_var_t* var;
  hgvs_tx_t* tx;

  // Position of variant in other coordinate systems (see note on
  // coordinate systems above)
  int ts_var_start;
  int ts_var_stop;

  int sp_var_start;
  int sp_var_stop;

  int cd_var_start;
  int cd_var_stop;

  // Overlap a single exon:
  hgvs_exon_t* exon;

  // Intronic, between two flanking exons
  hgvs_exon_t* flanking_donar_exon; // splice donar (5' end)
  hgvs_exon_t* flanking_acceptor_exon; // splice acceptor (3' end)

  // Varaint reference and alternate (on coding strand)
  const char* c_ref;
  const char* c_alt;

  // Only set when is_snp. The seq triplet and AA char for ref/alt codon
  char* seq_ref_codon;
  char* seq_alt_codon;
  char aa_ref_codon;
  char aa_alt_codon;

  // Only set when !is_snp
  char* c_alt_tx; // Full alternate transcript, compared to tx->coding_tx;

  // When !is_snp, compute full transcript amino acid sequence with and
  // without variant mutation.
  char* aa_ref_tx;
  char* aa_alt_tx;

} hgvs_vartx_t;


/*
  const char* chr;
  int start;
  int stop;
  qint64 id() const                           { return _id; }
  const CString& chr() const                  { return _chr; }
  const GAIntervalI32& bounds() const         { return _bounds; }
  int start() const                           { return _bounds.start; }
  int stop() const                            { return _bounds.stop; }
  
  const CString& geneName() const             { return _geneName; }
  const CString& transcriptName() const       { return _transcriptName; }
  const CString& LRGId() const                { return _LRGId; }
  int strand() const                          { return _strand; } //1 for forward, -1 for reverse
  bool isCoding() const                       { return _isCoding; }
  const GAIntervalI32& cds() const            { return _cds_chr; }
  int cdsStart() const                        { return _cds_chr.start; }
  int cdsStop() const                         { return _cds_chr.stop; }
  int tsStart() const                         { return exonAt(0).start_ts(); }
  int tsStop() const                          { return tsStart() + size_ts(); }
  int spStart() const                         { return exonAt(0).start_sp(); }
  int spStop() const                          { return spStart() + size_sp(); }
  int startSite_chr() const                   { return _startSite_chr; }
  int size_ts() const                         { return _bounds.sizeAbs(); }
  int size_sp() const                         { return _size_sp; }
  int size_cd() const                         { return _cds_sp.sizeAbs(); }
  const GAIntervalI32 &cds_sp() const         { return _cds_sp; }

} hgvs_tx;
*/

// Validate a transcript internal data, construct derived data structures
bool hgvs_tx_validate(hgvs_tx_t* tx_data, char* err_msg);

// Parse the provided genomic reference sequence read from tx_data->chr
// interval [tx_data->g_start to tx_data->g_stop] into encoded form.
bool hgvs_tx_encode(hgvs_tx_t* tx_data, const char* g_tx_seq, char* err_msg);

// Free data owned by a tx
void hgvs_free_tx(hgvs_tx_t* tx);

// Prepare a variant
bool hgvs_parse_var(hgvs_var_t* var, const char* chr, int pos, const char* ref, const char* alt,
                    char* err_msg);
// Free data owned by a var
void hgvs_free_var(hgvs_var_t* var);


// Compute transcript variant pair
bool hgvs_map_var_to_tx(hgvs_vartx_t* var_tx, const hgvs_tx_t* tx_data, const hgvs_var_t* var,
                        char* err_msg);
// Free data owned by a vartx
void hgvs_free_vartx(hgvs_vartx_t* var);

// Format g.
bool hgvs_g_dot(hgvs_var_t* var, char* hgvs_gdot, char* err_msg);

// Format c.
bool hgvs_c_dot(hgvs_vartx_t* var_tx, char* hgvs_pdot, char* err_msg);

// Format p.
bool hgvs_p_dot(hgvs_vartx_t* var_tx, char* hgvs_pdot, char* err_msg);

// Coordinate macros

// Interval i(start,stop) containment value o
// Interval i contains o iff every value contained by o is also contained by i
#define CONTAINS_OPEN(i_start, i_stop, value) (i_start < value && value < i_stop)

// Half open interval intersection test
// Interval i intersects o iff there exists at least one value v
//   such that i.containsHalfOpen(v) and o.containsHalfOpen(v)
// Adjacent intervals do not intersect
#define INTERSECTS_HALF_OPEN(i_start, i_stop, o_start, o_stop) (i_start < o_stop && o_start < i_stop)
#endif
