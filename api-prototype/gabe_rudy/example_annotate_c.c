#include "hgvs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// ----
// Example of a client of the HGVS lib.
//
// Read `annotate_variant_stdout_tsv` which gives a feel of how this lib
// would be used to take a variant, look up transcripts (not implemented)
// and use the HGVS lib to not only format c. and p. but also derive SO
// terms from the same computed attributes about a Variant+Transcript
// pair.
// ----

typedef struct tx_itr {
  void* ptr; // some implementation specific data
} tx_itr;

// Takes 0-based half-open interval and returns iterator that will yield
// transcript records that overlap the interval by subsequent calls to
// query_iter_next
tx_itr* _query_overlapping_transcripts(const char* chr, int g_start, int g_stop);
bool _query_iter_next(tx_itr* itr, hgvs_tx_t* tx_data);

// Needs to be implemented to read genomic reference for 0-based half-open interval
const char* read_genomic_reference(const char* chr, int g_start, int g_stop);

// Takes the hgvs_var_tx structure and utalizes is to compute the
// appropriate Sequence Ontology term for a variant-transcript pair.
const char* _get_ontology(hgvs_vartx_t* var_tx);

// Example annotating a variant, printing to stdout a TSV line
// representing the annotations. Input is VCF columns 1-4.
//
// Each line contains:
// Chr Ref Alt HGVS g.
bool annotate_variant_stdout_tsv(char* chr, int pos, char* ref, char* alts)
{
  // Split alts, as we are only annotating a single variant
  char err_msg[HGVS_MAX_ERR_MSG];
  char hgvs_gdot[HGVS_MAX_OUTPUT_LEN];
  char hgvs_cdot[HGVS_MAX_OUTPUT_LEN];
  char hgvs_pdot[HGVS_MAX_OUTPUT_LEN];
  char* alt;

  while ((alt = strsep(&alts, ",")) != NULL) {
    // Parse variant, does some minimal normalization (remove common prefix/suffix)
    hgvs_var_t var;
    if (!hgvs_parse_var(&var, chr, pos, ref, alt, err_msg)) {
      fprintf(stderr, "Error while parsing %s:%d %s/%s: %s", chr, pos, ref, alt, err_msg);
      return false;
    }

    // Construct HGVS g. - no transcript required
    if (!hgvs_g_dot(&var, hgvs_gdot, err_msg)) {
      fprintf(stderr, "Error while constructing g. %s:%d %s/%s: %s", chr, pos, ref, alt, err_msg);
      return false;
    }

    // This lib is not responsible for looking up proper NCBI accession
    // for a (chr_name, genomic_assembly_version), but that *should* be
    // used here to produce the full HGVS g. instead of the 'chr' field.
    const char* g_accession = var.chr;

    // Look up transcripts that overlap this position
    int tx_overlaps = 0;
    tx_itr* iter = _query_overlapping_transcripts(var.chr, var.g_start, var.g_stop);
    hgvs_tx_t tx_data;
    while (_query_iter_next(iter, &tx_data)) {
      // For each transcript now in memory in a record form, we can start
      // to figure out how this variant interacts with that transcript.
      hgvs_tx_validate(&tx_data, err_msg);
      if(!tx_data.is_valid){
        fprintf(stderr, "Invalid tx record: %s", err_msg);
        continue; //Bad tx record data
      }

      if(!tx_data.is_coding){
        continue; // We don't report on non-coding transcripts
      }

      const char* tx_seq = read_genomic_reference(tx_data.chr, tx_data.g_start, tx_data.g_stop);
      hgvs_tx_encode(&tx_data, tx_seq, err_msg);
      if(!tx_data.is_complete){
        continue; // We don't report on incomplete transcripts
      }

      // Map variant to transcript features
      hgvs_vartx_t var_tx;
      hgvs_map_var_to_tx(&var_tx,&tx_data, &var, err_msg);
      if(!var_tx.is_valid) {
        fprintf(stderr, "Unable to map var to tx: %s", err_msg);
        continue;
      }
      tx_overlaps++;

      // For "c.", mostly need to translate ref/alt to tx strand and get to tx start/stop
      // For intronic variants we need to find the nearest exon
      if (!hgvs_c_dot(&var_tx, hgvs_cdot, err_msg)) {
        fprintf(stderr, "Error while constructing c. %s:%d %s/%s: %s", chr, pos, ref, alt, err_msg);
        return false;
      }

      // Compute ful reference and alternate amino-acid sequence for this
      // transcript and variant, describe the delta between them
      // appropriately, using some SO like classifications
      // (i.e. stop-gains are not describe as "fs" variants.
      if (!hgvs_p_dot(&var_tx, hgvs_pdot, err_msg)) {
        fprintf(stderr, "Error while constructing c. %s:%d %s/%s: %s", chr, pos, ref, alt, err_msg);
        return false;
      }

      const char* ontology = _get_ontology(&var_tx);

      printf("%s\t%s\t%s:%s\t%s\t%s\t%s\t%s", var.chr, var.g_ref, chr, g_accession, hgvs_gdot, hgvs_cdot, hgvs_pdot, ontology);
      hgvs_free_vartx(&var_tx);
    }

    if (tx_overlaps == 0) {
      const char* ontology = "Intergenic";
      printf("%s\t%s\t%s:%s\t\t\t%s\t%s", var.chr, var.g_ref, chr, g_accession, hgvs_gdot, ontology);
    }
    // Cleanup
    free(iter);
    hgvs_free_var(&var);
  }

  return true;
}

// Takes 0-based half open interval and returns iterator that will yield
// transcript records that overlap the interval by subsequent calls to
// query_iter_next
tx_itr* _query_overlapping_transcripts(const char* chr, int g_start, int g_stop)
{
  // This function would query a transcript data source backend (file
  // based, or otherwise) and create an tx_itr capable of yielding
  // all the transcript records that overlap this genomic interval.
  //
  // Note, it's not the responsibility of the library to:
  //
  // * Curate or define an appropriate set of transcripts for a genomic
  //  build
  //
  // * Guarantee the matching of the genomic build of the variant and
  //   transcript source
  //
  // * Filter transcripts down to only those appropriate to fully
  //   annotate (i.e. filter out XM_ or pseudo-gene transcripts).
  //
  // * Pick the "best" transcript for a gene to report

  tx_itr* iter = calloc(sizeof(tx_itr), 1);
  return iter;
}

bool _query_iter_next(tx_itr* itr, hgvs_tx_t* tx_data)
{
  // ...
  return false;
}

#define ACCEPTOR_DISTANCE 2 // Could be a param

// Example of some representative classifier functions


bool _is_stop_gain(hgvs_vartx_t* var_tx)
{
  if(!var_tx->exon) return false;
  if(!var_tx->tx->is_coding) return false;

  hgvs_var_t* var = var_tx->var;

  if(var->is_snp){
    return var_tx->aa_ref_codon != '*' && var_tx->aa_alt_codon != '*';
  } else {
    // TODO: Port more Golden Helix code to form out API

    // we need to localize the search to the codon we are interested
    //in otherwise a frameshift variant will almost always get called
    //as stop gain since the fs will introduce 1 (or many) stop
    //codons downstream
    /*
    GAIntervalI32 var_cd = _tx->cdFromChr(_iv);
    int var_start_cd = var_cd.start;
    int var_stop_cd = var_cd.stop;
    //clamp to exon
    if(var_start_cd == INT_MISSING)
      var_start_cd = _ex.start_cd();
    if(var_stop_cd == INT_MISSING)
      var_stop_cd = _ex.stop_cd();
    //adjust to codon boundaries
    if(var_start_cd == var_stop_cd)
      var_stop_cd += 1;
    while(var_start_cd % 3 != 0)
      var_start_cd--;
    while(var_stop_cd % 3 != 0)
      var_stop_cd++;
    //make interval for ease of comparing
    GAIntervalI32 var_iv_cd = GAIntervalI32(var_start_cd, var_stop_cd);
    var_iv_cd.normalize();

    int ref_stop_idx = _tx->getSeq().proteinSeq.indexOf("*");
    int alt_stop_idx = _altTxSeq.proteinSeq.indexOf("*");

    int sizecood = _altTxSeq.proteinSeq.size()-1;
    int sizecood2 = ref_stop_idx + _getTXVarDelta()/3;

    return alt_stop_idx != -1 &&
      var_iv_cd.containsHalfOpen(alt_stop_idx*3) && // localize to codon containing variant
      alt_stop_idx < ref_stop_idx + _getTXVarDelta()/3;
    */
    return true;
  }
}

bool _is_splice_acceptor(hgvs_vartx_t* var_tx)
{
  // flanking_exon_donar and flanking_exon_acceptor are set when
  // projecing a variant to a transcript that is contained by a intron
  // between two exons (does not overlap either exon).
  if(!var_tx->flanking_acceptor_exon) return false;

  int acceptor_region_stop = var_tx->flanking_acceptor_exon->ts_start;
  int acceptor_region_start = acceptor_region_stop - ACCEPTOR_DISTANCE;

  if(var_tx->var->ref_length == 0 &&
     CONTAINS_OPEN(acceptor_region_start, acceptor_region_stop, var_tx->var->g_start))
    return true;
  if(var_tx->var->ref_length > 0 &&
     INTERSECTS_HALF_OPEN(acceptor_region_start, acceptor_region_stop, var_tx->var->g_start, var_tx->var->g_stop))
    return true;

  return false;
}

const char* _get_ontology(hgvs_vartx_t* var_tx)
{
  // 
  if(_is_exon_deleted(var_tx)){
    return "exon_loss_variant";
  } else if(_is_initiator_codon(var_tx)){
    return "initiator_codon_variant"; // start loss
  } else if(_is_stop_gain(var_tx)){
    return "stop_gained";
  } else if(_is_stop_lost(var_tx)){
    return "stop_lost";
  } else if(_is_start_gain(var_tx)){
    return "5_prime_UTR_premature_start_codon_gain_variant";
  } else if(_is_splice_acceptor(var_tx)) {
    return "splice_acceptor_variant";
  } else if(_is_splice_donor(var_tx)) {
    return "splice_donor_variant";
  } else if(_is_stop_retained(var_tx)) {
    return "stop_retained_variant";
  } else if(_is_inframe_deletion(var_tx)) {
    return "inframe_deletion";
  } else if(_is_inframe_insertion(var_tx)) {
    return "inframe_insertion";
  } else if(_is_disruptive_inframe_deletion(var_tx)) {
    return "disruptive_inframe_deletion";
  } else if(_is_disruptive_inframe_insertion(var_tx)) {
    return "disruptive_inframe_insertion";
  } else if(_is_frameshift(var_tx)) {
    return "frameshift_variant";
  } else if(_is_synonymous(var_tx)) {
    return "synonymous_variant";
  } else if(_is_missense(var_tx)) {
    return "missense_variant";
  } else if(_is_splice_region(var_tx)) {
    return "splice_region_variant";
  } else if(_is_non_coding_exon(var_tx)) {
    return "non_coding_exon_variant";
  } else if(_is_intron(var_tx)) {
    return "intron_variant";  // crucial that donor/acceptor region vars have already been found
  } else if(_is_3_prime_UTR(var_tx)) {
    return "3_prime_UTR_variant";
  } else if(_is_5_prime_UTR(var_tx)) {
    return "5_prime_UTR_variant";
  } else {
    return "unkown"; // Should not happen
  }
}

