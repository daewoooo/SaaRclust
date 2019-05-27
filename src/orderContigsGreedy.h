#ifndef _contiBAIT_greedyContigs_
#define _contiBAIT_greedyContigs_

RcppExport  SEXP orderContigsGreedy ( 
  							SEXP SStrandMatrix// Strand state for all contigs across all libraries. Should be a numeric matrix, with WW=1, WC=2, CC=3, NA=0. Should be ordered by quality.
              						
);
float score_order(std::vector<int> order_vector);
#endif