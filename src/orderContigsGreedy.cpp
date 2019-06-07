//Greedy ordering based on contig quality.
#include <Rcpp.h>
#include <limits.h>
#include "orderContigsGreedy.h" 
using namespace Rcpp;

int nRows;
int nCols;
IntegerMatrix strandMatrix;

int quick_score_cell(int *order, int pos, int size, int cell, int numSCE) {
  int elm = strandMatrix(order[size], cell);
  
  int k = pos - 1;
  // skip backward over 0 contigs to find prev_state at pos
  while(k >= 0 && strandMatrix(order[k], cell) == 0) --k;
  int p0 = k < 0 ? 0 : strandMatrix(order[k], cell);
  
  // if new contig is 0 or same as prev_state then score unchanged
  if(elm == 0 || p0 == elm) { return numSCE; }
  
  // add score for new contig
  if(p0 != 0) ++numSCE;
  if((p0 == 1 && elm == 2) || (p0 == 2 && elm == 1)) numSCE += 5;
  
  int p1 = elm, j = pos;
  
  // skip 0 contigs; they have no contribution to score but let changed prev_state pass through	
  while(j < size && strandMatrix(order[j], cell) == 0) { ++j; }
  
  if(j < size) {
    int s = strandMatrix(order[j], cell);
    
    // remove old contribution from shifted contig
    if(p0 != s) {
      if(p0 != 0) --numSCE;
      if((p0 == 1 && s == 2) || (p0 == 2 && s == 1)) numSCE -= 5;
    }
    
    // calculate new contribution from shifted contig
    if(p1 != s) {
      ++numSCE; // we know p1 != 0;
      if((p1 == 1 && s == 2) || (p1 == 2 && s == 1)) numSCE += 5;
    }
  }
  
  return numSCE;  
}

/*****************************************************************************
* Function to greedily order contigs based on their strand state using contig quality 
* 
****************************************************************************/

//@Param SStrandMatrix  Strand state for all contigs across all libraries
RcppExport  SEXP orderContigsGreedy(
    SEXP SStrandMatrix                
)
{
  IntegerMatrix sMatrix(SStrandMatrix);
  strandMatrix = sMatrix; 
  nRows = strandMatrix.nrow();
  nCols = strandMatrix.ncol();
  
  int *order = new int[nRows * nCols];
  int *score = new int[nCols];
  for(int cell = 0; cell < nCols; ++cell) {
    score[cell] = 0;
  }
  
  order[0] = 0;
  int best_score = 0;
  for(int i = 1; i < nRows; i++) {
    order[i] = i;
    best_score = INT_MAX;

    int best_index = 0;
    for(int j = 0; j <= i; j++) {
      int total = 0;      
      for(int cell = 0; cell < nCols; ++cell) {
        total += quick_score_cell(order, j, i, cell, score[cell]);
      }
      
      if(total < best_score) {
        best_index = j;
        best_score = total;
      }
    }
    
    // update cell scores
    for(int cell = 0; cell < nCols; ++cell) {
      score[cell] = quick_score_cell(order, best_index, i, cell, score[cell]);
    }
    
    // rearrange order to best new order
    for(int k = i; k > best_index; --k) {
      order[k] = order[k-1];
    }
    order[best_index] = i;
  }
  
  IntegerVector result(nRows);
  for(int i = 0; i < nRows; i++) {
    result[i] = order[i] + 1;
  }
  
  delete[] order;
  delete[] score;
  
  return Rcpp::List::create(Rcpp::Named("order") = result, Rcpp::Named("score") = best_score);
}