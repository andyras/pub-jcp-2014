#include "plots.hpp"

//#define DEBUG_PLOT

using namespace std;

void makePlots(map<const string, bool> &outs, struct PARAMETERS * p) {
  // populations in subsystems
  if (isOutput(outs, "populations.plt")) {
    plotPopulations("populations.plt", p);
  }

  // probabilities in k states
  if (isOutput(outs, "kprobs.plt") && (p->Nk > 1)) {
    plotKProbs("kprobs.plt", p);
  }
  else {
    if (p->Nk < 2) {
      cerr << "WARNING: <= 1 k states, not making kprobs.plt" << endl;
    }
  }

  // probabilities in c states
  if (isOutput(outs, "cprobs.plt") && (p->Nc > 1)) {
    plotCProbs(p);
  }
  else {
    if (p->Nc < 2) {
      cerr << "WARNING: <= 1 c states, not making cprobs.plt" << endl;
    }
  }

  // calculated Fermi level
  if (isOutput(outs, "muFromPops.plt") && (p->Nk > 1)) {
    if (!isOutput(outs, "muFromPops.out")) {
      std::cout << "WARNING: making muFromPops.plt although muFromPops.out does not exist." << std::endl;
    }
    plotMuFromPops("muFromPops.plt", p);
  }

  // density matrix in time
  if (isOutput(outs, "dmt_z.plt")) {
    plotDMt_z("dmt_z.plt", p);
  }

  // populations in k states as a movie
  if (isOutput(outs, "kprobs_movie.plt") && (p->Nk > 1)) {
    plotKProbsMovie("kprobs_movie.plt", p);
  }
  else {
    if (p->Nk < 2) {
      cerr << "WARNING: <= 1 k states, not making kprobs_movie.plt" << endl;
    }
  }

  // populations in c states as a movie
  if (isOutput(outs, "cprobs_movie.plt") && (p->Nc > 1)) {
    plotCProbsMovie("cprobs_movie.plt", p);
  }
  else {
    if (p->Nc < 2) {
      cerr << "WARNING: <= 1 c states, not making cprobs_movie.plt" << endl;
    }
  }

  return;
}

/* Plots the populations in different subsystems over time,
 * separately and on the same plot. */
void plotPopulations(char * fileName, struct PARAMETERS * p) {
  ofstream o(fileName);

  o << "#!/usr/bin/env gnuplot" << endl;
  o << endl;
  o << "reset" << endl;
  o << endl;
  o << "set style data lines" << endl;
  o << "set style line 1 lt 1 lc rgb 'red'" << endl;
  o << "set style line 2 lt 1 lc rgb 'dark-green'" << endl;
  o << "set style line 3 lt 1 lc rgb 'blue'" << endl;
  o << "set style increment user" << endl;
  o << endl;
  o << "set terminal pdfcairo enhanced dashed size 5,5 lw 2 font 'Arial-Bold,12'" << endl;
  o << "set output 'figures/populations.pdf'" << endl;
  o << endl;
  o << "stats './outs/tkprob.out' u ($1/41.3414):2 nooutput" << endl;
  o << "set xrange [STATS_min_x:STATS_max_x]" << endl;
  o << endl;
  o << "set rmargin 3" << endl;
  o << "set tmargin 0" << endl;
  o << "set bmargin 3" << endl;
  o << "set tics scale 0 nomirror" << endl;
  o << "unset xtics" << endl;
  o << "set format y '%.2e'" << endl;
  o << endl;
  o << "set multiplot" << endl;
  if (p->bridge_on) {
    o << "set size 1,0.31" << endl;
  }
  else {
    o << "set size 1,0.47" << endl;
  }
  o << endl;
  if (p->bridge_on) {
    o << "set origin 0,0.62" << endl;
  }
  else {
    o << "set origin 0,0.48" << endl;
  }
  o << "set ylabel 'Population (a.u.)'" << endl;
  o << "set title 'Bulk Population vs. Time'" << endl;
  o << "plot './outs/tkprob.out' u ($1/41.3414):2 lt 1 notitle" << endl;
  o << endl;
  if (p->bridge_on) {
    o << "set origin 0,0.32" << endl;
    o << "set title 'Bridge Population vs. Time'" << endl;
    o << "plot './outs/bprobs.out' u ($1/41.3414):2 lt 3 notitle";
    if (p->Nb > 1) {
      for (int ii = 0; ii < (p->Nb - 1); ii++) {
	o << ", \\" << endl << "     '' u ($1/41.3414):" << (ii+3) << " lt " << (ii+4) << " notitle";
      }
    }
    o << endl << endl;
  }
  o << "set origin 0,0.01" << endl;
  o << "set xlabel 'Time (fs)'" << endl;
  o << "set xtics scale 0" << endl;
  o << "set title 'QD Population vs. Time'" << endl;
  o << "plot './outs/tcprob.out' u ($1/41.3414):2 lt 2 notitle" << endl;
  o << endl;
  o << "unset multiplot" << endl;
  o << endl;
  o << "set terminal pdfcairo enhanced size 5,3 lw 2 font 'Arial-Bold,12' dashed rounded" << endl;
  o << "set output 'figures/populationsTogether.pdf'" << endl;
  o << endl;
  o << "set key tmargin right" << endl;
  o << "set size 1,1" << endl;
  o << "set tmargin 4" << endl;
  o << "set ytics 1 format '%.0f'" << endl;
  o << "unset style line" << endl;
  o << "set style line 1 lc rgb 'red'" << endl;
  o << "set style line 2 lc rgb 'dark-green'" << endl;
  o << "set style line 3 lc rgb 'blue'" << endl;
  o << endl;
  o << "set title 'Subsystem Populations vs. Time' offset 0,0.6" << endl;
  o << "plot './outs/tkprob.out' u ($1/41.3414):2 lt 1 lw 2 title 'Bulk', \\" << endl;
  if (p->bridge_on) {
    o<< "'./outs/tbprob.out' u ($1/41.3414):2 lt 3 lw 2 title 'Bridge', \\" << endl;
  }
  o<< "'./outs/tcprob.out' u ($1/41.3414):2 lt 2 lw 2 title 'QD'" << endl;

  return;
}

void plotKProbs(char * fileName, struct PARAMETERS * p) {
#ifdef DEBUG_PLOT
  cout << "\nMaking kprobs.plt" << endl;
#endif
  ofstream o(fileName);

  o << "#!/usr/bin/env gnuplot" << endl;
  o << endl;
  o << "set terminal pdfcairo size 5,3 font 'Arial-Bold,12'" << endl;
  o << endl;
  o << "set output '/dev/null'" << endl;
  o << "plot '<cut -d \" \" -f 2- ./outs/kprobs.out' u ($2*"
    << p->tout << "/" << p->numOutputSteps << "):($1*(" 
    << p->kBandTop << "-" << p->kBandEdge << ")/(";
  o << p->Nk << "-1)):3 matrix with image" << endl;
  o << endl;
  o << "set output 'figures/kprobs.pdf'" << endl;
  o << "set title 'Electron probability density in bulk conduction band'" << endl;
  o << "unset key " << endl;
  o << "set border 0" << endl;
  o << "set tics scale 0" << endl;
  o << "set ylabel 'Energy above band edge (a.u.)'" << endl;
  o << "set xlabel 'Time (a.u.)'" << endl;
  o << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]" << endl;
  o << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]" << endl;
  o << endl;
  o << "set palette defined(\\" << endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << endl;
  o << "1      0.7059  0.0157  0.1490\\" << endl;
  o << ")" << endl;
  o << endl;
  o << "replot" << endl;

  return;
}

/* Makes a gnuplot file to plot the QD populations over time */
void plotCProbs(struct PARAMETERS * p) {
#ifdef DEBUG_PLOT
  cout << "\nMaking cprobs.plt" << endl;
#endif
  ofstream o("cprobs.plt");

  o << "#!/usr/bin/env gnuplot" << endl;
  o << endl;
  o << "set terminal pdfcairo size 5,3 font 'Arial-Bold,12'" << endl;
  o << endl;
  o << "set output '/dev/null'" << endl;
  o << "plot '<cut -d \" \" -f 2- ./outs/cprobs.out' u ($2*"
    << p->tout << "/" << p->numOutputSteps << "):($1):3 matrix with image" << endl;
  o << endl;
  o << "set output 'figures/cprobs.pdf'" << endl;
  o << "set title 'Electron probability density in QD'" << endl;
  o << "unset key " << endl;
  o << "set border 0" << endl;
  o << "set tics scale 0" << endl;
  o << "set ylabel 'State (index above band edge)'" << endl;
  o << "set xlabel 'Time (a.u.)'" << endl;
  o << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]" << endl;
  o << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]" << endl;
  o << endl;
  o << "set palette defined(\\" << endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << endl;
  o << "1      0.7059  0.0157  0.1490\\" << endl;
  o << ")" << endl;
  o << endl;
  o << "replot" << endl;

  return;
}

void plotMuFromPops(char * fileName, struct PARAMETERS * p) {
#ifdef DEBUG_PLOT
  cout << "\nMaking muFromPops.plt" << endl;
#endif
  ofstream o(fileName);

  o << "#!/usr/bin/env gnuplot" << endl;
  o << endl;
  o << "set terminal pdfcairo size 5,3 font 'Arial-Bold,12'" << endl;
  o << endl;
  o << "set output '/dev/null'" << endl;
  o << "plot '<cut -d \" \" -f 2- ./outs/muFromPops.out' u ($2*"
    << p->tout << "/" << p->numOutputSteps << "):($1*(" 
    << p->kBandTop << "-" << p->kBandEdge << ")/(";
  o << p->Nk << "-1)):3 matrix with image" << endl;
  o << endl;
  o << "set output 'figures/muFromPops.pdf'" << endl;
  o << "set title 'Calculated Fermi level in bulk conduction band'" << endl;
  o << "unset key " << endl;
  o << "set border 0" << endl;
  o << "set tics scale 0" << endl;
  o << "set ylabel 'Energy above band edge (a.u.)'" << endl;
  o << "set xlabel 'Time (a.u.)'" << endl;
  o << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]" << endl;
  o << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]" << endl;
  o << endl;
  o << "set palette defined(\\" << endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << endl;
  o << "1      0.7059  0.0157  0.1490\\" << endl;
  o << ")" << endl;
  o << endl;
  o << "replot" << endl;

  return;
}

/* Plots the density matrix in time as a series of .png files */
void plotDMt_z(char * fileName, struct PARAMETERS * p) {
  ofstream o(fileName);

  o << "#!/usr/bin/env gnuplot" << endl;
  o << endl;
  o << "![ -d img ] && rm -rf img" << endl;
  o << "!mkdir -p img" << endl;
  o << endl;
  o << "set palette defined(\\" << endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << endl;
  o << "1      0.7059  0.0157  0.1490\\" << endl;
  o << ")" << endl;
  o << endl;
  o << "set terminal pngcairo font 'Arial-Bold,12'" << endl;
  o << "set output '/dev/null'" << endl;
  o << "plot 'outs/dmt_z.out' index 0 matrix with image" << endl;
  o << "set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]" << endl;
  o << "set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]" << endl;
  o << endl;
  o << "set tics out scale 0.5 nomirror" << endl;
  o << "set size square" << endl;
  o << endl;
  o << "do for [ii=0:" << p->numOutputSteps << "] {" << endl;
  o << "set output sprintf(\"img/dm%.5d.png\", ii)" << endl;
  o << "set title sprintf(\"Density Matrix at t = %f fs\", 1.0*ii*"
    << p->tout << "/" << p->numOutputSteps << ")" << endl;
  o << "plot 'outs/dmt_z.out' index ii matrix with image" << endl;
  o << "}" << endl;
  o << endl;
  o << "!mencoder -really-quiet -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$opt -mf type=png:fps=25 -nosound -o /dev/null \"mf://img/*.png\"" << endl;
  o << "!mencoder -really-quiet -ovc lavc -lavcopts vcodec=msmpeg4v2:vpass=1:$opt -mf type=png:fps=25 -nosound -o figures/dmt_z.avi \"mf://img/*.png\"" << endl;
  o << "!rm -f divx2pass.log" << endl;
  o << "!rm -rf img" << endl;
  return;
}

/* Plots the populations in the bulk states over time as a movie */
void plotKProbsMovie(char * fileName, struct PARAMETERS * p) {
  ofstream o(fileName);

  o << "#!/bin/bash" << endl;
  o << endl;
  o << "if ! command -v mencoder &> /dev/null; then" << endl;
  o << "  echo \"ERROR [$0]: need mencoder to run this script\"" << endl;
  o << "  exit" << endl;
  o << "fi" << endl;
  o << "" << endl;
  o << "if ! command -v transpose &> /dev/null; then" << endl;
  o << "  echo \"ERROR [$0]: need 'transpose' program to run this script\"" << endl;
  o << "  exit" << endl;
  o << "fi" << endl;
  o << "" << endl;
  o << "if ! command -v parallel &> /dev/null; then" << endl;
  o << "  echo \"ERROR [$0]: need 'parallel' command to run this script\"" << endl;
  o << "  exit" << endl;
  o << "fi" << endl;
  o << "" << endl;
  o << "function plotNumber {" << endl;
  o << "gnuplot << EOF" << endl;
  o << "set terminal pngcairo font 'Arial-Bold,12'" << endl;
  o << "set output sprintf(\"img/%.5d.png\", ${1})" << endl;
  o << endl;
  o << "set palette defined(\\" << endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << endl;
  o << "1      0.7059  0.0157  0.1490\\" << endl;
  o << ")" << endl;
  o << endl;
  o << "set tics scale 0" << endl;
  o << "unset key" << endl;
  o << "idx=${1}+1" << endl;
  o << "set output sprintf(\"img/%.5d.png\", ${1})" << endl;
  o << "set multiplot layout 2,1" << endl;
  o << endl;
  o << "set xrange [0:" << p->kBandWidth << "]" << endl;
  o << "set yrange [0:*]" << endl;
  o << "set xlabel 'Energy above band edge (a.u.)'" << endl;
  o << "set ylabel 'Population in state'" << endl;
  o << "set title sprintf(\"Bulk Populations at t = %.6f a.u.\", ${1}*"
    << p->tout << "/" << p->numOutputSteps << ")" << endl;
  o << "plot 'outs/kprobs_t.out' u (\\$0*" << p->kBandWidth << "/" << (p->Nk - 1)
    << "):idx every ::1 with filledcurves x1" << endl;
  o << endl;
  o << "set xrange [0:" << p->tout << "]" << endl;
  o << "set yrange [0.00:" << p->kBandWidth << "]" << endl;
  o << "set xlabel 'Time (a.u.)'" << endl;
  o << "set ylabel 'Energy above band edge (a.u.)'" << endl;
  o << "set title 'State Populations vs. Energy and Time'" << endl;
  o << "plot '<cut -d \" \" -f 2- outs/kprobs.out' u (\\$2*" << p->tout << "/"
    << p->numOutputSteps << "):(\\$1*" << p->kBandWidth << "/" << (p->Nk - 1) << "):3 matrix with image" << endl;
  o << "unset multiplot" << endl;
  o << "EOF" << endl;
  o << "}" << endl;
  o << endl;
  o << "export -f plotNumber	# do this so parallel can find the plotNumber command" << endl;
  o << endl;
  o << "rm -rf img" << endl;
  o << "mkdir img" << endl;
  o << "transpose -o _t outs/kprobs.out" << endl;
  o << endl;
  o << "parallel --gnu plotNumber ::: {0.." << p->numOutputSteps << "}" << endl;
  o << endl;
  o << "mencoder -really-quiet \"mf://img/*png\" -mf type=png:fps=18 -ovc lavc -o kprobs_movie.avi" << endl;
  o << endl;
  o << "rm -rf img" << endl;
  return;
}

/* Plots the populations in the QD states over time as a movie */
void plotCProbsMovie(char * fileName, struct PARAMETERS * p) {
  ofstream o(fileName);

  int Ic = p->Ic;
  int Nc = p->Nc;

  o << "#!/bin/bash" << endl;
  o << endl;
  o << "if ! command -v mencoder &> /dev/null; then" << endl;
  o << "  echo \"ERROR [$0]: need mencoder to run this script\"" << endl;
  o << "  exit" << endl;
  o << "fi" << endl;
  o << "" << endl;
  o << "if ! command -v transpose &> /dev/null; then" << endl;
  o << "  echo \"ERROR [$0]: need 'transpose' program to run this script\"" << endl;
  o << "  exit" << endl;
  o << "fi" << endl;
  o << "" << endl;
  o << "if ! command -v parallel &> /dev/null; then" << endl;
  o << "  echo \"ERROR [$0]: need 'parallel' command to run this script\"" << endl;
  o << "  exit" << endl;
  o << "fi" << endl;
  o << "" << endl;
  o << "function plotNumber {" << endl;
  o << "gnuplot << EOF" << endl;
  o << "set terminal pngcairo font 'Arial-Bold,12'" << endl;
  o << "set output sprintf(\"img/%.5d.png\", ${1})" << endl;
  o << endl;
  o << "set palette defined(\\" << endl;
  o << "0.0    0.8667  0.8667  0.8667,\\" << endl;
  o << "0.0625 0.8980  0.8471  0.8196,\\" << endl;
  o << "0.125  0.9255  0.8275  0.7725,\\" << endl;
  o << "0.1875 0.9451  0.8000  0.7255,\\" << endl;
  o << "0.25   0.9608  0.7686  0.6784,\\" << endl;
  o << "0.3125 0.9686  0.7333  0.6275,\\" << endl;
  o << "0.375  0.9686  0.6941  0.5804,\\" << endl;
  o << "0.4375 0.9686  0.6510  0.5294,\\" << endl;
  o << "0.5    0.9569  0.6039  0.4824,\\" << endl;
  o << "0.5625 0.9451  0.5529  0.4353,\\" << endl;
  o << "0.625  0.9255  0.4980  0.3882,\\" << endl;
  o << "0.6875 0.8980  0.4392  0.3451,\\" << endl;
  o << "0.75   0.8706  0.3765  0.3020,\\" << endl;
  o << "0.8125 0.8353  0.3137  0.2588,\\" << endl;
  o << "0.875  0.7961  0.2431  0.2196,\\" << endl;
  o << "0.9375 0.7529  0.1569  0.1843,\\" << endl;
  o << "1      0.7059  0.0157  0.1490\\" << endl;
  o << ")" << endl;
  o << endl;
  o << "set tics scale 0" << endl;
  o << "unset key" << endl;
  o << "idx=${1}+1" << endl;
  o << "set output sprintf(\"img/%.5d.png\", ${1})" << endl;
  o << "set multiplot layout 2,1" << endl;
  o << endl;
  o << "set xrange [" << p->energies[Ic] << ":" << p->energies[Ic + Nc - 1] << "]" << endl;
  o << "set yrange [0:*]" << endl;
  o << "set xlabel 'Energy (a.u.)'" << endl;
  o << "set ylabel 'Population in state'" << endl;
  o << "set title sprintf(\"QD Populations at t = %.6f a.u.\", ${1}*"
    << p->tout << "/" << p->numOutputSteps << ")" << endl;
  o << "plot 'outs/cprobs_t.out' u (\\$0*" << (p->energies[Ic+Nc-1] - p->energies[Ic]) << "/" << (p->Nc - 1)
    << "):idx every ::1 with filledcurves x1" << endl;
  o << endl;
  o << "set xrange [0:" << p->tout << "]" << endl;
  o << "set yrange [0.00:" << (p->energies[Ic+Nc-1] - p->energies[Ic]) << "]" << endl;
  o << "set xlabel 'Time (a.u.)'" << endl;
  o << "set ylabel 'Energy above band edge (a.u.)'" << endl;
  o << "set title 'State Populations vs. Energy and Time'" << endl;
  o << "plot '<cut -d \" \" -f 2- outs/cprobs.out' u (\\$2*" << p->tout << "/"
    << p->numOutputSteps << "):(\\$1*" << (p->energies[Ic+Nc-1] - p->energies[Ic]) << "/" << (p->Nc - 1) << "):3 matrix with image" << endl;
  o << "unset multiplot" << endl;
  o << "EOF" << endl;
  o << "}" << endl;
  o << endl;
  o << "export -f plotNumber	# do this so parallel can find the plotNumber command" << endl;
  o << endl;
  o << "rm -rf img" << endl;
  o << "mkdir img" << endl;
  o << "transpose -o _t outs/cprobs.out" << endl;
  o << endl;
  o << "parallel --gnu plotNumber ::: {0.." << p->numOutputSteps << "}" << endl;
  o << endl;
  o << "mencoder -really-quiet \"mf://img/*png\" -mf type=png:fps=18 -ovc lavc -o cprobs_movie.avi" << endl;
  o << endl;
  o << "rm -rf img" << endl;
  return;
}
