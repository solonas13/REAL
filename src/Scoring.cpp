/**
    REAL: An efficient REad ALigner for next generation sequencing reads.
    Copyright (C) 2010 German Tischler, Solon Pissis, Kimon Frousios

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#if defined(HAVE_CONFIG_H)
#include "real_config.hpp"
#endif

#include "Scoring.hpp"
#include <limits>
#include <iostream>

/* Intitialize constant Q_PRB[i] = pow(.1, i/10.0) */
const double Scoring::Q_PRB[65] = {
   1.0000000, 0.7943282, 0.6309573, 0.5011872, 0.3981072, 0.3162278, 0.2511886, 0.1995262, 0.1584893, 0.1258925,
   0.1000000, 0.0794328, 0.0630957, 0.0501187, 0.0398107, 0.0316228, 0.0251189, 0.0199526, 0.0158489, 0.0125893,
   0.0100000, 0.0079433, 0.0063096, 0.0050119, 0.0039811, 0.0031623, 0.0025119, 0.0019953, 0.0015849, 0.0012589,
   0.0010000, 0.0007943, 0.0006310, 0.0005012, 0.0003981, 0.0003162, 0.0002512, 0.0001995, 0.0001585, 0.0001259,
   0.0001000, 0.0000794, 0.0000631, 0.0000501, 0.0000398, 0.0000316, 0.0000251, 0.0000200, 0.0000158, 0.0000126,
   0.0000100, 0.0000079, 0.0000063, 0.0000050, 0.0000040, 0.0000032, 0.0000025, 0.0000020, 0.0000016, 0.0000013,
   0.0000010, 0.0000008, 0.0000006, 0.0000005, 0.0000004
};

/* 
 * Empty Constructor
 * Builds an odds-ratios table based on the default values.
 */
Scoring::Scoring()
{
   init( DFLT_SIMILARITY, DFLT_GC, DFLT_TRANS, DFLT_ERR, DFLT_GCMUT_BIAS );
}


/*
 * Constructor
 * Builds an odds-ratios matrix as dictated by the supplied parameters
 */
Scoring::Scoring(double similarity, double gcContent, double transitRate, double errorRate, double gcMutBias)
{
   init( similarity, gcContent, transitRate, errorRate, gcMutBias );
}


/*
 * Initialize class
 */
void Scoring::init(double similarity, double gcContent, double transitRate, double errorRate, double gcMutBias)
{
   double psum = 0;
   double transit = transitRate * (1 - similarity);
   double transver = (1 - transitRate) * (1 - similarity);
   
   this->bgFreq[0] = (1 - gcContent) / 2;   /* background frequencies of A and T */
   this->bgFreq[3] = (1 - gcContent) / 2;
   this->bgFreq[1] = gcContent / 2;         /* background frequencies of G and C */
   this->bgFreq[2] = gcContent / 2;
   
   /* Adjust gcMutBias taking the gcContent into account */
   gcMutBias = gcMutBias * (1-gcContent) / gcContent;

   /* Transition probabilities */
   /*AG*/ this->oddsRatio[0][2] = transit / (gcMutBias + 1) / (1-gcContent);
   /*TC*/ this->oddsRatio[3][1] = transit / (gcMutBias + 1) / (1-gcContent);
   /*GA*/ this->oddsRatio[2][0] = transit / (gcMutBias + 1) / gcContent * gcMutBias;
   /*CT*/ this->oddsRatio[1][3] = transit / (gcMutBias + 1) / gcContent * gcMutBias;
   /* Transversion probabilities */
   /*AC*/ this->oddsRatio[0][1] = transver / 2 / (gcMutBias + 1) / (1-gcContent);
   /*TG*/ this->oddsRatio[3][2] = transver / 2 / (gcMutBias + 1) / (1-gcContent);
   /*AT*/ this->oddsRatio[0][3] = transver / 2 / (gcMutBias + 1) / (1-gcContent);
   /*TA*/ this->oddsRatio[3][0] = transver / 2 / (gcMutBias + 1) / (1-gcContent);
   /*CA*/ this->oddsRatio[1][0] = transver / 2 / (gcMutBias + 1) / gcContent * gcMutBias;
   /*GT*/ this->oddsRatio[2][3] = transver / 2 / (gcMutBias + 1) / gcContent * gcMutBias;
   /*CG*/ this->oddsRatio[1][2] = transver / 2 / (gcMutBias + 1) / gcContent * gcMutBias;
   /*GC*/ this->oddsRatio[2][1] = transver / 2 / (gcMutBias + 1) / gcContent * gcMutBias;
   /* Conservation probabilities */
   /*AA*/ this->oddsRatio[0][0] = 1 - this->oddsRatio[0][1] - this->oddsRatio[0][2] - this->oddsRatio[0][3];
   /*TT*/ this->oddsRatio[3][3] = 1 - this->oddsRatio[3][0] - this->oddsRatio[3][1] - this->oddsRatio[3][2];
   /*GG*/ this->oddsRatio[2][2] = 1 - this->oddsRatio[2][0] - this->oddsRatio[2][1] - this->oddsRatio[2][3];
   /*CC*/ this->oddsRatio[1][1] = 1 - this->oddsRatio[1][0] - this->oddsRatio[1][2] - this->oddsRatio[1][3];
      
   /* Odds Ratios */
   for (int x=0; x<4; x++)
   {
      for (int y=0; y<4; y++)
      {
         this->oddsRatio[x][y] *= 1 - errorRate;
#if defined(DEBUG_SCORES)
         std::cerr << oddsRatio[x][y] << "   ";
#endif
         psum += oddsRatio[x][y];
         
         this->oddsRatio[x][y] /= this->bgFreq[y];
      }
#if defined(DEBUG_SCORES)
      std::cerr << std::endl;
#endif
   }
#if defined(DEBUG_SCORES)
   std::cerr << psum << std::endl << std::endl;
#endif

#if defined(DEBUG_SCORES)   
   for (int x=0; x<4; x++)
   {
      for (int y=0; y<4; y++)
      {
         std::cerr << oddsRatio[x][y] << "   ";
      }
      std::cerr << std::endl;
   }
   std::cerr << std::endl;
#endif

   LL = AutoArray<double>(4*4*64);
   for ( unsigned int c0 = 0; c0 < 4; ++c0 )
      for ( unsigned int c1 = 0; c1 < 4; ++c1 )
         for ( unsigned int q = 0; q < 64; ++q )
            LL [ (c0<<(2+6)) | (c1<<6) | q ] = getScore(c0,c1,q);
}



/*
 * Calculate score of two aligned bases if probabilities for each base call are known
 */
double Scoring::getScore(char refBase, double aPrb, double cPrb, double gPrb, double tPrb)
{
   double prb[4] = {aPrb,cPrb,gPrb,tPrb};
   double xyR = 0;
   for (int y=0; y<4; y++)
   {
      xyR += this->oddsRatio[static_cast<int>(refBase)][y] * prb[y] ;
   }
   return std::log(xyR)/std::log(2.0);
}


/*
 * Calculate score of two aligned bases, if only the probability of the optimal base-call is known (fastQ quality)
 */
double Scoring::getScore(char refBase, char readBase, int quality)
{
//    double xyR = 0;
//    for (int y=0; y<4; y++)
//    {
//       double prb;
//       if ( readBase == y )             /* optimal base call, weighted by the confidence in the call */
//          prb = 1 - Q_PRB[quality] ;
//       else                             /* suboptimal base calls, 
//                                         * the call-error probability is equally distributed among the 3 suboptimal calls */
//          prb =  Q_PRB[quality] / 3;
      
//       xyR += this->oddsRatio[static_cast<int>(refBase)][y] * prb;
//    }
//    return std::log(xyR)/std::log(2.0);
   return std::log( this->oddsRatio[static_cast<int>(refBase)][static_cast<int>(readBase)] )/std::log(2.0) * ( 1 - Q_PRB[quality] );
}


/*
* Calculate score of two aligned bases, using ONLY the substitution matrix, without any quality data (fastA - basically assume 100% certainty of base-call).
 */
double Scoring::getScore(char refBase, char readBase)
{
   return std::log( oddsRatio[static_cast<int>(refBase)][static_cast<int>(readBase)] ) / std::log(2.0);
}

std::ostream & operator<<(std::ostream & out, Scoring const & scoring)
{
   out << "Scoring(\n";
   out << " oddsRatio(\n";
   for ( unsigned int i = 0; i < 4; ++i )
   {
      out << "  ";
      for ( unsigned int j = 0; j < 4; ++j )
         out << scoring.oddsRatio[i][j] << (((j+1)<4)?"\t":"");
      out << "\n";
   }
   out << " )\n";
   out << " RawLogScoreTable(\n";
   for ( unsigned int i = 0; i < 4; ++i )
      for ( unsigned int j = 0; j < 4; ++j )
         for ( unsigned int q = 0; q < 63; ++q )
            out << " refBase=" << toollib::remapChar(i) << " readBase=" << toollib::remapChar(j) << " quality=" << q << " entry=" << scoring.getRawLogScoreTable(i,j,q) << "\n";
   out << " )\n";
   out << ")\n";
   return out;
}

const double Scoring::DFLT_SIMILARITY = 0.995;   /* default sequence similarity (0-1) */
const double Scoring::DFLT_ERR = 0.00;           /* default sequencing error rate per base [0-1) - Error-free = 0 */
const double Scoring::DFLT_TRANS = 0.71;         /* default transitions fraction of mutations (0-1) - Unbiased = 4/12 */
const double Scoring::DFLT_GC = 0.41;            /* default composition bias (0-1) - Unbiased = 0.5 */
const double Scoring::DFLT_GCMUT_BIAS = 2;       /* default mutability bias of G&C over A&T (>0) - Unbiased = GC/(1-GC) */

