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

#if ! defined(SCORING_HPP)
#define SCORING_HPP

#include <cmath>
#include <cassert>
#include <ostream>
#include "acgtnMap.hpp"
#include "AutoArray.hpp"

struct Scoring;
std::ostream & operator<<(std::ostream & out, Scoring const & scoring);

struct Scoring
{
      public:
      /* Scoring parameters START */
      static const double DFLT_SIMILARITY;   /* default sequence similarity (0-1) */
      static const double DFLT_ERR;           /* default sequencing error rate per base [0-1) - Error-free = 0 */
      static const double DFLT_TRANS;         /* default transitions fraction of mutations (0-1) - Unbiased = 4/12 */
      static const double DFLT_GC;            /* default composition bias (0-1) - Unbiased = 0.5 */
      static const double DFLT_GCMUT_BIAS;       /* default mutability bias of G&C over A&T (>0) - Unbiased = 1*/
      /* Scoring parameters END */

      private:

      double oddsRatio[4][4];  /* 4x4 table of base substitution odds ratios */
      double bgFreq[4];        /* background frequency of each base */


      void init(double similarity, double gcContent, double transitRate, double errorRate, double gcMutBias);   /* initialization */

      static const double Q_PRB[65];                /* look-up table for quality to error probability conversion */

      AutoArray<double> LL;

      friend std::ostream & operator<<(std::ostream & out, Scoring const & scoring);

      public:

      /* constructors */
      Scoring();
      Scoring(double similarity, double gcContent, double transitRate, double errorRate, double gcMutBias);

      /* (1) calculates the score for two bases when base-call quality is available.
       * (2) calculates the score for a position when the raw base-probabilities are available.
       * Bases are assumed to be in numeric format: A=0, C=1, G=2, T=3
       */
      double getScore(char refBase, double aPrb, double gPrb, double tPrb, double cPrb);
      double getScore(char refBase, char readBase, int quality);
      double getScore(char refBase, char readBase);

      inline double getRawLogScoreTable(unsigned int refBase, unsigned int readBase, int quality) const
      {
            return LL [ (refBase<<(2+6)) | (readBase<<6) | quality ];
      }
};

std::ostream & operator<<(std::ostream & out, Scoring const & scoring);
#endif
