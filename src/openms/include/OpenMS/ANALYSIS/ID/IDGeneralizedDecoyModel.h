// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_ID_IDGENERALIZEDDECOYMODEL_H
#define OPENMS_ANALYSIS_ID_IDGENERALIZEDDECOYMODEL_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <vector>


namespace OpenMS
{
  /**
    @brief Calculates an PEP from extended identifications

    One run containing extended search results can be used to annotate
    each of the peptide hits with a PEP.

    Also q-values can be reported instead.
    q-values are basically only adjusted p-values, also ranging from 0 to 1, with lower values being preferable.
    When looking at the list of hits ordered by q-values, then a hit with q-value of @em x means that there is an
    @em x*100 percent chance that all hits with a q-value <= @em x are a false positive hit.

    @todo implement KLR for traditional target and decoy analysis.
    @improvement report error properly for input with missing information

    @htmlinclude OpenMS_XXX.parameters

    @ingroup Utils
  */
  class OPENMS_DLLAPI IDGeneralizedDecoyModel :
    public DefaultParamHandler
  {
public:
    ///Default constructor
    IDGeneralizedDecoyModel();
    
    /**
      @brief Calculates the PEP of one run from an extended sequence db search

      @param id peptide identifications, containing target and (generalized) decoy hits
    */
    void apply(std::vector<PeptideIdentification>& id);

private:
    ///Not implemented
    IDGeneralizedDecoyModel(const IDGeneralizedDecoyModel&);

    ///Not implemented
    IDGeneralizedDecoyModel& operator=(const IDGeneralizedDecoyModel&);

    /// calculates the fdr stored into fdrs, given two vectors of scores
    void calculatePEPs_(Map<double, double>& score_to_fdr, std::vector<double>& target_scores, std::vector<double>& decoy_scores, bool q_value, bool higher_score_better, bool add_decoy_peptides, Size number_of_bins, double bandwidth);

    /// divide target and decoy scores into bins for model and decoy probability estimation, given two vectors of target and decoy scores
    void binData_(std::vector<double>& target, std::vector<double>& decoy, std::vector<double>& medians, std::vector<Int>& decoys, std::vector<Int>& binscores, Size number_of_bins);

    /// calculating PEPs stored into peps by using a mixture model, given also the bandwidth h of the kernel logistic regression and parameters alphas and beta from IRLS.
    void estimatePEPGeneralized_(std::vector<double>& peps, std::vector<double> target, std::vector<double> x, std::vector<Int> y, std::vector<Int> m, double h, bool higher_score_better);

    ///calculating Q-values from estimated PEPs
    void estimateQvaluesFromPEP_(std::vector<double>& q_values, std::vector<double>& peps);

    /// determining best bandwidth for kernel logistic regression.
    double crossValidation_(std::vector<double> x, std::vector<Int> y, std::vector<Int> m);
    /// calculating kernel function
    double kernel(double x1, double x2, double h);
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDGENERALIZEDDECOYMODEL_H
