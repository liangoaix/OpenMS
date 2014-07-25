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


#include <OpenMS/ANALYSIS/ID/IDGeneralizedDecoyModel.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <Eigen/Dense>
#include <numeric>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/StandardTypes.h>


#define ID_GENERALIZED_DECOY_MODEL_DEBUG
#undef  ID_GENERALIZED_DECOY_MODEL_DEBUG

double mymin(double a, double b)
{
  return a > b ? b : a;
}

using namespace std;
namespace OpenMS
{
  IDGeneralizedDecoyModel::IDGeneralizedDecoyModel() :
    DefaultParamHandler("IDGeneralizedDecoyModel")
  {
    defaults_.setValue("number_of_bins", 200, "Number of bins used for visualization.", ListUtils::create<String>("advanced"));
    defaults_.setValue("bandwidth_kernel", 0.0, "Bandwidth (double) of kernel for logistic regression procedure. Default 0 will implement cross-validation procedure to estimate best bandwidth", ListUtils::create<String>("advanced"));
    defaults_.setValue("output_name", "", "If output_plots is on, the output files will be saved in the following manner: <output_name>scores.txt for the scores and <output_name> which contains each step of the EM-algorithm e.g. output_name = /usr/home/OMSSA123 then /usr/home/OMSSA123_scores.txt, /usr/home/OMSSA123 will be written. If no directory is specified, e.g. instead of '/usr/home/OMSSA123' just OMSSA123, the files will be written into the working directory.", ListUtils::create<String>("advanced,output file"));
    defaults_.setValue("q_value", "false", "If 'true', q-values will be calculated based on the PEPs");
    defaults_.setValidStrings("q_value", ListUtils::create<String>("true,false"));
    defaults_.setValue("use_all_hits", "false", "If 'true' not only the first hit, but all are used (peptides only)");
    defaults_.setValidStrings("use_all_hits", ListUtils::create<String>("true,false"));
    defaults_.setValue("split_charge_variants", "false", "If set to 'true' charge variants are treated separately (for peptides of combined target/decoy searches only).");
    defaults_.setValidStrings("split_charge_variants", ListUtils::create<String>("true,false"));
    defaults_.setValue("treat_runs_separately", "false", "If set to 'true' different search runs are treated separately (for peptides of combined target/decoy searches only).");
    defaults_.setValidStrings("treat_runs_separately", ListUtils::create<String>("true,false"));
    defaults_.setValue("decoy_string", "_rev", "String which is appended at the accession of the protein to indicate that it is a decoy protein (for proteins only).");
    defaults_.setValue("add_decoy_peptides", "true", "If set to true, decoy peptides will be written to output file, too. The q-value is set to the closest target score.");
    defaults_.setValidStrings("add_decoy_peptides", ListUtils::create<String>("true,false"));
    defaults_.setValue("generalization_type", "inflated_MT", "Select different generalized methods");
    defaults_.setValidStrings("generalization_type", ListUtils::create<String>("inflated_MT,nonEnzymatic,PTM"));
    defaultsToParam_();
  }

  void IDGeneralizedDecoyModel::apply(vector<PeptideIdentification>& ids)
  {
    bool q_value = param_.getValue("q_value").toBool();
    bool use_all_hits = param_.getValue("use_all_hits").toBool();
    bool treat_runs_separately = param_.getValue("treat_runs_separately").toBool(); 
    bool split_charge_variants = param_.getValue("split_charge_variants").toBool();
    bool add_decoy_peptides = param_.getValue("add_decoy_peptides").toBool();
    Size number_of_bins = param_.getValue("number_of_bins");
    double bandwidth = param_.getValue("bandwidth_kernel");
    String generalization_type = param_.getValue("generalization_type"); 
#ifdef ID_GENERALIZED_DECOY_MODEL_DEBUG
    cerr << "Parameters: q_value=" << q_value << ", use_all_hits=" << use_all_hits << ", treat_runs_separately=" << treat_runs_separately << ", split_charge_variants=" << split_charge_variants << ", number_of_bins=" << number_of_bins << ", bandwidth" << bandwidth << ", generalization_type = " << generalization_type << endl;
#endif

    if (ids.empty())
    {
      LOG_WARN << "No peptide identifications given to IDGeneralizedDecoyModel! No calculation performed.\n";
      return;
    }

    // first search for all identifiers and charge variants
    set<String> identifiers;
    set<SignedSize> charge_variants;
    for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
    {
      identifiers.insert(it->getIdentifier());
      it->assignRanks();

      if (!use_all_hits)
      {
        vector<PeptideHit> hits = it->getHits();
        hits.resize(1);
        it->setHits(hits);
      }

      for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
      {
        charge_variants.insert(pit->getCharge());
      }
    }

#ifdef ID_GENERALIZED_DECOY_MODEL_DEBUG
    cerr << "#id-runs: " << identifiers.size() << " ";
    for (set<String>::const_iterator it = identifiers.begin(); it != identifiers.end(); ++it)
    {
      cerr << "," << *it;
    }
    cerr << endl;
    cerr << "#of charge states: " << charge_variants.size() << " ";
    for (set<SignedSize>::const_iterator it = charge_variants.begin(); it != charge_variants.end(); ++it)
    {
      cerr << "," << *it;
    }
    cerr << endl;
#endif

    for (set<SignedSize>::const_iterator zit = charge_variants.begin(); zit != charge_variants.end(); ++zit)
    {
#ifdef GENERALIZED_DECOY_MODEL_DEBUG
      cerr << "Charge variant=" << *zit << endl;
#endif
      // for all identifiers
      for (set<String>::const_iterator iit = identifiers.begin(); iit != identifiers.end(); ++iit)
      {
        if (!treat_runs_separately && iit != identifiers.begin())
        {
          continue;
        }

#ifdef GENERALIZED_DECOY_MODEL_DEBUG
        cerr << "Id-run: " << *iit << endl;
#endif
        // get the scores of all peptide hits
        vector<double> target_scores, decoy_scores;
        for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
        {
          // if runs should be treated separately, the identifiers must be the same
          if (treat_runs_separately && it->getIdentifier() != *iit)
          {
            continue;
          }

          for (Size i = 0; i < it->getHits().size(); ++i)
          {
            if (split_charge_variants && it->getHits()[i].getCharge() != *zit)
            {
              continue;
            }
            if (generalization_type == "nonEnzymatic")
            {
              if (it->getHits()[i].getAABefore() != 'R' && it->getHits()[i].getAABefore() != 'K')
              {
                decoy_scores.push_back(it->getHits()[i].getScore());
              }
              else
              {
                target_scores.push_back(it->getHits()[i].getScore());
              }
            }
            else if (generalization_type == "PTM")
            {
              String str = it->getHits()[i].getSequence().toString();
              if (str.find("Phospho") || str.find("Oxidation"))
              {
                decoy_scores.push_back(it->getHits()[i].getScore());
              }
              else
              {
                target_scores.push_back(it->getHits()[i].getScore());
              }
            }
            else
            {
              double cal_mass = (it->getHits()[i].getSequence().getMonoWeight() + (double)(it->getHits()[i].getCharge()) * Constants::PROTON_MASS_U) / (double)(it->getHits()[i].getCharge());
              double exp_mass = (double)(it->getMetaValue("MZ"));

              if (fabs(cal_mass - exp_mass) / exp_mass * 1000000 > 10)
              {
                decoy_scores.push_back(it->getHits()[i].getScore());
              }
              else
              {
                target_scores.push_back(it->getHits()[i].getScore());
              }
            }
          }
        }

#ifdef GENERALIZED_DECOY_MODEL_DEBUG
        cerr << "#target-scores=" << target_scores.size() << ", #decoy-scores=" << decoy_scores.size() << endl;
#endif

        // check decoy scores
        if (decoy_scores.empty())
        {
          String error_string = "IDGeneralizedDecoyModel: #decoy sequences is zero! Setting all target sequences to q-value/FDR 0! ";
          if (split_charge_variants || treat_runs_separately)
          {
            error_string += "(";
            if (split_charge_variants)
            {
              error_string += "charge_variant=" + String(*zit) + " ";
            }
            if (treat_runs_separately)
            {
              error_string += "run-id=" + *iit;
            }
            error_string += ")";
          }
          LOG_ERROR << error_string << std::endl;
        }

        // check target scores
        if (target_scores.empty())
        {
          String error_string = "IDGeneralizedDecoyModel: #target sequences is zero! Ignoring. ";
          if (split_charge_variants || treat_runs_separately)
          {
            error_string += "(";
            if (split_charge_variants)
            {
              error_string += "charge_variant=" + String(*zit) + " ";
            }
            if (treat_runs_separately)
            {
              error_string += "run-id=" + *iit;
            }
            error_string += ")";
          }
          LOG_ERROR << error_string << std::endl;
        }

        if (target_scores.empty() || decoy_scores.empty())
        {
          // no remove the the relevant entries, or put 'pseudo-scores' in
          for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
          {
            // if runs should be treated separately, the identifiers must be the same
            if (treat_runs_separately && it->getIdentifier() != *iit)
            {
              continue;
            }

            vector<PeptideHit> hits(it->getHits()), new_hits;
            for (Size i = 0; i < hits.size(); ++i)
            {
              if (split_charge_variants && hits[i].getCharge() != *zit)
              {
                new_hits.push_back(hits[i]);
                continue;
              }

              if (!hits[i].metaValueExists("target_decoy"))
              {
                LOG_FATAL_ERROR << "Meta value 'target_decoy' does not exists, reindex the idXML file with 'PeptideIndexer' first (run-id='" << it->getIdentifier() << ", rank=" << i + 1 << " of " << hits.size() << ")!" << endl;
                throw Exception::MissingInformation(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Meta value 'target_decoy' does not exist!");
              }

              String target_decoy(hits[i].getMetaValue("target_decoy"));
              if (target_decoy == "target")
              {
                // if it is a target hit, there are now decoys, fdr/q-value should be zero then
                new_hits.push_back(hits[i]);
                String score_type = it->getScoreType() + "_score";
                new_hits.back().setMetaValue(score_type, new_hits.back().getScore());
                new_hits.back().setScore(0);
              }
              else
              {
                if (target_decoy != "decoy" && target_decoy != "target+decoy")
                {
                  LOG_FATAL_ERROR << "Unknown value of meta value 'target_decoy': '" << target_decoy << "'!" << endl;
                }
              }
            }
            it->setHits(new_hits);
          }
          continue;
        }

        // calculate peps
        bool higher_score_better(ids.begin()->isHigherScoreBetter());
        Map<double, double> score_to_pep;
        calculatePEPs_(score_to_pep, target_scores, decoy_scores, q_value, higher_score_better, add_decoy_peptides, number_of_bins, bandwidth);
        String new_score_type;
        if (q_value)
        {
          new_score_type = "q-value";
        }
        else
        {
          new_score_type = "PEP";
        }

        // annotate peps
        for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
        {
          // if runs should be treated separately, the identifiers must be the same
          if (treat_runs_separately && it->getIdentifier() != *iit)
          {
            continue;
          }

          // String score_type = it->getScoreType() + "_score";
          vector<PeptideHit> hits;
          hits.metaRegistry().registerName(new_score_type, "calculated score of hit");
          for (vector<PeptideHit>::const_iterator pit = it->getHits().begin(); pit != it->getHits().end(); ++pit)
          {
            PeptideHit hit = *pit;

            if (split_charge_variants && pit->getCharge() != *zit)
            {
              hits.push_back(*pit);
              continue;
            }
            hit.setMetaValue(new_score_type, score_to_pep[pit->getScore()]);
            hit.setScore(pit->getScore());
            hits.push_back(hit);
          }
          it->setHits(hits);
        } // end: assign peps to peptide hits
        
      } // end loop of identifiers.
      if (!split_charge_variants)
      {
        break;
      }
    } // end loop of charge_variants.

    return;
  }

  void IDGeneralizedDecoyModel::calculatePEPs_(Map<double, double>& score_to_pep,
                                             vector<double>& target_scores, vector<double>& decoy_scores,
                                             bool q_value, bool higher_score_better, bool add_decoy_peptides, Size number_of_bins, double bandwidth)
  {
    vector<double> medians;
    vector<Int> decoysizes(number_of_bins);
    vector<Int> binscores(number_of_bins);
    binData_(target_scores, decoy_scores, medians, decoysizes, binscores, number_of_bins);

    double h = bandwidth;
    if (h == 0) { h = crossValidation_(medians, decoysizes, binscores); }
    if (h <= 0 || h >= 50) { h = 10; }
    vector<double> peps;
    estimatePEPGeneralized_(peps, target_scores, medians, decoysizes, binscores, h, higher_score_better);
    if (!q_value)
    {
      for (Size i = 0; i != target_scores.size(); ++i)
      {
        score_to_pep[target_scores[i]] = peps[i];
      }
    }
    else
    {
      vector<double> q_values;
      estimateQvaluesFromPEP_(q_values, peps);
      for (Size i = 0; i != target_scores.size(); ++i)
      {
        score_to_pep[target_scores[i]] = q_values[i];
      }
    }

    // assign q-value of decoy_score to closest target_score
    if (add_decoy_peptides)
    {
      estimatePEPGeneralized_(peps, decoy_scores, medians, decoysizes, binscores, h, higher_score_better);
      for (Size i = 0; i != decoy_scores.size(); ++i)
      {
        score_to_pep[decoy_scores[i]] = peps[i];
      }
    }
  }

  void IDGeneralizedDecoyModel::estimatePEPGeneralized_(vector<double>& peps, vector<double> target,
                                                      vector<double> x, std::vector<Int> y, std::vector<Int> m, double h, bool higher_score_better)
  {
    if (target[0] < target.back() || higher_score_better)
    {
      reverse(target.begin(), target.end()); //target need to be descending
    }
    Size n = x.size();
    double beta;
    Eigen::VectorXd g(n), ones(n), z(n), alphas(n);
    Eigen::MatrixXd matrix_k(n, n), identity(n, n);
    for (Size i = 0; i < n; i++)
    {
      ones(i) = 1;
      for (Size j = 0; j < n; j++)
      {
        matrix_k(i, j) = kernel(x[i], x[j], h);
      }
    }
    alphas.setZero();
    identity.setIdentity();
    Eigen::MatrixXd matrix_w = identity;
    double beta_old = 0.0;
    double step;
    Size iter = 0;
    double step_epsilon = 0.0001;
    do
    {
      g = matrix_k * alphas + beta_old * ones;
      double sigma, mu, p;
      for (Size ix = 0; ix < n; ix++)
      {
        assert(isfinite(g(ix)));
        p = 1 / (1 + exp(-g(ix)));
        mu = m[ix] * p;
        if (mu == 0) {sigma = 1.0; }else{sigma = mu * (1 - p); }
        z(ix) = g(ix) + (((double)y[ix]) - mu) / sigma;
        matrix_w(ix, ix) = 1 / sigma;
      }
      Eigen::MatrixXd matrix_m = matrix_k + matrix_w;
      Eigen::MatrixXd L = matrix_m.llt().matrixL();
      Eigen::VectorXd xi = (L * L.transpose()).colPivHouseholderQr().solve(ones); // (L*L.transpose()).inverse()*ones;
      Eigen::VectorXd zeta = (L * L.transpose()).colPivHouseholderQr().solve(z); //(L*L.transpose()).inverse()*z;
      beta = ones.dot(zeta) / ones.dot(xi);
      alphas = zeta - beta * xi;
      step = (beta_old - beta) * (beta_old - beta);
      beta_old = beta;
    }
    while ((step > step_epsilon || step < 0.0) && (++iter < 100));
// alphas and betas have been estimated...

    peps.clear();
    for (std::vector<double>::iterator it = target.begin(); it < target.end(); ++it)
    {
      double gx;
      Eigen::VectorXd ker(n);
      for (Size i = 0; i < n; i++)
      {
        ker(i) = kernel(*it, x[i], h);
      }
      gx = beta + alphas.dot(ker);
      peps.push_back(1 / (1 + exp(-gx)));
    }

    double p1 = *min_element(peps.begin(), peps.end()) - 1 / (double)(target.size());
    double p0 = *max_element(peps.begin(), peps.end());
    double top = 1;
    bool crap = false;
    vector<double>::iterator pep = peps.begin();

    for (; pep != peps.end(); ++pep)
    {
      if (crap)
      {
        *pep = top;
        continue;
      }
      *pep = (*pep - p1) / (p0 - p1);
      if (*pep >= top)
      {
        *pep = top;
        crap = true;
      }
    }
    partial_sum(peps.rbegin(), peps.rend(), peps.rbegin(), mymin);
  }

  double IDGeneralizedDecoyModel::kernel(double x1, double x2, double h)
  {
    double ker;
    double u = (x1 - x2) / h;
    ker = 1 / sqrt(2 * Constants::PI) * exp(-u * u / 2); //Guassian kernel. could be changed.
    return ker;
  }

  void IDGeneralizedDecoyModel::binData_(vector<double>& target, vector<double>& decoy, vector<double>& medians,
                                       vector<Int>& decoys, vector<Int>& binscores, Size number_of_bins)
  {
    sort(target.begin(), target.end());
    sort(decoy.begin(), decoy.end());

    double interval = (target.back() - target[0]) / number_of_bins;
    for (Size i = 0; i < number_of_bins - 1; ++i)
    {
      medians.push_back(target[0] + i * interval + interval / 2);
    }
    Size temp1 = 0; //number of target bin
    Size temp2 = 0; //number of decoy in this bin
    Int bin = 0;
    for (std::vector<double>::iterator it = target.begin(); it < target.end(); ++it)
    {
      double temp_divider = medians[bin] + interval / 2;
      if (*it < temp_divider && temp_divider - *it <= interval)
      {
        temp1 = (temp1 + 1);
      }
      else if (*it > temp_divider + interval)
      {
        binscores[bin] = temp1;
        ++bin;
        temp1 = 0;
        do
        {
          temp_divider += interval;
          binscores[bin] = temp1;
          ++bin;
        }
        while (*it > temp_divider + interval);
      }
      else if (*it >= temp_divider && *it <= temp_divider + interval)
      {
        binscores[bin] = temp1;
        ++bin;
        temp1 = 0;
        temp_divider += interval;
      }
    }
    bin = 0;
    for (std::vector<double>::iterator it = decoy.begin(); it < decoy.end(); ++it)
    {
//*it = *it + fabs(smallest_score) + 0.001;
      double temp_divider = medians[bin] + interval / 2;
      if (*it < temp_divider && temp_divider - *it <= interval)
      {
        temp2 = (temp2 + 1);
      }
      else if (*it > temp_divider + interval)
      {
        binscores[bin] += temp2;
        decoys[bin] = temp2;
        ++bin;
        temp2 = 0;
        do
        {
          temp_divider += interval;
          decoys[bin] = temp2;
          ++bin;
        }
        while (*it > temp_divider + interval);
      }
      else if (*it >= temp_divider && *it <= temp_divider + interval)
      {
        binscores[bin] += temp2;
        decoys[bin] = temp2;
        ++bin;
        temp2 = 0;
        temp_divider += interval;
      }
    }
  }

  void IDGeneralizedDecoyModel::estimateQvaluesFromPEP_(vector<double>& q_values, vector<double>& peps)
  {
    Int num = 0;
    double pep_avg = 0.0;
    if (peps.size())
    {
      q_values.clear();
    }
    for (vector<double>::const_iterator myP = peps.begin(); myP != peps.end(); ++myP, ++num)
    {
      pep_avg += *myP;
      q_values.push_back((pep_avg / (double)(num)));
    }
    partial_sum(q_values.rbegin(), q_values.rend(), q_values.rbegin(), mymin);
  }

  double IDGeneralizedDecoyModel::crossValidation_(vector<double> x, vector<Int> y, vector<Int> m)
  {
    vector<double> cross_validation;
    bool higher_score_better = true;
    vector<double> targets = x;
    vector<double> q_values, peps;
    for (Size hh = 1; hh < 50; ++hh)
    {
      double h = hh;
      estimatePEPGeneralized_(peps, targets, x, y, m, h, higher_score_better);
      estimateQvaluesFromPEP_(q_values, peps);
      Size num = 0;
      for (vector<double>::const_iterator q = q_values.begin(); q != q_values.end(); ++q)
      {
        if (*q < 0.03) {++num; }
      }
      double idrate = (double)(num) / (double)(q_values.size());
      cross_validation.push_back(idrate); //or calculate root mean square error.
/*
vector<double> Q_FDR;
for (Size i = 0; i != targets.size(); ++i)
{
Size closest_idx = 0;
for (Size j = 0; j != target_scores.size(); ++j)
{
if (fabs(targets[i] - target_scores[j]) < fabs(targets[i] - target_scores[closest_idx]))
{
closest_idx = j;
}
}
Q_FDR.push_back( score_to_fdr[target_scores[closest_idx]]);
}
double Deltarms= Math::meanSquareError( q_values.begin(), q_values.end(), Q_FDR.begin(), Q_FDR.end());
*/
    }
    double best_h = (double)(max_element(cross_validation.begin(), cross_validation.end()) - cross_validation.begin());
    return best_h;
  }

} // namespace OpenMS
