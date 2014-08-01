// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Test separate target and Generalized-decoy scores... 
#include <iostream>
#include <fstream>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <Eigen/Dense>
#include <numeric> 
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/ID/FalseDiscoveryRate.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

//#include <OpenMS/MATH/STATISTICS/PosteriorErrorProbabilityModel.h>
//#include <OpenMS/FORMAT/FileHandler.h>

using namespace Eigen;
using namespace OpenMS;
using namespace std;

double mymin(double a, double b);
double kernel(double x1, double x2,double h);

void Bindata(  vector<double> &target,  vector<double> &decoy,  vector<double> &medians,
				vector<Int> &decoys, vector<Int> &binscores, UInt number_of_bins);			
void ReadGeneralizedDecoy(vector<double> &target,vector<double> &decoy,vector<PeptideIdentification> pep_ids);	
void estimateQvaluesFromPEP(vector<double> &Qvalues,vector<double> &peps);

double estimateDeltaRMS(vector<double> Qvalues,vector<double> target, vector<double> decoy);



double mymin(double a, double b) {
  return a > b ? b : a;
}
double kernel(double x1, double x2, double h) {
        double ker;
        double u = (x1-x2)/h;
        ker = 1/sqrt(2* Constants::PI )*exp(-u*u/2);//Guassian kernel. could be changed.
        return ker;
}

void Bindata(  vector<double> &target,  vector<double> &decoy,  vector<double> &medians,
				vector<Int> &decoys, vector<Int> &binscores, UInt number_of_bins)
{
  sort(target.begin(), target.end());
  sort(decoy.begin(), decoy.end());
  double interval = (target.back() - target[0]) / number_of_bins;
  for(UInt i =0; i<number_of_bins-1; ++i){
	  medians.push_back(target[0]  + i*interval+ interval / 2) ; 
  }
  DPosition<2> temp;
  temp[0] = 0;//number of target bin
  temp[1] = 0;//number of decoy in this bin
  Int bin = 0;
  
  for (std::vector<double>::iterator it = target.begin(); it < target.end(); ++it)
  {
    double temp_divider = medians[bin]+interval/2;
    if (*it <temp_divider && temp_divider -*it <= interval)
    {
          temp[0] = (temp[0] + 1);
    }
    else if(*it >temp_divider + interval )
    {
        binscores[bin]= temp[0];
        ++bin;
        temp[0] = 0;
        do {
			temp_divider += interval;
            binscores[bin]= temp[0];
            ++bin;
            }while(*it > temp_divider+interval);
    }
    else if(*it >=temp_divider  && *it <=temp_divider + interval )
    {
		binscores[bin]= temp[0];
		++bin;
        temp[0] = 0;
        temp_divider += interval;
	 }
   }
   bin = 0;
   for (std::vector<double>::iterator it = decoy.begin(); it < decoy.end(); ++it)
   {
	//*it = *it + fabs(smallest_score) + 0.001;
	double temp_divider = medians[bin]+interval/2;
    if (*it <temp_divider && temp_divider -*it <= interval)
    {
          temp[0] = (temp[0] + 1);
    }
    else if(*it >temp_divider + interval )
    {
        binscores[bin]+= temp[0];
        decoys[bin] = temp[0];
        ++bin;
        temp[0] = 0;
        do {
			temp_divider += interval;
            decoys[bin] = temp[0];
            ++bin;
            }while(*it > temp_divider+interval);
     }
     else if(*it >=temp_divider  && *it <=temp_divider + interval )
     {
		binscores[bin]+= temp[0];
		decoys[bin] = temp[0];
		++bin;
        temp[0] = 0;
        temp_divider += interval;
	 }
   }
}


void ReadTargetDecoy(vector<double> &target,vector<double> &decoy,vector<PeptideIdentification> pep_ids)
{
  vector<PeptideIdentification>::const_iterator pid;
  for (pid = pep_ids.begin(); pid != pep_ids.end(); ++pid)
	{
		PeptideIdentification pit=*pid;
		for (vector<PeptideHit>::const_iterator hit_it = pit.getHits().begin();
         hit_it != pit.getHits().end(); ++hit_it)
		{			
			if(hit_it->getMetaValue("target-decoy") !="decoy" ) 
			{
				decoy.push_back(hit_it->getScore());
			}
			else 
			{
				target.push_back(hit_it->getScore());
			}			
		}
	}
}



void estimatePEPGeneralized(vector<double> &peps, vector<double> target, VectorXd &alphas, 
								double &beta ,vector<double> x,double &h,UInt N)
{
  if(target[0]<target.back()){
	reverse(target.begin(), target.end());  //target need to be descending
  }
  peps.clear();
  for (std::vector<double>::iterator it = target.begin(); it < target.end(); ++it)
  {  
        double gx;
        VectorXd ker(N);
        for(UInt i =0; i<N; i++){
                ker(i) =  kernel(*it,x[i],h);
        }
        gx = beta+alphas.dot(ker);
        peps.push_back(1/(1+exp(-gx)));
  }

  double p1=*min_element(peps.begin(), peps.end())-1/(double)(target.size());
  double p0=*max_element(peps.begin(), peps.end());
  double top = 1;
  bool crap = false;
  vector<double>::iterator pep = peps.begin();
  
  for (; pep != peps.end(); ++pep) 
  {
    if (crap) {
		*pep = top;
		continue;
    }
    *pep = (*pep-p1)/(p0-p1);
    if (*pep >= top) {
      *pep = top;
      crap = true;
    }
  }
  partial_sum(peps.rbegin(), peps.rend(), peps.rbegin(), mymin);
}

void estimateQvaluesFromPEP(vector<double> &Qvalues,vector<double> &peps)
{
	Int num = 0;
	double pep_avg = 0.0;
	if (peps.size()){
		Qvalues.clear();
	}
	for (vector<double>::const_iterator myP = peps.begin(); myP != peps.end(); ++myP, ++num) {
		pep_avg += *myP;
		Qvalues.push_back( (pep_avg / (double)(num)));
	}
	partial_sum(Qvalues.rbegin(), Qvalues.rend(), Qvalues.rbegin(), mymin);
}

double getIDrate(vector<double> Qvalues){ 
	UInt num = 0;
	for (vector<double>::const_iterator q = Qvalues.begin(); q != Qvalues.end(); ++q) {
		if(*q<0.05){++num;}
	}
	double idrate = (double)(num)/(double)(Qvalues.size());
	return idrate;
}

double estimateDeltaRMS(vector<double> Qvalues,vector<double> target, vector<double> decoy)
{
	vector<double>::iterator score = target.begin(); //ascending order...
	vector<double>::iterator scoreD = decoy.begin(); 
	vector<double> Q_FDR;
 	for (;score!=target.end(); ++score)
	{
		do{ ++scoreD;}while((*scoreD < *score) && (scoreD!=decoy.end()));	
		Q_FDR.push_back(  (double)(distance(score,target.end()) / distance(scoreD,decoy.end())) );
	}
	partial_sum(Q_FDR.begin(), Q_FDR.end(), Q_FDR.end(), mymin);
	vector<double>::const_iterator q1 = Qvalues.begin(); //ascending order...
	vector<double>::const_iterator q2 = Q_FDR.begin();
	double Deltarms= Math::meanSquareError( Qvalues.begin(), Qvalues.end(), Q_FDR.begin(), Q_FDR.end());
	return Deltarms;
}

int main(int argc, const char** argv)
{
  if (argc < 2) return 1;
  Param param;
  
  param.setValue("file_name", argv[1]);
  param.setValue("file_out", argv[2]);

  vector<ProteinIdentification> prot_ids;
  vector<PeptideIdentification> pep_ids;
  vector<FASTAFile::FASTAEntry> prot_tp;

  IdXMLFile().load((String)(param.getValue("file_name")), prot_ids, pep_ids);
  FASTAFile().load("/home/mi/liangoaix/storage/ground-truth/Sigma49/ms_ups2.fasta.TP", prot_tp);
  
  vector<String> protein_tp;
  for (Size i=0; i < prot_tp.size(); ++i)
  {
    protein_tp.push_back(prot_tp[i].description.substr(0,6));
  }

  map<String, double> fdr_ids;
  if (argc == 4)
  {
    vector<ProteinIdentification> prot_ids_fdr;
    vector<PeptideIdentification> pep_ids_fdr;
    
    IdXMLFile().load((String)(argv[3]), prot_ids_fdr, pep_ids_fdr);
    for (vector<PeptideIdentification>::iterator pid = pep_ids_fdr.begin(); pid != pep_ids_fdr.end(); ++pid)
    {
      for (vector<PeptideHit>::iterator hit_it = pid->getHits().begin(); hit_it != pid->getHits().end(); ++hit_it)
      {	
        PeptideHit hit = *hit_it;		
        String unique_id = String(pid->getMZ()) + "_" + String(pid->getRT()) +"_" + String(hit.getCharge())+"_" + hit.getAABefore();
        fdr_ids.insert( pair<String,double>( unique_id,hit.getScore()) );
      }
    }
  }
  
  //pep_ids[0].getReferencingHits(prot_ids[0].getHits()[0].getAccession(), hits); //
  //int num = 0;
  for (vector<PeptideIdentification>::iterator id = pep_ids.begin(); id != pep_ids.end(); ++id)
  {
    vector<PeptideHit> hits; 
    for (vector<PeptideHit>::iterator hit_it = id->getHits().begin(); hit_it != id->getHits().end(); ++hit_it)
    {	
      PeptideHit hit = *hit_it;		
      StringList protein_refs = hit.getProteinAccessions();
	  bool his_is_true_positive = false;
	  for (Size i = 0; i<protein_refs.size();++i)
	  {
        StringList sp = ListUtils::create<String>(protein_refs[i], '|');  
        if(std::find(protein_tp.begin(), protein_tp.end(), sp[1]) != protein_tp.end())
        {
	      his_is_true_positive = true;
	    }
	  }
	  if(his_is_true_positive)
	  {
        hit.setMetaValue("TP","true");
      }
      else
      {
        hit.setMetaValue("TP","false");
      }
      String unique_id = String(id->getMZ()) + "_" + String(id->getRT()) +"_" + String(hit.getCharge())+"_" + hit.getAABefore();
      hit.setMetaValue("UniqueID", unique_id);
      hits.push_back(hit);
      /*if(! his_is_true_positive)
      {
        if (hit.getMetaValue("target-decoy") == "target")
        { num++;}
	  }*/
    }
    id->setHits(hits);
    id->assignRanks();
  }

  //IdXMLFile().store((String)(param.getValue("file_out")), prot_ids, pep_ids);
  
//---input target and traditional decoy----estimate FDR-------------------------------------

//---------------------------------------------------------------
//--------------------------results-writing-----------------------------------

  ofstream f(argv[2]);
  f.precision(writtenDigits<double>(0.0));
  f << "ScanIndex\tUniqueID\tQ-Value\tScore\tE-value\tTarget-Decoy\tTruePositive\tSequence\tProteinLists\tFDR_qvalue\n";

  for (Size i = 0; i< pep_ids.size(); ++i )
  { 
	vector<PeptideHit> hits = pep_ids[i].getHits();
    for (vector<PeptideHit>::iterator hit_it = hits.begin(); hit_it != hits.end(); ++hit_it)
    {	
      PeptideHit hit = *hit_it;
      f << hit.getMetaValue("RT_index")<<"\t";
      f << hit.getMetaValue("UniqueID")<<"\t";
      f << hit.getScore() <<"\t"	;
      f << hit.getMetaValue("XTandem_score")<<"\t";
      f << hit.getMetaValue("E-Value") << "\t";
      f << hit.getMetaValue("target-decoy") <<"\t";
      f << hit.getMetaValue("TP")<<"\t";
      f << hit.getSequence() <<"\t" ;
      String proteins;
	  for (Size j = 0; j<hit.getProteinAccessions().size();++j)
	  {
        StringList sp = ListUtils::create<String>(hit.getProteinAccessions()[j], '|');  
        proteins += sp[1]+",";
	  }
      f << proteins << "\t";
      if(argc == 4)
      {
        if (fdr_ids[hit.getMetaValue("UniqueID")])
        {
		  f << fdr_ids[hit.getMetaValue("UniqueID")];
		} else
		{
	      f<< "null";
        }
      }  
      f<< "\n";
    }
  }
  
  f.close();

  return 0;  
} //end of main
