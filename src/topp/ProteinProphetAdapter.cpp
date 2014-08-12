// --------------------------------------------------------------------------
// OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
// * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
// * Neither the name of any author or any participating institution
// may be used to endorse or promote products derived from this software
// without specific prior written permission.
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
#include <OpenMS/FORMAT/PepXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/ProtXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <fstream>

#include <QtCore/QFile>
#include <QtCore/QProcess>
#include <QDir>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------
/**
  @page TOPP_ProteinProphetAdapter ProteinProphetAdapter
  @brief Computes a protein identification based on the number of identified peptides.
<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ ProteinProphetAdapter \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter (or other ID engines)</td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=3> @ref TOPP_PeptideIndexer </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter </td>
    </tr>
  </table>
</CENTER>
  @experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!
  
  TPP PeptideProphet and ProteinProphet should be installed properly before running this tool.
  
  Some information about the supported input types:
  @ref OpenMS::PepXMLFile "pepXML"
  @ref OpenMS::IdXMLFile "idXML"
  
  <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_ProteinProphetAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_ProteinProphetAdapter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPProteinProphetAdapter :
  public TOPPBase
{
public:
  TOPPProteinProphetAdapter() :
    TOPPBase("ProteinProphetAdapter", "TPP adapter, a protein inference tool.", false)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file");
    setValidFormats_("in", ListUtils::create<String>("idXML,pepXML"));
    registerOutputFile_("out", "<file>", "", "Output file");
    setValidFormats_("out", ListUtils::create<String>("idXML,pepXML,protXML"));
    addEmptyLine_();
    registerStringOption_("out_type", "<type>", "idXML", "Output format", false);
    setValidStrings_("out_type", ListUtils::create<String>("idXML,pepXML,protXML"));
    registerInputFile_("mz_file", "<file>", "", "Experiment data file", false);
    setValidFormats_("mz_file", ListUtils::create<String>("mzML"));
    registerInputFile_("database", "<file>", "", "FASTA file or pro file. Non-existing relative file-names are looked up via'OpenMS.ini:id_db_dir'", false);
    registerStringOption_("database_type", "<database type>", "AA", "Specify 'AA' for amino acid, 'NA' for nucleic acid (default 'AA')", false);
    setValidStrings_("database_type", ListUtils::create<String>("AA,NA"));
    addEmptyLine_();
    registerFlag_("proteinprophet_off", "Only PeptideProphet will run and ProteinProphet will be disabled; output can be pepXML or idXML containing peptide prophet probability.");
    registerInputFile_("tpp_executable", "<executable>", "/tpp/bin/xinteract",
                       "TPP bin directory e.g. '/usr/local/tpp/bin/xinteract'", true, false, ListUtils::create<String>("skipexists"));
    registerInputFile_("default_input_file", "<file>", "", "Default parameters input file, if not given default parameters are used", false);
    
    registerFlag_("refinement", "Discard search results with PeptideProphet probabilities below 0.05.");
    registerIntOption_("num_extra_interation", "<num>", 20, "Number of extra PeptideProphet interations; default <num>=20", false);
    registerIntOption_("ignore_charge", "<num>", -1, "Ignore charge <num>+", false);
    registerStringOption_("decoy_prefix", "<tag>", "", "Use decoy hits to pin down the negative distribution; the decoy protein names must begin with <tag> (whitespace is not allowed). e.g. decoy", false);
    registerIntOption_("conservative_level", "<num>", 0, "Specify how conservative the model is to be in number of standard deviations from negative mean to allow positive model to cover, higher is more conservative.", false);
    registerStringOption_("precursor_error_units", "<unit>", "dalton", "Specify the precusor error unit for the accurate mass model.", false);
    setValidStrings_("precursor_error_units", ListUtils::create<String>("dalton,ppm"));
    registerStringOption_("experiment_label", "<tag>", "", "Used to commonly label all spectra belonging to one experiment (required by iProphet). ", false);
    registerIntOption_("minimun_pep_length", "<num>", 7, "Minimum peptide length considered in the analysis (default 7).", false);
    registerDoubleOption_("filter_result", "<float>", 0.05, "Filter results below PeptideProphet probability <value> (default 0.05).", false);
    // -mw [calculate protein molecular weights]" ;
    registerStringOption_("fragment_type", "<unit>", "MONO", "Calculate monoisotopic/average peptide masses during conversion to pepXML.", false);
    setValidStrings_("fragment_type", ListUtils::create<String>("MONO,AVE"));
    registerStringOption_("sample_enzyme", "<enzyme>", "eT", "Specify sample enzyme: -eT = Trypsin, -eS = StrictTrypsin, -eC = Chymotrypsin, -eR = RalphTrypsin, -eA = AspN, -eD = Trypsin/CNBr, -eE = Elastase, -eN = Nonspecific or None, -eG, -eB, -eM, -e3, -eK, -eL, -eP, details see the TPP documentation.", false);
    setValidStrings_("sample_enzyme", ListUtils::create<String>("eT,eS,eC,eR,eA,eD,eE,eN,eG,eB,eM,e3,eK,eL,eP,"));
    
    // (generaloptions) (-Oprophetoptions) (-Xxpressoptions) (-Aasapoptions)
    addEmptyLine_();
    registerStringOption_("iprophet_option", "<option>", "", "iProphet options, run iProphet on the PeptideProphet result. Starting with 'i', details see the TPP documentation. ", false);
    registerStringOption_("ptmprophet_option", "<option>", "", "PTMProphet options, details see the TPP documentation. ", false);
    registerStringOption_("peptideprophet_option", "<option>", "", "PeptideProphet options, details see the TPP documentation. ", false);
    registerStringOption_("xpress_option", "<option>", "", "Run XPRESS analysis with any specified options that follow the 'X', e.g. 'X-nC'. Details see the TPP documentation. ", false);
    registerStringOption_("asapratio_option", "<option>", "", "Run ASAPRatio analysis with any specified options that follow the 'A', e.g. 'A-lDE-S'. details see the TPP documentation. ", false);
    registerStringOption_("proteinprophet_option", "<option>", "", "ProteinProphet options. details see the TPP documentation. ", false);
    registerFlag_("refreshparser_off", "RefreshParser will be disabled (by -nR in TPP).");
  }

  ExitCodes main_(int, const char**)
  {
// path to the log file
    String logfile(getStringOption_("log"));
    StringList parameters;
    String inputfile_name;
    String outputfile_name;
    inputfile_name = getStringOption_("in");
    writeDebug_(String("Input file: ") + inputfile_name, 1);
    if (inputfile_name == "")
    {
      writeLog_("No input file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    outputfile_name = getStringOption_("out");
    writeDebug_(String("Output file: ") + outputfile_name, 1);
    if (outputfile_name == "")
    {
      writeLog_("No output file specified. Aborting!");
      printUsage_();
      return ILLEGAL_PARAMETERS;
    }
    
    String temp_directory = QDir::toNativeSeparators((File::getTempDirectory() + "/" + File::getUniqueName() + "/").toQString()); // body for the tmp files
    {
      QDir d;
      d.mkpath(temp_directory.toQString());
    }
    String xinteract_output_filename(temp_directory + "xinteract_output_file.pep.xml");
    String xinteract_input_filename;
    // checking if input type is pepxml
    if (FileHandler().getType(inputfile_name) == FileTypes::PEPXML)
    {
      xinteract_input_filename = inputfile_name;
    }
    else
    {
      xinteract_input_filename = temp_directory + "tpp_input_file.pep.xml";
      // reading idXML input
      vector<ProteinIdentification> protein_ids_in;
      vector<PeptideIdentification> peptide_ids_in;
      IdXMLFile().load(inputfile_name, protein_ids_in, peptide_ids_in);
      String mz_file = getStringOption_("mz_file");
      String mz_dir = "";
      String base_name = "";
      if (!mz_file.empty())
      {
#ifdef OPENMS_WINDOWSPLATFORM
        StringList sp = ListUtils::create<String>(mz_file, '\\'); 
        for (Size i = 0; i != sp.size() - 1; ++i)
        {
          mz_dir += sp[i] + "\\";
        }
        if (mz_dir.hasSubstring(" ")) // contain space
        {
          parameters << "-a" + String("\"\\\"") + String(mz_dir) + String("\\\"\"");
        }
        else
        {
          parameters << "-a" + mz_dir; // xinteract only needs the path to the mz file.
        }
        base_name = sp.back(); 
#else
        StringList sp = ListUtils::create<String>(mz_file, '/'); 
        for (Size i = 0; i != sp.size() - 1; ++i)
        {
          mz_dir += sp[i] + "/";
        }
        if (mz_dir.hasSubstring(" "))
        {
          parameters << "-a" + String("\"") + String(mz_dir) + String("\"");
        }
        else
        {
          parameters << "-a" + mz_dir; // xinteract only needs the path to the mz file.
        }
        base_name = sp.back();
#endif
      }
      PepXMLFile().store(xinteract_input_filename, protein_ids_in, peptide_ids_in, mz_file, base_name, false);
    }
    
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    parameters << "-N" + xinteract_output_filename; // store the output temp file
    // checking database
    String db_name = getStringOption_("database");
    if (!db_name.empty())
    {
      bool db_name_contains_space = false;
      if (db_name.hasSubstring(" "))
      {
        db_name_contains_space = true;
      }
      if (db_name_contains_space)
      {
#ifdef OPENMS_WINDOWSPLATFORM
// Windows: use doubly escaped double quotes (and do a system call instead of QProcess later)
        parameters << "-D" + String("\"\\\"") + String(db_name) + String("\\\"\"");
#else
        parameters << "-D" + String("\"") + String(db_name) + String("\"");
#endif
      }
      else
      {
        parameters << "-D" + String(db_name);
      }
    }
    parameters << "-T" + getStringOption_("database_type");
    if (!getFlag_("refinement"))
    {
      parameters << "-p0";
    }
    parameters << "-x" + String(getIntOption_("num_extra_interation"));
    if (getIntOption_("ignore_charge") >= 0)
    {
      parameters << "-I" + String(getIntOption_("ignore_charge"));
    }
    if (!getStringOption_("decoy_prefix").empty())
    {
      parameters << "-d" + getStringOption_("decoy_prefix");
    }
    parameters << "-c" + String(getIntOption_("conservative_level"));
    if (getStringOption_("precursor_error_units") == "ppm")
    {
      parameters << "-PPM";
    }
    if (!getStringOption_("experiment_label").empty())
    {
      parameters << "-E" + getStringOption_("experiment_label");
    }
    parameters << "-l" + String(getIntOption_("minimun_pep_length"));
    parameters << "-p" + String(getDoubleOption_("filter_result"));
    if (!getStringOption_("fragment_type").empty())
    {
      parameters << "-" + getStringOption_("fragment_type");
    }
    if (!getStringOption_("sample_enzyme").empty())
    {
      parameters << "-" + getStringOption_("sample_enzyme");
    }
    if (!getStringOption_("iprophet_option").empty())
    {
      parameters << "-" + getStringOption_("iprophet_option");
    }
    if (!getStringOption_("ptmprophet_option").empty())
    {
      parameters << "-" + getStringOption_("ptmprophet_option");
    }
    if (!getStringOption_("peptideprophet_option").empty())
    {
      parameters << "-" + getStringOption_("peptideprophet_option");
    }
    if (!getStringOption_("xpress_option").empty())
    {
      parameters << "-" + getStringOption_("xpress_option");
    }
    if (!getStringOption_("asapratio_option").empty())
    {
      parameters << "-" + getStringOption_("asapratio_option");
    }
    if (getFlag_("refreshparser_off"))
    {
      parameters << "-nR";
    }
    //-------------------------------------------------------------
    // run PeptideProphet: xinteract
    //-------------------------------------------------------------
    // TODO: TPP can support more than one peptideID input.
    ProgressLogger pl;
    pl.setLogType(ProgressLogger::CMD);
    pl.startProgress(0, 1, "running xinteract...");
    
    String xinteract_executable(getStringOption_("tpp_executable"));
    QStringList qparam;
    qparam << xinteract_input_filename.toQString();
    for (Size ip = 0; ip < parameters.size(); ++ip)
    {
      qparam << parameters[ip].toQString();
    }
    Int status = 0;
#ifdef OPENMS_WINDOWSPLATFORM
    if (db_name_contains_space)
    {
      String call_string = xinteract_executable + " " + ListUtils::concatenate(parameters, " ");
      writeDebug_(call_string, 5);
      status = system(call_string.c_str());
    }
    else
    {
      writeDebug_(xinteract_executable + " " + ListUtils::concatenate(parameters, " "), 5);
      status = QProcess::execute(xinteract_executable.toQString(), qparam);
    }
#else
    writeDebug_(xinteract_executable + " " + ListUtils::concatenate(parameters, " "), 5);
    status = QProcess::execute(xinteract_executable.toQString(), qparam);
#endif
    if (status != 0)
    {
      writeLog_("TPP problem. Aborting! Calling command was: '" + xinteract_executable + " \"" + xinteract_input_filename + "\"'.\nDoes the TPP executable exist?");
// clean temporary files
      if (this->debug_level_ < 2)
      {
        File::removeDirRecursively(temp_directory);
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory << "'" << std::endl;
      }
      else
      {
        LOG_WARN << "Keeping the temporary files at '" << temp_directory << "'. Set debug level to <2 to remove them." << std::endl;
      }
      return EXTERNAL_PROGRAM_ERROR;
    }
    StringList files;
    File::fileList(temp_directory, "xinteract_output_file*.pep.xml", files);
    if (files.size() != 1)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, xinteract_output_filename);
    }
    pl.endProgress();
    //-------------------------------------------------------------
    // run ProteinProphet afterwards
    //-------------------------------------------------------------
    // "usage:\tProteinProphet <interact_pepxml_file1> [<interact_pepxml_file2>[....]] <output_protxml_file> (ICAT) (GLYC) (XPRESS) (ASAP_PROPHET) (FPKM) (NONSP) (ACCURACY) (ASAP) (PROTLEN) (NOPROTLEN) (IPROPHET) (NORMPROTLEN) (GROUPWTS) (INSTANCES) (REFRESH) (DELUDE) (NOOCCAM) (NOPLOT) (PROTMW)
    if (!getFlag_("proteinprophet_off"))
    {
      pl.startProgress(0, 1, "Running ProteinProphet...");
#ifdef OPENMS_WINDOWSPLATFORM
      StringList exe_path = ListUtils::create<String>(xinteract_executable, '\\'); // TODO: different platform
      String pp_executable;
      for (Size i = 0; i != exe_path.size() - 1; ++i)
      {
        pp_executable += exe_path[i] + '\\';
      }
#else
      StringList exe_path = ListUtils::create<String>(xinteract_executable, '/'); // TODO: different platform
      String pp_executable;
      for (Size i = 0; i != exe_path.size() - 1; ++i)
      {
        pp_executable += exe_path[i] + '/';
      }
#endif
      pp_executable += "ProteinProphet";
      String proteinprophet_input_filename(xinteract_output_filename); // TODO: more than one interact_pepxml files as input.
      String proteinprophet_output_filename(temp_directory + "proteinprophet_output_file.prot.xml");
      QStringList qparam_pp;
      qparam_pp << proteinprophet_input_filename.toQString();
      qparam_pp << proteinprophet_output_filename.toQString();
      if (!getStringOption_("proteinprophet_option").empty())
      {
        qparam_pp << String(getStringOption_("proteinprophet_option")).toQString();
      }
      Int new_status = 0;
      writeDebug_(pp_executable + " ", 5);
      status = QProcess::execute(pp_executable.toQString(), qparam_pp);
      if (new_status != 0)
      {
        writeLog_("ProteinProphet problem. Aborting! Calling command was: '" + pp_executable + " \"" + proteinprophet_input_filename + "\"'.\nIs ProteinProphet executable in the same folder as xinteract.");
        if (this->debug_level_ < 2)
        {
          File::removeDirRecursively(temp_directory);
          LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory << "'" << std::endl;
        }
        else
        {
          LOG_WARN << "Keeping the temporary files at '" << temp_directory << "'. Set debug level to <2 to remove them." << std::endl;
        }
        return EXTERNAL_PROGRAM_ERROR;
      }
      StringList files_pp;
      File::fileList(temp_directory, "proteinprophet_output_file*.prot.xml", files_pp);
      if (files_pp.size() != 1)
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, proteinprophet_output_filename);
      }
      if (getStringOption_("out_type") == "idXML")
      {
        vector<ProteinIdentification> protein_ids;
        vector<PeptideIdentification> peptide_ids;
        protein_ids.resize(1);
        peptide_ids.resize(1);
        ProtXMLFile().load(proteinprophet_output_filename, protein_ids[0], peptide_ids[0]);
        protein_ids[0].setSearchEngine("XTandem"); // TODO: add other search engines.
        IdXMLFile().store(outputfile_name, protein_ids, peptide_ids);
      }
      else
      {
        QFile::copy(proteinprophet_output_filename.toQString(), outputfile_name.toQString());
      }
      if (this->debug_level_ < 2)
      {
        File::removeDirRecursively(temp_directory);
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory << "'" << std::endl;
      }
      else
      {
        LOG_WARN << "Keeping the temporary files at '" << temp_directory << "'. Set debug level to <2 to remove them." << std::endl;
      }
      pl.endProgress();
    }
    else
    {
      //-------------------------------------------------------------
      // only storing PeptideProphet result.
      //-------------------------------------------------------------
      if (getStringOption_("out_type") == "idXML")
      {
        vector<ProteinIdentification> protein_ids;
        ProteinIdentification protein_id;
        vector<PeptideIdentification> peptide_ids;
        PepXMLFile().load(xinteract_output_filename, protein_ids, peptide_ids);
        protein_id.setSearchEngine("XTandem"); // TODO: add other search engines.
        protein_ids.push_back(protein_id);
        IdXMLFile().store(outputfile_name, protein_ids, peptide_ids);
      }
      else
      {
        QFile::copy(xinteract_output_filename.toQString(), outputfile_name.toQString());
      }
      // Deletion of temporary files
      if (this->debug_level_ < 2)
      {
        File::removeDirRecursively(temp_directory);
        LOG_WARN << "Set debug level to >=2 to keep the temporary files at '" << temp_directory << "'" << std::endl;
      }
      else
      {
        LOG_WARN << "Keeping the temporary files at '" << temp_directory << "'. Set debug level to <2 to remove them." << std::endl;
      }
    }
    return EXECUTION_OK;
  }

};
int main(int argc, const char** argv)
{
  TOPPProteinProphetAdapter tool;
  return tool.main(argc, argv);
}

/// @endcond