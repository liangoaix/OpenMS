// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(FeatureGroupingAlgorithmQT, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmQT* ptr = 0;
FeatureGroupingAlgorithmQT* nullPointer = 0;
START_SECTION((FeatureGroupingAlgorithmQT()))
	ptr = new FeatureGroupingAlgorithmQT();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmQT()))
	delete ptr;
END_SECTION

START_SECTION((static FeatureGroupingAlgorithm* create()))
	FeatureGroupingAlgorithm* ptr2 = 0;
  FeatureGroupingAlgorithm* base_NullPointer = 0;
	ptr2 = FeatureGroupingAlgorithmQT::create();
  TEST_NOT_EQUAL(ptr2, base_NullPointer)
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(FeatureGroupingAlgorithmQT::getProductName(), "unlabeled_qt")
END_SECTION

START_SECTION((virtual void group(const std::vector< FeatureMap >& maps, ConsensusMap& out)))
	// This is tested extensively in TEST/TOPP
	NOT_TESTABLE;
END_SECTION

START_SECTION((virtual void group(const std::vector<ConsensusMap>& maps, ConsensusMap& out)))
	// This is tested extensively in TEST/TOPP
	NOT_TESTABLE;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



