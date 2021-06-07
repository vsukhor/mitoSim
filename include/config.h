/* =============================================================================
   Copyright (C) 2015 Valerii Sukhorukov & Michael Meyer-Hermann,
   Helmholtz Center for Infection Research (Braunschweig, Germany).
   All Rights Reserved.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
================================================================================
*/

/**
* @file config.h
* @brief Configuration parameters, import from the configuration file.
* @author Valerii Sukhorukov.
*/

#ifndef MITOSIM_CONFIG_H
#define MITOSIM_CONFIG_H

#include <filesystem>
#include <string>

#include "utils/common/msgr.h"
#include "utils/config/reader.h"

namespace mitosim {

using utils::common::bools;
using utils::common::huge;
using utils::common::onehuge;
using utils::common::zerohuge;

/**
 * @brief Encapsulates and reads configuration parameters.
 * @details Reads from a file and stores initial configuration
 * parameters for the hole modeling session.
 */
template<typename realT>
class Config {

public:

    using Msgr = utils::common::Msgr;
    using szt = utils::common::szt;
    using ulong = utils::common::ulong;

private:

    const std::filesystem::path file;  ///< Configuration file name.
    utils::config::Reader read;  ///< Generic reader of configuraion files.

public:

    const std::string fnameSuffix;  ///< Application-specific suffix of the configuration file.
    const std::string runName;      ///< Current run index as a string.
    const std::filesystem::path workingDirOut;  ///< Output directory.

    const realT timeTotal;      ///< Total time to simulate.
    const szt   logFrequency;   ///< (iterations): Every logfreq iteration steps log is updated.
    const szt   saveFrequency;  ///< (iterations): Every saveFrequency steps save_mitos is executed.

    const realT edgeLength;  ///< (micrometers) Edge Length.
    const szt   mtmassini;   ///< Total chondriome length, edges.
    const szt   segmassini;  ///< Segment lengths at the beginning of simulation.

    // FISSION
    const bool  use_fission;   ///< A flag for the activation of the fission reaction.
    const realT rate_fission;  ///< Basic probability of breaking up a junction.

    // FUSION
    const bool  use_11_fusion;   ///< A flag for the activation of the fussion_11 reaction.
    const realT fusion_rate_11;  ///< Probability of a free end to bind to another free end.
    const bool  use_12_fusion;   ///< A flag for the activation of the fusion_12 reaction.
    const realT fusion_rate_12;  ///< Probability of a free end to bind to a side.
    const bool  use_1L_fusion;   ///< A flag for the activation of the fusion_1L reaction.
    const realT fusion_rate_1L;  ///< Probability of a free end to bind to a separate cycle junction node.

    /**
     * @brief Constructor.
     * @par workingDirOut Directory for the output.
     * @par configSuffix Application-specific suffix of the configuration file.
     * @par runName Run index.
     * @par msgr Output message processor.
     */
    explicit Config(
            const std::filesystem::path& workingDirOut,
            const std::string& fnameSuffix,
            const std::string& runName,
            Msgr& msgr
            )
        : file {workingDirOut / (std::string("config")+"_"+fnameSuffix+".txt")}
        , read {file, &msgr}
        , fnameSuffix {fnameSuffix}
        , runName {runName}
        , workingDirOut {workingDirOut}
        , timeTotal     {read("timeTotal", zerohuge<realT>)}
        , logFrequency  {read("logFrequency", onehuge<szt>)}
        , saveFrequency {read("saveFrequency", onehuge<szt>)}

        , edgeLength {read("edgeLength", zerohuge<realT>)}
        , mtmassini  {read("mtmassini", onehuge<szt>)}
        , segmassini {read("segmassini", onehuge<szt>)}

    // FISSION
        , use_fission  {read("use_fission", bools)}
        , rate_fission {read("rate_fission", zerohuge<realT>)}

    // FUSION
        , use_11_fusion  {read("use_11_fusion", bools)}
        , fusion_rate_11 {read("fusion_rate_11", zerohuge<realT>)}
        , use_12_fusion  {read("use_12_fusion", bools)}
        , fusion_rate_12 {read("fusion_rate_12", zerohuge<realT>)}
        , use_1L_fusion  {read("use_1L_fusion", bools)}
        , fusion_rate_1L {read("fusion_rate_1L", zerohuge<realT>)}
    {}
};

}  // namespace mitosim

#endif  // MITOSIM_CONFIG_H
