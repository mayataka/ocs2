/******************************************************************************
Copyright (c) 2020, Farbod Farshidian. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#pragma once

#include <memory>

#include <ocs2_core/loopshaping/LoopshapingDefinition.h>

#include "ocs2_oc/synchronized_module/ModeScheduleManager.h"

namespace ocs2 {

class LoopshapingModeScheduleManager : public ModeScheduleManager {
 public:
  /**
   * Constructor.
   * @param [in] modeScheduleManagerPtr: A shared pointer to the original ModeScheduleManager.
   * @param [in] loopshapingDefinitionPtr: A shared pointer to the loopshaping definition.
   */
  LoopshapingModeScheduleManager(std::shared_ptr<ModeScheduleManager> modeScheduleManagerPtr,
                                 std::shared_ptr<LoopshapingDefinition> loopshapingDefinitionPtr);

  /** Destructor */
  ~LoopshapingModeScheduleManager() override = default;

 private:
  void preSolverRunImpl(scalar_t initTime, scalar_t finalTime, const vector_t& currentState,
                        const CostDesiredTrajectories& costDesiredTrajectory, ModeSchedule& modeSchedule) override;

  std::shared_ptr<ocs2::ModeScheduleManager> modeScheduleManagerPtr_;
  std::shared_ptr<ocs2::LoopshapingDefinition> loopshapingDefinitionPtr_;
};

}  // namespace ocs2
