#include <ocs2_ipm/ipm/SlackDualData.h>

namespace ocs2 {

void setBarrier(size_t nc, scalar_t barrier, SlackDualData& slackDual) {
  slackDual.slack.resize(nc);
  slackDual.dual.resize(nc);
  setBarrier(barrier, slackDual);
}


void setBarrier(scalar_t barrier, SlackDualData& slackDual) {
  slackDual.barrier = barrier;
  slackDual.slack.fill(std::sqrt(barrier));
  slackDual.dual.fill(std::sqrt(barrier));
}


std::string checkSize(int constraintDim, const SlackDualData& slackDual, const std::string& dataName) {
  std::stringstream errorDescription;

  if (slackDual.slack.size() != constraintDim) {
    errorDescription << dataName << ".slack.size() != " << constraintDim << "\n";
  }
  if (slackDual.dual.size() != constraintDim) {
    errorDescription << dataName << ".dual.size() != " << constraintDim << "\n";
  }
  return errorDescription.str();
}


std::string checkPositive(const SlackDualData& slackDual, const std::string& dataName) {
  std::stringstream errorDescription;
  if (slackDual.barrier <= 0.0) {
    errorDescription << dataName << ".barrier = " << slackDual.barrier << "\n";
  }
  if (slackDual.slack.size() > 0) {
    if (slackDual.slack.minCoeff() <= 0.0) {
      errorDescription << dataName << ".slack = " << slackDual.slack.transpose() << "\n";
    }
  }
  if (slackDual.dual.size() > 0) {
    if (slackDual.dual.minCoeff() <= 0.0) {
      errorDescription << dataName << ".dual = " << slackDual.dual.transpose() << "\n";
    }
  }
  return errorDescription.str();
}

}  // namespace ocs2
