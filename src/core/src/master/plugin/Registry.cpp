#include "master/plugin/Registry.hpp"
#include <stdexcept>

using namespace core::master::plugin;

std::unique_ptr<IProgram> Registry::make_program(const std::string& key, const KV& kv, const core::master::RunContext& rc) const {
  auto it = programs_.find(key);
  if (it == programs_.end()) throw std::runtime_error("No program factory for key: " + key);
  return it->second(kv, rc);
}

std::shared_ptr<IAction> Registry::make_action(const std::string& key, const KV& kv, const core::master::RunContext& rc) const {
  auto it = actions_.find(key);
  if (it == actions_.end()) throw std::runtime_error("No action factory for key: " + key);
  return it->second(kv, rc);
}

std::shared_ptr<IGlobal> Registry::make_global(const std::string& key, const KV& kv, const core::master::RunContext& rc) const {
  auto it = globals_.find(key);
  if (it == globals_.end()) throw std::runtime_error("No global factory for key: " + key);
  return it->second(kv, rc);
}
