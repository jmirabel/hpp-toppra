#ifndef HPP_TOPPRA_SERIALIZATION_HH
#define HPP_TOPPRA_SERIALIZATION_HH

#include <fstream>
#include <hpp/core/path.hh>
#include <hpp/pinocchio/serialization.hh>

namespace hpp {
namespace core {
namespace parser {
template <class Archive>
void serializePath(
    DevicePtr_t& device,
    PathPtr_t& path,
    const std::string& filename)
{
  typename std::conditional<Archive::is_saving::value, std::ofstream,
                            std::ifstream>::type fs(filename);
  Archive ar(fs);
  if (device)
    ar.insert(device->name(), device.get());
  ar.initialize();
  //ar& hpp::serialization::make_nvp("device", device);
  ar& hpp::serialization::make_nvp("path", path);
}
/// \}
}  // namespace parser
}  // namespace core
}  // namespace hpp
#endif  // HPP_CORE_PARSER_ROADMAP_HH
