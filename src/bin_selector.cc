/* Taxonomic include/exclude filtering for Kraken 2 rank bins. */
#include "bin_selector.h"

#include <err.h>
#include <sysexits.h>
#include <unordered_map>

namespace kraken2 {

BinSelector::BinSelector(const Taxonomy &taxonomy,
                         const std::vector<uint64_t> &include_external_taxids,
                         const std::vector<uint64_t> &exclude_external_taxids) {
  const size_t n = taxonomy.node_count();
  allowed_.assign(n, 0);

  std::unordered_map<uint64_t, taxid_t> ext_to_int;
  ext_to_int.reserve(n);
  for (taxid_t node = 1; node < n; ++node)
    ext_to_int[taxonomy.nodes()[node].external_id] = node;

  std::vector<uint8_t> include_root(n, 0), exclude_root(n, 0);
  for (uint64_t ext : include_external_taxids) {
    auto it = ext_to_int.find(ext);
    if (it == ext_to_int.end())
      errx(EX_DATAERR, "--bin-include-clade taxid %llu is not present in taxonomy",
           static_cast<unsigned long long>(ext));
    include_root[it->second] = 1;
  }
  for (uint64_t ext : exclude_external_taxids) {
    auto it = ext_to_int.find(ext);
    if (it == ext_to_int.end())
      errx(EX_DATAERR, "--bin-exclude-clade taxid %llu is not present in taxonomy",
           static_cast<unsigned long long>(ext));
    exclude_root[it->second] = 1;
  }

  const bool has_includes = !include_external_taxids.empty();
  std::vector<uint8_t> in_include(n, 0), in_exclude(n, 0);

  for (taxid_t node = 1; node < n; ++node) {
    const taxid_t parent = taxonomy.nodes()[node].parent_id;
    const bool parent_include = parent < n && in_include[parent];
    const bool parent_exclude = parent < n && in_exclude[parent];
    in_include[node] = include_root[node] || parent_include;
    in_exclude[node] = exclude_root[node] || parent_exclude;
    allowed_[node] = ((!has_includes || in_include[node]) && !in_exclude[node]) ? 1 : 0;
  }
}

bool BinSelector::Allows(taxid_t taxon) const {
  return taxon != 0 && taxon < allowed_.size() && allowed_[taxon] != 0;
}

}  // namespace kraken2
