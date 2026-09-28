/* Taxonomic include/exclude filtering for Kraken 2 rank bins. */
#ifndef KRAKEN2_BIN_SELECTOR_H_
#define KRAKEN2_BIN_SELECTOR_H_

#include <cstdint>
#include <vector>

#include "kraken2_data.h"
#include "taxonomy.h"

namespace kraken2 {

class BinSelector {
 public:
  BinSelector(const Taxonomy &taxonomy,
              const std::vector<uint64_t> &include_external_taxids,
              const std::vector<uint64_t> &exclude_external_taxids);

  // taxon is an internal Kraken taxonomy ID after rank resolution.
  bool Allows(taxid_t taxon) const;

 private:
  std::vector<uint8_t> allowed_;
};

}  // namespace kraken2

#endif  // KRAKEN2_BIN_SELECTOR_H_
