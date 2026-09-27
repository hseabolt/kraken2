/*
 * Rank-aware taxonomic bin resolution for Kraken 2.
 *
 * This file is an optional output-layer extension and does not alter
 * Kraken 2 database formats or classification semantics.
 */
#ifndef KRAKEN2_RANK_BIN_RESOLVER_H_
#define KRAKEN2_RANK_BIN_RESOLVER_H_

#include <string>
#include <vector>

#include "kraken2_data.h"
#include "taxonomy.h"

namespace kraken2 {

class RankBinResolver {
 public:
  RankBinResolver(const Taxonomy &taxonomy, const std::string &target_rank);

  // Return the nearest ancestor (including the called node itself) whose rank
  // matches the requested target rank.  Returns 0 when the Kraken call is
  // above the requested rank or no matching ancestor exists.
  taxid_t Resolve(taxid_t classified_taxon) const;

  const std::string &target_rank() const { return target_rank_; }

  static std::string NormalizeRank(const std::string &rank);

 private:
  const Taxonomy &taxonomy_;
  std::string target_rank_;
  std::vector<taxid_t> rank_bin_map_;
};

}  // namespace kraken2

#endif  // KRAKEN2_RANK_BIN_RESOLVER_H_
