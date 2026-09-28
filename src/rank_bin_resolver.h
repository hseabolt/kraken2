/*
 * Rank-aware taxonomic bin resolution for Kraken 2.
 *
 * Optional output-layer extension.  Does not alter Kraken database formats
 * or classification semantics.
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
  // matches target_rank.  Return 0 when no such ancestor exists.  In the
  // special target rank "taxid", return the classified node unchanged.
  taxid_t Resolve(taxid_t classified_taxon) const;

  const std::string &target_rank() const { return target_rank_; }
  bool exact_taxid_mode() const { return exact_taxid_mode_; }

  static std::string NormalizeRank(const std::string &rank);

 private:
  const Taxonomy &taxonomy_;
  std::string target_rank_;
  bool exact_taxid_mode_ = false;
  std::vector<taxid_t> rank_bin_map_;
};

}  // namespace kraken2

#endif  // KRAKEN2_RANK_BIN_RESOLVER_H_
