/*
 * Rank-aware taxonomic bin resolution for Kraken 2.
 */
#include "rank_bin_resolver.h"

#include <algorithm>
#include <cctype>
#include <err.h>
#include <sysexits.h>

namespace kraken2 {

std::string RankBinResolver::NormalizeRank(const std::string &rank) {
  std::string normalized;
  normalized.reserve(rank.size());

  bool last_was_space = false;
  for (char c : rank) {
    unsigned char uc = static_cast<unsigned char>(c);
    char out = static_cast<char>(std::tolower(uc));
    if (out == '_' || out == '-')
      out = ' ';

    if (std::isspace(static_cast<unsigned char>(out))) {
      if (!normalized.empty() && !last_was_space) {
        normalized.push_back(' ');
        last_was_space = true;
      }
    } else {
      normalized.push_back(out);
      last_was_space = false;
    }
  }

  while (!normalized.empty() && normalized.back() == ' ')
    normalized.pop_back();

  // NCBI taxonomy uses "superkingdom" for the rank Kraken reports as domain.
  if (normalized == "domain")
    normalized = "superkingdom";

  return normalized;
}

RankBinResolver::RankBinResolver(const Taxonomy &taxonomy,
                                 const std::string &target_rank)
    : taxonomy_(taxonomy), target_rank_(NormalizeRank(target_rank)) {
  if (target_rank_.empty())
    errx(EX_USAGE, "bin rank cannot be empty");

  rank_bin_map_.assign(taxonomy_.node_count(), 0);
  bool rank_seen = false;

  // Kraken's TaxonomyNode contract requires parent_id to be lower-numbered
  // than the child.  This lets us propagate the nearest target-rank ancestor
  // in a single O(N) pass and resolve every read in O(1).
  for (taxid_t node_id = 1; node_id < taxonomy_.node_count(); ++node_id) {
    const TaxonomyNode &node = taxonomy_.nodes()[node_id];
    const char *rank_ptr = taxonomy_.rank_data() + node.rank_offset;

    if (rank_ptr != nullptr && target_rank_ == rank_ptr) {
      rank_bin_map_[node_id] = node_id;
      rank_seen = true;
    } else if (node.parent_id < rank_bin_map_.size()) {
      rank_bin_map_[node_id] = rank_bin_map_[node.parent_id];
    }
  }

  if (!rank_seen) {
    errx(EX_DATAERR,
         "requested bin rank '%s' is not present in the loaded taxonomy",
         target_rank_.c_str());
  }
}

taxid_t RankBinResolver::Resolve(taxid_t classified_taxon) const {
  if (classified_taxon == 0 || classified_taxon >= rank_bin_map_.size())
    return 0;
  return rank_bin_map_[classified_taxon];
}

}  // namespace kraken2
