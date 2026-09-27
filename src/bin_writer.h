/*
 * Streaming taxonomic FASTA/FASTQ bin writer for Kraken 2.
 */
#ifndef KRAKEN2_BIN_WRITER_H_
#define KRAKEN2_BIN_WRITER_H_

#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "kraken2_data.h"
#include "seqreader.h"
#include "taxonomy.h"

namespace kraken2 {

// Data accumulated by a classification worker for one taxonomic bin in one
// input block.  taxon is Kraken's internal taxonomy ID.  taxon == 0 is
// reserved here for a classified read that cannot be resolved to the requested
// bin rank (e.g. a family-level call when --bin-rank genus is requested).
struct BinBlock {
  taxid_t taxon = 0;
  std::string mate1_data;
  std::string mate2_data;
  uint64_t fragments = 0;
  uint64_t reads = 0;
  uint64_t bases = 0;
};

class BinWriter {
 public:
  BinWriter(const Taxonomy &taxonomy,
            const std::string &output_directory,
            const std::string &target_rank,
            size_t max_open_files = 128);
  ~BinWriter();

  // Idempotent.  All inputs in one run must use the same sequence format and
  // pairing mode.
  void Initialize(SequenceFormat format, bool paired);

  // Called only from Kraken's ordered output section.
  void WriteBlocks(const std::vector<BinBlock> &blocks);

  // Flush/close bin files and write bins.tsv.  Safe to call more than once.
  void Finalize();

 private:
  struct OpenStream {
    std::ofstream *stream = nullptr;
    uint64_t last_used = 0;
  };

  struct BinStats {
    uint64_t fragments = 0;
    uint64_t reads = 0;
    uint64_t bases = 0;
    std::string filename1;
    std::string filename2;
  };

  const Taxonomy &taxonomy_;
  std::string output_directory_;
  std::string target_rank_;
  size_t max_open_files_;
  bool initialized_ = false;
  bool finalized_ = false;
  bool paired_ = false;
  SequenceFormat format_ = FORMAT_AUTO_DETECT;
  std::string extension_;
  uint64_t access_tick_ = 0;
  size_t open_file_count_ = 0;

  std::map<std::string, OpenStream> open_streams_;
  std::set<std::string> initialized_paths_;
  std::map<taxid_t, BinStats> stats_;

  void EnsureOutputDirectory();
  std::ofstream &GetStream(const std::string &pathname);
  void CloseLeastRecentlyUsed();
  void CloseAll();

  std::string BinStem(taxid_t internal_taxid) const;
  std::string PathFor(taxid_t internal_taxid, int mate) const;
  std::string ScientificName(taxid_t internal_taxid) const;
  uint64_t ExternalTaxid(taxid_t internal_taxid) const;
  static std::string SanitizeFilename(const std::string &value);
  static std::string JoinPath(const std::string &dir, const std::string &name);
};

}  // namespace kraken2

#endif  // KRAKEN2_BIN_WRITER_H_
