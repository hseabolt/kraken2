/* Streaming taxonomic FASTA/FASTQ bin writer for Kraken 2. */
#include "bin_writer.h"

#include <algorithm>
#include <cerrno>
#include <cctype>
#include <climits>
#include <cstdio>
#include <cstring>
#include <err.h>
#include <sys/stat.h>
#include <sysexits.h>

namespace kraken2 {

BinWriter::BinWriter(const Taxonomy &taxonomy,
                     const std::string &output_directory,
                     const std::string &target_rank,
                     size_t max_open_files,
                     const std::string &compression,
                     int compression_level,
                     uint64_t min_fragments,
                     bool write_unresolved)
    : taxonomy_(taxonomy),
      output_directory_(output_directory),
      target_rank_(target_rank),
      max_open_files_(max_open_files),
      compression_(compression),
      compression_level_(compression_level),
      min_fragments_(min_fragments),
      write_unresolved_(write_unresolved) {
  if (output_directory_.empty())
    errx(EX_USAGE, "bin output directory cannot be empty");
  if (target_rank_.empty())
    errx(EX_USAGE, "bin rank cannot be empty");
  if (max_open_files_ < 2)
    errx(EX_USAGE, "max open bin files must be at least 2");
  if (compression_ != "none" && compression_ != "gzip")
    errx(EX_USAGE, "bin compression must be 'none' or 'gzip'");
  if (compression_level_ < 1 || compression_level_ > 9)
    errx(EX_USAGE, "bin compression level must be in [1,9]");
  if (min_fragments_ < 1)
    errx(EX_USAGE, "minimum bin fragments must be at least 1");
  EnsureOutputDirectory();
}

BinWriter::~BinWriter() {
  CloseAll();
}

void BinWriter::Initialize(SequenceFormat format, bool paired) {
  if (format != FORMAT_FASTA && format != FORMAT_FASTQ)
    errx(EX_DATAERR, "unable to determine input format for taxonomic bins");

  if (!initialized_) {
    initialized_ = true;
    format_ = format;
    paired_ = paired;
    extension_ = (format == FORMAT_FASTQ) ? ".fastq" : ".fasta";
    if (compression_ == "gzip")
      extension_ += ".gz";
    return;
  }

  if (format_ != format)
    errx(EX_DATAERR, "all inputs in a binned Kraken run must use the same FASTA/FASTQ format");
  if (paired_ != paired)
    errx(EX_DATAERR, "pairing mode changed during a binned Kraken run");
}

void BinWriter::WriteBlocks(const std::vector<BinBlock> &blocks) {
  if (!initialized_)
    errx(EX_SOFTWARE, "BinWriter used before initialization");

  for (const BinBlock &block : blocks) {
    if (block.fragments == 0)
      continue;
    if (block.taxon == 0 && !write_unresolved_)
      continue;

    BinStats &stats = stats_[block.taxon];
    if (stats.filename1.empty()) {
      stats.filename1 = PathFor(block.taxon, paired_ ? 1 : 0);
      if (paired_)
        stats.filename2 = PathFor(block.taxon, 2);
    }

    if (!block.mate1_data.empty())
      WriteData(stats.filename1, block.mate1_data);
    if (paired_ && !block.mate2_data.empty())
      WriteData(stats.filename2, block.mate2_data);

    stats.fragments += block.fragments;
    stats.reads += block.reads;
    stats.bases += block.bases;
  }
}

void BinWriter::Finalize() {
  if (finalized_)
    return;
  finalized_ = true;
  CloseAll();

  const std::string manifest_path = JoinPath(output_directory_, "bins.tsv");
  std::ofstream manifest(manifest_path.c_str(), std::ios::out | std::ios::trunc);
  if (!manifest)
    errx(EX_CANTCREAT, "unable to create bin manifest %s: %s",
         manifest_path.c_str(), std::strerror(errno));

  manifest << "bin_taxid\trank\tscientific_name\tfragments\treads\tbases"
           << "\tsequence_file_1\tsequence_file_2\n";

  for (const auto &kv : stats_) {
    const taxid_t internal_taxid = kv.first;
    const BinStats &stats = kv.second;
    if (stats.fragments < min_fragments_) {
      if (!stats.filename1.empty())
        std::remove(stats.filename1.c_str());
      if (!stats.filename2.empty())
        std::remove(stats.filename2.c_str());
      continue;
    }

    manifest << ExternalTaxid(internal_taxid) << '\t'
             << target_rank_ << '\t'
             << ScientificName(internal_taxid) << '\t'
             << stats.fragments << '\t'
             << stats.reads << '\t'
             << stats.bases << '\t'
             << stats.filename1 << '\t'
             << stats.filename2 << '\n';
  }
}

void BinWriter::EnsureOutputDirectory() {
  struct stat sb;
  if (stat(output_directory_.c_str(), &sb) == 0) {
    if (!S_ISDIR(sb.st_mode))
      errx(EX_CANTCREAT, "bin output path exists but is not a directory: %s",
           output_directory_.c_str());
    return;
  }
  if (errno != ENOENT)
    errx(EX_CANTCREAT, "unable to inspect bin output directory %s: %s",
         output_directory_.c_str(), std::strerror(errno));

  std::string current;
  if (!output_directory_.empty() && output_directory_[0] == '/')
    current = "/";

  size_t start = 0;
  while (start < output_directory_.size()) {
    size_t slash = output_directory_.find('/', start);
    std::string part = output_directory_.substr(
        start, slash == std::string::npos ? std::string::npos : slash - start);
    start = (slash == std::string::npos) ? output_directory_.size() : slash + 1;
    if (part.empty() || part == ".")
      continue;
    if (!current.empty() && current.back() != '/')
      current.push_back('/');
    current += part;
    if (stat(current.c_str(), &sb) == 0) {
      if (!S_ISDIR(sb.st_mode))
        errx(EX_CANTCREAT, "path component is not a directory: %s", current.c_str());
      continue;
    }
    if (mkdir(current.c_str(), 0775) != 0 && errno != EEXIST)
      errx(EX_CANTCREAT, "unable to create bin output directory %s: %s",
           current.c_str(), std::strerror(errno));
  }
}

BinWriter::OpenStream &BinWriter::GetStream(const std::string &pathname) {
  OpenStream &entry = open_streams_[pathname];
  entry.last_used = ++access_tick_;
  if (entry.plain != nullptr || entry.gzip != nullptr)
    return entry;

  if (open_file_count_ >= max_open_files_)
    CloseLeastRecentlyUsed();

  const bool append = initialized_paths_.count(pathname) != 0;
  if (compression_ == "gzip") {
    entry.gzip = gzopen(pathname.c_str(), append ? "ab" : "wb");
    if (entry.gzip == nullptr)
      errx(EX_CANTCREAT, "unable to open gzip bin output %s", pathname.c_str());
    gzsetparams(entry.gzip, compression_level_, Z_DEFAULT_STRATEGY);
    gzbuffer(entry.gzip, 1U << 20);
  } else {
    std::ios_base::openmode mode = std::ios::out | std::ios::binary;
    mode |= append ? std::ios::app : std::ios::trunc;
    entry.plain = new std::ofstream(pathname.c_str(), mode);
    if (!*entry.plain) {
      delete entry.plain;
      entry.plain = nullptr;
      errx(EX_CANTCREAT, "unable to open bin output %s: %s",
           pathname.c_str(), std::strerror(errno));
    }
  }

  initialized_paths_.insert(pathname);
  ++open_file_count_;
  return entry;
}

void BinWriter::WriteData(const std::string &pathname, const std::string &data) {
  if (data.empty())
    return;
  OpenStream &entry = GetStream(pathname);
  if (entry.plain != nullptr) {
    entry.plain->write(data.data(), static_cast<std::streamsize>(data.size()));
    if (!*entry.plain)
      errx(EX_IOERR, "failed writing bin output %s", pathname.c_str());
    return;
  }

  size_t offset = 0;
  while (offset < data.size()) {
    const size_t remain = data.size() - offset;
    const unsigned int chunk = static_cast<unsigned int>(
        std::min<size_t>(remain, static_cast<size_t>(INT_MAX)));
    const int written = gzwrite(entry.gzip, data.data() + offset, chunk);
    if (written <= 0) {
      int zerr = Z_OK;
      const char *msg = gzerror(entry.gzip, &zerr);
      errx(EX_IOERR, "failed writing gzip bin output %s: %s",
           pathname.c_str(), msg ? msg : "zlib error");
    }
    offset += static_cast<size_t>(written);
  }
}

void BinWriter::CloseStream(OpenStream &entry) {
  if (entry.plain != nullptr) {
    entry.plain->flush();
    entry.plain->close();
    delete entry.plain;
    entry.plain = nullptr;
    --open_file_count_;
  } else if (entry.gzip != nullptr) {
    if (gzclose(entry.gzip) != Z_OK)
      errx(EX_IOERR, "failed closing gzip bin output");
    entry.gzip = nullptr;
    --open_file_count_;
  }
}

void BinWriter::CloseLeastRecentlyUsed() {
  std::map<std::string, OpenStream>::iterator victim = open_streams_.end();
  for (auto it = open_streams_.begin(); it != open_streams_.end(); ++it) {
    if (it->second.plain == nullptr && it->second.gzip == nullptr)
      continue;
    if (victim == open_streams_.end() ||
        it->second.last_used < victim->second.last_used)
      victim = it;
  }
  if (victim != open_streams_.end())
    CloseStream(victim->second);
}

void BinWriter::CloseAll() {
  for (auto &kv : open_streams_) {
    if (kv.second.plain != nullptr || kv.second.gzip != nullptr)
      CloseStream(kv.second);
  }
  open_file_count_ = 0;
}

std::string BinWriter::BinStem(taxid_t internal_taxid) const {
  if (internal_taxid == 0)
    return "unresolved_" + SanitizeFilename(target_rank_);
  return "taxid_" + std::to_string(ExternalTaxid(internal_taxid)) + "_" +
         SanitizeFilename(ScientificName(internal_taxid));
}

std::string BinWriter::PathFor(taxid_t internal_taxid, int mate) const {
  std::string filename = BinStem(internal_taxid);
  if (paired_ && mate > 0)
    filename += (mate == 1) ? "_R1" : "_R2";
  filename += extension_;
  return JoinPath(output_directory_, filename);
}

std::string BinWriter::ScientificName(taxid_t internal_taxid) const {
  if (internal_taxid == 0)
    return "rank unresolved";
  if (internal_taxid >= taxonomy_.node_count())
    return "unknown";
  const TaxonomyNode &node = taxonomy_.nodes()[internal_taxid];
  const char *name = taxonomy_.name_data() + node.name_offset;
  return name ? std::string(name) : std::string("unknown");
}

uint64_t BinWriter::ExternalTaxid(taxid_t internal_taxid) const {
  if (internal_taxid == 0 || internal_taxid >= taxonomy_.node_count())
    return 0;
  return taxonomy_.nodes()[internal_taxid].external_id;
}

std::string BinWriter::SanitizeFilename(const std::string &value) {
  std::string result;
  result.reserve(value.size());
  bool last_underscore = false;
  for (char c : value) {
    unsigned char uc = static_cast<unsigned char>(c);
    if (std::isalnum(uc) || c == '.' || c == '-') {
      result.push_back(c);
      last_underscore = false;
    } else if (!last_underscore) {
      result.push_back('_');
      last_underscore = true;
    }
  }
  while (!result.empty() && result.back() == '_')
    result.pop_back();
  if (result.empty())
    result = "unnamed";
  return result;
}

std::string BinWriter::JoinPath(const std::string &dir, const std::string &name) {
  if (dir.empty() || dir == ".")
    return name;
  if (dir.back() == '/')
    return dir + name;
  return dir + "/" + name;
}

}  // namespace kraken2
