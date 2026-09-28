/*
 * Post-classification taxonomic read binning for Kraken 2.
 *
 * This utility is primarily used by `k2 classify` after multi-database
 * classification has been merged.  It intentionally consumes the final
 * standard Kraken output plus the original sequences rather than changing
 * Kraken's classification/report formats.
 */
#include <err.h>
#include <getopt.h>
#include <sysexits.h>

#include <cerrno>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "bin_selector.h"
#include "bin_writer.h"
#include "rank_bin_resolver.h"
#include "seqreader.h"
#include "taxonomy.h"

using namespace kraken2;
using std::ifstream;
using std::map;
using std::string;
using std::vector;

namespace {

static const size_t DEFAULT_FLUSH_FRAGMENTS = 10000;

struct Options {
  string taxonomy_filename;
  string kraken_filename;
  string output_directory;
  string target_rank;
  size_t max_open_files = 128;
  string compression = "none";
  int compression_level = 1;
  uint64_t min_fragments = 1;
  bool write_unresolved = true;
  bool paired = false;
  bool interleaved = false;
  vector<uint64_t> include_taxids;
  vector<uint64_t> exclude_taxids;
  vector<string> filenames;
};

void usage(int exit_code = EX_USAGE) {
  std::cerr
      << "Usage: bin_reads -t TAXO.K2D -k KRAKEN.OUT -B DIR -r RANK [options] <sequence files>\n\n"
      << "Required:\n"
      << "  -t FILE   Kraken taxonomy file (taxo.k2d or merged taxonomy)\n"
      << "  -k FILE   Final standard Kraken classification output\n"
      << "  -B DIR    Output directory for taxonomic sequence bins\n"
      << "  -r RANK   Target rank; arbitrary taxonomy rank or special value 'taxid'\n\n"
      << "Sequence mode:\n"
      << "  -P        Paired-end input in two files\n"
      << "  -S        Interleaved paired-end input in one file\n\n"
      << "Binning controls:\n"
      << "  -L INT    Maximum simultaneously open bin files (default: 128)\n"
      << "  -Z STR    Bin compression: none or gzip (default: none)\n"
      << "  -l INT    gzip compression level 1-9 (default: 1)\n"
      << "  -f INT    Remove bins with fewer than INT fragments (default: 1)\n"
      << "  -I TAXID  Include only bins in this clade; repeatable\n"
      << "  -E TAXID  Exclude bins in this clade; repeatable; exclusion wins\n"
      << "  -W STR    Unresolved-at-rank policy: write or drop (default: write)\n"
      << "  -h        Show this help\n\n"
      << "Input files supplied directly to bin_reads must be uncompressed. k2 handles\n"
      << "gzip/bzip2/xz inputs by feeding decompressed FIFOs to this utility.\n";
  std::exit(exit_code);
}

uint64_t parse_uint64(const char *text, const char *label, uint64_t min_value = 0) {
  if (text == nullptr || *text == '\0')
    errx(EX_USAGE, "%s requires an integer", label);
  if (*text == '-')
    errx(EX_USAGE, "invalid %s: %s", label, text);
  char *end = nullptr;
  errno = 0;
  unsigned long long value = std::strtoull(text, &end, 10);
  if (errno != 0 || end == text || *end != '\0' || value < min_value)
    errx(EX_USAGE, "invalid %s: %s", label, text);
  return static_cast<uint64_t>(value);
}

Options parse_options(int argc, char **argv) {
  Options opts;
  int c;
  while ((c = getopt(argc, argv, "h?t:k:B:r:L:Z:l:f:I:E:W:PS")) != -1) {
    switch (c) {
      case 'h': case '?':
        usage(c == 'h' ? 0 : EX_USAGE);
        break;
      case 't': opts.taxonomy_filename = optarg; break;
      case 'k': opts.kraken_filename = optarg; break;
      case 'B': opts.output_directory = optarg; break;
      case 'r': opts.target_rank = optarg; break;
      case 'L': opts.max_open_files = static_cast<size_t>(parse_uint64(optarg, "-L", 2)); break;
      case 'Z': opts.compression = optarg; break;
      case 'l': opts.compression_level = static_cast<int>(parse_uint64(optarg, "-l", 1)); break;
      case 'f': opts.min_fragments = parse_uint64(optarg, "-f", 1); break;
      case 'I': opts.include_taxids.push_back(parse_uint64(optarg, "-I", 1)); break;
      case 'E': opts.exclude_taxids.push_back(parse_uint64(optarg, "-E", 1)); break;
      case 'W': {
        string mode(optarg);
        if (mode == "write") opts.write_unresolved = true;
        else if (mode == "drop") opts.write_unresolved = false;
        else errx(EX_USAGE, "-W must be 'write' or 'drop'");
        break;
      }
      case 'P': opts.paired = true; break;
      case 'S': opts.interleaved = true; break;
      default: usage();
    }
  }

  if (opts.taxonomy_filename.empty() || opts.kraken_filename.empty() ||
      opts.output_directory.empty() || opts.target_rank.empty())
    usage();
  if (opts.compression != "none" && opts.compression != "gzip")
    errx(EX_USAGE, "-Z/--bin-compression must be 'none' or 'gzip'");
  if (opts.compression_level < 1 || opts.compression_level > 9)
    errx(EX_USAGE, "-l/--bin-compression-level must be in [1,9]");
  if (opts.paired && opts.interleaved)
    errx(EX_USAGE, "-P and -S are mutually exclusive");

  for (int i = optind; i < argc; ++i)
    opts.filenames.push_back(argv[i]);

  if (opts.filenames.empty())
    errx(EX_USAGE, "at least one sequence file is required");
  if (opts.paired && opts.filenames.size() % 2 != 0)
    errx(EX_USAGE, "paired mode requires an even number of sequence files");
  return opts;
}

string normalized_read_id(const string &id) {
  string out = id;
  if (out.size() > 2 && out[out.size() - 2] == '/' &&
      (out.back() == '1' || out.back() == '2'))
    out.resize(out.size() - 2);
  return out;
}

bool ids_match(const string &a, const string &b) {
  return a == b || normalized_read_id(a) == normalized_read_id(b);
}

uint64_t parse_external_taxid(const string &field, uint64_t line_no) {
  bool all_digits = !field.empty();
  for (char c : field) {
    if (!std::isdigit(static_cast<unsigned char>(c))) {
      all_digits = false;
      break;
    }
  }
  if (all_digits)
    return static_cast<uint64_t>(std::stoull(field));

  // --use-names outputs the third field as "scientific name (taxid N)".
  const string marker = "(taxid ";
  const size_t pos = field.rfind(marker);
  if (pos == string::npos || field.back() != ')')
    errx(EX_DATAERR, "unable to parse taxid from Kraken output line %llu: %s",
         static_cast<unsigned long long>(line_no), field.c_str());
  const size_t begin = pos + marker.size();
  const string digits = field.substr(begin, field.size() - begin - 1);
  if (digits.empty())
    errx(EX_DATAERR, "empty taxid in Kraken output line %llu",
         static_cast<unsigned long long>(line_no));
  for (char c : digits) {
    if (!std::isdigit(static_cast<unsigned char>(c)))
      errx(EX_DATAERR, "invalid taxid in Kraken output line %llu: %s",
           static_cast<unsigned long long>(line_no), field.c_str());
  }
  return static_cast<uint64_t>(std::stoull(digits));
}

struct KrakenRecord {
  bool classified = false;
  string read_id;
  uint64_t external_taxid = 0;
};

KrakenRecord parse_kraken_record(const string &line, uint64_t line_no) {
  const size_t t1 = line.find('\t');
  const size_t t2 = t1 == string::npos ? string::npos : line.find('\t', t1 + 1);
  const size_t t3 = t2 == string::npos ? string::npos : line.find('\t', t2 + 1);
  if (t1 == string::npos || t2 == string::npos || t3 == string::npos)
    errx(EX_DATAERR, "malformed Kraken output line %llu",
         static_cast<unsigned long long>(line_no));

  KrakenRecord rec;
  const string status = line.substr(0, t1);
  if (status != "C" && status != "U")
    errx(EX_DATAERR, "invalid Kraken status on line %llu: %s",
         static_cast<unsigned long long>(line_no), status.c_str());
  rec.classified = (status == "C");
  rec.read_id = line.substr(t1 + 1, t2 - t1 - 1);
  const string taxid_field = line.substr(t2 + 1, t3 - t2 - 1);
  rec.external_taxid = rec.classified ? parse_external_taxid(taxid_field, line_no) : 0;
  return rec;
}

void flush_buffers(BinWriter &writer, map<taxid_t, BinBlock> &buffers) {
  if (buffers.empty())
    return;
  vector<BinBlock> blocks;
  blocks.reserve(buffers.size());
  for (auto &kv : buffers)
    blocks.push_back(std::move(kv.second));
  writer.WriteBlocks(blocks);
  buffers.clear();
}

void append_fragment(map<taxid_t, BinBlock> &buffers,
                     taxid_t target,
                     Sequence &seq1,
                     Sequence *seq2) {
  BinBlock &block = buffers[target];
  block.taxon = target;
  block.mate1_data += seq1.to_string();
  block.fragments++;
  block.reads++;
  block.bases += seq1.seq.size();
  if (seq2 != nullptr) {
    block.mate2_data += seq2->to_string();
    block.reads++;
    block.bases += seq2->seq.size();
  }
}

}  // namespace

int main(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);
  Options opts = parse_options(argc, argv);

  Taxonomy taxonomy(opts.taxonomy_filename, true);
  taxonomy.GenerateExternalToInternalIDMap();
  RankBinResolver resolver(taxonomy, opts.target_rank);
  BinSelector selector(taxonomy, opts.include_taxids, opts.exclude_taxids);
  BinWriter writer(taxonomy, opts.output_directory, resolver.target_rank(),
                   opts.max_open_files, opts.compression,
                   opts.compression_level, opts.min_fragments,
                   opts.write_unresolved);

  ifstream kraken(opts.kraken_filename.c_str());
  if (!kraken)
    err(EX_NOINPUT, "unable to open Kraken output %s", opts.kraken_filename.c_str());

  map<taxid_t, BinBlock> buffers;
  size_t buffered_fragments = 0;
  uint64_t line_no = 0;
  string line;

  auto process_fragment = [&](Sequence &seq1, Sequence *mate2) {
    if (!std::getline(kraken, line))
      errx(EX_DATAERR, "sequence input contains more fragments than Kraken output");
    ++line_no;
    if (line.empty())
      errx(EX_DATAERR, "blank line in Kraken output at line %llu",
           static_cast<unsigned long long>(line_no));

    const KrakenRecord rec = parse_kraken_record(line, line_no);
    if (!ids_match(rec.read_id, seq1.header))
      errx(EX_DATAERR,
           "classification/read order mismatch at Kraken line %llu: expected '%s', observed '%s'",
           static_cast<unsigned long long>(line_no), rec.read_id.c_str(), seq1.header.c_str());
    if (mate2 != nullptr && !ids_match(seq1.header, mate2->header))
      errx(EX_DATAERR,
           "paired read identifier mismatch at Kraken line %llu: '%s' vs '%s'",
           static_cast<unsigned long long>(line_no), seq1.header.c_str(), mate2->header.c_str());

    writer.Initialize(seq1.format, mate2 != nullptr);
    if (mate2 != nullptr && mate2->format != seq1.format)
      errx(EX_DATAERR, "mates use different sequence formats at Kraken line %llu",
           static_cast<unsigned long long>(line_no));

    if (rec.classified) {
      const taxid_t internal = taxonomy.GetInternalID(rec.external_taxid);
      if (internal == 0 && rec.external_taxid != 0)
        errx(EX_DATAERR,
             "classified taxid %llu from Kraken line %llu is absent from the supplied taxonomy",
             static_cast<unsigned long long>(rec.external_taxid),
             static_cast<unsigned long long>(line_no));
      const taxid_t target = resolver.Resolve(internal);
      if (target == 0) {
        if (writer.write_unresolved())
          append_fragment(buffers, 0, seq1, mate2);
      } else if (selector.Allows(target)) {
        append_fragment(buffers, target, seq1, mate2);
      }
    }

    if (++buffered_fragments >= DEFAULT_FLUSH_FRAGMENTS) {
      flush_buffers(writer, buffers);
      buffered_fragments = 0;
    }
  };

  if (opts.paired) {
    for (size_t file_idx = 0; file_idx < opts.filenames.size(); file_idx += 2) {
      BatchSequenceReader reader1(opts.filenames[file_idx].c_str());
      BatchSequenceReader reader2(opts.filenames[file_idx + 1].c_str());
      Sequence seq1, seq2;
      while (reader1.NextSequence(seq1)) {
        if (!reader2.NextSequence(seq2))
          errx(EX_DATAERR, "mate 2 input ended before mate 1 input");
        process_fragment(seq1, &seq2);
      }
      if (reader2.NextSequence(seq2))
        errx(EX_DATAERR, "mate 2 input contains more reads than mate 1 input");
    }
  } else if (opts.interleaved) {
    for (const string &filename : opts.filenames) {
      BatchSequenceReader reader(filename.c_str());
      Sequence seq1, seq2;
      while (reader.NextSequence(seq1)) {
        if (!reader.NextSequence(seq2))
          errx(EX_DATAERR, "interleaved input contains an incomplete final pair");
        process_fragment(seq1, &seq2);
      }
    }
  } else {
    for (const string &filename : opts.filenames) {
      BatchSequenceReader reader(filename.c_str());
      Sequence seq1;
      while (reader.NextSequence(seq1))
        process_fragment(seq1, nullptr);
    }
  }

  if (std::getline(kraken, line))
    errx(EX_DATAERR, "Kraken output contains more fragments than sequence input (first extra line %llu)",
         static_cast<unsigned long long>(line_no + 1));

  flush_buffers(writer, buffers);
  writer.Finalize();
  return 0;
}
