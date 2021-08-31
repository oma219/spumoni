/* matching_statistics - Computes the matching statistics from BWT and Thresholds
    Copyright (C) 2020 Massimiliano Rossi
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*!
   \file matching_statistics.cpp
   \brief matching_statistics.cpp Computes the matching statistics from BWT and Thresholds.
   \author Massimiliano Rossi
   \date 13/07/2020
*/

extern "C" {
#include <xerrors.h>
}

#include <iostream>

#define VERBOSE

#include <common.hpp>

#include <sdsl/io.hpp>

#include <spumoni.hpp>

#include <malloc_count.h>

////////////////////////////////////////////////////////////////////////////////
/// kseq extra
////////////////////////////////////////////////////////////////////////////////

static inline size_t ks_tell(kseq_t *seq)
{
  return gztell(seq->f->f) - seq->f->end + seq->f->begin;
}

void copy_kstring_t(kstring_t &l, kstring_t &r)
{
  l.l = r.l;
  l.m = r.m;
  l.s = (char *)malloc(l.m);
  for (size_t i = 0; i < r.m; ++i)
    l.s[i] = r.s[i];
}
void copy_kseq_t(kseq_t *l, kseq_t *r)
{
  copy_kstring_t(l->name, r->name);
  copy_kstring_t(l->comment, r->comment);
  copy_kstring_t(l->seq, r->seq);
  copy_kstring_t(l->qual, r->qual);
  l->last_char = r->last_char;
}

////////////////////////////////////////////////////////////////////////////////
/// Parallel computation
////////////////////////////////////////////////////////////////////////////////

// This should be done using buffering.
size_t next_start_fastq(gzFile fp)
{
  int c;
  // Special case when we arr at the beginning of the file.
  if ((gztell(fp) == 0) && ((c = gzgetc(fp)) != EOF) && c == '@')
    return 0;

  // Strart from the previous character
  gzseek(fp, -1, SEEK_CUR);

  std::vector<std::pair<int, size_t>> window;
  // Find the first new line
  for (size_t i = 0; i < 4; ++i)
  {
    while (((c = gzgetc(fp)) != EOF) && (c != (int)'\n'))
    {
    }
    if (c == EOF)
      return gztell(fp);
    if ((c = gzgetc(fp)) == EOF)
      return gztell(fp);
    window.push_back(std::make_pair(c, gztell(fp) - 1));
  }

  for (size_t i = 0; i < 2; ++i)
  {
    if (window[i].first == '@' && window[i + 2].first == '+')
      return window[i].second;
    if (window[i].first == '+' && window[i + 2].first == '@')
      return window[i + 2].second;
  }

  return gztell(fp);
}

// test if the file is gzipped
static inline bool is_gzipped(std::string filename)
{
  FILE *fp = fopen(filename.c_str(), "rb");
  if(fp == NULL) error("Opening file " + filename);
  int byte1 = 0, byte2 = 0;
  fread(&byte1, sizeof(char), 1, fp);
  fread(&byte2, sizeof(char), 1, fp);
  fclose(fp);
  return (byte1 == 0x1f && byte2 == 0x8b);
}

// Return the length of the file
// Assumes that the file is not compressed
static inline size_t get_file_size(std::string filename)
{
  if (is_gzipped(filename))
  {
    std::cerr << "The input is gzipped!" << std::endl;
    return -1;
  }
  FILE *fp = fopen(filename.c_str(), "r");
  fseek(fp, 0L, SEEK_END);
  size_t size = ftell(fp);
  fclose(fp);
  return size;
}

std::vector<size_t> split_fastq(std::string filename, size_t n_threads)
{
  //Precondition: the file is not gzipped
  // scan file for start positions and execute threads
  size_t size = get_file_size(filename);

  gzFile fp = gzopen(filename.c_str(), "r");
  if (fp == Z_NULL)
  {
    throw new std::runtime_error("Cannot open input file " + filename);
  }

  std::vector<size_t> starts(n_threads + 1);
  for (int i = 0; i < n_threads + 1; ++i)
  {
    size_t start = (size_t)((size * i) / n_threads);
    gzseek(fp, start, SEEK_SET);
    starts[i] = next_start_fastq(fp);
  }
  gzclose(fp);
  return starts;
}

class ms_t
{
public:

  ms_t(std::string filename)
  {
    std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
    std::string filename_ms = filename + ms.get_file_extension();

    ifstream fs_ms(filename_ms);
    ms.load(fs_ms);
    fs_ms.close();

    std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();
    verbose("Matching Statistics Index Construction Complete");
    
    /* Commented these out since they are printed in main() */
    //verbose("Memory peak: ", malloc_count_peak());
    //verbose("Elapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());
}

  // Destructor
  ~ms_t() 
  {
      // NtD
  }

  // The outfile has the following format. The first size_t integer store the
  // length l of the query. Then the following l size_t integers stores the
  // pointers of the matching statistics, and the following l size_t integers
  // stores the lengths of the mathcing statistics.
  //void matching_statistics(kseq_t *read, FILE* out) 
  void matching_statistics(const char* read, size_t read_length, FILE* out) 
  {
    auto lengths = ms.query(read, read_length);

    size_t q_length = lengths.size();
    fwrite(&q_length, sizeof(size_t), 1,out);
    fwrite(lengths.data(), sizeof(size_t),q_length,out);
  }

protected:
  ms_pointers<> ms;
  size_t n = 0;
};



char complement(char n)
{
  switch (n)
  {
  case 'A':
    return 'T';
  case 'T':
    return 'A';
  case 'G':
    return 'C';
  case 'C':
    return 'G';
  default:
    return n;
  }
}

typedef struct{
  // Parameters
  ms_t *ms;
  std::string pattern_filename;
  std::string out_filename;
  size_t start;
  size_t end;
  size_t wk_id;
} mt_param;

void *mt_ms_worker(void *param)
{
  mt_param *p = (mt_param*) param;
  size_t n_reads = 0;
  size_t n_aligned_reads = 0;

  FILE *out_fd;
  gzFile fp;

  if ((out_fd = fopen(p->out_filename.c_str(), "w")) == nullptr)
    error("open() file " + p->out_filename + " failed");

  if ((fp = gzopen(p->pattern_filename.c_str(), "r")) == Z_NULL)
    error("open() file " + p->pattern_filename + " failed");

  gzseek(fp, p->start, SEEK_SET);

  kseq_t rev;
  int l;

  kseq_t *seq = kseq_init(fp);
  while ((ks_tell(seq) < p->end) && ((l = kseq_read(seq)) >= 0)) {
    std::string curr_read = std::string(seq->seq.s);
    transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); //Make sure all characters are upper-case

    p->ms->matching_statistics(curr_read.c_str(), seq->seq.l, out_fd);
  }

  kseq_destroy(seq);
  gzclose(fp);
  fclose(out_fd);

  return NULL;
}

void mt_ms(ms_t *ms, std::string pattern_filename, std::string out_filename, size_t n_threads)
{
  pthread_t t[n_threads] = {0};
  mt_param params[n_threads];
  std::vector<size_t> starts = split_fastq(pattern_filename, n_threads);
  for(size_t i = 0; i < n_threads; ++i)
  {
    params[i].ms = ms;
    params[i].pattern_filename = pattern_filename;
    params[i].out_filename = out_filename + "_" + std::to_string(i) + ".ms.tmp.out";
    params[i].start = starts[i];
    params[i].end = starts[i+1];
    params[i].wk_id = i;
    xpthread_create(&t[i], NULL, &mt_ms_worker, &params[i], __LINE__, __FILE__);
  }

  for(size_t i = 0; i < n_threads; ++i)
  {
    xpthread_join(t[i],NULL,__LINE__,__FILE__);
  }

  sleep(5);


  return;
}


////////////////////////////////////////////////////////////////////////////////
/// Single Thread
////////////////////////////////////////////////////////////////////////////////

size_t st_ms(ms_t *ms, std::string pattern_filename, std::string out_filename)
{
  size_t n_reads = 0;
  size_t n_aligned_reads = 0;
  kseq_t rev;
  int l;
  FILE *out_fd;

  out_filename += "_0.ms.tmp.out";

  if ((out_fd = fopen(out_filename.c_str(), "w")) == nullptr)
    error("open() file " + out_filename + " failed");

  gzFile fp = gzopen(pattern_filename.c_str(), "r");
  kseq_t* seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    std::string curr_read = std::string(seq->seq.s);
    transform(curr_read.begin(), curr_read.end(), curr_read.begin(), ::toupper); //Make sure all characters are upper-case

    ms->matching_statistics(curr_read.c_str(), seq->seq.l, out_fd);
  }

  kseq_destroy(seq);
  gzclose(fp);
  fclose(out_fd);

  return n_aligned_reads;
}

size_t st_ms_general(ms_t *ms, std::string pattern_filename, std::string out_filename)
{
  // Check if file is mistakenly labeled as non-fasta file
  std::string file_ext = pattern_filename.substr(pattern_filename.find_last_of(".") + 1);
  if (file_ext == "fa" || file_ext == "fasta") {
    error("The file extension for the patterns suggests it is a fasta file. Please run with -f option for correct results.");
  }

  size_t n_reads = 0;
  size_t n_aligned_reads = 0;
  kseq_t rev;
  int l;
  FILE *out_fd;
  out_filename += "_0.ms.tmp.out";

  if ((out_fd = fopen(out_filename.c_str(), "w")) == nullptr)
    error("open() file " + out_filename + " failed");

  std::ifstream input_fd (pattern_filename, std::ifstream::in | std::ifstream::binary);
  char ch = input_fd.get();
  char* buf = new char [2]; // To just hold separator character
  std::string read = "";

  while (input_fd.good()) {
    if (ch == '\x01') {ms->matching_statistics(read.c_str(), read.size(), out_fd); input_fd.read(buf, 2); read = "";}
    else {read += ch;}
    ch = input_fd.get();
  }
  fclose(out_fd);

  return n_aligned_reads;
}

typedef std::pair<std::string, std::vector<uint8_t>> pattern_t;

int main(int argc, char *const argv[]){

  // Parse arguments (struct info can be found in common.hpp)
  Args args;
  parseArgs(argc, argv, args);


  // Loads the RLEBWT and Thresholds ------------------------------------------------------------------------
  verbose("Construction of the matching statistics data structure");
  std::chrono::high_resolution_clock::time_point t_insert_start = std::chrono::high_resolution_clock::now();
  ms_t ms(args.filename);
  std::chrono::high_resolution_clock::time_point t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("\tMemory peak: ", malloc_count_peak());
  verbose("\tElapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());


  // Load patterns and generate PMLs -------------------------------------------------------------------------
  verbose("Processing patterns");
  t_insert_start = std::chrono::high_resolution_clock::now();
  
  std::string base_name = basename(args.filename.data());
  std::string out_filename = args.patterns + "_" + base_name;

  if (is_gzipped(args.patterns)) {
    verbose("The input is gzipped - forcing single threaded pseudo matching lengths.");
    args.th = 1;
  }

  // Determine approach to parse pattern files (based on whether it is general text or fasta format)
  if (args.is_fasta) {
    if(args.th == 1) {st_ms(&ms, args.patterns, out_filename);}
    else {mt_ms(&ms, args.patterns, out_filename, args.th);}
  }
  else {
        if(args.th == 1) {st_ms_general(&ms,args.patterns,out_filename);}
    else {error("Multi-threading not implemented yet for general-text querying.\n"); std::exit(1);}
  }

  t_insert_end = std::chrono::high_resolution_clock::now();

  verbose("\tMemory peak: ", malloc_count_peak());
  verbose("\tElapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  // Writing out the PMLs -------------------------------------------------------------------------
  verbose("Printing plain output");
  t_insert_start = std::chrono::high_resolution_clock::now();
  
  std::ofstream f_lengths(out_filename + ".pseudo_lengths");

  if (!f_lengths.is_open())
    error("open() file " + std::string(out_filename) + ".lengths failed");

  size_t n_seq = 0;
  for(size_t i = 0; i < args.th; ++i)
  {
    std::string tmp_filename = out_filename + "_" + std::to_string(i) + ".ms.tmp.out";
    FILE *in_fd;

    if ((in_fd = fopen(tmp_filename.c_str(), "r")) == nullptr)
      error("open() file " + tmp_filename + " failed");

    size_t length = 0;
    size_t m = 100; // Reserved size for pointers and lengths
    size_t *mem = (size_t*) malloc(m * sizeof(size_t));
    while(!feof(in_fd) and fread(&length,sizeof(size_t), 1,in_fd) > 0)
    {
      if( m < length)
      {
        // Resize lengths and pointers
        m = length;
        mem = (size_t*) realloc(mem, m * sizeof(size_t));
      }
        
      if ((fread(mem, sizeof(size_t), length, in_fd)) != length)
        error("fread() file " + std::string(tmp_filename) + " failed");

      f_lengths << ">" + std::to_string(n_seq) << endl;
      for(size_t i = 0; i < length; ++i)
        f_lengths  << mem[i] << " ";
      f_lengths << endl;

      n_seq ++;
    }
    fclose(in_fd);
  }

  f_lengths.close();

  t_insert_end = std::chrono::high_resolution_clock::now();
  
  auto mem_peak = malloc_count_peak();
  verbose("\tMemory peak: ", malloc_count_peak());
  verbose("\tElapsed time (s): ", std::chrono::duration<double, std::ratio<1>>(t_insert_end - t_insert_start).count());

  //auto mem_peak = malloc_count_peak();
  //verbose("Memory peak: ", malloc_count_peak());

  size_t space = 0;
  if (args.memo) {}
  if (args.store) {}

  if (args.csv)
    std::cerr << csv(args.filename.c_str(), time, space, mem_peak) << std::endl;

  return 0;
}
