#pragma once

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <vector>
#include <memory>
#include <string>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <ctime>
#include <chrono>
#include <omp.h>
#include <cstdlib>

#define MAX_LINES 40000000

#define GAP_COST -2
#define MISMATCH_COST -3

#define SPACER1 "GTGGCCGATGTTTCGCATCGGCGTACGACT"
#define SPACER2 "ATCCACGTGCTTGAGAGGCCAGAGCATTCG"

#define READS_FILE "mouseDNA/mouseDNA"
#define RESULTS_FILE "results/results"
#define BARCODES_FILE "cellBarcodes"

std::string FIND = "NNNNNNNNNNNNNNNNNN" + std::string(SPACER1) + std::string("NNNNNNNN") + std::string(SPACER2) + std::string("NNNNNNNN");

int SCORES[5][5] = {
  {5, MISMATCH_COST, MISMATCH_COST, MISMATCH_COST, 5},
  {MISMATCH_COST, 5, MISMATCH_COST, MISMATCH_COST, 5},
  {MISMATCH_COST, MISMATCH_COST, 5, MISMATCH_COST, 5},
  {MISMATCH_COST, MISMATCH_COST, MISMATCH_COST, 5, 5},
  {5, 5, 5, 5, 5}
};

std::unordered_map<char, int> INDEXES = {
  {'A', 0},
  {'C', 1},
  {'G', 2},
  {'T', 3},
  {'N', 4}
};

std::unordered_map<char, char> COMPLIMENTS = {
  {'A', 'T'},
  {'C', 'G'},
  {'G', 'C'},
  {'T', 'A'},
  {'N', 'N'}
};

struct ReadInData {
  std::vector<std::pair<std::string, std::string>> reads;
  std::vector<std::vector<std::string>> barcodes;
};

struct AlignmentData {
  int **table;
  int maxI;
  int maxJ;
  int maxValue;
};

struct TracebackData {
  std::string aligned1;
  std::string aligned2;
  int startI;
  int startJ;
};

struct CollapseData {
  std::unordered_map<std::string, std::vector<std::string>> barcodeToUMIMap;
};

struct UMINode {
  std::string umi;
  std::vector<std::string> similars;
};

struct UMINodeCompare {
    bool operator() (const std::shared_ptr<UMINode> lhs, const std::shared_ptr<UMINode> rhs) const {
      if (lhs->similars.size() != rhs->similars.size()) {
        return lhs->similars.size() > rhs->similars.size();
      } else {
        return lhs->umi < rhs->umi;
      }
    }
};

struct BaseCount {
  int a, c, g, t;
  BaseCount(std::string &seq) {
    a = 0;
    c = 0;
    g = 0;
    t = 0;
    for (unsigned int i = 0; i < seq.length(); i++) {
      if (seq.at(i) == 'A') {
        a++;
      } else if (seq.at(i) == 'C') {
        c++;
      } else if (seq.at(i) == 'G') {
        g++;
      } else {
        t++;
      }
    }
  }
  ~BaseCount() {}
  bool isSimilar(BaseCount other) {
    return abs(a - other.a) + abs(c - other.c) + abs(g - other.g) + abs(t - other.t) < 4;
  }
};

void FindBarcodes(int index);

void CollapseUMIs();

std::shared_ptr<ReadInData> ReadFiles(int index);

std::shared_ptr<CollapseData> ReadResultsFiles();

int GetValue(char c1, char c2);

std::shared_ptr<AlignmentData> Align(const std::string &seq1, const std::string &seq2);

std::shared_ptr<TracebackData> Traceback(const std::string &seq1, const std::string &seq2, std::shared_ptr<AlignmentData> alignData);

std::pair<std::string, std::string> GetAlignmentString(std::shared_ptr<TracebackData> traceData);

std::string ReverseCompliment(const std::string &seq);

void FreeTable(std::shared_ptr<AlignmentData> alignData, unsigned int length);
