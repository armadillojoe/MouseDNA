#include "MouseDNA.h"

int main(int argc, char **argv) {
  if (argc == 2) {
    if (std::string(argv[1]) == "collapse") {
      CollapseUMIs();
    } else {
      int index = std::stoi(argv[1]);
      std::ifstream reads(READS_FILE + std::to_string(index));
      if (reads.good()) {
        reads.close();
        FindBarcodes(index);
      } else {
        reads.close();
        std::cout << "Invalid File Number" << std::endl;
      }
    }
  } else {
    std::cout << "Usage: ./MouseDNA <File Index Number>" << std::endl;
  }
}

void FindBarcodes(int index) {
  std::shared_ptr<ReadInData> readData = ReadFiles(index);
  auto start = std::chrono::system_clock::now();
  #pragma omp parallel for
  for (auto p = readData->reads.begin(); p < readData->reads.end(); p++) {
    std::string seq1 = FIND;
    std::string seq2 = p->first.substr(0, 200);
    std::string seq2Reverse = ReverseCompliment(p->first).substr(0, 200);
    std::string correctSeq = seq2;
    std::shared_ptr<AlignmentData> alignData = Align(seq1, seq2);
    std::shared_ptr<AlignmentData> correctAlignData = alignData;
    std::shared_ptr<AlignmentData> alignDataReverse = Align(seq1, seq2Reverse);
    if (alignDataReverse->maxValue > alignData->maxValue) {
      correctSeq = seq2Reverse;
      correctAlignData = alignDataReverse;
    }
    std::shared_ptr<TracebackData> traceData = Traceback(seq1, correctSeq, correctAlignData);
    std::pair<std::string, std::string> codesAndAlignment = GetAlignmentString(traceData);
    if (correctAlignData->maxValue > 370) {
      std::cout << codesAndAlignment.second << codesAndAlignment.first << std::endl;
    }
    FreeTable(alignData, seq1.length() + 1);
    FreeTable(alignDataReverse, seq1.length() + 1);
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed Time: " << elapsed_seconds.count() << std::endl;
}

void CollapseUMIs() {
  auto start = std::chrono::system_clock::now();
  std::shared_ptr<CollapseData> collapseData = ReadResultsFiles();
  for (auto p : collapseData->barcodeToUMIMap) {
    std::cout << p.first + ":" + std::to_string(p.second.size()) + ":";
    std::unordered_set<std::string> umis;
    std::vector<std::string> uniqueUmis;
    for (unsigned int i = 0; i < p.second.size(); i++) {
      umis.insert(p.second[i]);
    }
    std::copy(umis.begin(), umis.end(), std::back_inserter(uniqueUmis));
    std::unordered_map<std::string, std::vector<std::string>> graph;
    for (unsigned int i = 0; i < uniqueUmis.size(); i++) {
      std::vector<std::string> vec;
      graph[uniqueUmis[i]] = vec;
    }
    omp_lock_t writelock;
    omp_init_lock(&writelock);
    #pragma omp parallel for
    for (unsigned int i = 0; i < uniqueUmis.size(); i++) {
      int maxValue = 0;
      std::string maxUMI = "";
      for (unsigned int j = i + 1; j < uniqueUmis.size(); j++) {
        std::string seq1 = uniqueUmis[i];
        std::string seq2 = uniqueUmis[j];
        std::shared_ptr<BaseCount> c1 = std::make_shared<BaseCount>(seq1);
        std::shared_ptr<BaseCount> c2 = std::make_shared<BaseCount>(seq2);
        if (c1->isSimilar(*c2)) {
          std::shared_ptr<AlignmentData> alignData = Align(seq1, seq2);
          std::shared_ptr<TracebackData> traceData = Traceback(seq1, seq2, alignData);
          std::pair<std::string, std::string> codesAndAlignment = GetAlignmentString(traceData);
          if (alignData->maxValue > 40 && alignData->maxValue > maxValue) {
            maxValue = alignData->maxValue;
            maxUMI = seq2;
          }
          FreeTable(alignData, seq1.length() + 1);
        }
      }
      if (maxUMI != "") {
        omp_set_lock(&writelock);
        graph[uniqueUmis[i]].push_back(maxUMI);
        graph[maxUMI].push_back(uniqueUmis[i]);
        omp_unset_lock(&writelock);
      }
    }
    omp_destroy_lock(&writelock);
    std::set<std::shared_ptr<UMINode>, UMINodeCompare> nodes;
    for (auto p : graph) {
      std::shared_ptr<UMINode> newNode = std::make_shared<UMINode>();
      newNode->umi = p.first;
      newNode->similars = p.second;
      nodes.insert(newNode);
    }
    std::unordered_set<std::string> duplicates;
    std::set<std::string> reals;
    for (auto n : nodes) {
      if (duplicates.count(n->umi) == 0) {
        reals.insert(n->umi);
        for (auto &u : n->similars) {
          duplicates.insert(u);
        }
      }
    }
    std::cout << reals.size() << std::endl;
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  std::cout << "Elapsed Time: " << elapsed_seconds.count() << std::endl;
}

void FreeTable(std::shared_ptr<AlignmentData> alignData, unsigned int length) {
  for (unsigned int i = 0; i < length; i++) {
    delete[] alignData->table[i];
  }
  delete[] alignData->table;
}

std::string ReverseCompliment(const std::string &seq) {
  std::string rev = "";
  for (unsigned int i = 0; i < seq.length(); i++) {
    rev = COMPLIMENTS[seq.at(i)] + rev;
  }
  return rev;
}

std::pair<std::string, std::string> GetAlignmentString(std::shared_ptr<TracebackData> traceData) {
  unsigned int maxLineLength = 700;
  std::string output = "";
  std::string commons = "";
  std::string codes = "";

  for (unsigned int i = 0; i < traceData->aligned1.length(); i++) {
    if (traceData->aligned1.at(i) == '-' || traceData->aligned2.at(i) == '-') {
      commons += " ";
    } else if (traceData->aligned1.at(i) == traceData->aligned2.at(i)) {
      commons += traceData->aligned1.at(i);
    } else if (GetValue(traceData->aligned1.at(i), traceData->aligned2.at(i)) > 0) {
      codes += traceData->aligned2[i];
      commons += "|";
    } else {
      commons += " ";
    }
  }

  unsigned int rows = commons.length() / maxLineLength;

  for (unsigned int i = 0; i < rows + 1; i++) {
    output += std::to_string(traceData->startJ) + "\t";
    for (unsigned int j = i * maxLineLength; j < (i + 1) * maxLineLength; j++) {
      if (j < traceData->aligned1.length()) {
        char curr = traceData->aligned1.at(j);
        if (curr == '-') {
          traceData->startJ++;
        }
        output += curr;
      }
    }
    output += "\n\t";
    for (unsigned int j = i * maxLineLength; j < (i + 1) * maxLineLength; j++) {
      if (j < commons.length()) {
        output += commons.at(j);
      }
    }
    output += "\n" + std::to_string(traceData->startI) + "\t";
    for (unsigned int j = i * maxLineLength; j < (i + 1) * maxLineLength; j++) {
      if (j < traceData->aligned2.length()) {
        char curr = traceData->aligned2.at(j);
        if (curr == '-') {
          traceData->startI++;
        }
        output += curr;
      }
    }
    output += "\n";
  }

  std::pair<std::string, std::string> ret;
  ret.first = codes;
  ret.second = output;
  return ret;
}

std::shared_ptr<TracebackData> Traceback(const std::string &seq1, const std::string &seq2, std::shared_ptr<AlignmentData> alignData) {
  std::shared_ptr<TracebackData> data = std::make_shared<TracebackData>();

  data->aligned1 = "";
  data->aligned2 = "";

  int currI = alignData->maxI;
  int currJ = alignData->maxJ;

  while (alignData->table[currJ][currI] > 0) {
    int left = alignData->table[currJ][currI - 1] + GAP_COST;
    int diagonal = alignData->table[currJ - 1][currI - 1] + GetValue(seq1.at(currJ - 1), seq2.at(currI - 1));
    int top = alignData->table[currJ - 1][currI] + GAP_COST;

    if (top >= left && top >= diagonal) {
      data->aligned2 = "-" + data->aligned2;
      data->aligned1 = seq1.at(currJ - 1) + data->aligned1;
      currJ--;
    } else if (left >= top && left >= diagonal) {
      data->aligned1 = "-" + data->aligned1;
      data->aligned2 = seq2.at(currI - 1) + data->aligned2;
      currI--;
    } else if (diagonal >= top && diagonal >= left) {
      data->aligned1 = seq1.at(currJ - 1) + data->aligned1;
      data->aligned2 = seq2.at(currI - 1) + data->aligned2;
      currI--;
      currJ--;
    }
  }

  data->startI = currI + 1;
  data->startJ = currJ + 1;

  return data;
}

std::shared_ptr<AlignmentData> Align(const std::string &seq1, const std::string &seq2) {
  std::shared_ptr<AlignmentData> data = std::make_shared<AlignmentData>();

  int x = seq1.length() + 1;
  int y = seq2.length() + 1;
  data->table = new int*[x]();
  for (int i = 0; i < x; i++) {
    data->table[i] = new int[y]();
  }

  for (int j = 1; j < x; j++) {
    for (int i = 1; i < y; i++) {
      int left = data->table[j][i - 1] + GAP_COST;
      int diagonal = data->table[j - 1][i - 1] + GetValue(seq1.at(j - 1), seq2.at(i - 1));
      int top = data->table[j - 1][i] + GAP_COST;
      int choices[] = {left, diagonal, top, 0};
      int max = *std::max_element(choices, choices + 4);
      if (max > data->maxValue) {
        data->maxValue = max;
        data->maxI = i;
        data->maxJ = j;
      }
      data->table[j][i] = max;
    }
  }

  return data;
}

int GetValue(char c1, char c2) {
  return SCORES[INDEXES[c1]][INDEXES[c2]];
}

std::shared_ptr<ReadInData> ReadFiles(int index) {
  std::shared_ptr<ReadInData> data = std::make_shared<ReadInData>();

  std::ifstream reads(READS_FILE + std::to_string(index));
  std::string line;
  std::string lastLine;
  int counter = 0;

  while (std::getline(reads, line)) {
    if (counter > MAX_LINES) {
      break;
    }
    if (counter % 4 == 0) {
      lastLine = line;
    }
    if (counter % 4 == 1) {
      data->reads.push_back(std::pair<std::string, std::string>(line, lastLine));
    }
    counter++;
  }

  reads.close();

  std::ifstream barcodes(BARCODES_FILE);

  while (std::getline(barcodes, line)) {
    std::vector<std::string> parts;
    boost::split(parts, line, [](char c){
      return c == ',';
    });
    data->barcodes.push_back(parts);
  }

  barcodes.close();

  return data;
}

std::shared_ptr<CollapseData> ReadResultsFiles() {
  int NUM_FILES = 4;
  std::shared_ptr<CollapseData> data = std::make_shared<CollapseData>();
  std::unordered_map<std::string, std::string> readToBarcodeMap;

  for (int i = 0; i < NUM_FILES; i++) {
    std::ifstream reads(RESULTS_FILE + std::to_string(i) + std::string(".txt"));
    std::string line;
    std::string lastLine;
    int counter = 0;
    while (std::getline(reads, line)) {
      if (counter > 0) {
        if (counter % 2 == 0) {
          lastLine = line;
        }
        if (counter % 2 == 1) {
          if (lastLine != "None" && lastLine.length() > 8) {
            readToBarcodeMap[line] = lastLine;
          }
        }
      }
      counter++;
    }
    reads.close();
  }

  std::vector<std::pair<std::string, std::string>> idReadsMap;

  for (int i = 0; i < NUM_FILES; i++) {
    std::ifstream reads(READS_FILE + std::to_string(i));
    std::string line;
    std::string lastLine;
    int counter = 0;
    while (std::getline(reads, line)) {
      if (counter % 4 == 0) {
        lastLine = line;
      }
      if (counter % 4 == 1) {
        idReadsMap.push_back(std::pair<std::string, std::string>(lastLine, line));
      }
      counter++;
    }
    reads.close();
  }

  omp_lock_t writelock;
  omp_init_lock(&writelock);
  #pragma omp parallel for
  for (auto it = idReadsMap.begin(); it < idReadsMap.end(); it++) {
    if (readToBarcodeMap.find(it->first) != readToBarcodeMap.end()) {
      std::string barcode = readToBarcodeMap[it->first];
      std::string find = "NNNNNNNNNN";
      find += barcode.substr(0, 8);
      find += SPACER1;
      find += barcode.substr(8, 8);
      find += SPACER2;
      find += barcode.substr(16, 8);
      std::string seq2 = it->second;
      std::string seq2Reverse = ReverseCompliment(it->second);
      std::string correctSeq = seq2;
      std::shared_ptr<AlignmentData> alignData = Align(find, seq2);
      std::shared_ptr<AlignmentData> correctAlignData = alignData;
      std::shared_ptr<AlignmentData> alignDataReverse = Align(find, seq2Reverse);
      if (alignDataReverse->maxValue > alignData->maxValue) {
        correctSeq = seq2Reverse;
        correctAlignData = alignDataReverse;
      }
      std::shared_ptr<TracebackData> traceData = Traceback(find, correctSeq, correctAlignData);
      std::pair<std::string, std::string> codesAndAlignment = GetAlignmentString(traceData);
      omp_set_lock(&writelock);
      data->barcodeToUMIMap[barcode.substr(0, 24)].push_back(codesAndAlignment.first);
      omp_unset_lock(&writelock);
      FreeTable(alignData, find.length() + 1);
      FreeTable(alignDataReverse, find.length() + 1);
    }
  }
  omp_destroy_lock(&writelock);

  return data;
}
