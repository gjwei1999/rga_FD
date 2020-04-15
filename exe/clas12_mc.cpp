#include "clas12_mc.hpp"
#include <future>
#include <thread>

int main(int argc, char** argv) {
  // Need this to make sure root doesn't break
  ROOT::EnableThreadSafety();

  int NUM_THREADS = 4;
  if (getenv("NUM_THREADS") != NULL) NUM_THREADS = atoi(getenv("NUM_THREADS"));
  if (NUM_THREADS > argc - NUM_THREADS) NUM_THREADS = 1;

  // Make a vector of vectors of strings the size of the number of threads
  std::vector<std::vector<std::string>> infilenames(NUM_THREADS);
  // Get the output file name
  std::string outfilename;

  if (argc >= 2) {
    // First argument is the output file
    outfilename = argv[1];
    // All other files are split evently by the under of threads
    for (int i = 2; i < argc; i++) infilenames[i % NUM_THREADS].push_back(argv[i]);
  } else {
    return 1;
  }

  // Make a set of threads (Futures are special threads which return a value)
  std::future<size_t> threads[NUM_THREADS];

  // Define events to be used to get Hz later
  size_t events = 0;

  // Make your histograms object as a shared pointer that all the threads will have
  auto hists = std::make_shared<Histogram>(outfilename);

  // Start timer
  auto start = std::chrono::high_resolution_clock::now();
  // For each thread
  for (size_t i = 0; i < NUM_THREADS; i++) {
    // Set the thread to run a task A-Syncroisly
    // The function we run is the first argument (run_files)
    // The functions areruments are all the remaining arguments
    threads[i] = std::async(run_files, infilenames.at(i), hists, i);
  }

  // For each thread
  for (size_t i = 0; i < NUM_THREADS; i++) {
    // Get the information from the thread in this case how many events each thread actually computed
    events += threads[i].get();
  }

  // Timer and Hz calculator functions that print at the end
  std::cout.imbue(std::locale(""));  // Puts commas in
  std::chrono::duration<double> elapsed_full = (std::chrono::high_resolution_clock::now() - start);
  std::cout << RED << elapsed_full.count() << " sec" << DEF << std::endl;
  std::cout << BOLDYELLOW << events / elapsed_full.count() << " Hz" << DEF << std::endl;
  return 0;
}
