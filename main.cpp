#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include <vector>
#include <map>
#include <set>

#include "RNAlandmap.h"

extern "C" {
  #include "RNAlandmap_cmdline.h"
}

using namespace std;

int main(int argc, char **argv)
{
  // clock
  time_t run_time = clock();

  // parse arguments
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info) != 0) {
    fprintf(stderr, "Argument parsing problem.\n");
    exit(EXIT_FAILURE);
  }

  // check arguments
  if (args_info.threshold_given && args_info.threshold_arg<0) {
    args_info.threshold_arg=0;
    fprintf(stderr, "Warning: --threshold(-n) set to negative -- clipping to 0.");
  }

  // output;
  vector<saddle> output;

  // map number -> number of struct (for minima)
  map<int, int> map_num;

  // set options
  SetOpt(args_info);

  // main function
  DoTheJob(map_num, output);

  // print output to stdout
  PrintOutput(map_num, output);

  // print dot file
  bool print_dot = true;
  if (print_dot) {
    PrintDot(args_info, map_num, output);
  }

  // free stuff
  FreeStuff();
  cmdline_parser_free(&args_info);

  // time?
  fprintf(stderr, "run_time: %.2f secs.\n", (clock() - run_time)/(double)CLOCKS_PER_SEC);

  return 0;
}
