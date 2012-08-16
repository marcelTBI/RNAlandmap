#include <vector>
#include <map>
#include <set>

#include "move_set.h"

extern "C" {
  #include "RNAlandmap_cmdline.h"
}

using namespace std;

// structure on output
struct saddle {
  short *str;
  int num;
  int energy;

  // connects
  set<set<int> > locmin_conn;
  set<int> saddle_conn;
};

// main functions
void SetOpt(bool debug, bool noLP, bool shifts);  // sets options and inits some things
int DoTheJob(map<int, int> &map_num, vector<saddle> &output);  // main job
void FreeStuff();

// print functions
void PrintDot(gengetopt_args_info &args_info, map<int, int> &map_num, vector<saddle> &output);
void PrintOutput(map<int, int> &map_num, vector<saddle> &output);
