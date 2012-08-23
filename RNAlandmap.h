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

// structure with all info about each structure;
struct struct_info {
  set<int> LM_nums;
  set<int> saddle_nums;
  int debug_num;
};

// energy and structure in short* format
struct en_struct {
  int energy;
  short *str;
};

// main functions
void SetOpt(gengetopt_args_info &args_info);  // sets options
int DoTheJob(map<int, int> &map_num, vector<saddle> &output);  // main job
void FreeStuff();

// print functions
void PrintDot(gengetopt_args_info &args_info, map<int, int> &map_num, vector<saddle> &output);
void PrintOutput(map<int, int> &map_num, vector<saddle> &output);
