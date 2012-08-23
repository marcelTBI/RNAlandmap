#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include "RNAlandmap.h"
#include "hash_util.h"

extern "C" {
  #include "fold.h"
}

// ============= GLOBAL ARRAYS ================= (used also in move_set.cpp)

// map structure - info about the structure
unordered_map<short*, struct_info, hash_fncts, hash_eq> struct_map (HASHSIZE);

// collecting sets of minima/saddles
set<set<int> > numbers;
set<int> saddles;

// ============= LOCAL ARRAYS ==================

// set of degen structures - for quick skipping
set<short*, setcomp> struct_set;

bool DEBUGG;
//map<short*, int, setcomp> debug_map; // map struct to their number
unordered_map<int, en_struct> invert_map; // map num to struct in debug_map/struct_map

// global variables for landmap
int gl_energy;
short *gl_str;
bool gl_direct;
int gl_threshold;

// global set of degeneracy:
set<short*, setcomp> degen_set;
bool degen_lm;
bool gl_degen;

// global options
encoded *enc;
degen deg;
char *seq;

// union-find set array
vector<int> parent;
unsigned int num_unions = 0;

// ===================== LOCAL FUNCTIONS ====================

// and union-find set functions
int find(int x) {
  if (x != parent[x] && parent[x] != parent[parent[x]])
    parent[x] = find(parent[x]);
  return parent[x];
}

void union_set(int x, int y) {
  int u, v;
  u = find(x);
  v = find(y);
  if (u != v) {
    parent[u] = v;
    num_unions++;
  }
}

bool connected_all() {
  return (num_unions == parent.size()-1);
}

bool joint(int x, int y) {
  return find(x) == find(y);
}

void enlarge_parent() {
  parent.push_back(parent.size());
}

// some checkings
inline bool isStruct(char *p)
{
  // check first two chars - should be enough
  if (strlen(p)<2) return false;
  if ((p[0]=='.' || p[0]=='(' || p[0]==')') && (p[1]=='.' || p[1]=='(' || p[1]==')')) return true;
  else return false;
}

inline bool isSeq(char *p)
{
  if (strlen(p)<2) return false;
  // check first two chars - should be enough
  switch (p[0]){
    case 'A':
    case 'C':
    case 'G':
    case 'T':
    case 'U':
    case 'a':
    case 'c':
    case 'g':
    case 't':
    case 'u': switch (p[1]){
      case 'A':
      case 'C':
      case 'G':
      case 'T':
      case 'U':
      case 'a':
      case 'c':
      case 'g':
      case 't':
      case 'u': return true;
    }
    default : return false;
  }
}

// function to do with every structure (always return false - to continue searching)
bool landmap (short *str, int energy)
{
  if (DEBUGG) fprintf (stderr, "tryin %s", pt_to_str(str).c_str());
  if (DEBUGG) fprintf (stderr, " %d\n", struct_map.count(str)>0 ? struct_map[str].debug_num: -1);

  if (energy > gl_energy) return false;
  if (energy == gl_energy) {
    // already done?
    if (struct_map.count(str)>0) {
      return false;
    }

    // find degeneracy (searching all neighbours with equal energy)
      // do only once per struct
    if (!gl_degen) {
      gl_degen = true;
      short *last = enc->pt;
      enc->pt = str;
      degen_set = find_equal_energy(*enc, energy, deg, gl_direct, degen_lm);
      enc->pt = last;
    }
  }

  // lower energy and we don't have it? (threshold or filtering)
  if (struct_map.count(str)==0) {
    //if (DEBUGG) fprintf(stderr, "MISSING!\n");
    return false;
  }

  set<int>::iterator set_it;
  struct_info si = struct_map[str];

  //fprintf (stderr, "doing %s", pt_to_str(str).c_str());
  if (DEBUGG) fprintf(stderr, "%2d: ", si.debug_num);
  if (DEBUGG) for (set_it=si.LM_nums.begin(); set_it!=si.LM_nums.end(); set_it++) fprintf(stderr, " %d", *set_it);
  if (DEBUGG) fprintf(stderr, "\n");


  // collect number sets (collect minima information)
  if (si.LM_nums.size()>0) numbers.insert(si.LM_nums);

  // collect saddle information
  saddles.insert(si.saddle_nums.begin(), si.saddle_nums.end());

  return false;
}

// ================== GLOBAL FUNCTIONS ==================

int DoTheJob(map<int, int> &map_num, vector<saddle> &output)
{
    // read stdin for sequence
  seq = my_getline(stdin);
  strtok(seq, " ");
  if (!isSeq(seq)) {
    if (seq) free(seq);
    return -2;
  }

  // create encoded structure
  enc=encode_seq(seq);

  // read structures and do the stuff...
  char *str;
  int number = 0; //must be from 0 - parents array
  int cnt = 0;
  while ((str = my_getline(stdin)) && (gl_threshold > number || !connected_all())) {  // stopping condition
    float en;
    int energy;
    bool have_en = false;
    bool have_str = false;


    char *p = strtok(str, " ");
    while (p) {
      if (isStruct(p)) {
        encode_str(enc, p);
        have_str = true;
      } else {
        if (!have_en && sscanf(p, "%f", &en)==1) {
          energy = (int) (en*100+(en>0.0?0.5:-0.5));
          have_en = true;
        }
      }
      p = strtok(NULL, " ");
    }

    if (!have_str) {
      free(str);
      break;
    }

    // if degeneracy have found this struct before: (keep the set small for efficiency)
    set<short*, setcomp>::iterator it;
    if ((it = struct_set.find(enc->pt))!=struct_set.end()) {
      struct_set.erase(it);
      if (DEBUGG) fprintf(stderr, "skip: %s %4d\n", pt_to_str(enc->pt).c_str(), cnt);
      if (struct_map.count(enc->pt)==0) fprintf(stderr, "WRONG: %s\n", pt_to_str(enc->pt).c_str());
      struct_map[enc->pt].debug_num = cnt; // we can complete the info
      en_struct &es = invert_map[cnt]; // and also invert_map
      es.str = enc->pt;
      es.energy = gl_energy;
      free(str);
      cnt++;
      continue;
    }

    //------- now we can assume that we have no degeneracy... or first in degenerated set

    // somehow energy_of_move doesnt work without first calling energy_of_structure
    if (1 || !have_en) energy = energy_of_structure_pt(seq, enc->pt, enc->s0, enc->s1, 0);

    // init
    numbers.clear();
    saddles.clear();
    degen_set.clear();
    degen_lm = false;
    gl_energy = energy;
    gl_str = allocopy(enc->pt);
    gl_degen = false;

    // find neighbours and do landmap function on them
    int dunno;
    if (DEBUGG) fprintf(stderr, "---- %d ----\n base %s %6.2f\n", cnt, pt_to_str(gl_str).c_str(), gl_energy/100.0);
    browse_neighs(*enc, energy, deg, dunno);

    // debrief
      // every time:
    if (DEBUGG && struct_map.count(gl_str)>0) fprintf(stderr, "WRONG: %s\n", pt_to_str(gl_str).c_str());
    struct_info &si = struct_map[gl_str];
    si.debug_num = cnt;   // create map entries with empty sets
    //saddle_map[gl_str];
    //debug_map[gl_str] = cnt;
    en_struct &es = invert_map[cnt];
    es.energy = energy;
    es.str = gl_str;
      // number the structure in struct_map and saddle_map
    for (set<set<int> >::iterator it=numbers.begin(); it!=numbers.end(); it++) si.LM_nums.insert(it->begin(), it->end());
    si.saddle_nums.insert(saddles.begin(), saddles.end());

      // minimum
    if (numbers.size() == 0) {
      // must be under threshold to add a new minimum
      if (gl_threshold > number) {
        enlarge_parent();
        map_num[number]=cnt;
        si.LM_nums.insert(number);
        number++;
      }
    } else {
      // maybe saddle point - more numbers joining
      if (numbers.size()>1) {
        // saddle point?
        saddle sad;
        sad.str = gl_str;
        sad.num = cnt;
        sad.energy = energy;

        vector<std::pair<int, int> > to_union;

        // refine numbers:
        set<int> tmp;
        do {
          tmp.clear();
          for (set<set<int> >::iterator it=numbers.begin(); it!=numbers.end(); it++) {
            set<set<int> >::iterator it2=it;
            it2++;
            for (; it2!=numbers.end(); it2++) {
              int first = *(it->begin());
              int second = *(it2->begin());
              if (joint(first, second)) {
                tmp = *it;
                tmp.insert(it2->begin(), it2->end());
                numbers.erase(it);
                numbers.erase(it2);
                break;
              }
            }
            if (tmp.size()!=0) break;
          }
          if (tmp.size()!=0) {
            numbers.insert(tmp);
          }
        } while (tmp.size()!=0);

        // crawl through sets of numbers, check if they are joint yet
        for (set<set<int> >::iterator it=numbers.begin(); it!=numbers.end(); it++) {
          set<set<int> >::iterator it2=it;
          it2++;
          for (; it2!=numbers.end(); it2++) {
            int first = *(it->begin());
            int second = *(it2->begin());
            // saddle point (they are not joint yet)
            if (!joint(first, second)) {
              to_union.push_back(make_pair(first, second));

              // found saddle - fill info
              sad.saddle_conn = si.saddle_nums; // saddle information
              sad.locmin_conn.insert(*it);      // add locmin info
              sad.locmin_conn.insert(*it2);
            }
          }
        }
        // number the saddle and store it
        if (to_union.size()>0) {
          si.saddle_nums.insert(cnt);
          output.push_back(sad);
        }
        // union things
        for (unsigned int i=0; i< to_union.size(); i++) {
          union_set(to_union[i].first, to_union[i].second);
        }
      }
    }

      // degeneracy:
    if (degen_set.size()>0) {
      set<short*, setcomp>::iterator degen_it;
      degen_it = degen_set.find(gl_str);
      free(*degen_it);
      degen_set.erase(degen_it);
      for (set<short*, setcomp>::iterator it=degen_set.begin(); it!=degen_set.end(); it++) {
        // assign, that we have done this structure
        if (DEBUGG) fprintf(stderr, "  deg %s %4d\n", pt_to_str(*it).c_str(), cnt);
        //if (struct_map.find(*it)!=struct_map.end()) fprintf(stderr, "wrong: %s == \n %4d: %s\n", pt_to_str(struct_map.find(*it)->first).c_str(), cnt, pt_to_str(*it).c_str());
        /*unordered_map<short*, struct_info, hash_fncts, hash_eq>::iterator mit;
        mit = struct_map.find(*it);
        if (mit!=struct_map.end() && DEBUGG) {
          fprintf(stderr, "wrong: %s == \n %4d: %s\n", pt_to_str(struct_map.find(*it)->first).c_str(), cnt, pt_to_str(*it).c_str());

        }*/
        struct_info &si2 = struct_map[*it];
        si2.LM_nums = si.LM_nums; // allocated memory is now on struct_map to free!!
        si2.saddle_nums = si.saddle_nums;
        si2.debug_num = -1;  // we dont know its number

        // for quicker access
        struct_set.insert(*it);
        //invert_map[?] = *it;
      }
    }
    free(str);
    cnt++;
  }

  if (str) free(str);

  if (DEBUGG) fprintf(stderr, "\n");

  // print debug saddle info
  if (DEBUGG) {
    fprintf(stderr, "assigned saddles:\n");
    for (unordered_map<short*, struct_info, hash_fncts, hash_eq>::iterator it=struct_map.begin(); it!=struct_map.end(); it++) {
      if (it->second.saddle_nums.size()==0) continue;
      fprintf(stderr, "%d:", it->second.debug_num);
      for (set<int>::iterator it2=it->second.saddle_nums.begin(); it2!=it->second.saddle_nums.end(); it2++) {
        fprintf(stderr, " %d", *it2);
      }
      fprintf(stderr, "\n");
    }
  }

  // remove minima, that we are uncertain of
  int i = map_num.size()-1;
  while (i>=0 && invert_map[map_num[i]].energy==gl_energy) {
    map_num.erase(i);
    i--;
  }

  return map_num.size();
}


void SetOpt(gengetopt_args_info &args_info)
{
  // create options:
  options *opt = (options*) malloc(sizeof(options));
  opt->verbose_lvl = 0;
  opt->noLP = args_info.noLP_flag;
  opt->shift = args_info.shift_flag;
  opt->first = false;
  opt->EOM = true;
  opt->f_point = landmap;

  deg.opt = opt;

  DEBUGG = args_info.debug_flag;

  gl_direct = args_info.direct_flag;

  gl_threshold = args_info.threshold_arg;
  if (gl_threshold==0) gl_threshold = INT_MAX;
}


void FreeStuff()
{
  unordered_map<short*, struct_info, hash_fncts, hash_eq>::iterator sm_it;
  for (sm_it = struct_map.begin(); sm_it!=struct_map.end(); sm_it++) free(sm_it->first);
  struct_map.clear();
  invert_map.clear();
  free_encode(enc);
  free(deg.opt);
  free(seq);
}

// print output
void PrintOutput(map<int, int> &map_num, vector<saddle> &output)
{
  // print minima
  for (unsigned int i=0; i<map_num.size(); i++) {
    int num = map_num[i];
    printf("%5d %s %6.2f\n", num, pt_to_str(invert_map[num].str).c_str(), invert_map[num].energy/100.0);
  }

  printf("\n");

  // print saddles
  for (unsigned int i=0; i<output.size(); i++) {
    printf("%5d %s %6.2f (", output[i].num, pt_to_str(output[i].str).c_str(), output[i].energy/100.0);
    // print connections
    for (set<set<int> >::iterator j=output[i].locmin_conn.begin(); j!=output[i].locmin_conn.end(); j++) {
      printf("[");
      for (set<int>::iterator k=j->begin(); k!=j->end(); k++) {
        if (k!=j->begin()) printf(" ");
        printf("%d", map_num[*k]);
      }
      printf("]");
    }
    printf(") {");
    // print saddle connections
    for (set<int>::iterator j=output[i].saddle_conn.begin(); j!=output[i].saddle_conn.end(); j++) {
      if (j!=output[i].saddle_conn.begin()) printf(" ");
      printf("%d", *j);
    }
    printf("}\n");
  }
}

// print dot file (according to the options in args_info)
void PrintDot(gengetopt_args_info &args_info, map<int, int> &map_num, vector<saddle> &output) {
  FILE *dot;
  int number = map_num.size();
  dot = fopen(args_info.name_dot_arg, "w");
  if (dot) {
    if (args_info.landmark_flag) {  // print landmark
      fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
      //nodes:
      for (int i=0; i<number; i++) {
        // print label also with number in barrier tree?
        if (1) fprintf(dot, "\"%d\" [label=\"%d(%d)\"]\n", map_num[i], map_num[i], i+1);
        else fprintf(dot, "\"%d\" [label=\"%d\"]\n", map_num[i], map_num[i]);
      }
      // edges (fuck me...)
      for (unsigned int i=0; i<output.size(); i++) {
        for (set<set<int> >::iterator it=output[i].locmin_conn.begin(); it!=output[i].locmin_conn.end(); it++) {
          set<set<int> >::iterator it2=it;
          for (it2++; it2!=output[i].locmin_conn.end(); it2++) {
            for (set<int>::iterator in=it->begin(); in!=it->end(); in++) {
              for (set<int>::iterator in2=it2->begin(); in2!=it2->end(); in2++) {
                fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n", map_num[*in], map_num[*in2], output[i].energy/100.0);
              }
            }
          }
        }
      }
      fprintf(dot, "}\n");
    } else {  // print landmap

      float color = 0.7;

      fprintf(dot, "Graph G {\n\tnode [width=0.1, height=0.1, shape=circle];\n");
      //nodes (minima):
      fprintf(dot, "/*nodes - minima*/\n");
      for (int i=0; i<number; i++) {
        // print label also with number in barrier tree?
        if (1) fprintf(dot, "\"%d\" [label=\"%d(%d)\"]\n", map_num[i], map_num[i], i+1);
        else fprintf(dot, "\"%d\" [label=\"%d\"]\n", map_num[i], map_num[i]);
      }
      //nodes (saddles):
      fprintf(dot, "/*nodes - saddles*/\n");
      for (unsigned int i=0; i<output.size(); i++) {
        fprintf(dot, "\"S%d\" [label=\"S%d\", color=\"0.0 0.0 %.2f\", fontcolor=\"0.0 0.0 %.2f\"]\n", output[i].num, output[i].num, color, color);
      }

      //edges (landmark = minima to minima)
      fprintf(dot, "/*edges - minima to minima (landmark)*/\n");
      for (unsigned int i=0; i<output.size(); i++) {
        for (set<set<int> >::iterator it=output[i].locmin_conn.begin(); it!=output[i].locmin_conn.end(); it++) {
          set<set<int> >::iterator it2=it;
          for (it2++; it2!=output[i].locmin_conn.end(); it2++) {
            for (set<int>::iterator in=it->begin(); in!=it->end(); in++) {
              for (set<int>::iterator in2=it2->begin(); in2!=it2->end(); in2++) {
                fprintf(dot, "\"%d\" -- \"%d\" [label=\"%.2f\"]\n", map_num[*in], map_num[*in2], output[i].energy/100.0);
              }
            }
          }
        }
      }
      // edges (saddles to minima)
      fprintf(dot, "/*edges - saddles to minima*/\n");
      for (unsigned int i=0; i<output.size(); i++) {
        for (set<set<int> >::iterator it=output[i].locmin_conn.begin(); it!=output[i].locmin_conn.end(); it++) {
          for (set<int>::iterator in=it->begin(); in!=it->end(); in++) {
            fprintf(dot, "\"S%d\" -- \"%d\" [color=\"0.0 0.0 %.2f\", fontcolor=\"0.0 0.0 %.2f\"]\n", output[i].num, map_num[*in], color, color);
          }
        }
      }
      // edges(saddles to saddles)
      fprintf(dot, "/*edges - saddles to saddles*/\n");
      for (unsigned int i=0; i<output.size(); i++) {
        for (set<int>::iterator it=output[i].saddle_conn.begin(); it!=output[i].saddle_conn.end(); it++) {
          fprintf(dot, "\"S%d\" -- \"S%d\" [color=\"0.0 0.0 %.2f\", fontcolor=\"0.0 0.0 %.2f\"]\n", *it, output[i].num, color, color);
        }
      }

      fprintf(dot, "}\n");
    }
    fclose(dot);

    // start neato:
    if (args_info.print_graph_flag) {
      char syst[200];
      sprintf(syst, "%s -Tps < %s > %s", (args_info.dot_flag ? "dot" : "neato"), args_info.name_dot_arg, args_info.name_graph_arg);
      system(syst);
    }
  }
}
