#include <TROOT.h>
#include <argp.h>
#include <iostream>

#ifndef _PARSER_
#define _PARSER_

static int parse_opt(int key, char *arg, struct argp_state *state);

// -- Global variables
struct argp_option options[] = {{"overwrite", 'w', 0, 0, "Overwrite results"},
                                {"debug", 'g', 0, 0, "Enable Debug Info"},
                                {"input-file", 'i', "FILE", 0, "Input file"},
                                {"output-file", 'o', "FILE", 0, "Output file"},
                                {"nTrees", 'n', "INT", 0, "Number of trees"},
                                {"maxDepth", 'm', "INT", 0, "Max depth"},
                                {"nCuts", 'c', "INT", 0, "Number of cuts"},
                                {0}};
struct argp argp = {options, parse_opt, 0, 0};

struct arguments {
  const char* infilename;
  const char* outfilename;
  bool flag_input;
  bool flag_output;
  bool flag_overwrite;
  bool flag_debug;
  Int_t nTrees;
  Int_t maxDepth;
  Int_t nCuts;
};

static int parse_opt(int key, char *arg, struct argp_state *state) {

  struct arguments *a = (struct arguments *)state->input;

  switch (key) {
  case 'w': {
    a->flag_overwrite = kTRUE;
    break;
  }
  case 'g': {
    a->flag_debug = kTRUE;
    break;
  }
  case 'i': {
    a->flag_input = kTRUE;
    a->infilename = arg;
    break;
  }
  case 'o': {
    a->flag_output = kTRUE;
    a->outfilename = arg;
    break;
  }
  case 'n': {
    a->nTrees = atoi(arg);
    break;
  }
  case 'm': {
    a->maxDepth = atoi(arg);
    break;
  }
  case 'c': {
    a->nCuts = atoi(arg);
    break;
  }
  case ARGP_KEY_ARG: {
    argp_failure(state, 1, 0, "No arguments requested");
    break;
  }
  case ARGP_KEY_INIT: {
    a->infilename     = NULL;
    a->outfilename    = NULL;
    // a->recotype       = "STD";
    a->flag_overwrite = kFALSE;
    a->flag_debug = kFALSE;
    a->flag_input = kFALSE;
    a->flag_output = kFALSE;
    a->nTrees = 300;
    a->maxDepth = 9;
    a->nCuts = 50;
    break;
  }
  case ARGP_KEY_END: {
    std::cout << std::endl;
    break;
  }
  }

  return 0;
}

#endif
