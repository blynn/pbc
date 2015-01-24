// Adapted from: XGetopt.h Version 1.2 by Hans Dietrich

#pragma once

extern int optind, opterr;
extern char *optarg;

int getopt(int argc, char *argv[], char *optstring);