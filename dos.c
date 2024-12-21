#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <getopt.h>

#define sizeof(x)    (ptrdiff_t) sizeof(x)
#define countof(a)   (sizeof(a) / sizeof(*(a)))
#define lengthof(s)  (countof(s) - 1)

typedef struct {
  char *data;
  ptrdiff_t len;
} str;

typedef struct {
  str head;
  str tail;
} strpair;

#define S(s)          (str) { s, lengthof(s) }
#define STR(s, l)     (str) { s, l }
#define STRPAIR(h, t) (strpair) { h, t }

typedef struct {
  double *data;
  ptrdiff_t len;
} fltlist;


static str readcontent(str path) {
  char *fcontent = NULL;
  int fsize = 0;
  FILE *fp;

  static char buf[256];
  if (path.len >= 256) return (str) {0};
  memcpy(buf, path.data, path.len);
  buf[path.len] = '\0';

  fp = fopen(buf, "r");
  if(fp) {
    fseek(fp, 0, SEEK_END);
    fsize = ftell(fp);
    rewind(fp);

    fcontent = malloc(sizeof(char) * fsize);
    fread(fcontent, 1, fsize, fp);

    fclose(fp);
  }

  return STR(fcontent, fsize);
}

static strpair chop(str s, char delim) {
  if (s.len == 0 || s.data == NULL) {
    return (strpair) {0};
  }

  char *delim_pos = memchr(s.data, delim, s.len);

  if (delim_pos == NULL) {
    return STRPAIR(STR(s.data, s.len), (str) {0});
  }

  ptrdiff_t head_len = delim_pos - s.data;
  ptrdiff_t tail_len = s.len - head_len - 1;

  return STRPAIR(STR(s.data, head_len), STR(delim_pos + 1, tail_len));
}

static double to_double(str s) {
    static char buf[256];
    if (s.len >= 256) return 0.0/0.0;
    memcpy(buf, s.data, s.len);
    buf[s.len] = '\0';
    return strtod(buf, NULL);
}

static ptrdiff_t count_lines(str content) {
  ptrdiff_t count = 0;
  str remaining = content;

  while (remaining.len) {
    remaining = chop(remaining, '\n').tail;
    count++;
  }

  return count;
}

static fltlist getitems(str path) {
  str content = readcontent(path);
  ptrdiff_t num_lines = count_lines(content);

  fltlist lst = {0};
  lst.data = malloc(sizeof(double) * num_lines);
  lst.len = 0;

  strpair pair = {0};
  pair.tail = content;

  while (pair.tail.len) {
    pair = chop(pair.tail, '\n');
    lst.data[lst.len++] = to_double(pair.head);
  }

  return lst;
}

static void find_min_max(fltlist lst, double *min, double *max) {
  if (lst.len == 0) {
    *min = 0.0/0.0;
    *max = 0.0/0.0;
    return;
  }

  *min = lst.data[0];
  *max = lst.data[0];

  for (ptrdiff_t i = 1; i < lst.len; i++) {
    if (lst.data[i] < *min) *min = lst.data[i];
    if (lst.data[i] > *max) *max = lst.data[i];
  }
}

#define PI 3.14159265358979323844

void print_dos(fltlist lst, double sigma, double padding_factor, ptrdiff_t num_points) {
  double min_val;
  double max_val;
  find_min_max(lst, &min_val, &max_val);

  double padding = (max_val - min_val) * padding_factor;
  min_val -= padding;
  max_val += padding;

  double dE = (max_val - min_val) / (double)(num_points - 1);
#if defined(LORENTZ_KERNEL)
  double norm = 1.0 / PI / sigma;
#else
  double norm = 1.0 / sqrt(2.0 * PI) / sigma;
#endif

  double integral    = 0;
  double* DOS_values = malloc(sizeof(double) * num_points);

#pragma omp parallel for
  for (ptrdiff_t i = 0; i < num_points; i++) {
    double E = min_val + (double)i * dE;
    double rho = 0.0;

    for (ptrdiff_t j = 0; j < lst.len; j++) {
      double delta_E = (E - lst.data[j]) / sigma;
#if defined(LORENTZ_KERNEL)
      rho += 1.0 / (1.0 + delta_E * delta_E);
#else
      rho += exp(-0.5 * delta_E * delta_E);
#endif
    }
    DOS_values[i] = rho * norm;
    integral     += rho * norm * dE;
  }

   for (ptrdiff_t i = 0; i < num_points; ++i) {
    double E = min_val + (double)i * dE;
    printf("%lf %lf\n", E, DOS_values[i] / integral);
  }
}

#define SIGMA 0.05
#define PADDING_FACTOR 0.1
#define NUM_POINTS 1000

static void usage(const char *progname, int status) {
  fprintf(status ? stderr : stdout,
          "Usage: %s [OPTION]... FILE\n"
          "Calculate density of states from eigenvalues in FILE.\n"
          "\n"
          "Options:\n"
          "  -s, --sigma=VALUE      set the broadening (default: %.2f)\n"
          "  -p, --padding=VALUE    set the padding (default: %.2f)\n"
          "  -n, --npoints=NUM      set number of points (default: %d)\n"
          "  -h, --help             display this help and exit\n"
          "\n"
          "Report bugs to: <https://github.com/keyehzy/dos.c/issues>\n",
          progname, SIGMA, PADDING_FACTOR, NUM_POINTS);
  exit(status);
}

int main(int argc, char **argv) {
  static struct option long_options[] = {
    {"sigma",   required_argument, 0, 's'},
    {"padding", required_argument, 0, 'p'},
    {"npoints", required_argument, 0, 'n'},
    {"help",    no_argument,       0, 'h'},
    {0, 0, 0, 0}
  };

  double sigma   = SIGMA;
  double padding = PADDING_FACTOR;
  int npoints    = NUM_POINTS;

  while (1) {
    int option_index = 0;
    int c = getopt_long(argc, argv, "s:p:n:h", long_options, &option_index);

    if (c == -1) break;

    switch (c) {
    case 's':
      sigma = strtod(optarg, NULL);
      if (sigma <= 0) {
        fprintf(stderr, "%s: sigma must be positive\n", argv[0]);
        return 1;
      }
      break;

    case 'p':
      padding = strtod(optarg, NULL);
      if (padding <= 0) {
        fprintf(stderr, "%s: padding must be positive\n", argv[0]);
        return 1;
      }
      break;

    case 'n':
      npoints = strtol(optarg, NULL, 10);
      if (npoints <= 0) {
        fprintf(stderr, "%s: number of points must be positive\n", argv[0]);
        return 1;
      }
      break;

    case 'h':
      usage(argv[0], 0);
      break;

    case '?':
      usage(argv[0], 1);
      break;
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "%s: missing input file\n", argv[0]);
    fprintf(stderr, "Try '%s --help' for more information.\n", argv[0]);
    return 1;
  }

  str path = STR(argv[optind], strlen(argv[optind]));
  fltlist items = getitems(path);
  print_dos(items, sigma, padding, npoints);
  return 0;
}
