#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>


#define sizeof(x)    (ptrdiff_t) sizeof(x)
#define countof(a)   (sizeof(a) / sizeof(*(a)))
#define lengthof(s)  (countof(s) - 1)


typedef struct {
  char *beg;
  char *end;
} arena;

arena newarena(ptrdiff_t cap) {
  arena a = {0};
  a.beg = malloc(cap);
  a.end = a.beg ? a.beg+cap : 0;
  return a;
}

__attribute((malloc, alloc_size(2, 4), alloc_align(3)))
void *alloc(arena *a, ptrdiff_t size, ptrdiff_t align, ptrdiff_t count)
{
    ptrdiff_t padding = -(uintptr_t)a->beg & (align - 1);
    ptrdiff_t available = a->end - a->beg - padding;
    if (available < 0 || count > available / size) {
        abort();
    }
    void *p = a->beg + padding;
    a->beg += padding + count*size;
    return memset(p, 0, count*size);
}

#define NEW(a, t, n)  (t *) alloc(a, sizeof(t), _Alignof(t), n)


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
  ptrdiff_t cap;
} fltlist;


#define MAX(x, y) (((x) > (y)) ? (x) : (y))

void list_push(fltlist *lst, double item, arena *a) {
  if (lst->len + 1 >= lst->cap) {
    ptrdiff_t newcap = MAX(1, 2 * lst->cap);
    double   *oldptr = lst->data;
    double   *newptr = NEW(a, double, newcap);
    memcpy(newptr, oldptr, sizeof(double) * lst->len);
    lst->data = newptr;
    lst->cap  = newcap;
  }
  lst->data[lst->len++] = item;
}


static str readcontent(str path, arena *a) {
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

    fcontent = NEW(a, char, fsize);
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

fltlist getitems(str path, arena *perm) {
  strpair pair = {0};
  pair.tail = readcontent(path, perm);

  fltlist lst = {0};

  while (pair.tail.len) {
    pair = chop(pair.tail, '\n');
    list_push(&lst, to_double(pair.head), perm);
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

#define NUM_POINTS 1000
#define SIGMA 0.05
#define PADDING_FACTOR 0.1
#define PI 3.14159265358979323844

void print_dos(fltlist lst, arena *a) {
  double min_val;
  double max_val;
  find_min_max(lst, &min_val, &max_val);

  double padding = (max_val - min_val) * PADDING_FACTOR;
  min_val -= padding;
  max_val += padding;

  double dE       = (max_val - min_val) / (double)(NUM_POINTS - 1);
#if defined(LORENTZ_KERNEL)
  double norm     = 1.0 / PI / SIGMA;
#else
  double norm     = 1.0 / sqrt(2.0 * PI) / SIGMA;
#endif
  double integral = 0.0;

  double *E_values   = NEW(a, double, NUM_POINTS);
  double *DOS_values = NEW(a, double, NUM_POINTS);

  for (ptrdiff_t i = 0; i < NUM_POINTS; i++) {
    double E = min_val + (double)i * dE;
    double rho = 0.0;

    for (ptrdiff_t j = 0; j < lst.len; j++) {
      double delta_E = (E - lst.data[j]) / SIGMA;
#if defined(LORENTZ_KERNEL)
      rho += 1.0 / (1.0 + delta_E * delta_E);
#else
      rho += exp(-0.5 * delta_E * delta_E);
#endif
    }

    E_values[i] = E;
    DOS_values[i] = rho * norm;
    integral += rho * norm * dE;
  }

  for (ptrdiff_t i = 0; i < NUM_POINTS; ++i) {
    printf("%lf %lf\n", E_values[i], DOS_values[i] / integral);
  }
}

int main(int argc, char **argv) {
  if (argc != 2) {
    fprintf(stderr,
            "Usage: %s [OPTION]... FILE\n"
            "Calculate density of states from eigenvalues in FILE.\n",
            argv[0]);
    return 1;
  }

  arena perm = newarena(1<<15);
  str path = STR(argv[1], strlen(argv[1]));
  fltlist items = getitems(path, &perm);
  print_dos(items, &perm);
  return 0;
}
