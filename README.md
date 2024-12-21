# dos.c

Given a list of eigenvalues, computes the density of states.

```bash
$ ./dos --help
Usage: ./dos [OPTION]... FILE
Calculate density of states from eigenvalues in FILE.

Options:
  -s, --sigma=VALUE      set the broadening (default: 0.05)
  -p, --padding=VALUE    set the padding (default: 0.10)
  -n, --npoints=NUM      set number of points (default: 1000)
  -h, --help             display this help and exit

Report bugs to: <https://github.com/keyehzy/dos.c/issues>
```

## Example

Example using the eigenvalues of a [simple 1d system](test.c):

```bash
make
./test > test.dat
./dos test.dat > output
```

# License

[MIT License](LICENSE).