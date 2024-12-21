# dos.c

Given a list of eigenvalues, computes the density of states.

## Example

Example using the eigenvalues of a [simple 1d system](test.c):

```bash
make
./test > test.dat
./dos test.dat > output
```

## Options

Options can be changed using preprocessor directives in the [source code](dos.c) itself.

# License

[MIT License](LICENSE).