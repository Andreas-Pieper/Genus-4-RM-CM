#define  _GNU_SOURCE
#include <assert.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>

#include "arb.h"
#include "acb.h"
#include "acb_mat.h"
#include "acb_theta.h"
#define LOG102 0.301029995663981195213738894725
//arb_print := func<elt | ReplaceCharacter(Sprintf("%o +/- %.*o", elt, 3, Abs(elt)*10^(1-Precision(Parent(elt)))), "E", "e")>;


bool check_usage(int argc, char* argv[]) {
    if(argc < 5) {
        fprintf(stderr, "expected 4 arguments:\n\t%s genus precision input_file output_file\n", argv[0]);
        fprintf(stderr, "got %d:\n\t", argc - 1);
        for(size_t i=0; i<argc; i++) {
            fprintf(stderr, "%s ", argv[i]);
        }
        fprintf(stderr, "\n");
        return false;
    }
    return true;
}

bool arb_set_str_line(arb_t res, FILE *stream, const slong prec) {
    static char *line = NULL;
    static size_t len = 0;
    static ssize_t nread;
    nread = getline(&line, &len, stream);
    if(nread == -1)
        return false;

    // drop the end of line
    if( (nread > 0) && line[nread - 1] == '\n' ) {
        line[nread - 1] = '\0';
    }
    return 0 == arb_set_str(res, line, prec);
}

bool acb_set_str_twolines(acb_t res, FILE *stream, const slong prec) {
    return arb_set_str_line(acb_realref(res), stream, prec) && arb_set_str_line(acb_imagref(res), stream, prec);
}

bool read_input(acb_mat_t tau, acb_ptr z, const char filename[], const slong prec) {
    assert( acb_mat_nrows(tau) == acb_mat_ncols(tau) );
    slong g = acb_mat_ncols(tau);

    FILE *stream = fopen(filename, "r");

    bool res = true;
    if (stream == NULL) {
        fprintf(stderr, "read_input: couldn't open %s\n", filename);
        return false;
    }

    if (res) {
     for(size_t i=0; i < g; ++i) {
          if( !acb_set_str_twolines(z + i, stream, prec) ){
            res = false;
            break;
          }
      }
    }

    if (res) {
      for(size_t i=0; i < g; ++i) {
          for(size_t j=0; j < g; ++j) {
              if( !acb_set_str_twolines(acb_mat_entry(tau, i, j), stream, prec) ) {
                // break both loops
                res = false;
                i = j = g;
                break;
              }
          }
      }
    }
    fclose(stream);
    return res;
}

bool write_output(acb_ptr theta, slong len, const char filename[]) {
    FILE *stream = fopen(filename, "w");
    bool res = true;
    if (stream == NULL) {
        fprintf(stderr, "write_output: couldn't open %s\n", filename);
        res = false;
    }
    if (res) {
      for (size_t i=0; i < len; ++i) {
        slong re_digits = floor(arb_rel_accuracy_bits(acb_realref(&theta[i]))*LOG102);
        arb_fprintd(stream, acb_realref(&theta[i]), re_digits);
        flint_fprintf(stream, "\n");
        slong imag_digits = floor(arb_rel_accuracy_bits(acb_imagref(&theta[i]))*LOG102);
        arb_fprintd(stream, acb_imagref(&theta[i]), imag_digits);
        flint_fprintf(stream, "\n");
      }
    }
    fclose(stream);
    return res;
}

int main(int argc, char* argv[])
{
    if( !check_usage(argc, argv) ) {
        return EXIT_FAILURE;
    }
    // input variables
    slong g = atoi(argv[1]);
    slong prec = atoi(argv[2]);
    char* input_filename = argv[3];
    char* output_filename = argv[4];
    slong verbose = false;
    if(argc > 5)
        verbose = atoi(argv[5]);
    if( verbose ) {
        flint_printf("g = %ld\n", g);
        flint_printf("prec = %ld\n", prec);
    }


    acb_mat_t tau;
    acb_mat_init(tau, g, g);

    acb_ptr z = _acb_vec_init(g);

    slong nb = 1 << (2 * g);
    acb_ptr theta = _acb_vec_init(nb); // 2^2g

    if( !read_input(tau, z, input_filename, prec) ) {
        fprintf(stderr, "Error reading input file %s\n", input_filename);
        return EXIT_FAILURE;
    }
    if( verbose ) {
        flint_printf("z:\n");
        _acb_vec_printd(z, g, 10);
        flint_printf("Tau:\n");
        acb_mat_printd(tau, 10);
    }


    // Compute theta values
    acb_theta_all(theta, z, tau, false, prec);



    if( verbose ) {
        flint_printf("\nTheta values:\n");
        for (size_t i = 0; i < nb; ++i) {
            acb_printd(&theta[i], 10); flint_printf("\n");
        }
    }
    if( !write_output(theta, nb, output_filename) )
        return EXIT_FAILURE;

    // Clear and exit
    acb_mat_clear(tau);
    _acb_vec_clear(theta, nb);
    _acb_vec_clear(z, g);
    return EXIT_SUCCESS;
}


