/* example.c libabpoa usage example
   To compile: 
gcc -g example.c -I ./include -L ./lib -labpoa -lz -lm -o example
   or:
gcc -g example.c -I ./include ./lib/libabpoa.a -lz -lm -o example
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "include/abpoa.h"

// AaCcGgTtNn ... ==> 0,1,2,3,4 ...
// BbDdEeFf   ... ==> 5,6,7,8 ...
unsigned char _char26_table[256] = {
	 0,  1,  2,  3,   4,  5,  6,  7,   8,  9, 10, 11,  12, 13, 14, 15, 
	16, 17, 18, 19,  20, 21, 22, 23,  24, 25, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26,  0,  5,  1,   6,  7,  8,  2,   9, 10, 11, 12,  13, 14,  4, 15, 
	16, 17, 18, 19,   3, 20, 21, 22,  23, 24, 25, 26,  26, 26, 26, 26, 
	26,  0,  5,  1,   6,  7,  8,  2,   9, 10, 11, 12,  13, 14,  4, 15, 
	16, 17, 18, 19,   3, 20, 21, 22,  23, 24, 25, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26, 
	26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26,  26, 26, 26, 26
};

// 0/1/2/3/4=>ACGTN
// 5/6/7/8=>BDEF ...
const char _char256_table[256] = {
	'A', 'C', 'G', 'T',  'N', 'B', 'D', 'E',  'F', 'H', 'I',  'J', 'K', 'L', 'M', 'O',
	'P', 'Q', 'R', 'S',  'U', 'V', 'W', 'X',  'Y', 'Z', '*', '-',  '*', '*', '*', '*',
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', 'A', 'B', 'C',  'D', 'E', 'F', 'G',  'H', 'I', 'J', 'K',  'L', 'M', 'N', 'O', 
	'P', 'Q', 'R', 'S',  'T', 'U', 'V', 'W',  'X', 'Y', 'Z', '*',  '*', '*', '*', '*', 
	'*', 'A', 'B', 'C',  'D', 'E', 'F', 'G',  'H', 'I', 'J', 'K',  'L', 'M', 'N', 'O', 
	'P', 'Q', 'R', 'S',  'T', 'U', 'V', 'W',  'X', 'Y', 'Z', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*', 
	'*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*',  '*', '*', '*', '*'
};

int main(void) {
    int i, j, n_seqs = 10;
    char seqs[10][100] = {
        "CGTCAATCTATCGAAGCATACGCGGGCAGAGCCGAAGACCTCGGCAATCCA",
        "CCACGTCAATCTATCGAAGCATACGCGGCAGCCGAACTCGACCTCGGCAATCAC",
        "CGTCAATCTATCGAAGCATACGCGGCAGAGCCCGGAAGACCTCGGCAATCAC",
        "CGTCAATGCTAGTCGAAGCAGCTGCGGCAGAGCCGAAGACCTCGGCAATCAC",
        "CGTCAATCTATCGAAGCATTCTACGCGGCAGAGCCGACCTCGGCAATCAC",
        "CGTCAATCTAGAAGCATACGCGGCAAGAGCCGAAGACCTCGGCCAATCAC",
        "CGTCAATCTATCGGTAAAGCATACGCTCTGTAGCCGAAGACCTCGGCAATCAC",
        "CGTCAATCTATCTTCAAGCATACGCGGCAGAGCCGAAGACCTCGGCAATC",
        "CGTCAATGGATCGAGTACGCGGCAGAGCCGAAGACCTCGGCAATCAC",
        "CGTCAATCTAATCGAAGCATACGCGGCAGAGCCGTCTACCTCGGCAATCACGT"
        };

    // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
     
    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    abpt->out_cons = 1; // generate consensus sequence, set 0 to disable
    abpt->w = 6, abpt->k = 9; abpt->min_w = 10; // minimizer-based seeding and partition
    abpt->progressive_poa = 1;

    abpoa_post_set_para(abpt);

    // --- ADD FIRST SEQUENCE ---

    char seq[100] = "ACG"; 
    int seq_len = strlen(seq);
    uint8_t *bseq = (uint8_t*)malloc(sizeof(uint8_t) * seq_len); 
    for (i=0; i <seq_len; i++) {
        bseq[i] = _char26_table[(int)seq[i]];
    }

    abpoa_res_t res;
    abpoa_add_graph_alignment(ab, abpt, bseq, seq_len, NULL, res, 0, 1, 1);
    //abpoa_add_graph_sequence(ab, abpt, bseq, seq_len, NULL, res, 0, 1, 1);

    //abpoa_generate_gfa(ab, abpt, stdout);

    // --- ADD SECOND SEQUENCE ---
    // TODO: use abpoa_add_subgraph_alignment, same as add_graph_alignment but requires ids

    char seq2[100] = "ATGCTTTTT";
    int seq_len2 = strlen(seq2);
    uint8_t *bseq2 = (uint8_t*)malloc(sizeof(uint8_t) * seq_len2); 
    for (i=0; i <seq_len2; i++) {
        bseq2[i] = _char26_table[(int)seq2[i]];
        printf("%c", bseq2[i]);
    }

    abpoa_res_t res2;
    abpoa_add_graph_alignment(ab, abpt, bseq2, seq_len2, NULL, res2, 1, 2, 1);
    //abpoa_add_subgraph_alignment(ab, abpt, seq_len, seq_len+seq_len2, bseq2, seq_len2, NULL, res2, 1, 2, 1);

    //abpoa_generate_gfa(ab, abpt, stdout);
    
    free(bseq); free(bseq2);
    abpoa_free(ab); abpoa_free_para(abpt); 
    return 0;
}
