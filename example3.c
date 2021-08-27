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

    int ids[100] = {0};
    
    // --- ADD FIRST SEQUENCE ---
    char seq[100] = "ACG"; 
    int seq_len = strlen(seq);
    uint8_t *bseq = (uint8_t*)malloc(sizeof(uint8_t) * seq_len); 
    for (i=0; i <seq_len; i++) {
        bseq[i] = _char26_table[(int)seq[i]];
        ids[i] = abpoa_add_graph_node(ab->abg, bseq[i]);
    }
    
    int last_node_id = ABPOA_SRC_NODE_ID;
    int curr_node_id = 0;

    for (i=0; i <seq_len; i++) {
        curr_node_id = ids[i];
        abpoa_add_graph_edge(ab->abg, last_node_id, curr_node_id, 0, 1, 0, 0, 0);
        last_node_id = curr_node_id;
    }

    // --- ADD SECOND SEQUENCE ---
    char seq2[100] = "TTT"; 
    int seq_len2 = strlen(seq2);
    uint8_t *bseq2 = (uint8_t*)malloc(sizeof(uint8_t) * seq_len2); 
    for (i=0; i <seq_len2; i++) {
        bseq2[i] = _char26_table[(int)seq2[i]];
        ids[seq_len+i] = abpoa_add_graph_node(ab->abg, bseq2[i]);
    }
    for (i=0; i <seq_len2; i++) {
        curr_node_id = ids[seq_len+i];
        abpoa_add_graph_edge(ab->abg, last_node_id, curr_node_id, 0, 1, 0, 0, 0);
        last_node_id = curr_node_id;
    }
    abpoa_add_graph_edge(ab->abg, last_node_id, ABPOA_SINK_NODE_ID, 0, 1, 0, 0, 0);

    abpoa_generate_gfa(ab, abpt, stdout);
    
    // --- NOW TRY ALIGNING A SEQUENCE ---
    char new_seq[100] = "AAACTGGTAT";
    int new_seq_len = strlen(new_seq);
    uint8_t *new_bseq = (uint8_t*)malloc(sizeof(uint8_t) * new_seq_len);
    for (i=0; i<new_seq_len; i++) {
        new_bseq[i] = _char26_table[(int)new_seq[i]];
    }

    abpoa_res_t *res;
    abpoa_align_sequence_to_graph(ab, abpt, new_seq, new_seq_len, res);


    /*
    // --- ADD THIRD SEQUENCE ---
    char seq3[100] = "AAACC"; 
    int seq_len3 = strlen(seq3);
    uint8_t *bseq3 = (uint8_t*)malloc(sizeof(uint8_t) * seq_len3); 
    for (i=0; i <seq_len3; i++) {
        bseq3[i] = _char26_table[(int)seq3[i]];
        ids[seq_len+seq_len2+i] = abpoa_add_graph_node(ab->abg, bseq3[i]);
    }
    last_node_id = ids[seq_len-1];
    for (i=0; i <seq_len3; i++) {
        curr_node_id = ids[seq_len+seq_len2+i];
        abpoa_add_graph_edge(ab->abg, last_node_id, curr_node_id, 0, 1, 0, 0, 0);
        last_node_id = curr_node_id;
    }
    abpoa_add_graph_edge(ab->abg, last_node_id, ABPOA_SINK_NODE_ID, 0, 1, 0, 0, 0);
    */
    
    /*
    // ----- Now align the reads -----
    int n_seqs_aln = 10;
    char seqs_aln[10][100] = {
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

    // collect sequence length, trasform ACGT to 0123
    int *seq_lens = (int*)malloc(sizeof(int) * n_seqs_aln);
    uint8_t **bseqs = (uint8_t**)malloc(sizeof(uint8_t*) * n_seqs_aln);
    for (i = 0; i < n_seqs_aln; ++i) {
        seq_lens[i] = strlen(seqs_aln[i]);
        bseqs[i] = (uint8_t*)malloc(sizeof(uint8_t) * seq_lens[i]);
        for (j = 0; j < seq_lens[i]; ++j)
            bseqs[i][j] = _char26_table[(int)seqs_aln[i][j]];
    }
    */


    /*
    // 2. variables to store result
    uint8_t **cons_seq; int **cons_cov, *cons_l, cons_n=0;
    uint8_t **msa_seq; int msa_l=0;

    // perform abpoa-msa
    ab->abs->n_seq = 0; // To re-use ab, n_seq needs to be set as 0
    abpoa_msa(ab, abpt, n_seqs_aln, NULL, seq_lens, bseqs, NULL, &cons_seq, &cons_cov, &cons_l, &cons_n, &msa_seq, &msa_l);
    fprintf(stdout, ">Multiple_sequence_alignment\n");
    for (i = 0; i < n_seqs_aln; ++i) {
        for (j = 0; j < msa_l; ++j) {
            fprintf(stdout, "%c", _char256_table[msa_seq[i][j]]);
        }
        fprintf(stdout, "\n");
    }
    
    //abpoa_generate_gfa(ab, abpt, stdout);
    */
    
    free(bseq); free(bseq2); //free(seq_lens); free(bseqs);
    free(new_bseq);
    abpoa_free(ab); abpoa_free_para(abpt); 
    return 0;
}
