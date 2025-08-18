#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include <stdint.h>
#include <unistd.h>
#include <getopt.h>  // For command line parsing
#include <glib.h>
#include <math.h>
#include <omp.h>
#include <limits.h>		/* for CHAR_BIT */
#define BARCODE_LENGTH 16
#define BARCODE_CODE_LENGTH 4
#define MAX_BARCODES 30000
#define MAX_BARCODE_MISMATCHES 3
#define MIN_POSTERIOR 0.975
#define MAX_FEATURE_N 2
#define MAX_BARCODE_N 1

#define LINE_STORAGE_CHUNK_SIZE 100000 * 4 * LINE_LENGTH
#define UNMATCHED_BARCODE_STORAGE_CHUNK_SIZE 100000 * sizeof(unmatched_barcodes)    
#define MAX_UNMATCHED_MALLOCS 100
#define MAX_BARCODE_MALLOCS 100
#define MAX_LINE_MALLOCS 1000
#define LINE_LENGTH 1024
#define MAX_FILES 128

#define LOOKUP_STRING "AAAAAAACAAAGAAATAACAAACCAACGAACTAAGAAAGCAAGGAAGTAATAAATCAATGAATTACAAACACACAGACATACCAACCCACCGACCTACGAACGCACGGACGTACTAACTCACTGACTTAGAAAGACAGAGAGATAGCAAGCCAGCGAGCTAGGAAGGCAGGGAGGTAGTAAGTCAGTGAGTTATAAATACATAGATATATCAATCCATCGATCTATGAATGCATGGATGTATTAATTCATTGATTTCAAACAACCAAGCAATCACACACCCACGCACTCAGACAGCCAGGCAGTCATACATCCATGCATTCCAACCACCCAGCCATCCCACCCCCCCGCCCTCCGACCGCCCGGCCGTCCTACCTCCCTGCCTTCGAACGACCGAGCGATCGCACGCCCGCGCGCTCGGACGGCCGGGCGGTCGTACGTCCGTGCGTTCTAACTACCTAGCTATCTCACTCCCTCGCTCTCTGACTGCCTGGCTGTCTTACTTCCTTGCTTTGAAAGAACGAAGGAATGACAGACCGACGGACTGAGAGAGCGAGGGAGTGATAGATCGATGGATTGCAAGCACGCAGGCATGCCAGCCCGCCGGCCTGCGAGCGCGCGGGCGTGCTAGCTCGCTGGCTTGGAAGGACGGAGGGATGGCAGGCCGGCGGGCTGGGAGGGCGGGGGGGTGGTAGGTCGGTGGGTTGTAAGTACGTAGGTATGTCAGTCCGTCGGTCTGTGAGTGCGTGGGTGTGTTAGTTCGTTGGTTTTAAATAACTAAGTAATTACATACCTACGTACTTAGATAGCTAGGTAGTTATATATCTATGTATTTCAATCACTCAGTCATTCCATCCCTCCGTCCTTCGATCGCTCGGTCGTTCTATCTCTCTGTCTTTGAATGACTGAGTGATTGCATGCCTGCGTGCTTGGATGGCTGGGTGGTTGTATGTCTGTGTGTTTTAATTACTTAGTTATTTCATTCCTTCGTTCTTTGATTGCTTGGTTGTTTTATTTCTTTGTTTT"
typedef struct _barcode_type{
    unsigned char sequence_code[BARCODE_CODE_LENGTH];
} barcode_type;

typedef struct _barcode_count{
    unsigned char sequence_code[BARCODE_CODE_LENGTH];
    uint16_t counts;
} barcode_count;

typedef struct _barcode_storage{
    barcode_type *storage[MAX_BARCODE_MALLOCS];
    barcode_type *next_entry;
    barcode_type *last_entry;
    int nMallocs;
} barcode_storage;

typedef struct _line_storage{
    char *storage[MAX_LINE_MALLOCS];
    char *next_entry;
    char *last_entry;
    int nMallocs;
}line_storage;

typedef struct _unmatched_barcodes{
    unsigned char barcode[BARCODE_CODE_LENGTH];
    unsigned char closest_barcodes[MAX_BARCODE_MISMATCHES][BARCODE_CODE_LENGTH];
    unsigned char number_of_closest_barcodes;
    unsigned char qscores[MAX_BARCODE_MISMATCHES];
    char *lines;
    unsigned char file_index;
    struct _unmatched_barcodes *next;
} unmatched_barcodes;

typedef struct _unmatched_barcodes_storage{
    unmatched_barcodes *storage[MAX_UNMATCHED_MALLOCS];
    unmatched_barcodes *next_entry;
    unmatched_barcodes *last_entry;
    int nMallocs;
} unmatched_barcodes_storage;

typedef struct _unmatched_barcode_list{
    unmatched_barcodes *first_entry;
    unmatched_barcodes *last_entry;
} unmatched_barcode_list;


int check_barcode(char *barcode, unsigned char *code);
barcode_count*  readWhiteList(char *whitelist_filename,GHashTable *hash);
void string2code(char *string, int sequence_length, unsigned char *code);
void string2all_codes(char *string, unsigned char codes[4][LINE_LENGTH/2+1]);
int find_matching_codes(unsigned char *a, unsigned char *b, size_t length);
void code2string(unsigned char *code, char *string, int code_length);
int fuzzy_matching_codes(unsigned char *a, unsigned char *b, size_t length, int maxHammingDistance);
char check_sequence(char *sequence, int sequence_length);
void allocate_unmatched_barcodes_storage(unmatched_barcodes_storage *storage_pool);
int find_closest_barcodes(unsigned char *code, barcode_type *corrected_codes, unsigned char *indices);
int find_variant_match(unsigned char *code, int sequence_index, unsigned char *corrected_bases );
int checkSequenceAndCorrectForN(char *line, char *corrected_lines[], char *buffer,int sequence_length, int maxN);
int generate_sequences(char* string, int string_length, int *indices, char **output, int index, int number_of_indices, char *output_indices[], int number_of_outputs);
char* printToBuffer(char *string, int sequence_length, char *buffer);
char* allocate_line_storage(line_storage *line_storage_pool);
int update_barcode_counts(unsigned char *code, char *correction_needed);

int checkAndCorrectBarcode(char **lines, int maxN, unsigned char file_index);
unmatched_barcodes* add_unmatched_barcode_store_lines(unsigned char *barcode, barcode_type *corrected_codes,  char **lines, unsigned char *qscores, int number_of_variants, unsigned char file_index);
unmatched_barcodes* allocate_unmatched_barcodes();
unsigned char* find_best_posterior_match (unmatched_barcodes *entry, double min_posterior);

