#include <stdio.h>
#include <string.h>

int main(int argc, char** argv) {
    long start = 0L, end = 0L;
    int prev = -1;
    char chrom[100] = {0}; 
    int chrom_len = 0;  
    int is_chrom_set = 0; 

    for (;;) {
        int c = fgetc(stdin);
        switch (c) {
            case EOF:
            case '>': {
                if (is_chrom_set) {
                    chrom[chrom_len] = '\0';
                }
                if (c == EOF) return 0;
                chrom_len = 0;
                while ((c = fgetc(stdin)) != EOF && c != '\n') {
                    if (chrom_len < sizeof(chrom) - 1) {
                        chrom[chrom_len++] = c;
                    }
                }
                chrom[chrom_len] = '\0';
                is_chrom_set = 1;
                start = 0L;
                prev = -1;
                break;
            }
            case '\n':
                break;
            case 'g':
            case 'G': {
                if (prev == 'C') {
                    end = start + 1;
                    if (is_chrom_set) {
                        printf("%s\t%ld\t%ld\n", chrom, start - 1, end);
                    }
                }
                prev = 'G';
                start++;
                break;
            }
            case 'c':
            case 'C': {
                if (prev == 'G') {
                    end = start + 1;
                    if (is_chrom_set) {
                        printf("%s\t%ld\t%ld\n", chrom, start - 1, end);
                    }
                }
                prev = 'C';
                start++;
                break;
            }
            default: {
                prev = c;
                start++;
                break;
            }
        }
    }
    return 0;
}
