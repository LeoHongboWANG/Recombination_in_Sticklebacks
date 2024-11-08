#!/usr/bin/env python3
import sys

def main():
    if len(sys.argv) != 3:
        sys.stderr.write(f"Usage: {sys.argv[0]} <gtf> <id.txt(transcript_id\tgene_id)> > out.gtf\n")
        sys.exit(1)
    gtf_file, id_file = sys.argv[1], sys.argv[2]
    transcript_ids = set()
    gene_ids = set()
    with open(id_file) as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                transcript_ids.add(parts[0])
                gene_ids.add(parts[1])
    with open(gtf_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            attrs = parts[8]
            gene_id = None
            transcript_id = None
            for attr in attrs.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_id '):
                    gene_id = attr.split('"')[1]
                elif attr.startswith('transcript_id '):
                    transcript_id = attr.split('"')[1]
            if (parts[2] == 'gene' and gene_id in gene_ids) or (transcript_id in transcript_ids):
                print(line, end='')

if __name__ == '__main__':
    main()
