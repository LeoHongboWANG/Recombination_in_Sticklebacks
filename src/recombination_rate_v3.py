import pandas as pd
from argparse import ArgumentParser

def calculate_recombination_rates(linkage_map_file, window_size, stride):
    linkage_map = pd.read_csv(linkage_map_file, sep="\t", header=None, names=["chr", "pos_bp", "pos_cM"])
    recombination_data = []

    for chrom in linkage_map["chr"].unique():
        chrom_data = linkage_map[linkage_map["chr"] == chrom]
        start_bp = 0
        
        while start_bp <= chrom_data["pos_bp"].max():
            end_bp = start_bp + window_size
            interval_data = chrom_data[(chrom_data["pos_bp"] >= start_bp) & (chrom_data["pos_bp"] < end_bp)]
            
            if not interval_data.empty:
                min_cM = interval_data["pos_cM"].min()
                max_cM = interval_data["pos_cM"].max()
                delta_cM = max_cM - min_cM
                rho = (delta_cM / (window_size / 1e6))  # Recombination rate per megabase pair
            else:
                rho = 0
            
            recombination_data.append([chrom, start_bp, end_bp, rho])
            start_bp += stride

    recombination_rates = pd.DataFrame(recombination_data, columns=["chr", "start_bp", "end_bp", "rho"])
    return recombination_rates

def main():
    parser = ArgumentParser(description="Calculate recombination rates given a specified window size and stride.")
    parser.add_argument("linkage_map_file", type=str, help="Linkage map file containing chromosome, position in base pairs, and genetic distance in centiMorgans.")
    parser.add_argument("window_size", type=int, help="Window size in base pairs (e.g., 1000000 for 1Mb windows).")
    parser.add_argument("stride", type=int, help="Stride length in base pairs (e.g., 500000 for a stride of 500kb).")
    parser.add_argument("-o", "--output", type=str, default="recombination_rates.txt", help="Output file to save the recombination rates.")
    args = parser.parse_args()

    recombination_rates = calculate_recombination_rates(args.linkage_map_file, args.window_size, args.stride)
    recombination_rates.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
    main()
