import pandas as pd
import numpy as np
from argparse import ArgumentParser
from scipy import stats

def calculate_recombination_rate_linear(window_data, window_size):
    slope, intercept, r_value, p_value, std_err = stats.linregress(window_data["pos_bp"], window_data["pos_cM"])
    rate = slope * 1e6
    return max(0, rate) 

def calculate_recombination_rate_maxmin(window_data, window_size):
    min_cM = window_data["pos_cM"].min()
    max_cM = window_data["pos_cM"].max()
    delta_cM = max_cM - min_cM
    rate = (delta_cM / (window_size / 1e6))
    return max(0, rate) 

def calculate_recombination_rates(linkage_map_file, window_size, shift, threshold, method):
    linkage_map = pd.read_csv(linkage_map_file, sep="\t", header=None, names=["chr", "pos_bp", "pos_cM"])
    linkage_map = linkage_map.sort_values(["chr", "pos_bp"])
    recombination_data = []
    negative_count = 0

    for chrom in linkage_map["chr"].unique():
        chrom_data = linkage_map[linkage_map["chr"] == chrom]
        
        start_bp = chrom_data["pos_bp"].min()
        end_bp = chrom_data["pos_bp"].max()
        
        positions = np.arange(start_bp + window_size/2, end_bp, shift)
        
        for pos in positions:
            window_start = int(pos - window_size/2)
            window_end = int(pos + window_size/2)
            
            window_data = chrom_data[(chrom_data["pos_bp"] >= window_start) & (chrom_data["pos_bp"] <= window_end)]
            
            if len(window_data) >= threshold:
                if method == 'linear':
                    rate = calculate_recombination_rate_linear(window_data, window_size)
                else:
                    rate = calculate_recombination_rate_maxmin(window_data, window_size)
                
                if rate < 0:
                    negative_count += 1
                    rate = 0 
            else:
                rate = np.nan
            
            recombination_data.append([chrom, window_start, window_end, rate])

    recombination_rates = pd.DataFrame(recombination_data, columns=["chr", "start_bp", "end_bp", "cM/Mb"])
    return recombination_rates, negative_count

def main():
    parser = ArgumentParser(description="Calculate recombination rates using a sliding window approach.")
    parser.add_argument("linkage_map_file", type=str, help="Linkage map file containing chromosome, position in base pairs, and genetic distance in centiMorgans.")
    parser.add_argument("--window_size", type=int, default=3000000, help="Window size in base pairs (default: 3000000).")
    parser.add_argument("--shift", type=int, default=500000, help="Shift between consecutive windows in base pairs (default: 500000).")
    parser.add_argument("--threshold", type=int, default=5, help="Minimum number of markers in a window (default: 5).")
    parser.add_argument("--method", type=str, choices=['linear', 'maxmin'], default='maxmin', help="Method to calculate recombination rate: 'linear' for linear regression, 'maxmin' for max-min method (default: maxmin).")
    parser.add_argument("-o", "--output_prefix", type=str, default="recombination_rates", help="Prefix for output files.")
    args = parser.parse_args()

    recombination_rates, negative_count = calculate_recombination_rates(
        args.linkage_map_file, args.window_size, args.shift, args.threshold, args.method
    )

    recombination_rates.to_csv(f"{args.output_prefix}_{args.method}.txt", sep="\t", index=False)

    print(f"Results saved to {args.output_prefix}_{args.method}.txt")
    print(f"Number of windows with originally negative recombination rates: {negative_count}")

if __name__ == "__main__":
    main()
