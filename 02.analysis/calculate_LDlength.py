import gzip
import sys

def process_file(filename):
    with gzip.open(filename, 'rt') as f:
        header = f.readline()  # Skip header line
        first_line = f.readline().strip().split('\t')
        max_r2 = float(first_line[1])
        
        half_r2 = max_r2 / 2
        dis2value = {}

        for line in f:
            parts = line.strip().split('\t')
            try:
                distance = int(parts[0])
                r2_value = float(parts[1])
                dis2value[distance] = r2_value
            except ValueError as e:
                print(f"Skipping line due to error: {line}")
                print(e)
                continue

    found = False
    for key in sorted(dis2value.keys()):
        next_key = key + 1
        if next_key in dis2value:
            current_value = dis2value[key]
            next_value = dis2value[next_key]
            if current_value >= half_r2 and next_value < half_r2:
                print(f"Processing {filename}")
                print(f"max LD: r2: {max_r2}")
                print(f"half LD: r2: {half_r2}\tLD length: {key}")
                found = True
                break
        else:
            print(f"Missing value for distance: {next_key}")
    
    if not found:
        print(f"Could not find a distance where r2 falls below half of max LD for {filename}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python calculate_LDlength.py <filename>")
        sys.exit(1)

    input_file = sys.argv[1]
    process_file(input_file)
