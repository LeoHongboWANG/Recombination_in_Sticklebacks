rule all:
    input:
        expand("LG{index}_order.mapped", index=range(1, 22))

rule extract_lg:
    input: "data_f.call"
    output:
        lg_file="LG{index}_f.call",
        nums_file="LG{index}_nums.txt"
    shell:
        """
        head -n 6 {input} > {output.lg_file}
        grep -w LG{wildcards.index} {input} >> {output.lg_file}
        grep -w LG{wildcards.index} {input} | wc -l > {output.nums_file}
        """

rule create_phy:
    input: "LG{index}_nums.txt"
    output: "LG{index}.phy"
    shell:
        """
        a=$(cat {input})
        seq 1 $((a - 1)) > {output}
        """

rule order_markers:
    input:
        lg_file="LG{index}_f.call",
        phy_file="LG{index}.phy"
    output: "LG{index}.map"
    log: "LG{index}_order_markers.log"
    shell:
        """
        java -cp /scratch/project_2005070/Linkagemap/Lepmap3/bin/ OrderMarkers2 \
            numThreads=32 \
            data={input.lg_file} \
            evaluateOrder={input.phy_file} \
            outputPhasedData=1 \
            improveOrder=0 \
            sexAveraged=1 > {output} 2> {log}
        """

rule map_snps:
    input:
        lg_file="LG{index}_f.call",
        map_file="LG{index}.map"
    output: "LG{index}_order.mapped"
    shell:
        """
        awk '(NR>=7)' {input.lg_file} | cut -f 1,2 > LG{wildcards.index}.snps.txt
        awk -v FS="\\t" -v OFS="\\t" '(NR==FNR){{s[NR-1]=$0}}(NR!=FNR){{if ($1 in s) $1=s[$1];print}}' \
            LG{wildcards.index}.snps.txt {input.map_file} > {output}
        """
