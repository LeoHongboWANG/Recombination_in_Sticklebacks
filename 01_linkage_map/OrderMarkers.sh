
#!/bin/bash 
 
java -cp /scratch/project_2005070/Linkagemap/Lepmap3/bin/ ParentCall2  data=eco.txt vcfFile=eco.vcf removeNonInformative=1 XLimit=2 halfSibs=1  > data.call
java -cp /scratch/project_2005070/Linkagemap/Lepmap3/bin/ Filtering2 data=data.call removeNonInformative=1 MAFLimit=0.05 missingLimit=0.1 dataTolerance=0.001 > data_f.call



for i in {1..21}; do

    head -n 6 data_f.call > LG${i}_f.call
    grep -w LG${i} data_f.call >> LG${i}_f.call
    
    # Step 2
    grep -w LG${i} data_f.call | wc -l > LG${i}_nums.txt

    # Step 3
    while read -r a; do
        seq 1 $((a - 1)) > LG${i}.phy
    done < LG${i}_nums.txt

    # Step 4: Run the OrderMarkers2
    java -cp /Lepmap3/bin OrderMarkers2 \
        numThreads=32 \
        data=LG${i}_f.call \
        evaluateOrder=LG${i}.phy \
        useMorgan=1 \
        outputPhasedData=1 \
        improveOrder=0 > LG${i}_diff.map 

    # Step 5
    awk '(NR>=7)' LG${i}_f.call | cut -f 1,2 > LG${i}.snps.txt
    awk -v FS="\t" -v OFS="\t" '(NR==FNR){s[NR-1]=$0}(NR!=FNR){if ($1 in s) $1=s[$1];print}' \
        LG${i}.snps.txt LG${i}_diff.map > LG${i}_order_diff.mapped
done

