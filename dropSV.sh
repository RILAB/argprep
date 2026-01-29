#!/bin/bash


#DEFAULT VALUES
CUTOFF=9101264
SRC="."

#DISPLAY HELP
# Function to display the help message
display_help() {
    echo "Usage: $0 [OPTION]..."
    echo "This script removes large indels for later parsing with GATK."
    echo ""
    echo "Options:"
    echo "  -h, --help    Display this help message and exit"
    echo "  -c  --cutoff          INDEL size to remove (default: $CUTOFF)"
    echo "  -d  --directory          Source directory containing gVCF files (default: current directory"
    # Add more options as needed
}

# Check if no arguments are provided
if [ "$#" -eq 0 ]; then
    display_help
    exit 1 # Exit with an error status
fi


#GET NEW VALUES FROM COMMAND LINE
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      display_help
      exit 0
      ;;
    -d|--directory)
      SRC="$2"
      shift # past argument
      shift # past value
      ;;
    -c|--cutoff)
      CUTOFF="$2"
      shift # past argument
      shift # past value
      ;;
    -*|--*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

hash bcftools tabix 2>/dev/null
 if [ "$?" -ne 0 ]; then
    echo "Error: Please ensure bcftools and tabix are installed and available."
    exit 1
fi

cd $SRC || { echo "ERROR: Failed to change to $SRC. Please check directory exists"; exit 1; }
echo $(ls *.gvcf.gz | wc -l) "gVCF files found in" $SRC

mkdir -p cleangVCF
cd cleangVCF

for gvcf in $SRC/*.gvcf.gz; do
  bcftools view "$gvcf" | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' |
  awk -v cutoff=$CUTOFF -v file="$gvcf" '{d=length($3)-length($4); if(d<0)d=-d; if(d>cutoff) print file"\t"$1"\t"$2"\t"d}'
done | sort -k1,1 > super_indels.txt

cut -f1 super_indels.txt | sort -u > need_filter.txt

for gvcf in $SRC/*.gvcf.gz; do
    base=$(basename "$gvcf")
    echo "Processing $base"
    if grep -qx "$gvcf" need_filter.txt; then
		lines_before=$(zgrep -cv "^#" "$gvcf")
        awk -v f="$gvcf" '$1==f {print $2"\t"$3}' super_indels.txt > bad_sites.tmp
        bcftools view -T ^bad_sites.tmp -Oz -o "$base" "$gvcf"
        tabix -p vcf "$base"
        lines_after=$(zgrep -cv "^#" "$base")
        echo "  before: $lines_before, after: $lines_after, removed: $((lines_before - lines_after))"
    else
        cp "$gvcf" .		
		    cp "$gvcf.tbi" .
        echo "  no large indels found, file copied without changes."
    fi
done

cat super_indels.txt | wc -l | awk '{print "Total super large indels identified across all gVCFs: "$1}'
cat need_filter.txt | wc -l | awk '{print "Total gVCFs with super large indels removed: "$1}'

rm -f bad_sites.tmp need_filter.txt super_indels.txt
echo "All done. Cleaned gVCFs are in" $(pwd)