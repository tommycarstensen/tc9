bsub -q normal -M3000000 -R'select[mem>3000] rusage[mem=3000]' -P agv -o ped2vcf.out -e ped2vcf.err python ~/github/tc9/ped2vcf.py
