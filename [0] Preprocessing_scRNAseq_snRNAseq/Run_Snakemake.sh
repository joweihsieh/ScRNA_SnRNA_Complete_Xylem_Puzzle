#### NOTE
# The snakemake pipeline was previously established by Shang-Che Kuo.
# Publication: Tung, CC., Kuo, SC., Yang, CL. et al. Single-cell transcriptomics unveils xylem cell development and evolution. Genome Biol 24, 3 (2023). 
# https://doi.org/10.1186/s13059-022-02845-1
# https://github.com/Woodformation1136/SingleCell


####
#cd /home/f06b22037/SSD2/JW/1136project_SingleCell
#conda activate snakemake




nohup snakemake --use-conda -c 48 scRNATung_snRNABio2 > QZ_Tung_SDXsnRNA_2_snakemake.log &
nohup snakemake --use-conda -c 48 scRNATung_snRNABio3 > QZ_Tung_SDXsnRNA_3_snakemake.log &
nohup snakemake --use-conda -c 48 snRNA_Bio2_Bio3 > QZ_Tung_SDXsnRNA_6_snakemake.log &

nohup snakemake --use-conda -c 48 scRNATung_ProtosnRNABio1 > QZ_scRNATung_ProtosnRNABio1_snakemake.log &
nohup snakemake --use-conda -c 48 scRNATung_ProtosnRNABio2 > QZ_scRNATung_ProtosnRNABio2_snakemake.log &

nohup snakemake --use-conda -c 48 snRNABio2_ProtosnRNABio1 > QZ_snRNABio2_ProtosnRNABio1_snakemake.log &
nohup snakemake --use-conda -c 48 snRNABio3_ProtosnRNABio2 > QZ_snRNABio3_ProtosnRNABio2_snakemake.log &
nohup snakemake --use-conda -c 48 snRNABio3_ProtosnRNABio1 > QZ_snRNABio3_ProtosnRNABio1_snakemake.log &
nohup snakemake --use-conda -c 48 snRNABio2_ProtosnRNABio2 > QZ_snRNABio2_ProtosnRNABio2_snakemake.log &

nohup snakemake --use-conda -c 48 ProtosnRNABio1_ProtosnRNABio2 > QZ_ProtosnRNABio1_ProtosnRNABio2_snakemake.log &
