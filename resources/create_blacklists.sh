### Data downloaded with R package excluderanges (PMID: 37067481)
# T2T-CHM13v2.0
Rscript excluderanges.R AH107304 T2T-CHM13v2.0_excluded.bed
pigz T2T-CHM13v2.0_excluded.bed

# GRCm39 (mm39)
Rscript excluderanges.R AH107321 GRCm39_excluded.bed
sed -i 's/^chr//' GRCm39_excluded.bed
pigz GRCm39_excluded.bed