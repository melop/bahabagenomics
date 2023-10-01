#zcat ~/../../software/jbrowse2/bahaha_taipingensis/markers_XY.bed.gz | awk '{print $1"\t"$2"\t"$3}' > marker_XY.regions
#zcat ~/../../software/jbrowse2/bahaha_taipingensis/markers_ZW.bed.gz | awk '{print $1"\t"$2"\t"$3}' > marker_ZW.regions

bcftools view -R marker_XY.regions -s"hap.F0,hap.F1,hap.M0,hap.M1" ann.merged.vcf.gz | grep HIGH  > chr10_high_impact_XY.vcf
bcftools view -R marker_ZW.regions -s"hap.F0,hap.F1,hap.M0,hap.M1" ann.merged.vcf.gz | grep HIGH  > chr10_high_impact_ZW.vcf
