#!/bin/bash


# Author: @cb46


export PATH=/software/team118/phast/bin:$PATH


if [ -z $1 ]; then
	echo "Usage: ./neutral_model_ancestral_repeats.sh <reference.genome> <cactus.alignment>"
	exit -1
fi



REF=$1
HAL=$2


mkdir -p method_1 && mkdir -p method_1/mask_mrca && mkdir -p method_1/single_copy_regions


# Extract ancestral genome from the alignment representing the MRCA of all species
hal2fasta $HAL $REF > method_1/mask_mrca/$REF.fa

# Run windowmasker to find additional repeats in the MRCA genome: we run windowmasker with -dust false because we want to exclude low complexity regions
windowmasker -mk_counts -in method_1/mask_mrca/$REF.fa -out method_1/mask_mrca/$REF.counts
windowmasker -ustat method_1/mask_mrca/$REF.counts -in method_1/mask_mrca/$REF.fa -outfmt interval -out method_1/mask_mrca/$REF.wm -dust false

# Run RNAscan-SE to identify tRNAs
tRNAscan-SE -E -o method_1/mask_mrca/$REF.tRNA.txt method_1/mask_mrca/$REF.fa

# Exclude ancestral repeats found in tRNAs
python3 exclude_tRNAs.py --wm method_1/mask_mrca/$REF.sm.wm --tRNA method_1/mask_mrca/$REF.tRNA.txt

# Get target genomes of alignment
halStats --genomes $HAL | tr -s ' '  '\n' | grep -v '^Anc' > method_1/genomes.txt
targetGenomes=`cat method_1/genomes.txt | while read species;do echo $species;done  | paste -s -d, -`

# Check whether the candidate repeat sequences are single copy regions or not
for f in method_1/mask_mrca/*.bed
do
	cat $f | while read contig start end length
	do
		halSingleCopyRegionsExtract --refSequence $contig --start $start --length $length --targetGenomes $targetGenomes $HAL $REF >> method_1/single_copy_regions/$contig.singlecopy
	done
done

# Extract a MAF for each single copy ancestral contig
for f in method_1/single_copy_regions/*.singlecopy
do
	cat $f | while read contig start end
	do
		hal2maf --refGenome $REF --refSequence $contig --start $start --length `(expr $end - $start)` --onlyOrthologs --targetGenomes $targetGenomes --append $HAL method_1/single_copy_regions/$contig.singlecopy.maf
	done
done

# Prepare input alignment for GERP analysis
# This script will separate the ancestral contig based on the type of contig it aligns to in the other species, being this an autosome, W, or Z chromsome
python3 prepare_alignment.py - -maf method_1/single_copy_regions/ --refGenome $REF

# Run phyloFit to obtain the neutral model
mkdir -p method_1/neutral_model/
for f in ancestral_repeats_*.maf
do
	filename=$(basename $f .maf)
	cat $f | grep -v '#' | grep 'Anc00' | awk 'OFS="\t"{print $2,"Masking","AR",$3,$3+$4,".",".",".","type: Ancestral_repeat"}' > $filename.gff
	phyloFit --tree "((((apotomis_turbidana_gca905147355v1,hedya_salicella_gca905404275v1),notocelia_uddmanniana_gca905163555v1),((zeuzera_pyrina_gca907165235v1,zygaena_filipendulae_gca907165275v1),(((blastobasis_adustella_gca907269095v1,blastobasis_lacticolella_gca905147135v1),carposina_sasakii_gca014607495v1),((((((anthocharis_cardamines_gca905404175v1,(pieris_brassicae_gca905147105v1,pieris_napi_gca905475465v1)),(colias_croceus_gca905220415v1,zerene_cesonia_gca012273895v2)),leptidea_sinapis_gca905404315v1),((((((aricia_agestis_gca905147365v1,cyaniris_semiargus_gca905187585v1),(lysandra_bellargus_gca905333045v1,lysandra_coridon_gca905220515v1)),plebejus_argus_gca905404155v1),(celastrina_argiolus_gca905187575v1,glaucopsyche_alexis_gca905404095v1)),lycaena_phlaeas_gca905333005v1),(((((boloria_selene_gca905231865v2,fabriciana_adippe_gca905404265v1),limenitis_camilla_gca905147385v1),((((inachis_io_gca905147045v1,nymphalis_urticae_gca905147175v1),nymphalis_polychloros_gca905220585v1),(vanessa_atalanta_gca905147765v1,vanessa_cardui_gca905220365v1)),(melitaea_cinxia_gca905220565v1,mellicta_athalia_gca905220545v1))),((maniola_hyperantus_gca902806685v2,maniola_jurtina_gca905333055v1),pararge_aegeria_gca905163445v1)),danaus_plexippus_gca018135715v1))),(erynnis_tages_gca905147235v1,(hesperia_comma_gca905404135v1,ochlodes_sylvanus_gca905404295v1))),((cnaphalocrocis_medinalis_gca014851415v1,(endotricha_flammealis_gca905163395v1,ephestia_elutella_gca018467065v1)),((habrosyne_pyritoides_gca907165245v1,thyatira_batis_gca905147785v1),((((((biston_betularia_gca905404145v1,erannis_defoliaria_gca905404285v1),ectropis_grisescens_gca017562165v1),((crocallis_elinguaria_gca907269065v1,ennomos_fuscantarius_gca905220475v1),hylaea_fasciaria_gca905147375v1)),idaea_aversata_gca907269075v1),((bombyx_mori_gca014905235v2,((deilephila_porcellus_gca905220455v1,hemaris_fuciformis_gca907164795v1),((laothoe_populi_gca905220505v1,mimas_tiliae_gca905332985v1),manduca_sexta_gca014839805v1))),dendrolimus_punctatus_gca012273795v1)),((clostera_curtula_gca905475355v1,((notodonta_dromedarius_gca905147325v1,(pheosia_gnoma_gca905404115v1,pheosia_tremula_gca905333125v1)),phalera_bucephala_gca905147815v1)),(((euproctis_similis_gca905147225v1,lymantria_monacha_gca905163515v1),((hypena_proboscidalis_gca905147285v1,spilosoma_lubricipeda_gca905220595v1),(laspeyria_flexula_gca905147015v1,schrankia_costaestrigalis_gca905475405v1))),(((amphipyra_tragopoginis_gca905220435v1,craniophora_ligustri_gca905163465v1),(((((atethmia_centrago_gca905333075v2,cosmia_trapezina_gca905163495v1),(((noctua_fimbriata_gca905163415v1,noctua_pronuba_gca905220335v1),xestia_xanthographa_gca905147715v1),ochropleura_plecta_gca905475445v1)),phlogophora_meticulosa_gca905147745v1),((hecatera_dysodea_gca905332915v1,mamestra_brassicae_gca905163435v1),mythimna_impura_gca905147345v1)),(spodoptera_exigua_gca011316535v1,(spodoptera_frugiperda_gca012979215v2,spodoptera_litura_gca002706865v2)))),(abrostola_tripartita_gca905340225v1,(autographa_gamma_gca905146925v1,autographa_pulchrina_gca905475315v1)))))))))))),tinea_trinotella_gca905220615v1)" --subst-mod REV --features $filename.gff --do-cats AR $f
	mv phyloFit.AR.mod method_1/neutral_model/$filename.mod
done

