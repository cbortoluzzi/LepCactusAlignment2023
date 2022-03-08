#!/bin/bash



# Author: @cb46


export PATH=/software/team118/phast/bin:$PATH



if [ -z $1 ]; then
	echo "Usage: ./phyloFit.sh <reference.genome> <cactus.alignment>"
	exit -1
fi



REF=$1 # name of species to use as reference
HAL=$2 # name of cactus alignment in HAL format


mkdir -p neutral_model/$REF


# Tree topology that we will use in phyloFit
tree="((((apotomis_turbidana_gca905147355v1,hedya_salicella_gca905404275v1),notocelia_uddmanniana_gca905163555v1),((zeuzera_pyrina_gca907165235v1,zygaena_filipendulae_gca907165275v1),(((blastobasis_adustella_gca907269095v1,blastobasis_lacticolella_gca905147135v1),carposina_sasakii_gca014607495v1),((((((anthocharis_cardamines_gca905404175v1,(pieris_brassicae_gca905147105v1,pieris_napi_gca905475465v1)),(colias_croceus_gca905220415v1,zerene_cesonia_gca012273895v2)),leptidea_sinapis_gca905404315v1),((((((aricia_agestis_gca905147365v1,cyaniris_semiargus_gca905187585v1),(lysandra_bellargus_gca905333045v1,lysandra_coridon_gca905220515v1)),plebejus_argus_gca905404155v1),(celastrina_argiolus_gca905187575v1,glaucopsyche_alexis_gca905404095v1)),lycaena_phlaeas_gca905333005v1),(((((boloria_selene_gca905231865v2,fabriciana_adippe_gca905404265v1),limenitis_camilla_gca905147385v1),((((inachis_io_gca905147045v1,nymphalis_urticae_gca905147175v1),nymphalis_polychloros_gca905220585v1),(vanessa_atalanta_gca905147765v1,vanessa_cardui_gca905220365v1)),(melitaea_cinxia_gca905220565v1,mellicta_athalia_gca905220545v1))),((maniola_hyperantus_gca902806685v2,maniola_jurtina_gca905333055v1),pararge_aegeria_gca905163445v1)),danaus_plexippus_gca018135715v1))),(erynnis_tages_gca905147235v1,(hesperia_comma_gca905404135v1,ochlodes_sylvanus_gca905404295v1))),((cnaphalocrocis_medinalis_gca014851415v1,(endotricha_flammealis_gca905163395v1,ephestia_elutella_gca018467065v1)),((habrosyne_pyritoides_gca907165245v1,thyatira_batis_gca905147785v1),((((((biston_betularia_gca905404145v1,erannis_defoliaria_gca905404285v1),ectropis_grisescens_gca017562165v1),((crocallis_elinguaria_gca907269065v1,ennomos_fuscantarius_gca905220475v1),hylaea_fasciaria_gca905147375v1)),idaea_aversata_gca907269075v1),((bombyx_mori_gca014905235v2,((deilephila_porcellus_gca905220455v1,hemaris_fuciformis_gca907164795v1),((laothoe_populi_gca905220505v1,mimas_tiliae_gca905332985v1),manduca_sexta_gca014839805v1))),dendrolimus_punctatus_gca012273795v1)),((clostera_curtula_gca905475355v1,((notodonta_dromedarius_gca905147325v1,(pheosia_gnoma_gca905404115v1,pheosia_tremula_gca905333125v1)),phalera_bucephala_gca905147815v1)),(((euproctis_similis_gca905147225v1,lymantria_monacha_gca905163515v1),((hypena_proboscidalis_gca905147285v1,spilosoma_lubricipeda_gca905220595v1),(laspeyria_flexula_gca905147015v1,schrankia_costaestrigalis_gca905475405v1))),(((amphipyra_tragopoginis_gca905220435v1,craniophora_ligustri_gca905163465v1),(((((atethmia_centrago_gca905333075v2,cosmia_trapezina_gca905163495v1),(((noctua_fimbriata_gca905163415v1,noctua_pronuba_gca905220335v1),xestia_xanthographa_gca905147715v1),ochropleura_plecta_gca905475445v1)),phlogophora_meticulosa_gca905147745v1),((hecatera_dysodea_gca905332915v1,mamestra_brassicae_gca905163435v1),mythimna_impura_gca905147345v1)),(spodoptera_exigua_gca011316535v1,(spodoptera_frugiperda_gca012979215v2,spodoptera_litura_gca002706865v2)))),(abrostola_tripartita_gca905340225v1,(autographa_gamma_gca905146925v1,autographa_pulchrina_gca905475315v1)))))))))))),tinea_trinotella_gca905220615v1)"


# Extract 4d sites: since we have an alignment for each contig, we will estimate a nonconserved (neutral) model for each of them and then calculate the average nonconserved model
# A separate nonconserved model will be estimated for the W and Z chromosome
echo "Extract 4d sites from the alignment"
for f in sequences/$REF/*.maf
do
	chr=$(basename $f .maf | cut -f2 -d'.')
	# Obtain 4-fold generate sites
	msa_view sequences/$REF/$REF.$chr.maf --in-format MAF --4d --features annotation/$REF/$REF.$chr.CDS.gff3 --out-format SS > neutral_model/$REF/$REF.$chr.codons.ss
	msa_view neutral_model/$REF/$REF.$chr.codons.ss --in-format SS --out-format SS --tuple-size 1 > neutral_model/$REF/$REF.$chr.sites.ss
	# Estimate the nonconserved (neutral) model for each chromosome
	phyloFit --tree $tree --subst-mod REV --EM --msa-format SS neutral_model/$REF/$REF.$chr.sites.ss --out-root neutral_model/$REF/$REF.$chr.nonconserved-4d
done


# Obtain the average nonconserved model
#echo "Obtain the average nonconserved model"
#re='^[0-9]+$'
#autosomes=`cat sequences/$REF/$REF.bed | while read contig size;do if [[ $contig =~ $re ]]; then echo $contig; fi;done`
#mod=`for chr in $autosomes;do ls neutral_model/$REF/$REF.$chr.nonconserved-4d.mod;done`
#phyloBoot --read-mods $mod --output-average neutral_model/$REF/$REF.ave.nonconserved-4d.mod

#rm neutral_model/$REF/*.ss
