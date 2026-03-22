#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import subprocess
import time
start_time = time.time()

studies = {
	# "Alasoo_1_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_IFNg+Salmonella.all.tsv.gz",
	# "Alasoo_1_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/tx/Alasoo_2018_tx_macrophage_IFNg+Salmonella.all.tsv.gz",
	# "Alasoo_2_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_IFNg.all.tsv.gz",
	# "Alasoo_2_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/tx/Alasoo_2018_tx_macrophage_IFNg.all.tsv.gz",
	# "Alasoo_3_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_Salmonella.all.tsv.gz",
	# "Alasoo_3_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/tx/Alasoo_2018_tx_macrophage_Salmonella.all.tsv.gz",
	# "Alasoo_4_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_naive.all.tsv.gz",
	# "Alasoo_4_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/tx/Alasoo_2018_tx_macrophage_naive.all.tsv.gz",

	# "BLUEPRINT_1_ge": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/ge/BLUEPRINT_PE_ge_T-cell.all.tsv.gz",
	# "BLUEPRINT_1_iu": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/tx/BLUEPRINT_PE_tx_T-cell.all.tsv.gz",
	# "BLUEPRINT_2_ge": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/ge/BLUEPRINT_SE_ge_monocyte.all.tsv.gz",
	# "BLUEPRINT_2_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/tx/BLUEPRINT_SE_tx_monocyte.all.tsv.gz",
	# "BLUEPRINT_3_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/ge/BLUEPRINT_SE_ge_neutrophil.all.tsv.gz",
	# "BLUEPRINT_3_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/tx/BLUEPRINT_SE_tx_neutrophil.all.tsv.gz",
	
	# "Bossini_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Bossini-Castillo_2019/ge/Bossini-Castillo_2019_ge_Treg_naive.all.tsv.gz",
	# "Bossini_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Bossini-Castillo_2019/tx/Bossini-Castillo_2019_tx_Treg_naive.all.tsv.gz",
	
	# "Braineaq2_1_ge": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq2/ge/Braineac2_ge_putamen.all.tsv.gz",
	# "Braineaq2_1_iu": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq2/tx/Braineac2_tx_putamen.all.tsv.gz",
	# "Braineaq2_2_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq2/ge/Braineac2_ge_substantia_nigra.all.tsv.gz",
	# "Braineaq2_2_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq2/tx/Braineac2_tx_substantia_nigra.all.tsv.gz",
	
	# "CAP_1_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CAP/ge/CAP_ge_LCL_naive.all.tsv.gz",
	# "CAP_1_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CAP/tx/CAP_tx_LCL_naive.all.tsv.gz",
	# "CAP_2_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CAP/ge/CAP_ge_LCL_statin.all.tsv.gz",
	# "CAP_2_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CAP/tx/CAP_tx_LCL_statin.all.tsv.gz",
	
	# "CommonMind_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CommonMind/ge/CommonMind_ge_DLPFC_naive.all.tsv.gz",
	# "CommonMind_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CommonMind/tx/CommonMind_tx_DLPFC_naive.all.tsv.gz",
	
	# "FUSION_1_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/FUSION/ge/FUSION_ge_adipose_naive.all.tsv.gz",
	# "FUSION_1_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/FUSION/tx/FUSION_tx_adipose_naive.all.tsv.gz",
	# "FUSION_2_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/FUSION/ge/FUSION_ge_muscle_naive.all.tsv.gz",
	# "FUSION_2_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/FUSION/tx/FUSION_tx_muscle_naive.all.tsv.gz",
	
	# "GENCORD_1_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/ge/GENCORD_ge_LCL.all.tsv.gz",
	# "GENCORD_1_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/tx/GENCORD_tx_LCL.all.tsv.gz",
	# "GENCORD_2_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/ge/GENCORD_ge_T-cell.all.tsv.gz",
	# "GENCORD_2_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/tx/GENCORD_tx_T-cell.all.tsv.gz",
	# "GENCORD_3_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/ge/GENCORD_ge_fibroblast.all.tsv.gz",
	# "GENCORD_3_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/tx/GENCORD_tx_fibroblast.all.tsv.gz",
	
	# "GEUVADIS_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GEUVADIS/ge/GEUVADIS_ge_LCL.all.tsv.gz",
	# "GEUVADIS_tx":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GEUVADIS/tx/GEUVADIS_tx_LCL.all.tsv.gz",

	# "HipSci_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/HipSci/ge/HipSci_ge_iPSC.all.tsv.gz",
	# "HipSci_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/HipSci/tx/HipSci_tx_iPSC.all.tsv.gz",

	# "Lepik_ge":			"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Lepik_2017/ge/Lepik_2017_ge_blood.all.tsv.gz",
	# "Lepik_iu":			"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Lepik_2017/tx/Lepik_2017_tx_blood.all.tsv.gz",

	# "Peng_ge":			"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Peng_2018/ge/Peng_2018_ge_placenta_naive.all.tsv.gz",
	# "Peng_iu":			"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Peng_2018/tx/Peng_2018_tx_placenta_naive.all.tsv.gz",

	# "ROSMAP_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/ROSMAP/ge/ROSMAP_ge_brain_naive.all.tsv.gz",
	# "ROSMAP_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/ROSMAP/tx/ROSMAP_tx_brain_naive.all.tsv.gz",
	
	# "Schwartzentruber_ge":"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schwartzentruber_2018/ge/Schwartzentruber_2018_ge_sensory_neuron.all.tsv.gz",
	# "Schwartzentruber_iu":"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schwartzentruber_2018/tx/Schwartzentruber_2018_tx_sensory_neuron.all.tsv.gz",
	
	# "Steinberg_1_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Steinberg_2020/ge/Steinberg_2020_ge_high_grade_cartilage_naive.all.tsv.gz",
	# "Steinberg_1_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Steinberg_2020/tx/Steinberg_2020_tx_high_grade_cartilage_naive.all.tsv.gz",
	# "Steinberg_2_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Steinberg_2020/ge/Steinberg_2020_ge_low_grade_cartilage_naive.all.tsv.gz",
	# "Steinberg_2_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Steinberg_2020/tx/Steinberg_2020_tx_low_grade_cartilage_naive.all.tsv.gz",
	# "Steinberg_3_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Steinberg_2020/ge/Steinberg_2020_ge_synovium_naive.all.tsv.gz",
	# "Steinberg_3_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Steinberg_2020/tx/Steinberg_2020_tx_synovium_naive.all.tsv.gz",
	
	# "TwinsUK_1_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/TwinsUK/ge/TwinsUK_ge_LCL.all.tsv.gz",
	# "TwinsUK_1_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/TwinsUK/tx/TwinsUK_tx_LCL.all.tsv.gz",
	# "TwinsUK_2_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/TwinsUK/ge/TwinsUK_ge_blood.all.tsv.gz",
	# "TwinsUK_2_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/TwinsUK/tx/TwinsUK_tx_blood.all.tsv.gz",
	# "TwinsUK_3_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/TwinsUK/ge/TwinsUK_ge_fat.all.tsv.gz",
	# "TwinsUK_3_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/TwinsUK/tx/TwinsUK_tx_fat.all.tsv.gz",
	# "TwinsUK_4_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/TwinsUK/ge/TwinsUK_ge_skin.all.tsv.gz",
	# "TwinsUK_4_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/TwinsUK/tx/TwinsUK_tx_skin.all.tsv.gz",
	
	# "Young_ge":			"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Young_2019/ge/Young_2019_ge_microglia_naive.all.tsv.gz",
	"Young_iu":			"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Young_2019/tx/Young_2019_tx_microglia_naive.all.tsv.gz",
	
	"iPSCORE_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/iPSCORE/ge/iPSCORE_ge_iPSC.all.tsv.gz",
	"iPSCORE_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/iPSCORE/tx/iPSCORE_tx_iPSC.all.tsv.gz",

	"BrainSeq_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq/ge/BrainSeq_ge_brain.all.tsv.gz",
	"BrainSeq_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq/tx/BrainSeq_tx_brain.all.tsv.gz",

	"Nedelec_1_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Nedelec_2016/ge/Nedelec_2016_ge_macrophage_Listeria.all.tsv.gz",
	"Nedelec_1_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Nedelec_2016/tx/Nedelec_2016_tx_macrophage_Listeria.all.tsv.gz",
	"Nedelec_2_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Nedelec_2016/ge/Nedelec_2016_ge_macrophage_Salmonella.all.tsv.gz",
	"Nedelec_2_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Nedelec_2016/tx/Nedelec_2016_tx_macrophage_Salmonella.all.tsv.gz",
	"Nedelec_3_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Nedelec_2016/ge/Nedelec_2016_ge_macrophage_naive.all.tsv.gz",
	"Nedelec_3_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Nedelec_2016/tx/Nedelec_2016_tx_macrophage_naive.all.tsv.gz",
	
	"PhLiPS_1_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/PhLiPS/ge/PhLiPS_ge_HLC.all.tsv.gz",
	"PhLiPS_1_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/PhLiPS/tx/PhLiPS_tx_HLC.all.tsv.gz",
	"PhLiPS_2_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/PhLiPS/ge/PhLiPS_ge_iPSC.all.tsv.gz",
	"PhLiPS_2_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/PhLiPS/tx/PhLiPS_tx_iPSC.all.tsv.gz",

	"Quach_1_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/ge/Quach_2016_ge_monocyte_IAV.all.tsv.gz",
	"Quach_1_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/tx/Quach_2016_tx_monocyte_IAV.all.tsv.gz",
	"Quach_2_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/ge/Quach_2016_ge_monocyte_LPS.all.tsv.gz",
	"Quach_2_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/tx/Quach_2016_tx_monocyte_LPS.all.tsv.gz",
	"Quach_3_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/ge/Quach_2016_ge_monocyte_Pam3CSK4.all.tsv.gz",
	"Quach_3_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/tx/Quach_2016_tx_monocyte_Pam3CSK4.all.tsv.gz",
	"Quach_4_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/ge/Quach_2016_ge_monocyte_R848.all.tsv.gz",
	"Quach_4_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/tx/Quach_2016_tx_monocyte_R848.all.tsv.gz",
	"Quach_5_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/ge/Quach_2016_ge_monocyte_naive.all.tsv.gz",
	"Quach_5_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Quach_2016/tx/Quach_2016_tx_monocyte_naive.all.tsv.gz",

	"Schmiedel_1_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_B-cell_naive.all.tsv.gz",
	"Schmiedel_1_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_B-cell_naive.all.tsv.gz",
	"Schmiedel_2_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_CD4_T-cell_anti-CD3-CD28.all.tsv.gz",
	"Schmiedel_2_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_CD4_T-cell_anti-CD3-CD28.all.tsv.gz",
	"Schmiedel_3_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_CD4_T-cell_naive.all.tsv.gz",
	"Schmiedel_3_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_CD4_T-cell_naive.all.tsv.gz",
	"Schmiedel_4_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_CD8_T-cell_anti-CD3-CD28.all.tsv.gz",
	"Schmiedel_4_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_CD8_T-cell_anti-CD3-CD28.all.tsv.gz",
	"Schmiedel_5_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_CD8_T-cell_naive.all.tsv.gz",
	"Schmiedel_5_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_CD8_T-cell_naive.all.tsv.gz",
	"Schmiedel_6_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_NK-cell_naive.all.tsv.gz",
	"Schmiedel_6_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_NK-cell_naive.all.tsv.gz",
	"Schmiedel_7_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_Tfh_memory.all.tsv.gz",
	"Schmiedel_7_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_Tfh_memory.all.tsv.gz",
	"Schmiedel_8_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_Th1-17_memory.all.tsv.gz",
	"Schmiedel_8_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_Th1-17_memory.all.tsv.gz",
	"Schmiedel_9_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_Th17_memory.all.tsv.gz",
	"Schmiedel_9_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_Th17_memory.all.tsv.gz",
	"Schmiedel_10_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_Th1_memory.all.tsv.gz",
	"Schmiedel_10_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_Th1_memory.all.tsv.gz",
	"Schmiedel_11_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_Th2_memory.all.tsv.gz",
	"Schmiedel_11_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_Th2_memory.all.tsv.gz",
	"Schmiedel_12_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_Treg_memory.all.tsv.gz",
	"Schmiedel_12_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_Treg_memory.all.tsv.gz",
	"Schmiedel_13_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_Treg_naive.all.tsv.gz",
	"Schmiedel_13_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_Treg_naive.all.tsv.gz",
	"Schmiedel_14_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_monocyte_CD16_naive.all.tsv.gz",
	"Schmiedel_14_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_monocyte_CD16_naive.all.tsv.gz",
	"Schmiedel_15_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/ge/Schmiedel_2018_ge_monocyte_naive.all.tsv.gz",
	"Schmiedel_15_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Schmiedel_2018/tx/Schmiedel_2018_tx_monocyte_naive.all.tsv.gz",

	# "GTEx_1_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_LCL.all.tsv.gz",
	# "GTEx_1_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_LCL.all.tsv.gz",
	# "GTEx_2_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_adipose_subcutaneous.all.tsv.gz",
	# "GTEx_2_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_adipose_subcutaneous.all.tsv.gz",
	# "GTEx_3_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_adipose_visceral.all.tsv.gz",
	# "GTEx_3_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_adipose_visceral.all.tsv.gz",
	# "GTEx_4_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_adrenal_gland.all.tsv.gz",
	# "GTEx_4_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_adrenal_gland.all.tsv.gz",
	# "GTEx_5_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_artery_aorta.all.tsv.gz",
	# "GTEx_5_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_artery_aorta.all.tsv.gz",
	# "GTEx_6_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_artery_coronary.all.tsv.gz",
	# "GTEx_6_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_artery_coronary.all.tsv.gz",
	# "GTEx_7_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_artery_tibial.all.tsv.gz",
	# "GTEx_7_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_artery_tibial.all.tsv.gz",
	# "GTEx_8_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_blood.all.tsv.gz",
	# "GTEx_8_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_blood.all.tsv.gz",
	# "GTEx_9_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_amygdala.all.tsv.gz",
	# "GTEx_9_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_amygdala.all.tsv.gz",
	# "GTEx_10_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_anterior_cingulate_cortex.all.tsv.gz",
	# "GTEx_10_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_anterior_cingulate_cortex.all.tsv.gz",
	# "GTEx_11_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_caudate.all.tsv.gz",
	# "GTEx_11_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_caudate.all.tsv.gz",
	# "GTEx_12_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_cerebellar_hemisphere.all.tsv.gz",
	# "GTEx_12_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_cerebellar_hemisphere.all.tsv.gz",
	# "GTEx_13_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_cerebellum.all.tsv.gz",
	# "GTEx_13_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_cerebellum.all.tsv.gz",
	# "GTEx_14_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_cortex.all.tsv.gz",
	# "GTEx_14_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_cortex.all.tsv.gz",
	# "GTEx_15_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_frontal_cortex.all.tsv.gz",
	# "GTEx_15_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_frontal_cortex.all.tsv.gz",
	# "GTEx_16_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_hippocampus.all.tsv.gz",
	# "GTEx_16_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_hippocampus.all.tsv.gz",
	# "GTEx_17_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_hypothalamus.all.tsv.gz",
	# "GTEx_17_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_hypothalamus.all.tsv.gz",
	# "GTEx_18_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_nucleus_accumbens.all.tsv.gz",
	# "GTEx_18_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_nucleus_accumbens.all.tsv.gz",
	# "GTEx_19_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_putamen.all.tsv.gz",
	# "GTEx_19_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_putamen.all.tsv.gz",
	# "GTEx_20_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_spinal_cord.all.tsv.gz",
	# "GTEx_20_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_spinal_cord.all.tsv.gz",
	# "GTEx_21_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_brain_substantia_nigra.all.tsv.gz",
	# "GTEx_21_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_brain_substantia_nigra.all.tsv.gz",
	# "GTEx_22_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_breast.all.tsv.gz",
	# "GTEx_22_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_breast.all.tsv.gz",
	# "GTEx_23_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_colon_sigmoid.all.tsv.gz",
	# "GTEx_23_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_colon_sigmoid.all.tsv.gz",
	# "GTEx_24_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_colon_transverse.all.tsv.gz",
	# "GTEx_24_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_colon_transverse.all.tsv.gz",
	# "GTEx_25_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_esophagus_gej.all.tsv.gz",
	# "GTEx_25_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_esophagus_gej.all.tsv.gz",
	# "GTEx_26_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_esophagus_mucosa.all.tsv.gz",
	# "GTEx_26_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_esophagus_mucosa.all.tsv.gz",
	# "GTEx_27_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_esophagus_muscularis.all.tsv.gz",
	# "GTEx_27_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_esophagus_muscularis.all.tsv.gz",
	# "GTEx_28_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_fibroblast.all.tsv.gz",
	# "GTEx_28_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_fibroblast.all.tsv.gz",
	# "GTEx_29_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_heart_atrial_appendage.all.tsv.gz",
	# "GTEx_29_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_heart_atrial_appendage.all.tsv.gz",
	# "GTEx_30_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_heart_left_ventricle.all.tsv.gz",
	# "GTEx_30_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_heart_left_ventricle.all.tsv.gz",
	# "GTEx_31_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_kidney_cortex.all.tsv.gz",
	# "GTEx_31_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_kidney_cortex.all.tsv.gz",
	# "GTEx_32_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_liver.all.tsv.gz",
	# "GTEx_32_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_liver.all.tsv.gz",
	# "GTEx_33_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_lung.all.tsv.gz",
	# "GTEx_33_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_lung.all.tsv.gz",
	# "GTEx_34_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_minor_salivary_gland.all.tsv.gz",
	# "GTEx_34_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_minor_salivary_gland.all.tsv.gz",
	# "GTEx_35_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_muscle.all.tsv.gz",
	# "GTEx_35_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_muscle.all.tsv.gz",
	# "GTEx_36_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_nerve_tibial.all.tsv.gz",
	# "GTEx_36_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_nerve_tibial.all.tsv.gz",
	# "GTEx_37_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_ovary.all.tsv.gz",
	# "GTEx_37_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_ovary.all.tsv.gz",
	# "GTEx_38_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_pancreas.all.tsv.gz",
	# "GTEx_38_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_pancreas.all.tsv.gz",
	# "GTEx_39_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_pituitary.all.tsv.gz",
	# "GTEx_39_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_pituitary.all.tsv.gz",
	# "GTEx_40_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_prostate.all.tsv.gz",
	# "GTEx_40_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_prostate.all.tsv.gz",
	# "GTEx_41_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_skin_not_sun_exposed.all.tsv.gz",
	# "GTEx_41_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_skin_not_sun_exposed.all.tsv.gz",
	# "GTEx_42_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_skin_sun_exposed.all.tsv.gz",
	# "GTEx_42_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_skin_sun_exposed.all.tsv.gz",
	# "GTEx_43_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_small_intestine.all.tsv.gz",
	# "GTEx_43_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_small_intestine.all.tsv.gz",
	# "GTEx_44_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_spleen.all.tsv.gz",
	# "GTEx_44_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_spleen.all.tsv.gz",
	# "GTEx_45_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_stomach.all.tsv.gz",
	# "GTEx_45_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_stomach.all.tsv.gz",
	# "GTEx_46_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_testis.all.tsv.gz",
	# "GTEx_46_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_testis.all.tsv.gz",
	# "GTEx_47_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_thyroid.all.tsv.gz",
	# "GTEx_47_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_thyroid.all.tsv.gz",
	# "GTEx_48_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_uterus.all.tsv.gz",
	# "GTEx_48_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_uterus.all.tsv.gz",
	# "GTEx_49_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/ge/GTEx_ge_vagina.all.tsv.gz",
	# "GTEx_49_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GTEx/tx/GTEx_tx_vagina.all.tsv.gz"
}

study_count = len(studies)
i = 1
for prefix, filepath in studies.items():
	print( "######################")
	print(f"Study {i}/{study_count}")
	print( "######################")
	subprocess.run(f"./00pipeline.sh cluster {filepath} {prefix}", shell=True, capture_output=False, check=True)
	i += 1

# runtime calculation
total_time = time.time() - start_time
print(f"\n\n{30*'#'}\nfinished running pipeline on {study_count} studies in {total_time/60:.2f} minutes!\n{30*'#'}")








