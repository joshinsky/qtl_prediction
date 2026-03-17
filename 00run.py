#!/net/mimer/mnt/tank/projects2/kvs_students/2026/jl_qtl_prediction/conda_env/josh_env/bin/python3

import subprocess
import time
start_time = time.time()

studies = {
	"Alasoo_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/ge/Alasoo_2018_ge_macrophage_IFNg+Salmonella.all.tsv.gz",
	"Alasoo_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Alasoo_2018/tx/Alasoo_2018_tx_macrophage_IFNg+Salmonella.all.tsv.gz",
	"BLUEPRINT_1_ge": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/ge/BLUEPRINT_PE_ge_T-cell.all.tsv.gz",
	"BLUEPRINT_1_iu": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/tx/BLUEPRINT_PE_tx_T-cell.all.tsv.gz",
	"BLUEPRINT_2_ge": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/ge/BLUEPRINT_SE_ge_monocyte.all.tsv.gz",
	"BLUEPRINT_2_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/tx/BLUEPRINT_SE_tx_monocyte.all.tsv.gz",
	"BLUEPRINT_3_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/ge/BLUEPRINT_SE_ge_neutrophil.all.tsv.gz",
	"BLUEPRINT_3_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BLUEPRINT/tx/BLUEPRINT_SE_tx_neutrophil.all.tsv.gz",
	"Bossini_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Bossini-Castillo_2019/ge/Bossini-Castillo_2019_ge_Treg_naive.all.tsv.gz",
	"Bossini_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/Bossini-Castillo_2019/tx/Bossini-Castillo_2019_tx_Treg_naive.all.tsv.gz",
	"BrainSeq_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq/ge/BrainSeq_ge_brain.all.tsv.gz",
	"BrainSeq_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq/tx/BrainSeq_tx_brain.all.tsv.gz",
	"Braineaq2_1_ge": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq2/ge/Braineac2_ge_putamen.all.tsv.gz",
	"Braineaq2_1_iu": 	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq2/tx/Braineac2_tx_putamen.all.tsv.gz",
	"Braineaq2_2_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq2/ge/Braineac2_ge_substantia_nigra.all.tsv.gz",
	"Braineaq2_2_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/BrainSeq2/tx/Braineac2_tx_substantia_nigra.all.tsv.gz",
	"CAP_1_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CAP/ge/CAP_ge_LCL_naive.all.tsv.gz",
	"CAP_1_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CAP/tx/CAP_tx_LCL_naive.all.tsv.gz",
	"CAP_2_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CAP/ge/CAP_ge_LCL_statin.all.tsv.gz",
	"CAP_2_iu": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CAP/tx/CAP_tx_LCL_statin.all.tsv.gz",
	"CommonMind_ge":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CommonMind/ge/CommonMind_ge_DLPFC_naive.all.tsv.gz",
	"commonMind_iu":	"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/CommonMind/tx/CommonMind_tx_DLPFC_naive.all.tsv.gz",
	"FUSION_1_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/FUSION/ge/FUSION_ge_adipose_naive.all.tsv.gz",
	"FUSION_1_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/FUSION/tx/FUSION_tx_adipose_naive.all.tsv.gz",
	"FUSION_2_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/FUSION/ge/FUSION_ge_muscle_naive.all.tsv.gz",
	"FUSION_2_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/FUSION/tx/FUSION_tx_muscle_naive.all.tsv.gz",
	"GENCORD_1_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/ge/GENCORD_ge_LCL.all.tsv.gz",
	"GENCORD_1_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/tx/GENCORD_tx_LCL.all.tsv.gz",
	"GENCORD_2_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/ge/GENCORD_ge_T-cell.all.tsv.gz",
	"GENCORD_2_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/tx/GENCORD_tx_T-cell.all.tsv.gz",
	"GENCORD_3_ge":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/ge/GENCORD_ge_fibroblast.all.tsv.gz",
	"GENCORD_3_iu":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GENCORD/tx/GENCORD_tx_fibroblast.all.tsv.gz",
	"GEUVADIS_ge": 		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GEUVADIS/ge/GEUVADIS_ge_LCL.all.tsv.gz",
	"GEUVADIS_tx":		"/home/local/databases/ebi_eqtl_catalogue/pub/databases/spot/eQTL/sumstats/GEUVADIS/tx/GEUVADIS_tx_LCL.all.tsv.gz"
}

for prefix, filepath in studies.items():
	subprocess.run(f"./00pipeline.sh cluster {filepath} {prefix}", shell=True, capture_output=False, check=True)


print(f"\n\n{30*'#'}\nfinished running pipeline on {len(studies)} studies!\n{30*'#'}")