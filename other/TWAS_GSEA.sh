#!/bin/sh
#SBATCH --mem=50G
#SBATCH -t 5:0:0
#SBATCH -p shared

module load apps/R/3.6.0

TWAS_GSEA=~/PROGRAMS/opain/TWAS-GSEA/TWAS-GSEA.V1.2.R
FUSION_REP=~/data/EXTERNAL_DATA/FUSION
OUTPUT_DIR=~/results/results_TWAS_GSEA
GO_SET_FILE=~/data/EXTERNAL_DATA/MSigDB/c5.all.v7.0.symbols.gmt


for pheno in logWMHnorm_meta FA MD
do
        TWAS_DIR=~/results/results_postGWAS_${pheno}_15.10.2018/RESULTS_FUSION_GTEx_v7

        for tissue in Whole_Blood Artery_Aorta Artery_Tibial Artery_Coronary
        do
                FeaturePred_DIR=~/data/EXTERNAL_DATA/FUSION/Predicted_expression/GTEx_${tissue}
                
                count=1
                for file in ${TWAS_DIR}/*${tissue}*
                do
                        if (( $count == 1 ))
                        then
                                cat ${file} > ${OUTPUT_DIR}/TWASresults_${pheno}_${tissue}
                        else
                                sed '1d' ${file} >>  ${OUTPUT_DIR}/TWASresults_${pheno}_${tissue}
                        fi
                        count=$[$count +1]
                done

                Rscript ${TWAS_GSEA} \
                        --twas_results ${OUTPUT_DIR}/TWASresults_${pheno}_${tissue} \
                        --pos ${FUSION_REP}/GTEx_v7/${tissue}.P01.pos \
                        --gmt_file ${GO_SET_FILE} \
                        --expression_ref ${FeaturePred_DIR}/FeaturePredictions.csv \
                        --use_alt_id ID \
                        --qqplot F \
                        --output ${OUTPUT_DIR}/TWAS_GSEA_${pheno}_${tissue}

        done

        TWAS_DIR=~/results/results_postGWAS_${pheno}_15.10.2018/RESULTS_FUSION_Studies

        for tissue in CMC.BRAIN.RNASEQ YFS.BLOOD.RNAARR
        do

                FeaturePred_DIR=~/data/EXTERNAL_DATA/FUSION/Predicted_expression/${tissue}

                count=1
                for file in ${TWAS_DIR}/*${tissue}_chr*
                do
                        if (( $count == 1 ))
                        then
                        cat ${file} > ${OUTPUT_DIR}/TWASresults_${pheno}_${tissue}
                        else
                        sed '1d' ${file} >>  ${OUTPUT_DIR}/TWASresults_${pheno}_${tissue}
                        fi
                        count=$[$count +1]
                done

                Rscript ${TWAS_GSEA} \
                        --twas_results ${OUTPUT_DIR}/TWASresults_${pheno}_${tissue} \
                        --pos ${FUSION_REP}/${tissue}/${tissue}.pos \
                        --gmt_file ${GO_SET_FILE} \
                        --expression_ref ${FeaturePred_DIR}/FeaturePredictions.csv \
                        --use_alt_id ID \
                        --qqplot F \
                        --output ${OUTPUT_DIR}/TWAS_GSEA_${pheno}_${tissue}
        done
done
