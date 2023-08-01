library(plyr)
library(dplyr)
library(tibble)
library(stringr)
library(Seurat)
library(tidyr)
library(ggplot2)
library(data.table)
library(patchwork) # plot_annotation function
library(gridExtra)  # for 'grid.arrange' function
library(ggrepel)

## Cell type fraction. Look at line91: https://github.com/satijalab/seurat/issues/962
#library(dittoSeq)  # bioconductor - it requires "nloptr" CRAN package     
########################################################################################
dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(dir); print(dir)  # "/Users/lees18/Library/CloudStorage/OneDrive-UPMC/S80d_BRCAERpos_scRNAseq_Website/11_OutgridLinedotplot_Correl_PathwayActivity_GeneExp_ColorTAMFract"
SeuratObjFile<-"/Users/lees18/Library/CloudStorage/OneDrive-UPMC/S80_BRCA_scRNAseq_Macrophage_S100A/10_1f_IntegratedSeuratObj_SCTMergeHarmony_OnlyPreTreat/SeuratObj_HarmonyGSE176EGA6608_CellTypeMacroTcell_RNAScaled_Res1.5PC30KP15_20230708.rds"   # I need a RNA scaled object. 
MacroFractionFile<-"/Users/lees18/Library/CloudStorage/OneDrive-UPMC/S80_BRCA_scRNAseq_Macrophage_S100A/10_1f_IntegratedSeuratObj_SCTMergeHarmony_OnlyPreTreat/MacrophageFractionHighMidLow_ERpos29s_FromWholecellRes1.5_20230708.txt"   
#CytoskeletonGeneFile <- "/Users/lees18/OneDrive - UPMC/S80_Macrophage_ImmunoTherapy/5_1d_CytoskeletonGeneList/CytoskeletonGene_118g.txt"
PathwayActivityFile <- "/Users/lees18/Library/CloudStorage/OneDrive-UPMC/S80_BRCA_scRNAseq_Macrophage_S100A/13_1_PROGENy_scRNAseq_WholeCell/PROGENy_PathwayActivity_GSE176EGA6608_WholeCellScaledData_Top100Perm1000_20230708.rds";
MeanSum <- "Mean"
########################################################################################
####################################################################
## Step1. Read macrophage fraction file - Macro-Rich vs Macro-Poor    # Read CellPhoneDB metadata file
####################################################################
MacrophageFractHighMidLow_Sort <- fread(MacroFractionFile, header=TRUE, stringsAsFactors=FALSE); dim(MacrophageFractHighMidLow_Sort)  #  29  15    This is only ERpos
colnames(MacrophageFractHighMidLow_Sort)[which(colnames(MacrophageFractHighMidLow_Sort)=="V5")] <- "NotAvail"
colnames(MacrophageFractHighMidLow_Sort)[1] <- "CaseID"; head(MacrophageFractHighMidLow_Sort)
#       CaseID  CD4Tcell Mesenchymal  CD8Tcell  NotAvail   Myeloid CancerEpithelial Endothelial Proliferating NormalEpithelial      Bcell Plasmablast MacrophageRatio
# 1:    AH0319 23.762376   12.360268 19.770042  1.181731 15.777707        1.9163207    0.990099     4.0881508       13.3823060 5.49345257   1.2775471       MacroRich
# 2: BIOKEY_12 27.832903    7.164120 29.450288 15.796785  5.308593        0.1389165    3.383608     3.5225243        0.0000000 4.06826751   3.3339948        MacroMid


####################################################################
## Step2. Read pySCENIC TF activity file - M
####################################################################
PathwayActivity <- readRDS(PathwayActivityFile);PathwayActivity[1:3,1:3]
PathwayActivity <- PathwayActivity %>% data.frame %>% tibble::rownames_to_column("CellID")
colnames(PathwayActivity)[2:ncol(PathwayActivity)] <- paste0(colnames(PathwayActivity)[2:ncol(PathwayActivity)], "_PathwayActivity") 
dim(PathwayActivity); PathwayActivity[1:3,1:3] # 90532cells    15 (14 pathways)
#                             CellID Androgen_PathwayActivity EGFR_PathwayActivity
# 1 GSE176_CID3586_AAGACCTCAGCATGAG                -0.850            -0.994
# 2 GSE176_CID3586_AAGGTTCGTAGTACCT                -0.562            -0.872

## Some Pathway activities are all NA ??  
table(is.na(PathwayActivity)) # All FALSE
ColumnMean <- colMeans(PathwayActivity[,2:ncol(PathwayActivity)])
table(is.na(ColumnMean)) # FALSE 14

PathwayActivity_NoNACol <- PathwayActivity    # %>% data.frame %>% dplyr::select(- names(which(is.na(ColumnMean)))) ; dim(PathwayActivity_NoNACol)  #  90532  357
####################################################################
## Step3. Read seurat object; get metadata 
####################################################################
seurat_object_ERpos_NewMeta <- readRDS(file=SeuratObjFile) # this is already normalized and scaled 
seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608 <- gsub("Monocyte/Macrophage", "Monocyte", seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608)
seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608 <- gsub("cDC|pDC|DC",   "DendriticCell", seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608)
seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608 <- gsub("Tregulatory",   "Treg", seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608)
seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608 <- gsub("CancerEpithelial",  "EpithelialCancer", seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608)

table(seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608)
# Bcell CancerEpithelial         CD4Tcell         CD8Tcell    DendriticCell      Endothelial       Fibroblast       Macrophage         Mastcell         Monocyte    Myofibroblast           NKcell      Tregulatory 
# 4365            28423            13201            10516              486             7446             9606             2994              565             1252             6100             4552             1026 
seurat_object_ERpos_NewMeta <- SetIdent(seurat_object_ERpos_NewMeta, value=seurat_object_ERpos_NewMeta$CellTypeMacroTcell_GSE176EGA6608)
table(seurat_object_ERpos_NewMeta@active.ident)  # my mistake on "Messenchymal"
# Endothelial       Fibroblast            Bcell    Myofibroblast CancerEpithelial         CD8Tcell           NKcell         CD4Tcell         Mastcell       Macrophage         Monocyte    DendriticCell      Tregulatory 
# 7446             9606             4365             6100            28423            10516             4552            13201              565             2994             1252              486             1026 
unique(seurat_object_ERpos_NewMeta$orig.ident)  # 29 CaseID levels
####################################################################
## Step4. Read cytokine interaction genes / Secretome genes
####################################################################
# CytoskeletonGene <- fread(CytoskeletonGeneFile, header=TRUE, stringsAsFactors=FALSE); CytoskeletonGene[1:3,]  # 265 genes
# colnames(CytoskeletonGene)[1] <- "Cytoskeleton"
# CytoskeletonGene$Cytoskeleton <- toupper(CytoskeletonGene$Cytoskeleton );CytoskeletonGene[1:3,]
####################################################################
## Step4. Replace @active.ident and extract each cell type  -  in Minor cell type
####################################################################
# table(seurat_object_ERpos_NewMeta$CellTypeMajor)
# table(seurat_object_ERpos_NewMeta$CellTypeMinor)
# seurat_object_ERpos_NewMeta <- SetIdent(seurat_object_ERpos_NewMeta, value=seurat_object_ERpos_NewMeta$CellTypeMinor)
# table(seurat_object_ERpos_NewMeta@active.ident)

#Idents(seurat_object_ERpos_NewMeta)  # this calls 'active.ident'
#CellType <- names(table(Idents(seurat_object_ERpos_NewMeta)))

####################################################################
## Step4. Extract Metadata  
####################################################################
MetaData_ERpos <- seurat_object_ERpos_NewMeta@meta.data; dim(MetaData_ERpos) # 90532     19
rownames(MetaData_ERpos) <- gsub("-","_", rownames(MetaData_ERpos))
##MetaData_ERpos <- MetaData_ERpos %>% tibble::rownames_to_column("CellID") %>% dplyr::mutate(CaseID=MetaData_ERpos$orig.ident)   # Cell_ID and CaseID are already set up. 
#colnames(MetaData_ERpos)[which(colnames(MetaData_ERpos)=="Cell_ID")] <- "CellID"
MetaData_ERpos[,which(colnames(MetaData_ERpos)=="Cell_ID")] <- NULL   ## Delete "Cell_ID" column because it has some NAs.  "CellID" column also has some NAs. Fix it. 
MetaData_ERpos$CellID <- rownames(MetaData_ERpos)     
MetaData_ERpos$CaseID <- MetaData_ERpos$orig.ident
#MetaData_ERpos$CellTypeMacroTcell_GSE176EGA6608 <- gsub("DC_Macrophage", "Macrophage", MetaData_ERpos$CellTypeMacroTcell_GSE176EGA6608)     # ##  <<<<<<=======  ## change the cell type name: DC_Macrophage => Macrophage  
  
MetaData_ERpos[1:2,]        # I have "CaseID" column already
#                           orig.ident nCount_RNA nFeature_RNA  CaseID percent.mt CellTypeMajor     CellTypeMinor    CellTypeSubset RNA_snn_res.0.1 seurat_clusters CellTypeMacrophage                  Cell_ID DatasetID
# CID3586_AAGACCTCAGCATGAG    CID3586       4581         1689 CID3586   1.506221   Endothelial Endothelial_ACKR1 Endothelial_ACKR1               1               8        Endothelial CID3586_AAGACCTCAGCATGAG GSE176078
# CID3586_AAGGTTCGTAGTACCT    CID3586       1726          779 CID3586   5.793743   Endothelial Endothelial_ACKR1 Endothelial_ACKR1               1               8        Endothelial CID3586_AAGGTTCGTAGTACCT GSE176078

table(rownames(MetaData_ERpos)[80000:80050] == seurat_object_ERpos_NewMeta$Cell_ID[80000:80050])    # In "Cell_ID" column, there are "NA" values. 
#rownames(MetaData_ERpos) <- MetaData_ERpos$CellID
rownames(MetaData_ERpos)[60000:60050] 
#fwrite(MetaData_ERpos, file="Metadata_GSE176EGA6608_90532cell.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE )

################################################################################
## Step5. inner_join between MetaData_ERpos and MacrophageFractHighMidLow_Sort
################################################################################
###  rownames_to_column for "CaseID"
#MacrophageFractHighMidLow_CaseIDCol <- MacrophageFractHighMidLow_Sort %>% tibble::rownames_to_column("CaseID")

## inner_join between MetaData and MacroFraction
MetaData_ERpos_MacroPoorMidRich <- dplyr::inner_join(MetaData_ERpos, MacrophageFractHighMidLow_Sort) %>%   # by "CaseID"
        dplyr::select(-c(orig.ident, nCount_RNA, nFeature_RNA, percent.mt, nCount_SCT, nFeature_SCT, CellTypeMajor, CellTypeMinor,CellTypeSubset, seurat_clusters,
                         SCT_snn_res.1.5,  BC_type, cellType, CellTypeByMarker_GSE176EGA6608, DatasetID)); 
dim(MetaData_ERpos_MacroPoorMidRich) # 90532     16    
table(MetaData_ERpos_MacroPoorMidRich$MacrophageRatio)
# MacroMid MacroPoor MacroRich 
#  47200     26352     16980 

## Replace Myeloid with Macrophage in celltype_major, and make a new column, "celltype_Macrophage"
#MetaData_ERpos_MacroPoorMidRich$celltype_Macrophage <- ifelse(MetaData_ERpos_MacroPoorMidRich$CellTypeByMarker_GSE176EGA6608 == "Myeloid", "Myeloid", MetaData_ERpos_MacroPoorMidRich$CellTypeByMarker_GSE176EGA6608)

table(MetaData_ERpos_MacroPoorMidRich$CellTypeMacroTcell_GSE176EGA6608)
# Bcell    CancerEpithelial            CD4Tcell            CD8Tcell                 cDC                  DC         Endothelial          Fibroblast          Macrophage            Mastcell 
# 4365               28423               14227               10516                 266                 175                7446                9606                3340                 565 
# Monocyte/Macrophage       Myofibroblast              NKcell                 pDC 
#               906                6100                4552                  45
MetaData_ERpos_CellType_MacroPoorMidRich <- MetaData_ERpos_MacroPoorMidRich %>% dplyr::mutate(celltype_MacroFraction = paste0(MetaData_ERpos_MacroPoorMidRich$CellTypeMacroTcell_GSE176EGA6608, "_", MetaData_ERpos_MacroPoorMidRich$MacrophageRatio))
table(MetaData_ERpos_CellType_MacroPoorMidRich$celltype_MacroFraction); dim(MetaData_ERpos_CellType_MacroPoorMidRich)  # 90532     17
# BasalEpithelial_MacroMid  BasalEpithelial_MacroPoor  BasalEpithelial_MacroRich             Bcell_MacroMid            Bcell_MacroPoor            Bcell_MacroRich  CancerEpithelial_MacroMid CancerEpithelial_MacroPoor 
# 1603                       1078                        333                        406                       3138                        246                      17749                      11265 
# CancerEpithelial_MacroRich          CD4Tcell_MacroMid         CD4Tcell_MacroPoor         CD4Tcell_MacroRich       Endothelial_MacroMid      Endothelial_MacroPoor      Endothelial_MacroRich          Mastcell_MacroMid 
# 10318                       6423                      27204                       5291                       4138                       3694                        833                        129 
# Mastcell_MacroPoor         Mastcell_MacroRich       Mesenchymal_MacroMid      Mesenchymal_MacroPoor      Mesenchymal_MacroRich           Myeloid_MacroMid          Myeloid_MacroPoor          Myeloid_MacroRich 
# 555                        154                       5640                       5841                       3562                       2418                       3166                       2424 
# Plasmablast_MacroMid      Plasmablast_MacroPoor      Plasmablast_MacroRich     Proliferating_MacroMid    Proliferating_MacroPoor    Proliferating_MacroRich 
# 535                       1735                         87                        710                        790                        319 

MetaData_ERpos_CellType_MacroPoorMidRich[1:2,]
#     CaseID CellTypeByMarker_GSE176EGA6608 CellTypeMacrophage_GSE176EGA6608 BasalEpithelial   Bcell CancerEpithelial CD4Tcell Endothelial Macrophage  Mastcell Mesenchymal  Myeloid Plasmablast Proliferating MacrophageRatio
# 1 CID3586                    Endothelial                      Endothelial        5.582565 5.11316         6.269908 75.27242    1.626153   1.005868 0.1173512    2.112322 2.045264   0.3017603     0.5532272       MacroPoor
# 2 CID3586                    Endothelial                      Endothelial        5.582565 5.11316         6.269908 75.27242    1.626153   1.005868 0.1173512    2.112322 2.045264   0.3017603     0.5532272       MacroPoor

################################################################################
## Step6. calculate celltype fraction for each sample
################################################################################
## ===== ### Method1.   count the cells per cell type and per patients
CellTypeCount <- table(MetaData_ERpos_CellType_MacroPoorMidRich$CaseID, MetaData_ERpos_CellType_MacroPoorMidRich$CellTypeMacroTcell_GSE176EGA6608)  ## this is the code to get cell type fraction
## Combine cDC, pDC, and DC.   Change "Monocyte/Macrophage" to "Monocyte"
colnames(CellTypeCount)[colnames(CellTypeCount)=="Monocyte/Macrophage"] <- "Monocyte"
CellTypeCount_DF <- data.frame(matrix(CellTypeCount, nrow=nrow(CellTypeCount)))
rownames(CellTypeCount_DF) <- rownames(CellTypeCount); colnames(CellTypeCount_DF) <- colnames(CellTypeCount)
#CellTypeCount_DF <- CellTypeCount_DF %>% dplyr::mutate(DendriticCell = cDC + pDC + DC) %>% dplyr::select(-c(cDC, pDC, DC))
dim(CellTypeCount_DF); CellTypeCount_DF   # 29 12
#           BasalEpithelial Bcell CancerEpithelial CD4Tcell Endothelial Macrophage Mastcell Mesenchymal Myeloid Plasmablast Proliferating
# BIOKEY_12              18   398             1977     5548         344         34      150         370     503         311           146
# BIOKEY_17             157    13             2808      259         153        217       19        1580     127          15           312
# BIOKEY_18              87   390              649     1994         731         55       50         646     278          70            57

CellTypeRatio <- CellTypeCount_DF/rowSums(CellTypeCount_DF); class(CellTypeRatio)  # table
rowSums(CellTypeRatio)

CellTypeRatio_CaseID <- data.frame(rbind(CellTypeRatio)) %>% tibble::rownames_to_column("CaseID")
colnames(CellTypeRatio_CaseID) <- gsub("\\.","_", colnames(CellTypeRatio_CaseID))

## ===== ### Method2.  https://www.r-bloggers.com/2020/05/calculating-ratios-with-tidyverse/
#MetaData_ERpos_CellTypeRatio <- MetaData_ERpos_CellType_MacroPoorMidRich %>% group_by(CaseID, celltype_MacroFraction) %>% dplyr::summarise(N=n(), .groups="keep") %>%    
             # group_by(CaseID) %>%   dplyr::mutate(TotalCellCount = sum(N),  Ratio = round(N/TotalCellCount, 3))
##   `summarise()` has grouped output by 'CaseID'. You can override using the `.groups` argument.
#head(MetaData_ERpos_CellTypeRatio); class(MetaData_ERpos_CellTypeRatio)  # tbl    #  120  5
# A tibble: 6 × 5
# Groups:   CaseID [1]
#   CaseID  celltype_Macrophage     N TotalCellCount Ratio
#   <chr>   <chr>               <int>          <int> <dbl>
# 1 CID3586 B-cells               321           6178 0.052
# 2 CID3586 CAFs                  185           6178 0.03 
# 3 CID3586 Endothelial           157           6178 0.025

################################################################################
## Step7. Extract gene scaled expression data from seurat object and calcualte mean expression for cytokine genes. 
################################################################################
## Count Data 
# CountExp <- readRDS("/Users/lees18/Library/CloudStorage/OneDrive-UPMC/S80_Macrophage_ImmunoTherapy/10_1f_IntegratedSeuratObj_SCTMergeHarmony_OnlyPreTreat/CountData_MacroWhole_MainCellType_80390c19229g_Tp.rds")
# dim(CountExp)  # 80390c 19229g

## Get scaled expression matrix. 
#HighVarGene <- seurat_object_ERpos_NewMeta@assays$SCT@var.features;length(HighVarGene)  # 2000
#seurat_object_ERpos_Scaled <- ScaleData(seurat_object_ERpos_NewMeta)   # features=HighVarGene
ScaledExp_GSE176EGA6608 <- seurat_object_ERpos_NewMeta@assays[["RNA"]]@scale.data; 
colnames(ScaledExp_GSE176EGA6608) <- gsub("-","_", colnames(ScaledExp_GSE176EGA6608)); dim(ScaledExp_GSE176EGA6608); ScaledExp_GSE176EGA6608[1:3,1:3]  #  33186 90532
#               GSE176_CID3586_AAGACCTCAGCATGAG GSE176_CID3586_AAGGTTCGTAGTACCT GSE176_CID3586_ACCAGTAGTTGTGGCC
# RP11-34P13.7                     -0.02307934                     -0.02307934                     -0.02307934
# FO538757.3                       -0.02069071                     -0.02069071                     -0.02069071

## Checiing CellID
#table(rownames(CountExp) %in% colnames(ScaledExp_GSE176EGA6608)) # TRUE 80390. 
table(PathwayActivity$CellID %in% colnames(ScaledExp_GSE176EGA6608)) # TRUE 90532 

# InterestGene <- c("S100A11","S100A13","S100A8", "S100A10", "S100A4", "S100A5", "S100A3", "ITGA7","ACTG1", "ITGB3", "RAC2","DIAPH1","ITGAL", 
#                 "SOD2","IRF1","VEGFA","FLT1","CD40","TGFBR2","CCL17","CCL19","TNFRSF13C","TNFRSF25","DSTN","WASL","ARF1","S100P","MPRIP","GIT1","RAF1")
# length(InterestGene); # 30
# table(rownames(ScaledExp_GSE176EGA6608) %in% InterestGene) # TRUE: 30
# table(InterestGene %in% rownames(ScaledExp_GSE176EGA6608)); # TRUE 30


########## ========== #### Remove non-genes ####### ======== ###########
# ScaleExp_GSE176EGA6608_Subset30g <- subset(ScaledExp_GSE176EGA6608, rownames(ScaledExp_GSE176EGA6608) %in% "S100A11" )
ScaleExp_GSE176EGA6608_Subset30g <- subset(ScaledExp_GSE176EGA6608, ! grepl("^LOC|^LIN|\\.|^AC\\d|^AL\\d|^AP|\\-" ,  rownames(ScaledExp_GSE176EGA6608) ) ); dim(ScaleExp_GSE176EGA6608_Subset30g) # 18625 90532
rownames(ScaleExp_GSE176EGA6608_Subset30g) <- paste0(rownames(ScaleExp_GSE176EGA6608_Subset30g), "_GeneExp")
dim(ScaleExp_GSE176EGA6608_Subset30g); ScaleExp_GSE176EGA6608_Subset30g[1:3,1:3]  #  30 90532 # [1] 33186 90532   # aftr remove ^LOC, ^LIN and "." genes;  20277 90532


## Checking CellID
table(PathwayActivity$CellID %in% colnames(ScaleExp_GSE176EGA6608_Subset30g)) # TRUE 90532 

## ========= ######### ============== Filter out genes of low variance ========= ############## =============== ######
# vars <- apply(ScaledExp_CytoskeletonGene, 1, var)
# saveRDS(vars, "CalculateVariance_AcrossAllCell.rds")
# HighVarGene <- vars > quantile(vars, (nrow(ScaledExp_CytoskeletonGene) - 10000)/nrow(ScaledExp_CytoskeletonGene))
# table(HighVarGene) # FALSE 28374  TRUE 10000
# ScaledExp_HighVarGene <- ScaledExp_CytoskeletonGene[HighVarGene,]  # This takes long
# dim(ScaledExp_HighVarGene) # 10000 121784
# saveRDS(ScaledExp_HighVarGene, file="ScaledExp_HighVarGene_10000g.rds")
## ========= ######### ============== ############# =======  ######## ========= ############## =============== ######

ScaledExp_GSE176EGA6608_Tp <- ScaleExp_GSE176EGA6608_Subset30g  %>% t %>% data.frame 
#saveRDS(GSE176EGA6608_ScaledExp_Tp, file="ScaledExp_CytoskeletonGene_Tp_HighVar10kg.rds")
ScaledExp_GSE176EGA6608_Tp[1:2,1:5]  
#                                 TNFRSF25_GeneExp S100A10_GeneExp S100A11_GeneExp
# GSE176_CID3586_AAGACCTCAGCATGAG        -0.244133       1.4583165       0.6099090
# GSE176_CID3586_AAGGTTCGTAGTACCT        -0.244133       0.8460192       1.0453171

#rownames(GSE176EGA6608_ScaledExp_Tp) <- MetaData_ERpos$CellID   #  239794  253
## Let's remove genes that have low expression in 239,794 cells
#GeneExpVar <- apply(GSE176EGA6608_ScaledExp_Tp, 2, var)
#names(GeneExpVar)[GeneExpVar<0.005]  # "INHBC" "CCL1"  "IFNK"  "GH2"   "IFNA5" "IFNA2"
#ScaledExp_CytoskeletonGene_RmvGene <- GSE176EGA6608_ScaledExp_Tp %>% dplyr::select (-names(GeneExpVar)[GeneExpVar<0.005] ); dim(ScaledExp_CytoskeletonGene_RmvGene)   ## 121784 244 

ScaledExp_GSE176EGA6608_CellType <- ScaledExp_GSE176EGA6608_Tp %>% tibble::rownames_to_column("CellID") %>% dplyr::inner_join(MetaData_ERpos[,c("CellID", "CaseID", "CellTypeMacroTcell_GSE176EGA6608")]) # %>% # inner_join by "CellID"
                               # dplyr::select(-CellID)
table(ScaledExp_GSE176EGA6608_CellType$CaseID); dim(ScaledExp_GSE176EGA6608_CellType);   # 90532    33
# BIOKEY_12 BIOKEY_17 BIOKEY_18 BIOKEY_20 BIOKEY_21 BIOKEY_22 BIOKEY_24 BIOKEY_27 BIOKEY_29  BIOKEY_3 BIOKEY_30  BIOKEY_4  BIOKEY_5  BIOKEY_6  BIOKEY_7   CID3586   CID3941   CID3948   CID3963   CID4040   CID4066   CID4067 
# 9799      5660      5007      3221      6117      2329      3911      3870      1335      5162      6728      8190      3970      7370      2315      5965       566      2226      3399      2360      4571      3207 
# CID4290A   CID4398   CID4461   CID4463   CID4471  CID4530N   CID4535 
# 4518      4297       433      1055      7079      3732      3392 
#GSE176EGA6608_ScaledExp_CellType$CellTypeMacroTcell_GSE176EGA6608 <- gsub("Myeloid", "MonoMacro", ScaledExp_CytoskeletonGene_CellType$CellTypeMacroTcell_GSE176EGA6608)

####### =========  Let's keep the 90532 cells of PathwayActivity_NoNACol in ScaledExp data too. ========= ###########
ScaledExp_GSE176EGA6608_CellType_90532c <- ScaledExp_GSE176EGA6608_CellType %>% dplyr::filter(CellID %in% PathwayActivity$CellID); dim(ScaledExp_GSE176EGA6608_CellType_90532c) # 90532c 10g

## Checking CellID
table(ScaledExp_GSE176EGA6608_CellType_90532c$CellID %in% PathwayActivity_NoNACol$CellID)  # TRUE 90532
table(PathwayActivity_NoNACol$CellID %in% ScaledExp_GSE176EGA6608_CellType_90532c$CellID)  # TRUE 90532

## subset cytokine gene expression data by cell type
#MyCellTypeSubset <- names(table(MetaData_ERpos_MacroPoorMidRich$CellTypeMacroTcell_GSE176EGA6608)); print(MyCellTypeSubset); length(MyCellTypeSubset)
#MyCellType_NoNA <- MyCellTypeSubset[ -which(MyCellTypeSubset == "NA") ];print(MyCellType_NoNA)
# [1] "BasalEpithelial"  "Bcell"            "CancerEpithelial" "CD4Tcell"         "Endothelial"      "Macrophage"       "Mastcell"         "Mesenchymal"      "Myeloid"          "Plasmablast"      "Proliferating"   
#MyCellType_NoNA <- gsub("Myeloid","MonoMacro",MyCellType_NoNA )

MyCellTypeSubset <- unique(ScaledExp_GSE176EGA6608_CellType_90532c$CellTypeMacroTcell_GSE176EGA6608); print(MyCellTypeSubset)  # 13
# [1] "Endothelial"      "Fibroblast"       "Bcell"            "Myofibroblast"    "CancerEpithelial" "CD8Tcell"         "NKcell"           "CD4Tcell"         "Mastcell"         "Macrophage"       "Monocyte"         "DendriticCell"   

# 
# CellTypeNumb<-0;  All_CorrTestSummary_ByCellType <- data.frame
# for(EachCellTypeSubset in MyCellTypeSubset) {
#         # EachCellTypeSubset <- MyCellTypeSubset[6]
#         CellTypeNumb <- CellTypeNumb+1
#         print(paste0("currnt cell type: ", EachCellTypeSubset))
#         
#         MetaData_MacroPoorMidRich_CelltypeSubset <- ScaledExp_GSE176EGA6608_CellType_90532c %>% dplyr::filter(CellTypeMacroTcell_GSE176EGA6608==EachCellTypeSubset) # %>% dplyr::select(-PreOnTreat);
#         dim(MetaData_MacroPoorMidRich_CelltypeSubset); MetaData_MacroPoorMidRich_CelltypeSubset[1:2,]  # 7446   18 (30 genes, first col is CellID, last column is CaseID)
#         #                             CellID S100A11_GeneExp  CaseID CellTypeMacroTcell_GSE176EGA6608
#         # 1 GSE176_CID3586_AGGGATGTCGTTACAG       -1.526304 CID3586                           NKcell
#         # 2 GSE176_CID3586_CCGTACTTCCTTTCTC       -1.526304 CID3586                           NKcell
#         
#         ################################################################################
#         ## Step8. inner_join between MetaData_MacroPoorMidRich_CelltypeSubset and Pathway activity  by CellID
#         ################################################################################
#         ## Checking CellIDs between scaeled exp data and PathwayActivity_NoNACol data. 
#         table(MetaData_MacroPoorMidRich_CelltypeSubset$CellID %in% PathwayActivity_NoNACol$CellID) # TRUE 4552
#         table(PathwayActivity_NoNACol$CellID %in% MetaData_MacroPoorMidRich_CelltypeSubset$CellID) # FALSE: 85980  TRUE 4552
#         
#         CelTypeFract_PathwayActivity <- dplyr::inner_join(MetaData_MacroPoorMidRich_CelltypeSubset, PathwayActivity_NoNACol); dim(CelTypeFract_PathwayActivity) #  4552   18
#         colnames(CelTypeFract_PathwayActivity) <- gsub("\\.","",colnames(CelTypeFract_PathwayActivity))
#         CelTypeFract_PathwayActivity[1:2,]   ## first column is CaseID, 4th column is CellID
#         #     CaseID CellTypeMacroTcell_GSE176EGA6608                          CellID    Bcell CancerEpithelial TNFa_PathwayActivity Trail_PathwayActivity VEGF_PathwayActivity WNT_PathwayActivity
#         # 1 CID3586                      Endothelial GSE176_CID3586_AAGACCTCAGCATGAG 6.011107         9.506697                0.664                -0.004               -0.392              -0.746
#         # 2 CID3586                      Endothelial GSE176_CID3586_AAGGTTCGTAGTACCT 6.011107         9.506697               -0.984                 0.478               -0.018              -0.702
#         
#         ################################################################################
#         ## Step9. Calculate Mean Pathway activity in each CaseID. ========== #########
#         ################################################################################
#         MyCaseID <- unique(MetaData_ERpos_CellType_MacroPoorMidRich$CaseID); length(MyCaseID) #  29 
#         # All_CorrTestSummary_ByCaseID <- data.frame(matrix(ncol=1, nrow=length(MyCytokine))); CaseIDCount <- 0;
#         All_CaseMeanPathwayActivity <- data.frame();  CaseIDCount <-0;
#         for (EachCaseID in MyCaseID) { 
#                 # EachCaseID <- MyCaseID[1]
#                 print(paste0("Each CaseID : ", EachCaseID))  # S100A10_GeneExp
#                 CaseIDCount=CaseIDCount+1;
#                 print(paste0("CaseID count: ", CaseIDCount)) 
#                 
#                 #Subset by CaseID
#                 CelTypeFract_PathwayActivity_ByCaseID <- CelTypeFract_PathwayActivity %>% dplyr::filter(CaseID == EachCaseID); dim(CelTypeFract_PathwayActivity_ByCaseID) # 149  31
#                 
#                 ## Calcualte mean pathway activity per CaseID
#                 Celltype_MeanPathwayActivity <- data.frame(t(colMeans(CelTypeFract_PathwayActivity_ByCaseID[, c(2, 5:ncol(CelTypeFract_PathwayActivity_ByCaseID))  ] ) ))
#                 rownames(Celltype_MeanPathwayActivity) <- EachCaseID
#                 All_CaseMeanPathwayActivity <- rbind(All_CaseMeanPathwayActivity, Celltype_MeanPathwayActivity)
#         }
#         dim(All_CaseMeanPathwayActivity); All_CaseMeanPathwayActivity[1:3,1:3]
#         All_CaseMeanPathwayActivity_CaseID <- All_CaseMeanPathwayActivity %>% tibble::rownames_to_column("CaseID")
#         
#         ## Macrophage fraction will be Corrx
#         CorrX <- All_CaseMeanPathwayActivity_CaseID[, "S100A11_GeneExp"]; length(CorrX) # 29
#         
#         MyPathway <- grep("_PathwayActivity", colnames(All_CaseMeanPathwayActivity_CaseID), value=TRUE); length(MyPathway) # 14 
#         All_CorrTestSummary_Pathway <- data.frame(); SubLoopCount=0;
#         for(EachPathway in MyPathway) {
#                   # EachPathway <- MyPathway[1]
#                   print(paste0("Each Pathway : ", EachPathway))  # ALX4_PathwayActivity
#                   SubLoopCount <- SubLoopCount + 1; 
#                   print(paste0("current subloopnumb: ", SubLoopCount))
#                   CorrY <- All_CaseMeanPathwayActivity_CaseID[, EachPathway];
#                   
#                   ## Corr test
#                   CorrTest <- cor.test(CorrX, CorrY, method=c("spearman"), exact=FALSE)
#                   # plot(CorrX, CorrY)
#                   
#                   if(is.na(CorrTest$estimate)) {
#                     CorrTest$estimate <- 0
#                     CorrTest$p.value <- 1
#                   }
#                   
#                   CorrTestSummary_Pathway <- data.frame(t(c(CorrTest$estimate, CorrTest$p.value)))
#                   colnames(CorrTestSummary_Pathway) <- c(paste0(EachCellTypeSubset,"_Rho"),paste0(EachCellTypeSubset,"_pval"))
#                   rownames(CorrTestSummary_Pathway) <- EachPathway
#                   
#                   All_CorrTestSummary_Pathway <- rbind(All_CorrTestSummary_Pathway, CorrTestSummary_Pathway)
#                   
#         }  # end of For loop
#         
#         if(CellTypeNumb ==1) {
#           All_CorrTestSummary_ByCellType <- All_CorrTestSummary_Pathway
#         } else if (CellTypeNumb > 1) {
#           All_CorrTestSummary_ByCellType<-cbind(All_CorrTestSummary_ByCellType,  All_CorrTestSummary_Pathway)
#         }
#         
# } # end of the CellType FOR loop
# dim(All_CorrTestSummary_ByCellType) #  14 24 



################ ================= ################ ================= ################ ================= ################ ================= 
## Make correaltion outgrid dotplot between Pathway activity and S100 gene expression in cancer epithelial cell type. - Color dots per macrophage fraction 
################ ================= ################ ================= ################ ================= ################ ================= 
CellTypeNumb<-0; 
for(EachCellTypeSubset in MyCellTypeSubset) {
      # EachCellTypeSubset <- MyCellTypeSubset[6]  # 5 is  "CancerEpithelial"
      CellTypeNumb <- CellTypeNumb+1
      print(paste0("currnt cell type: ", EachCellTypeSubset))
      
      ## Subset by cell type.
      ScaledExp_GSE176EGA6608_CellTypeSubset <- ScaledExp_GSE176EGA6608_CellType_90532c %>% dplyr::filter(CellTypeMacroTcell_GSE176EGA6608==EachCellTypeSubset) %>% dplyr::select(-CellTypeMacroTcell_GSE176EGA6608); # 28423 3
      dim(ScaledExp_GSE176EGA6608_CellTypeSubset); ScaledExp_GSE176EGA6608_CellTypeSubset[1:3,1:3]  # 28423    32 (30 genes, first col is CellID, last column is CaseID)
      # dim(ScaledExp_GSE176EGA6608_CellTypeSubset); ScaledExp_GSE176EGA6608_CellTypeSubset[1:2,c(1:3,28:32)]  # 28423    32 (30 genes, first col is CellID, last column is CaseID)
      #                             CellID TNFRSF25_GeneExp S100A10_GeneExp DSTN_GeneExp CD40_GeneExp RAC2_GeneExp TNFRSF13C_GeneExp  CaseID
      # 1 GSE176_CID3586_CGAACATAGGTGATTA        -0.244133       -1.167087    -1.009784    6.8146516    3.2882694       -0.09155154 CID3586
      # 2 GSE176_CID3586_CATGGCGAGCTCCTCT        -0.244133       -1.167087    -1.009784   -0.3012185   -0.5856517       -0.09155154 CID3586
      
      ################################################################################
      ## Step8. inner_join between ScaledExp_GSE176EGA6608_CellTypeSubset and Pathway activity  by CellID
      ################################################################################
      ## Checking CellIDs between scaeled exp data and PathwayActivity_NoNACol data. 
      table(colnames(ScaledExp_GSE176EGA6608) %in% PathwayActivity_NoNACol$CellID) #   TRUE 90532   # ScaledExp_GSE176EGA6608 has whole 90532 cells. 
      table(ScaledExp_GSE176EGA6608_CellTypeSubset$CellID %in% PathwayActivity_NoNACol$CellID) # TRUE 28423
      table(PathwayActivity_NoNACol$CellID %in% ScaledExp_GSE176EGA6608_CellTypeSubset$CellID) # FALSE: 62109  TRUE 28423
      
      ScaleExp_PathwayActivity <- dplyr::inner_join(ScaledExp_GSE176EGA6608_CellTypeSubset, PathwayActivity_NoNACol); dim(ScaleExp_PathwayActivity) # 28423    17
      ScaleExp_PathwayActivity[1:2,c(1:10)]   ## first column is CellID, 32nd column is CaseID
      #                             CellID TNFRSF25_GeneExp S100A10_GeneExp S100A11_GeneExp S100A8_GeneExp RAC2_GeneExp TNFRSF13C_GeneExp  CaseID Androgen_PathwayActivity EGFR_PathwayActivity Estrogen_PathwayActivity
      # 1 GSE176_CID3586_CGAACATAGGTGATTA        -0.244133       -1.167087       -1.526304     -0.1257482    3.2882694       -0.09155154 CID3586                   -0.660               -0.968                   -0.974
      # 2 GSE176_CID3586_CATGGCGAGCTCCTCT        -0.244133       -1.167087       -1.526304     -0.1257482   -0.5856517       -0.09155154 CID3586                   -0.998                0.370                   -0.242
      
      ################################################################################
      ## Step9. inner_join between ScaledExp_PathwayActivity and MacrophageFraction by CellID
      ################################################################################
      ScaleExp_PathwayActivity_MacroFract <- dplyr::inner_join(ScaleExp_PathwayActivity, MetaData_ERpos_CellType_MacroPoorMidRich[, c("CellID","Macrophage")]); 
      dim(ScaleExp_PathwayActivity_MacroFract); ScaleExp_PathwayActivity_MacroFract[1:2,1:3] # 28423    13
      
      MyCytokine <- grep("_GeneExp", colnames(ScaledExp_GSE176EGA6608_CellTypeSubset), value=TRUE); length(MyCytokine)  #  18625
      #MyCytokine <- grep("S100A", colnames(ScaledExp_GSE176EGA6608_CellTypeSubset), value=TRUE); length(MyCytokine)  
      # [1] "S100A10_GeneExp" "S100A11_GeneExp" "S100A8_GeneExp"  "S100A5_GeneExp"  "S100A4_GeneExp"  "S100A3_GeneExp"  "S100A13_GeneExp"

      GeneExpCount=0;
      for(EachCytoGene in MyCytokine[7000:length(MyCytokine)]) {
            if(EachCytoGene=="NA_") next;  
            # EachCytoGene <- MyCytokine[147]; print(paste0("EachCyto Gene: ", EachCytoGene)) # [1] "S100A10_GeneExp" "S100A11_GeneExp" "S100A8_GeneExp"  "S100A5_GeneExp"  "S100A4_GeneExp"  "S100A3_GeneExp"  "S100A13_GeneExp"
            print(paste0("Each cytogene: ", EachCytoGene))  # S100A10_GeneExp
            GeneExpCount=GeneExpCount+1;
            print(paste0("Gene count: ", GeneExpCount)) 
            
            #CorrX <- ScaleExp_PathwayActivity_MacroFract[, EachCytoGene]; length(CorrX) # 552
            
            MyPathway <- grep("_PathwayActivity", colnames(ScaleExp_PathwayActivity_MacroFract), value=TRUE); length(MyPathway) # 
            #MyPathway <- MyPathway[!MyPathway %in% c( )]; print(MyPathway)
            
                  MyPlot_AllCellType <- list(); PathayNumber=0;
                  for(EachPathway in MyPathway) {
                        # EachPathway <- MyPathway[1]
                        print(paste0("Each Pathway : ", EachPathway))  # ALX4_PathwayActivity
                        PathayNumber <- PathayNumber + 1; 
                        print(paste0("current PathayNumber: ", PathayNumber))
    
                        ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= 
                        ## Step10. select only EachCytoGene and EachPathay columns, including 'CellID' and 'CaseID' column
                        ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= 
                        colnames(ScaleExp_PathwayActivity_MacroFract)[colnames(ScaleExp_PathwayActivity_MacroFract) == EachCytoGene] <- "CytoGeneSelect"
                        colnames(ScaleExp_PathwayActivity_MacroFract)[colnames(ScaleExp_PathwayActivity_MacroFract) == EachPathway] <- "PathwaySelect"
                        ScaleExpPathwayActivityMacroFract_Subset <- ScaleExp_PathwayActivity_MacroFract %>% dplyr::select(c(CellID,CaseID, CytoGeneSelect,PathwaySelect, Macrophage ))
                        
                        ## Strangely, ScaledExp data has -1.167087 in so many cells. I will delete these cells. 
                        MinScaledExp <- min(ScaleExpPathwayActivityMacroFract_Subset$CytoGeneSelect); print(MinScaledExp) # -1.167087
                        if(mean(ScaleExpPathwayActivityMacroFract_Subset$CytoGeneSelect) == min(ScaleExpPathwayActivityMacroFract_Subset$CytoGeneSelect)) {
                              ## Return column names in ScaleExp_PathwayActivity_MacroFract
                              colnames(ScaleExp_PathwayActivity_MacroFract)[colnames(ScaleExp_PathwayActivity_MacroFract) == "CytoGeneSelect"] <- EachCytoGene 
                              colnames(ScaleExp_PathwayActivity_MacroFract)[colnames(ScaleExp_PathwayActivity_MacroFract) == "PathwaySelect"] <- EachPathway
                              next; 
                        }
                        
                        ScaleExpPathwayActivityMacroFract_ExcludeMinExp <- ScaleExpPathwayActivityMacroFract_Subset %>% dplyr::filter(CytoGeneSelect>MinScaledExp)
                        dim(ScaleExpPathwayActivityMacroFract_ExcludeMinExp); min(ScaleExpPathwayActivityMacroFract_ExcludeMinExp$CytoGeneSelect) # 18336     5
                        ## Calculate mean of S100A gene expression and 
                        
                        ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= 
                        ## Step11. Calculate mean of Gene expression, pathway activity per CaseID
                        ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= ###### ======= 
                        MyCaseID <- unique(ScaleExpPathwayActivityMacroFract_ExcludeMinExp$CaseID); print(MyCaseID) # 29 CaseID
                        
                        All_CaseCytokineMeanExp<- data.frame()
                        for(EachCase in MyCaseID) {
                                # EachCase <- MyCaseID[1]
                                ScaleExpPathwayActivityMacroFract_Subset <- ScaleExpPathwayActivityMacroFract_ExcludeMinExp %>% dplyr::filter(CaseID==EachCase); dim(ScaleExpPathwayActivityMacroFract_Subset)
                                ScaledExp_MeanCytoExp <- data.frame(t(colMeans(ScaleExpPathwayActivityMacroFract_Subset[,3:ncol(ScaleExpPathwayActivityMacroFract_Subset)])))
                                rownames(ScaledExp_MeanCytoExp) <- EachCase
                                All_CaseCytokineMeanExp <- rbind(All_CaseCytokineMeanExp, ScaledExp_MeanCytoExp)
                        }
                        colnames(All_CaseCytokineMeanExp)[3] <- "TAM_Infiltration"
                        dim(All_CaseCytokineMeanExp); head(All_CaseCytokineMeanExp) # 29  3
                        
                        # CorrX <- ScaleExpPathwayActivityMacroFract_ExcludeMinExp[, "CytoGeneSelect"]; length(CorrX) # 552
                        # CorrY <- ScaleExpPathwayActivityMacroFract_ExcludeMinExp[, "PathwaySelect"]; length(CorrY) # 552
                        CorrX <- All_CaseCytokineMeanExp[, "CytoGeneSelect"]; length(CorrX) # 552
                        CorrY <- All_CaseCytokineMeanExp[, "PathwaySelect"]; length(CorrY) # 552
                        
                        if(length(CorrX) < 5 | length(CorrY) < 5) {
                              ## Return column names in ScaleExp_PathwayActivity_MacroFract
                              colnames(ScaleExp_PathwayActivity_MacroFract)[colnames(ScaleExp_PathwayActivity_MacroFract) == "CytoGeneSelect"] <- EachCytoGene 
                              colnames(ScaleExp_PathwayActivity_MacroFract)[colnames(ScaleExp_PathwayActivity_MacroFract) == "PathwaySelect"] <- EachPathway
                              next;
                        }
                        
                        
                        ## Corr test
                        CorrTest <- cor.test(CorrX, CorrY, method=c("spearman"), exact=FALSE); print(CorrTest)
                        # plot(CorrX, CorrY)
                        
                        if(is.na(CorrTest$estimate)) {
                          CorrTest$estimate <- 0
                          CorrTest$p.value <- 1
                        }
                    
                        #get intercept and slope value
                        reg<-lm(formula = CorrX ~ CorrY, data=All_CaseCytokineMeanExp)
                        coeff<-coefficients(reg); intercept<-coeff[1]; slope<- coeff[2]
                        
                        # ggplot gradient color lecture:  https://www.datanovia.com/en/blog/ggplot-gradient-color/
                        # sp <- ggplot(ScaleExpPathwayActivityMacroFract_ExcludeMinExp, aes(CytoGeneSelect, PathwaySelect))+
                        #   geom_point(aes(color = Macrophage))
                        # midMac <- mean(ScaleExpPathwayActivityMacroFract_ExcludeMinExp$Macrophage); print(midMac) # 4.576187
                        # sp +  scale_color_gradient2(midpoint = midMac, low = "blue", mid = "white",high = "red", space = "Lab" )
                        
                        #https://stackoverflow.com/questions/58698075/assign-specific-color-to-definite-value-in-bar-plot-using-scale-fill-gradientn
                        IPcol <- c("blue", "#66B2FF","#CCFFFF","#E0E0E0","#FFFFCC", "#feb24c", "#fd8d3c","#FF0000","red", "#CC0000")#, "#bd0026")
                         
    
                        MidMac <-  mean(All_CaseCytokineMeanExp$TAM_Infiltration); print(MidMac) # 4.048
                        MyPlot <- ggplot(All_CaseCytokineMeanExp, aes(CytoGeneSelect, PathwaySelect))+
                          geom_point(aes(color = TAM_Infiltration), alpha=1, size=2.2, show.legend=T) +  # increase 'alpha' to make the color stronger. 1 is maximum.
                          guides(color = guide_legend(override.aes = list(size = 2.5), reverse=TRUE)) + 
                           #scale_color_manual(values=c( "blue","#302626","red"))+                # c("#ECEBEB", "red", "blue", "#302626"))+  # orange, blue, black, white grey ("#ECEBEB")
                          #scale_color_gradient2(midpoint=MidMac, low="#3333FF", mid="#FFFFCC",high="red", space="Lab")+  # color strength is weak.
                          #scale_color_gradient(low="#3333FF", high="#FF0000", space="Lab")+ 
                         # scale_color_gradientn(colours = rainbow(5))+
                          scale_colour_gradientn(colours = IPcol)+
                          #xlim(c(-5.0, 11)) + ylim(c(-2.2, 2.0)) +
                          #geom_vline(xintercept=c(0.,0),lty=4,col="black",lwd=0.8) +
                          #geom_hline(yintercept=-log(0.2, 10), lty=4, col="black", lwd=0.8) +
                          labs(color="TAM infiltration (%)", x=paste0(gsub("_(.*)","", EachCytoGene)," exp"), y=paste0(EachPathway),
                               title=paste0(EachCellTypeSubset, "\nrho=", print(round(CorrTest$estimate,3)), ", pval=", print(round(CorrTest$p.value,3)) ) ) + theme_bw() +  # title="ILC ER+/HER2- (n=13) vs IDC ER+/HER2+ (n=27)"
                          theme(#legend.title="TAM infiltration (%)", #element_blank(),  #legend.position="none",
                                # legend.title=element_text(size=7), legend.text = element_text(size=8), 
                                legend.position="none",
                                legend.margin=margin(2,2,2,2), legend.box.margin=margin(0,0,0,0),
                                axis.text.x=element_text(colour="Black",size=7,angle=0,hjust=.5,vjust=.5,face="plain"),  #face="plain" or "bold"
                                axis.text.y=element_text(colour="Black",size=7,angle=0,hjust=.5,vjust=.5,face="plain"),
                                plot.title = element_text(size=8),
                                axis.title.x=element_text(colour="black",size=8,angle=0,hjust=.5,vjust=.5,face="plain"),
                                axis.title.y=element_text(colour="black",size=8,angle=90,hjust=.5,vjust=.5,face="plain"),
                                axis.line.y=element_line(color="black", linewidth=0.7),
                                axis.line.x=element_line(color="black", linewidth =0.7),
                                axis.ticks=element_line(colour="black",linewidth=1),
                                axis.ticks.length=unit(.10, "cm"), text=element_text(size=22),
                                panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
                          # geom_abline(intercept=intercept, slope=slope, color="black",linetype="dashed", size=1.5) 
                          stat_smooth(method = "glm", formula = y ~ x, geom = "smooth", color="grey", se=FALSE)
    
                        MyPlot 
                        #OutFile<- paste0("Scatterdot_Correlation_",EachCellType2,"_",EachCytokine, "Rho",print(round(CorrTest$estimate,3)),"_p",print(round(CorrTest$p.value,3)),"_W6H5.pdf")
                        #ggsave(MyPlot, file=OutFile, width=7, height=6)
                        MyPlot_AllCellType <- c(MyPlot_AllCellType, list(MyPlot))
                        
                        ## Return column names in ScaleExp_PathwayActivity_MacroFract
                        colnames(ScaleExp_PathwayActivity_MacroFract)[colnames(ScaleExp_PathwayActivity_MacroFract) == "CytoGeneSelect"] <- EachCytoGene 
                        colnames(ScaleExp_PathwayActivity_MacroFract)[colnames(ScaleExp_PathwayActivity_MacroFract) == "PathwaySelect"] <- EachPathway
                  }  # end of For loop
         
                  length(MyPlot_AllCellType)  # 14 
                  if(length(MyPlot_AllCellType)==0) {
                      next;
                  }
                  EachCytoGene <- gsub("Gene","", EachCytoGene); print(EachCytoGene)
                  OutBoxplot <- paste0("OutgridDotplot_Correl_", EachCytoGene, "_PathwayActivity_",EachCellTypeSubset, ".png"); print(OutBoxplot)
                  #pdf(OutBoxplot, width=8.7, height=10.2)
                  png(file=OutBoxplot, width=640, height=650, pointsize=10)
                      grid.arrange(MyPlot_AllCellType[[1]], MyPlot_AllCellType[[2]], MyPlot_AllCellType[[3]],MyPlot_AllCellType[[4]], 
                                   MyPlot_AllCellType[[5]], MyPlot_AllCellType[[6]], MyPlot_AllCellType[[7]],MyPlot_AllCellType[[8]], 
                                   MyPlot_AllCellType[[9]], MyPlot_AllCellType[[10]], MyPlot_AllCellType[[11]],MyPlot_AllCellType[[12]], 
                                   MyPlot_AllCellType[[13]], MyPlot_AllCellType[[14]],   nrow=4, top="" )
                  dev.off()

          }  # End of CytoGene FOR loop
}   

