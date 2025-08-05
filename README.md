# vemurafenib-transcriptome-analysis
Transcriptome analysis of BRAFV600E melanoma cells treated with vemurafenib (GSE42872)
## 🎯 Purpose
To identify differentially expressed genes and key signaling pathways in melanoma cells treated with vemurafenib, in order to understand its molecular impact and uncover potential therapeutic targets or biomarkers.

## 📁 Data Source
- **Dataset**: [GSE42872 – Gene Expression Omnibus (NCBI)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42872)
- **Platform**: GPL10558 (Illumina HumanHT-12 V4.0 expression beadchip)
- **Sample**: Vemurafenib-treated vs untreated A375 melanoma cells

## 🧪 Method Summary
- Raw data pre-processing and normalization  
- Differential expression analysis (limma)  
- KEGG pathway enrichment  
- Visualization (heatmap, barplots)

## 🔍 Main Pathways Identified
- MicroRNAs in cancer  
- Transcriptional misregulation in cancer  
- Prostate cancer  
- p53 signaling pathway  

### ⚙️ Notable Upregulated Genes
- `NUPR1`, `CD36`, `MOXD1`, `DCT`, `ST6GALNAC2`

These genes were significantly upregulated following vemurafenib treatment and are involved in key cancer-related pathways.

## 📊 Visualizations
- Differential expression heatmap  
- KEGG pathway barplots  
- Volcano plot 
