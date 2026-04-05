# CTSS: Cell-type purification strategies
## 1.	Introduction to CTPS
### 1.1	Overview of CTPS
The CTPS (Cell-Type-Purification Strategie) is designed to achieve high-purity cell sorting with minimal antibody reliance, by integrating flow cytometry data with single-cell transcriptomic data.. CTPS total contains four steps: pre-training, generation of training dataset, computational training and validation of fates.

![实验流程图](assets/CTPSworkflow.png)

### 1.2 How to use CTPS
CTPS can be downloaded locally and used in R. Prior to use, the required dependent packages need to be installed.
R (version ≥ 4.3.3)
GA (version ≥ 3.2.4)
geometry (version ≥ 0.5.2)
sp (version ≥ 2.1-3)
Rmalschains (version ≥ 0.5.2)
overlap (version ≥ 0.2-10)
scales (version ≥ 1.3.0)
flowCore (version ≥ 2.14.2)
XML (version ≥ 3.99-0.16.1)
plyr (version ≥ 1.8.9)
affy (version ≥ 1.80.0)
grid (version ≥ 4.3.3)
caret (version ≥ 7.0-1)
e1071 (version ≥ 1.7.6)
pROC (version ≥ 1.18.5)
ggplot (version ≥ 3.5.1)
