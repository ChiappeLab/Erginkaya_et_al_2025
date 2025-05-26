# **Guide to Setting Up FlyWire Analysis in RStudio**


## **1. Install R and RStudio**
Ensure you have:
- **R** (latest version) from [CRAN](https://cran.r-project.org/)
- **RStudio** from [Posit](https://posit.co/download/rstudio-desktop/)

---

## **2. Install Miniconda**
FlyWire's `fafbseg` package relies on Python dependencies, so we use Miniconda for package management.

### **Install Miniconda in R**
```r
install.packages("reticulate")
reticulate::install_miniconda(force = TRUE)
```
After installation, **restart RStudio**.

---

## **3. Find the Python Environment (`r-reticulate`)**
Check that R detects the environment:
```r
reticulate::conda_list()
```
You should see `r-reticulate` listed.

---

## **4. Install Required Python Packages**
Since some dependencies are not available in Conda, install them using `pip`.

### **Activate Conda in Terminal (CMD or Anaconda Prompt)**
```bash
conda activate r-reticulate
```

### **Install Dependencies**
```bash
pip install caveclient numpy cloud-volume
```

### Other option: Install using reticulate
```r
reticulate::py_install("cloud-volume", envname = "r-reticulate", pip = TRUE) 
```

### **Verify Installation**
```bash
python -c "import caveclient, cloudvolume, numpy; print('All dependencies installed correctly!')"
```

If you see: `All dependencies installed correctly!`, you're good to go!

---

## **5. Ensure R Uses the Correct Python**

```r
reticulate::use_condaenv("r-reticulate", required = TRUE)
reticulate::py_config()
```

You should see:
```
python: C:/Users/user/AppData/Local/r-miniconda/envs/r-reticulate/python.exe
```

Verify module availability:
```r
reticulate::py_module_available("caveclient")
reticulate::py_module_available("cloudvolume")
```
Both should return `TRUE`.

---

## **6. Install the R Packages**
For the full adult female brain (FAFB) dataset:
- Documentation can be found in: https://github.com/natverse/fafbseg

```r
install.packages("natmanager")
natmanager::install(pkgs = "fafbseg")
install.packages(c("arrow", "dplyr", "ggpubr", "ggplot2"))
```
For the Drosophila Male Adult Nerve Cord dataset (aka MANC) dataset:
- Documentation can be found in: https://github.com/natverse/malevnc
```r
install.packages("natmanager")
natmanager::install(pkgs="malevnc")
install.packages(c("arrow", "dplyr", "ggpubr", "ggplot2"))
```

---

## **7. Run Your Script**
```r
library(fafbseg)
bIPS.skid <- flywire_latestid('720575940622581173', version = 630)
```

If everything is set up correctly, this should work!

---

## **8. Troubleshooting**

### ❌ **Error: Python version conflict**
```r
.rs.restartR()
reticulate::use_python("C:/Users/user/AppData/Local/r-miniconda/envs/r-reticulate/python.exe", required = TRUE)
reticulate::py_config()
```

### ❌ **Error: `caveclient` or `cloudvolume` not found**
```r
reticulate::py_install(c("caveclient", "cloudvolume"), envname = "r-reticulate", method = "pip")
```

### ❌ **Error: `conda` is not recognized in CMD**
#### **Fix:** Add Conda to **System PATH**:
1. Press `Win + R`, type `sysdm.cpl`, and hit **Enter**.
2. Go to **Advanced** → **Environment Variables**.
3. Edit **Path** under **System Variables** and add:
   ```
   C:\Users\user\AppData\Local\r-miniconda
   C:\Users\user\AppData\Local\r-miniconda\condabin
   C:\Users\user\AppData\Local\r-miniconda\Scripts
   ```
4. Restart your computer.

---

## **9. Summary**
✅ **Steps:**
1. Install **Miniconda** and create the `r-reticulate` environment.
2. Install **Python dependencies** via `pip`.
3. Ensure **R uses the correct Python environment**.
4. Install **R packages** (`fafbseg`, `ggplot2`, `dplyr`, etc.).
5. Run your **FlyWire analysis script**.

🚀 **You’re now ready to use FlyWire in RStudio!** 🎉

---

## **10. Acknowledgements**
This work utilizes the following datasets and software:

### **malevnc: The Drosophila Male Adult Nerve Cord Dataset**
- Provides access to the **MANC** dataset, enabling analysis of neuronal structures in the male *Drosophila* nerve cord.
- Developed as part of the **Janelia FlyEM Project Team**, supporting connectomics research.
- Learn more: [malevnc GitHub](https://github.com/natverse/malevnc)
- **Citations**:
  - Bates et al. (2020), *eLife*, DOI: [10.7554/eLife.53350](https://doi.org/10.7554/eLife.53350)
  - Takemura et al. (2023), *bioRxiv*, DOI: [10.1101/2023.06.05.543757](https://doi.org/10.1101/2023.06.05.543757)
  - Marin et al. (2023), *bioRxiv*, DOI: [10.1101/2023.06.05.543407](https://doi.org/10.1101/2023.06.05.543407)
  - Cheong et al. (2023), *bioRxiv*, DOI: [10.1101/2023.06.07.543976](https://doi.org/10.1101/2023.06.07.543976)

### **fafbseg: Full Adult Female Brain (FAFB) Dataset**
- Provides tools for analyzing **FlyWire-segmented EM data**.
- Integrated with the **NeuroAnatomy Toolbox (natverse)**.
- Developed to support the **FlyWire** automated segmentation project.
- Learn more: [fafbseg GitHub](https://github.com/natverse/fafbseg)
- **Citations**:
  - Bates et al. (2020), *eLife*, DOI: [10.7554/eLife.53350](https://doi.org/10.7554/eLife.53350)
  - Dorkenwald et al. (2023), *bioRxiv*, DOI: [10.1101/2023.06.27.546656](https://doi.org/10.1101/2023.06.27.546656)
  - Schlegel et al. (2023), *bioRxiv*, DOI: [10.1101/2023.06.27.546055](https://doi.org/10.1101/2023.06.27.546055)

Development of the **natverse** including the **fafbseg** package has been supported by:
- **NIH BRAIN Initiative** (grant **1RF1MH120679-01**)
- **NSF/MRC Neuronex2** (NSF **2014862/MC_EX_MR/T046279/1**)
- **Medical Research Council** (MC_U105188491)

---