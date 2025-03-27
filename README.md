# Repository for "A Competitive Disinhibitory Network for Robust Optic Flow Processing in Drosophila"

**Authors:**  
Mert Erginkaya $^{1,2}$, Tomás Cruz $^{1,3}$, Margarida Brotas $^{1,4}$, André Marques $^{1}$, Kathrin Steck $^{5}$, Aljoscha Nern $^{6}$, Filipa Torrão $^{1}$, Nélia Varela $^{1}$, Davi D. Bock $^{7}$, Michael Reiser $^{6}$, and M. Eugenia Chiappe $^{1,8,*}$

$^{1}$ Champalimaud Neuroscience Programme, Champalimaud Centre for the Unknown, Lisbon, Portugal  
$^{2}$ Neurobiology and Genetics, Theodor-Boveri-Institute, Biocenter, Julius-Maximilians-University of Würzburg, Würzburg, Germany  
$^{3}$ Friedrich Miescher Institute for Biomedical Research, and Biozentrum, Department of Cell Biology, University of Basel, Basel, Switzerland  
$^{4}$ CEDOC, iNOVA4Health, NOVA Medical School, Universidade Nova de Lisboa, Lisbon, Portugal  
$^{5}$ Center for Teacher Education and School Research, Universität Leipzig, Leipzig, Germany  
$^{6}$ Janelia Research Campus, Howard Hughes Medical Institute, Ashburn, VA  
$^{7}$ University of Vermont Larner College of Medicine, Burlington, VT  
$^{8}$ Lead contact

$^{*}$ Correspondence should be addressed to [eugenia.chiappe@neuro.fchampalimaud.org](mailto:eugenia.chiappe@neuro.fchampalimaud.org)

**Published in:**  
*Journal Name, Year*

**DOI:** [10.xxxx/xxxxxxx](link_to_doi)

---

## Overview
This repository contains the accompanying code and raw data associated with the publication “A Competitive Disinhibitory Network for Robust Optic Flow Processing in Drosophila.” It includes MATLAB scripts for data analysis, as well as scripts to download large data files from Zenodo.

---

## 📁 Repository Structure

```
| data/
| ├── raw/                                       # 
| └── Preprocessed Data/                         # For MATLAB  
|
| results/
| ├── matlab_outputs/                            # Outputs from MATLAB script           
| |
| └── r_outputs/                                 # Outputs from R scripts
|      ├── diads_triads_analysis/                  # Outputs from the Diads and Triads analysis
|      ├── inputs_analysis/                        # Outputs from the Inputs analysis
|      ├── outputs_analysis/                       # Outputs from the Outputs analysis
|      └── rootpoint_analysis/                     # Outputs from the rootpoint analysis
|
| scripts/
| ├── matlab/                                    # Matlab scripts
| |    ├── CaImagingPopulationAnalysis240808.m     # Main analysis and plotting script
| |    └── helpers/                                # Helper functions for the main script
| |
| └── r/                                         # Folder Containing .R scripts used for figures 
|
| download_data_linux.sh                         # Script to download data on Linux/macOS
| download_data_windows.bat                      # Script to download data on Windows
| README.md
| License
```

---

## 🛠️ Instructions

### Step 1: Clone the Repository

You can clone the repository in two ways:

#### 1.1 Command Line

```bash
git clone https://github.com/ChiappeLab/Erginkaya_et_al_2025.git
cd Erginkaya_et_al_2025
```

#### 1.2 GitHub Desktop

1. Open GitHub Desktop.  
2. Click **File** → **Clone repository**.  
3. Select the repository from your list or paste the URL (above).  
4. Choose the local path where you want to clone it.

---

### Step 2: Download large Data

Large data files are hosted externally in Zenodo and can be downloaded using the provided scripts or directly from the Zenodo repository and added to the data folder from the following link:
- DOI: [10.5281/zenodo.14967806](10.5281/zenodo.14967806)

#### 2.1 Linux/macOS

1. Make sure you have `curl` or `wget` installed (for most Linux distributions, they are pre-installed).
2. Make sure you have `unzip` installed if you are downloading and extracting zipped folders.

   - **Debian/Ubuntu**:
     ```bash
     sudo apt-get update
     sudo apt-get install unzip -y
     ```
   - **macOS (Homebrew)**:
     ```bash
     brew install unzip
     ```

3. Navigate to the Repo directory and make the script executable and run it:
   ```bash
   chmod +x download_data_linux.sh
   ./download_data_linux.sh
   ```
   This will download and place the files into the appropriate folders (e.g., `data/Preprocessed Data/`).

#### 2.2 Windows

1. Double click the `download_data_windows.bat` press any key if prompted and wait for Download to finish

This will download and place the files into the appropriate folders (e.g., e.g., `data/Preprocessed Data/`).

---

### Step 3: Run the Main Script in MATLAB

Once the data is downloaded and placed into the correct folders, you can open MATLAB, add the entire repository folder to PATH and run the analysis script:

```matlab
run('scripts/matlab/CaImagingPopulationAnalysis240808.m')
```

This script will handle preprocessing, analysis, and plotting of results as described in the publication.

---

## 🐝 License

This work is licensed under the [MIT License](License).

---

## 📧 Contact

For any questions or issues, please contact:  
M. Eugenia Chiappe: [eugenia.chiappe@neuro.fchampalimaud.org](mailto:eugenia.chiappe@neuro.fchampalimaud.org)

